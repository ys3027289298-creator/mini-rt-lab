"""Minimal deterministic ray tracer used to produce the HDR target image.

Scene: a unit sphere on a checkerboard floor, one Lambert sky-light, and two
small spheres as specular highlights (soft mirror). Everything is analytic and
vectorised, so a frame is reproducible bit-for-bit on the same machine.
"""
import numpy as np


def _norm(v):
    n = np.linalg.norm(v, axis=-1)
    n = np.where(n > 0.0, n, 1.0)
    return v / n[..., None]


def _dot(a, b):
    return np.sum(a * b, axis=-1)


def _checker(p, scale=1.0):
    """3D checkerboard, returns 0/1 float."""
    c = np.floor(p * scale).astype(np.int64)
    return (np.sum(c, axis=-1) % 2).astype(np.float64)


def _sphere_hit(orig, dirs, c, r):
    """Return hit t (>0) and mask for an analytic sphere."""
    oc = orig - c
    b = np.sum(oc * dirs, axis=-1)
    cc = np.sum(oc * oc, axis=-1) - r * r
    disc = b * b - cc
    hit = disc > 0.0
    sq = np.sqrt(np.where(hit, disc, 0.0))
    t0 = -b - sq
    t1 = -b + sq
    t = np.where(t0 > 1e-4, t0, t1)
    hit &= t > 1e-4
    return t, hit


def _plane_hit(orig, dirs, y):
    """Floor plane y = -1, returns t and mask."""
    dy = dirs[..., 1]
    safe = np.where(np.abs(dy) < 1e-9, 1e-9, dy)
    t = (y - orig[..., 1]) / safe
    return t, (t > 1e-4) & (np.abs(dy) >= 1e-9)


def render(width=320, height=240, spp=96, seed=7, soft_shadow=True):
    """Render and return a float64 HxWx3 radiance image."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:height, 0:width].astype(np.float64)
    aspect = width / height

    # camera: look down -z toward the scene
    cam = np.array([0.0, 0.9, 4.2])
    forward = _norm(np.array([0.0, -0.15, -1.0]))
    right = _norm(np.cross(forward, np.array([0.0, 1.0, 0.0])))
    up = np.cross(right, forward)
    tan_f = np.tan(np.deg2rad(38.0) * 0.5)

    img = np.zeros((height, width, 3), dtype=np.float64)
    light = np.array([-0.35, 0.9, 0.55])
    light = light / np.linalg.norm(light)

    # sphere centres / radii
    spheres = [
        (np.array([0.0, -0.15, 0.0]), 0.85, np.array([0.85, 0.22, 0.18])),
        (np.array([1.45, -0.55, -0.9]), 0.45, np.array([0.20, 0.75, 0.30])),
        (np.array([-1.5, -0.6, -0.5]), 0.40, np.array([0.25, 0.35, 0.9])),
    ]

    # one jittered sample per (spp) pass per pixel, accumulated
    for _ in range(spp):
        jx = (rng.random((height, width)) - 0.5)
        jy = (rng.random((height, width)) - 0.5)
        px = (2.0 * (xx + 0.5 + jx) / width - 1.0) * tan_f * aspect
        py = (1.0 - 2.0 * (yy + 0.5 + jy) / height) * tan_f
        pix = px[..., None] * right + py[..., None] * up + forward
        dirs = _norm(pix)

        best_t = np.full((height, width), np.inf)
        hit_kind = np.zeros((height, width), dtype=np.int32)  # 0 none, 1 floor, 2 sphere
        hit_idx = np.full((height, width), -1, dtype=np.int32)

        t, m = _plane_hit(cam, dirs, -1.0)
        upd = m & (t < best_t)
        best_t = np.where(upd, t, best_t)
        hit_kind = np.where(upd, 1, hit_kind)

        for i, (c, r, _col) in enumerate(spheres):
            t, m = _sphere_hit(cam, dirs, c, r)
            upd = m & (t < best_t)
            best_t = np.where(upd, t, best_t)
            hit_kind = np.where(upd, 2, hit_kind)
            hit_idx = np.where(upd, i, hit_idx)

        safe_t = np.where(np.isfinite(best_t), best_t, 0.0)
        p = cam + dirs * safe_t[..., None]
        rad = np.zeros((height, width, 3), dtype=np.float64)

        # --- floor: checker albedo, N.L shading, soft shadow from spheres,
        #     plus a faint glossy reflection of the sky so it is not flat.
        floor = hit_kind == 1
        if np.any(floor):
            n = np.array([0.0, 1.0, 0.0])
            alb = 0.55 * _checker(p, 1.7)[..., None] + 0.12
            ndl = np.clip(_dot(np.broadcast_to(n, dirs.shape), light), 0.0, 1.0)
            rad_floor = alb * ndl[..., None]
            if soft_shadow:
                # sample 4 jittered light directions, shadow-ray to each sphere
                acc = np.zeros_like(rad_floor)
                lj = np.array([[-0.9, 0.0, -0.4], [0.9, 0.0, -0.2],
                               [-0.6, 0.0, 0.9], [0.7, 0.0, 0.7]])
                for k in range(lj.shape[0]):
                    ld = _norm(light + 0.12 * lj[k])
                    dirs_sh = np.broadcast_to(ld, dirs.shape)
                    vis = np.ones((height, width), dtype=np.float64)
                    for c, r, _col in spheres:
                        ts, ms = _sphere_hit(p, dirs_sh, c, r)
                        vis = np.where(ms, 0.0, vis)
                    ndl_k = np.clip(_dot(np.broadcast_to(n, dirs.shape), ld), 0.0, 1.0)
                    acc += alb * ndl_k[..., None] * vis[..., None]
                rad_floor = 0.65 * rad_floor + 0.35 * acc
            # glossy sky reflection, roughness blur via 3 jittered mirror rays
            refl = dirs - 2.0 * _dot(dirs, np.broadcast_to(n, dirs.shape))[..., None] * n
            refl = _norm(refl)
            sky = np.clip(refl[..., 1], 0.0, 1.0) * 0.18
            rad_floor = rad_floor + sky[..., None] * np.array([0.75, 0.82, 1.0])
            rad = np.where(floor[..., None], rad_floor, rad)

        # --- spheres: Lambert + specular highlight + rim
        for i, (c, r, col) in enumerate(spheres):
            m = (hit_kind == 2) & (hit_idx == i)
            if not np.any(m):
                continue
            n = _norm(p - c)
            ndl = np.clip(_dot(n, light), 0.0, 1.0)
            diff = col * (0.15 + 0.85 * ndl[..., None])
            hv = _norm(light + _norm(cam - p))
            spec = np.clip(_dot(n, hv), 0.0, 1.0) ** 48.0
            rim = (1.0 - np.clip(np.abs(_dot(n, -dirs)), 0.0, 1.0)) ** 3.0
            rad_s = diff + spec[..., None] * 0.5 + rim[..., None] * col * 0.35
            rad = np.where(m[..., None], rad_s, rad)

        # background: vertical gradient sky
        bg = hit_kind == 0
        if np.any(bg):
            tsky = np.clip(dirs[..., 1] * 0.5 + 0.5, 0.0, 1.0)
            sky_col = np.stack([0.55 + 0.35 * tsky,
                                0.68 + 0.25 * tsky,
                                0.95 - 0.05 * tsky], axis=-1)
            rad = np.where(bg[..., None], sky_col, rad)

        # slight exposure curve, keeps values in a benign HDR range
        img += rad
    img /= float(spp)
    return np.clip(img, 0.0, None)
