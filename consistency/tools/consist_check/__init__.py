"""Framebuffer IO and pixel-level diff classification.

The two C renderers emit identical .fbuf files:
  uint32 magic, width, height, channels(=5), then float32 pixels of
  [r, g, b, linear_view_depth, primitive_id].
"""

import numpy as np

CH = 5
MAGIC = 0x46425546
BG_ID = -1.0
BG_DEPTH = 1.0e30

# A pixel counts as "agreeing" only when it is indistinguishable above the
# float32 quantization floor. COLOR_EPS sits just above that floor (~1e-6
# here) and well below the genuine precision-class signal (~1e-2), so it
# separates noise from real disagreement without carving into any class.
COLOR_EPS = 5.0e-6
DEPTH_REL_EPS = 1.0e-6

OCCLUSION = "occlusion"
PRECISION = "precision"
BOUNDARY = "boundary"
UNCLASSIFIED = "unclassified"
CLASSES = (OCCLUSION, PRECISION, BOUNDARY, UNCLASSIFIED)


def read_fbuf(path):
    with open(path, "rb") as f:
        raw = np.fromfile(f, dtype=np.uint32, count=4)
        if len(raw) != 4 or int(raw[0]) != MAGIC:
            raise ValueError(f"{path}: bad framebuffer header")
        w, h, ch = int(raw[1]), int(raw[2]), int(raw[3])
        if ch != CH:
            raise ValueError(f"{path}: expected {CH} channels, got {ch}")
        pix = np.fromfile(f, dtype=np.float32)
    expected = w * h * CH
    if pix.size != expected:
        raise ValueError(f"{path}: truncated ({pix.size} floats, want {expected})")
    return pix.reshape(h, w, CH)


def channels(fb):
    rgb = fb[:, :, 0:3].astype(np.float64)
    depth = fb[:, :, 3].astype(np.float64)
    ids = fb[:, :, 4].astype(np.float64).round().astype(np.int64)
    return rgb, depth, ids


def _neighbor_id_sets(ids):
    """For each pixel, a boolean (H, W, max_id+1) volume marking which
    primitive ids occur in its 3x3 window. Background (-1) is dropped."""
    h, w = ids.shape
    max_id = int(ids.max())
    present = np.zeros((h, w, max(1, max_id + 1)), dtype=bool)
    hit = ids >= 0
    present[hit, ids[hit]] = True
    out = np.zeros_like(present)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            ys = slice(max(0, dy), h + min(0, dy))
            xs = slice(max(0, dx), w + min(0, dx))
            yd = slice(max(0, -dy), h + min(0, -dy))
            xd = slice(max(0, -dx), w + min(0, -dx))
            out[ys, xs] |= present[yd, xd]
    return out


def _background_neighbor(ids):
    """True where the 3x3 window contains a background (-1) pixel."""
    bg = (ids < 0).astype(np.uint8)
    out = np.zeros_like(bg)
    h, w = ids.shape
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            ys = slice(max(0, dy), h + min(0, dy))
            xs = slice(max(0, dx), w + min(0, dx))
            yd = slice(max(0, -dy), h + min(0, -dy))
            xd = slice(max(0, -dx), w + min(0, -dx))
            out[ys, xs] |= bg[yd, xd]
    return out.astype(bool)


def classify(rast_fb, trace_fb):
    """Return a label map and per-pixel diagnostic metrics.

    The buckets are mutually exclusive and decided by geometry, not by an
    arbitrary numeric cutoff. Any diff pixel matching no rule lands in
    `unclassified`, which the harness treats as a hard failure (the
    taxonomy must cover everything the two renderers can disagree on).

    precision : same primitive id in both frames; only the color/depth
                value drifts (float32 raster path vs float64 intersection).
    boundary  : the two frames pick different ids AND both contested ids
                coexist in this pixel's 3x3 neighborhood — a silhouette or
                shared-edge membership dispute, exactly one pixel wide.
    occlusion : the frames pick different ids but the contested ids do NOT
                coexist locally; visibility/order flips across an interior
                depth-crossing band rather than along a geometric edge.
    """
    if rast_fb.shape != trace_fb.shape:
        raise ValueError("framebuffer size mismatch")
    r_rgb, r_z, r_id = channels(rast_fb)
    t_rgb, t_z, t_id = channels(trace_fb)

    color_delta = np.abs(r_rgb - t_rgb).max(axis=2)
    z_scale = np.maximum(np.abs(r_z), np.abs(t_z))
    depth_delta = np.abs(r_z - t_z) / np.where(z_scale > 1e-6, z_scale, 1.0)

    hit = (r_id >= 0) | (t_id >= 0)
    id_diff = r_id != t_id
    color_bad = color_delta > COLOR_EPS
    depth_bad = (depth_delta > DEPTH_REL_EPS) & hit
    diff = hit & (id_diff | color_bad | depth_bad)

    # A boundary dispute means BOTH contested primitives are geometrically
    # present in the local 3x3 window — the two ids coexist locally, so the
    # disagreement is a one-pixel membership choice along a shared edge.
    r_nbhd = _neighbor_id_sets(r_id)
    t_nbhd = _neighbor_id_sets(t_id)
    h, w = r_id.shape
    rid_hit, tid_hit = r_id >= 0, t_id >= 0
    r_valid = np.where(rid_hit, r_id, 0)
    t_valid = np.where(tid_hit, t_id, 0)
    # raster's chosen id is also present in tracer's local window and vice versa
    cross_a = rid_hit & tid_hit & np.take_along_axis(
        t_nbhd, r_valid[:, :, None], axis=2)[:, :, 0]
    cross_b = rid_hit & tid_hit & np.take_along_axis(
        r_nbhd, t_valid[:, :, None], axis=2)[:, :, 0]
    # silhouette against background: one frame hits where the other misses,
    # and the hit frame has background in its local window (true edge pixel)
    bg_dispute = (rid_hit ^ tid_hit) & (
        _background_neighbor(r_id) | _background_neighbor(t_id))
    boundary = id_diff & ((cross_a & cross_b) | bg_dispute)

    labels = np.full(r_id.shape, "", dtype=object)
    labels[~diff] = "match"
    labels[diff & boundary] = BOUNDARY
    labels[diff & id_diff & ~boundary] = OCCLUSION
    same_id_prec = diff & ~id_diff & (color_bad | depth_bad)
    labels[same_id_prec] = PRECISION
    labels[diff & (labels == "")] = UNCLASSIFIED

    metrics = {
        "diff": diff,
        "id_diff": id_diff & hit,
        "color_delta": color_delta,
        "depth_delta": depth_delta,
        "boundary": boundary,
    }
    return labels, metrics
