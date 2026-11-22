"""Block-wise low-rank compression of a rendered HDR frame.

Pipeline
--------
1. Split each colour channel into BxB blocks (edge-padded with zeros).
2. Take the SVD of every block.  Truncating to rank r gives the Frobenius-optimal
   rank-r approximant; by Eckart-Young the squared error of *any* rank-r
   approximant of that block is bounded below by the sum of the discarded
   singular values squared:

       ||A - A_r||_F^2  >=  sum_{i>r} sigma_i^2                (lower bound)

   This bound is per block; summing over blocks gives the lower bound for the
   block-structured codec.  A measured residual that dips below it means the
   residual or the bound is computed inconsistently (the implementation is
   biased) -- `check_bounds` fails hard on that.
3. Reconstruct and measure the true square error (which includes factor
   quantisation), the theoretical bound, and the measured storage volume.

Rank 0 is the degenerate case: it cannot carry any image data, so it falls back
to lossless storage and the record is flagged (`fallback: true`).
"""
from __future__ import annotations

import struct
import numpy as np
from dataclasses import dataclass, field, asdict

import coding

MAGIC = b"LRK1"


@dataclass
class Block:
    row: int
    col: int
    h: int
    w: int


@dataclass
class ChannelData:
    index: int
    padded: np.ndarray          # (Hp, Wp) float64
    blocks: list[Block]
    U: np.ndarray               # (nblocks, bs, bs) float64
    s: np.ndarray               # (nblocks, bs) float64
    Vt: np.ndarray              # (nblocks, bs, bs) float64
    H: int = 0
    W: int = 0
    bs: int = 0
    rmax: int = 0


def _pad_to_blocks(img2d, bs):
    H, W = img2d.shape
    if H % bs or W % bs:
        # Zero padding would make the block matrices larger than the visible
        # region, so Eckart-Young's bound (a bound on the *padded* block) would
        # no longer lower-bound the error on the visible pixels. Refuse rather
        # than report a bound that is not rigorous.
        raise ValueError(
            f"image {H}x{W} is not divisible by block size {bs}; "
            "choose a block size that divides both dimensions"
        )
    return np.array(img2d, dtype=np.float64, copy=True), H // bs, W // bs


def prepare_channels(volume, bs):
    """Compute full per-block SVDs for every channel (done once, reused for all ranks)."""
    H, W, C = volume.shape
    chans = []
    for c in range(C):
        padded, nh, nw = _pad_to_blocks(volume[:, :, c], bs)
        blocks = []
        stack = np.empty((nh * nw, bs, bs), dtype=np.float64)
        for by in range(nh):
            for bx in range(nw):
                i = by * nw + bx
                stack[i] = padded[by * bs:(by + 1) * bs, bx * bs:(bx + 1) * bs]
                blocks.append(Block(by, bx, bs, bs))
        U, s, Vt = np.linalg.svd(stack, full_matrices=False)
        chans.append(ChannelData(c, padded, blocks, U, s, Vt,
                                 H=H, W=W, bs=bs, rmax=bs))
    return chans


def _quantize_factors(A, B, bits):
    """Fixed-point quantise the two factor tensors to `bits` signed levels.

    A = U_r * diag(s_r), B = Vt_r.  A's magnitude is bounded by s_max, B's
    entries lie in [-1, 1], so one per-channel scale each is enough (tiny header).
    Returns (Aq uint, Bq uint, scaleA, shift) with reconstruction
    A ~= (Aq - shift) * scaleA / qmax, B ~= (Bq - shift) / qmax.
    """
    qmax = float((1 << (bits - 1)) - 1)
    shift = 1 << (bits - 1)
    dtype = np.uint16 if bits <= 16 else np.uint32
    if bits == 32:
        # "32-bit factors" = native float32, lossless bit-exact storage
        Aq = np.ascontiguousarray(A, dtype=np.float32).view(np.uint32)
        Bq = np.ascontiguousarray(B, dtype=np.float32).view(np.uint32)
        return Aq, Bq, 1.0, 0, np.float32
    scaleA = float(np.max(np.abs(A))) if A.size else 1.0
    if scaleA <= 0.0:
        scaleA = 1.0
    Aq = np.rint(A / scaleA * qmax) + shift
    Bq = np.rint(B * qmax) + shift
    Aq = np.clip(Aq, 0, (1 << bits) - 1).astype(dtype)
    Bq = np.clip(Bq, 0, (1 << bits) - 1).astype(dtype)
    return Aq, Bq, scaleA, shift, dtype


def _dequantize_factors(Aq, Bq, scaleA, shift, bits, dtype):
    if bits == 32:
        return Aq.view(np.float32), Bq.view(np.float32)
    qmax = float((1 << (bits - 1)) - 1)
    A = (Aq.astype(np.float64) - shift) * (scaleA / qmax)
    B = (Bq.astype(np.float64) - shift) / qmax
    return A, B


def _header(nch, bs, ranks, scales, bits, mode):
    ranks_b = struct.pack("<%dH" % len(ranks), *ranks)
    scales_b = struct.pack("<%dd" % len(scales), *scales)
    head = MAGIC + struct.pack("<BHHB", mode, bs, nch, bits)
    return head + ranks_b + scales_b


def compress_blocks(volume, chans, rank, factor_bits=32):
    """Block-wise rank-r compression.  Returns (recon, stats).

    Rank 0 carries no image information, so it degenerates to the lossless path
    and the record is flagged (`fallback: true`) -- never silently stored as
    zeros.
    """
    if rank <= 0:
        recon, st = compress_lossless(volume, rank=0, reason="rank_zero")
        st["mode"] = "block"
        st["rank_effective"] = [0] * volume.shape[2]
        st["block_size"] = int(chans[0].bs)
        st["factor_bits"] = int(factor_bits)
        return recon, st

    H, W, C = volume.shape
    recon = np.zeros_like(volume, dtype=np.float64)
    payload = []          # (uint array, dtype) for the size measurement
    ranks = []
    scales = []
    bound = 0.0           # sum of discarded singular values squared
    sq_err = 0.0          # measured squared residual on valid pixels
    for ch in chans:
        r = int(min(rank, ch.rmax))
        ranks.append(r)
        # per-block factors
        A = ch.U[:, :, :r] * ch.s[:, None, :r] if r > 0 else np.zeros((len(ch.blocks), ch.bs, 0))
        B = ch.Vt[:, :r, :]
        Aq, Bq, scaleA, shift, dtype = _quantize_factors(A, B, factor_bits)
        payload.append((Aq, dtype))
        payload.append((Bq, dtype))
        scales.append(scaleA)
        # theoretical lower bound: discarded energy
        if r < ch.rmax:
            bound += float(np.sum(ch.s[:, r:] ** 2))
        # reconstruct
        if r == 0:
            blk_rec = np.zeros((len(ch.blocks), ch.bs, ch.bs))
        else:
            Ad, Bd = _dequantize_factors(Aq, Bq, scaleA, shift, factor_bits, dtype)
            blk_rec = Ad @ Bd
        # scatter + accumulate error on valid region only
        nh = ch.padded.shape[0] // ch.bs
        nw = ch.padded.shape[1] // ch.bs
        rec_pad = np.zeros_like(ch.padded)
        for i, bl in enumerate(ch.blocks):
            by, bx = divmod(i, nw)
            rec_pad[by * ch.bs:(by + 1) * ch.bs, bx * ch.bs:(bx + 1) * ch.bs] = blk_rec[i]
        recon[:, :, ch.index] = rec_pad[:H, :W]
        diff = ch.padded[:H, :W] - rec_pad[:H, :W]
        sq_err += float(np.sum(diff * diff))
    head = _header(C, chans[0].bs, ranks, scales, factor_bits, mode=1)
    bytes_est = len(head) + coding.coded_size(payload)
    return recon, {
        "mode": "block",
        "rank": int(rank),
        "rank_effective": ranks,
        "block_size": int(chans[0].bs),
        "factor_bits": int(factor_bits),
        "squared_error": sq_err,
        "bound_sq": bound,
        "bytes": int(bytes_est),
        "header_bytes": len(head),
        "fallback": False,
    }


def compress_global(volume, rank, factor_bits=32):
    """Whole-channel SVD (no blocking) -- comparison point, not the primary codec."""
    H, W, C = volume.shape
    recon = np.zeros_like(volume)
    payload = []
    ranks, scales = [], []
    bound = 0.0
    sq_err = 0.0
    for c in range(C):
        A = volume[:, :, c]
        U, s, Vt = np.linalg.svd(A, full_matrices=False)
        r = int(min(rank, len(s)))
        ranks.append(r)
        Fa = U[:, :r] * s[None, :r]
        Fb = Vt[:r, :]
        # float32 factors; global mode reports the idealised storage of k(H+W)
        # float32 numbers, which is what a real global codec would store.
        payload.append((np.ascontiguousarray(Fa, dtype=np.float32).view(np.uint32), np.uint32))
        payload.append((np.ascontiguousarray(Fb, dtype=np.float32).view(np.uint32), np.uint32))
        scales.append(1.0)
        bound += float(np.sum(s[r:] ** 2))
        rec = Fa.astype(np.float64) @ Fb.astype(np.float64)
        recon[:, :, c] = rec
        sq_err += float(np.sum((A - rec) ** 2))
    head = _header(C, 0, ranks, scales, 32, mode=2)
    bytes_est = len(head) + coding.coded_size(payload)
    return recon, {
        "mode": "global",
        "rank": int(rank),
        "rank_effective": ranks,
        "block_size": 0,
        "factor_bits": 32,
        "squared_error": sq_err,
        "bound_sq": bound,
        "bytes": int(bytes_est),
        "header_bytes": len(head),
        "fallback": False,
    }


def compress_lossless(volume, rank=0, reason="rank_zero"):
    """Lossless reference path; also the rank-0 fallback."""
    H, W, C = volume.shape
    payload = [(np.ascontiguousarray(volume[:, :, c], dtype=np.float32).view(np.uint32).reshape(-1),
                np.uint32) for c in range(C)]
    head = _header(C, 0, [rank] * C, [1.0] * C, 32, mode=0)
    coded = len(head) + coding.coded_size(payload)
    raw = int(volume.size * 4)
    return volume.copy(), {
        "mode": "lossless",
        "rank": int(rank),
        "rank_effective": [int(rank)] * C,
        "block_size": 0,
        "factor_bits": 32,
        "squared_error": 0.0,
        "bound_sq": float(np.sum(volume.astype(np.float64) ** 2)),
        "bytes": int(coded),
        "raw_bytes": raw,
        "header_bytes": len(head),
        "fallback": bool(rank == 0),
        "fallback_reason": reason if rank == 0 else "",
    }


def finish_stats(rec, stats, volume):
    """Attach PSNR / RMSE / compression ratios to a stats dict."""
    H, W, C = volume.shape
    n = H * W * C
    mse = stats["squared_error"] / n
    peak = float(np.max(volume))
    stats["rmse"] = float(np.sqrt(mse))
    stats["psnr_db"] = float(10.0 * np.log10(peak * peak / mse)) if mse > 0 else float("inf")
    stats["peak"] = peak
    raw = int(volume.size * 4)
    stats["raw_bytes"] = raw
    lossless = stats.get("lossless_bytes")
    stats["ratio_vs_raw"] = float(raw / stats["bytes"]) if stats["bytes"] > 0 else float("inf")
    stats["ratio_vs_lossless"] = (float(lossless / stats["bytes"])
                                  if lossless and stats["bytes"] > 0 else float("nan"))
    stats["bytes_per_pixel"] = float(stats["bytes"] / n)
    return stats


def check_bounds(records, rel_tol=1e-6):
    """Verify measured residual never undercuts the Eckart-Young lower bound.

    A residual below the bound by more than a tolerance would mean the
    implementation is biased (the bound omits a term, or blocks are
    double-counted).  The tolerance absorbs float round-off in the SVD and in
    the reconstruction; it is tight (1e-6 relative) so a real bias is caught.
    Returns (ok, worst_margin).
    """
    ok = True
    worst = float("inf")
    for rec in records:
        if rec.get("mode") == "lossless" or rec.get("squared_error", 0.0) == 0.0:
            continue
        scale = max(rec["bound_sq"], rec.get("energy_sq", 0.0), 1e-30)
        tol = rel_tol * scale
        margin = rec["squared_error"] - rec["bound_sq"]
        worst = min(worst, margin)
        rec["bound_rel_margin"] = margin / scale
        if margin < -tol:
            rec["biased"] = True
            ok = False
        else:
            rec.setdefault("biased", False)
    return ok, (0.0 if worst == float("inf") else worst)


def volume_energy(volume):
    return float(np.sum(volume.astype(np.float64) ** 2))
