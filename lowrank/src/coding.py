"""Lossless payload coder used to measure storage volume.

Each value is split into bit planes (bitshuffle), the planes are *bit-packed*
and then deflated. Bit-plane separation is what real codecs lean on: the
exponent/sign bits of a low-rank payload are highly redundant and compress
well, while the low mantissa bits are close to uniform and stay at ~raw size.
Packing the planes means incompressible data costs its true size (no 8x
byte-per-bit blowup), so the reported volume is honest rather than inflated.

Round-trip is verified on every measurement, so "lossless" is tested, not
asserted.
"""
import zlib
import numpy as np


def _n_planes(arr):
    return arr.size * np.dtype(arr.dtype).itemsize * 8


def _to_bits(arr):
    a = np.ascontiguousarray(arr)
    n = a.size
    nbytes = np.dtype(a.dtype).itemsize
    raw = np.frombuffer(a.tobytes(), dtype=np.uint8)          # (n*nbytes,)
    bits = np.unpackbits(raw, bitorder="little")              # (n*nbytes*8,)
    bits = bits.reshape(n, nbytes, 8)                         # (value, byte, bit)
    planes = np.transpose(bits, (2, 1, 0))                    # (bit, byte, value)
    planes = planes.reshape(8 * nbytes, n)                    # -> pack rows
    pad = (-n) % 8
    if pad:
        planes = np.pad(planes, ((0, 0), (0, pad)))
    return np.packbits(planes, axis=1, bitorder="little")     # (8*nbytes, ceil(n/8))


def _from_bits(packed, shape, dtype):
    n = int(np.prod(shape))
    nbytes = np.dtype(dtype).itemsize
    planes = np.frombuffer(np.ascontiguousarray(packed).tobytes(), dtype=np.uint8)
    planes = planes.reshape(8 * nbytes, -1)
    bits = np.unpackbits(planes, axis=1, bitorder="little")[:, :n]  # (8*nbytes, n)
    bits = bits.reshape(8, nbytes, n)
    vals = np.transpose(bits, (2, 1, 0)).reshape(-1)          # (value, byte, bit)
    raw = np.packbits(vals, bitorder="little")[: n * nbytes]
    return np.frombuffer(raw.tobytes(), dtype=dtype).reshape(shape)


def coded_size(arrays):
    """Sum of deflated packed-plane sizes (bytes); verifies lossless round-trip."""
    total = 0
    for arr, dtype in arrays:
        a = np.ascontiguousarray(np.asarray(arr, dtype=dtype))
        packed = _to_bits(a)
        comp = zlib.compress(packed.tobytes(), 6)
        back = _from_bits(zlib.decompress(comp), a.shape, dtype)
        if not np.array_equal(back, a):
            raise RuntimeError("lossless coder round-trip mismatch")
        total += len(comp)
    return total


def raw_size(arrays):
    return int(sum(np.asarray(a).size * np.dtype(d).itemsize for a, d in arrays))
