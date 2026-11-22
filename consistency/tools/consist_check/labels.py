"""Hand-labeled difference sets and recall checking.

Labels are resolution independent and support two entry shapes:

  {"type": "point", "u", "v", "cls"}
      a single inspected pixel, stored in normalized [0,1] coordinates.

  {"type": "roi", "u0", "v0", "u1", "v1", "cls"}
      a hand-drawn rectangle in normalized coordinates that the human
      asserts contains differences of one class. ROIs are the robust way
      to label edge/line features: integer re-gridding can shift such a
      feature by a pixel column, but it stays inside the rectangle.

At the base resolution a point label must match the exact pixel. At other
resolutions a point is allowed a small, scale-proportional neighborhood
(an integer re-sampling tolerance). An ROI only requires the asserted class
to be present inside its mapped rectangle, at every resolution. In both
cases a label with no matching pixel is a miss; the gate requires zero.
"""

import json

import numpy as np

from . import classify, read_fbuf

CLASSES_FOR_LABELS = ("occlusion", "precision", "boundary")


def load_labels(path):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    for entry in data["labels"]:
        if entry["cls"] not in CLASSES_FOR_LABELS:
            raise ValueError(f"{path}: unknown class {entry['cls']}")
        entry.setdefault("type", "point")
        if entry["type"] not in ("point", "roi"):
            raise ValueError(f"{path}: unknown label type {entry['type']}")
    return data


def _nx(u, w):
    return min(w - 1, max(0, int(round(u * (w - 1)))))


def _ny(v, h):
    return min(h - 1, max(0, int(round(v * (h - 1)))))


def recall_against_labels(rast_path, trace_path, label_path):
    """Check every hand label against the automatic classification.

    Returns (per_class_breakdown, misses). An empty `misses` list is the
    "zero missed classifications" gate.
    """
    labels = load_labels(label_path)
    base_w, _ = labels["base_resolution"]
    rast = read_fbuf(rast_path)
    trace = read_fbuf(trace_path)
    h, w = rast.shape[0], rast.shape[1]
    auto, _ = classify(rast, trace)

    if w == base_w:
        point_radius = 0
    else:
        point_radius = min(2, max(1, int(round(abs(w - base_w) / base_w))))

    misses = []
    seen = {cls: {"labeled": 0, "hit": 0} for cls in CLASSES_FOR_LABELS}
    for entry in labels["labels"]:
        cls = entry["cls"]
        seen[cls]["labeled"] += 1
        found = False
        where = None

        if entry["type"] == "point":
            x, y = _nx(entry["u"], w), _ny(entry["v"], h)
            r = point_radius
            y0, y1 = max(0, y - r), min(h, y + r + 1)
            x0, x1 = max(0, x - r), min(w, x + r + 1)
            window = auto[y0:y1, x0:x1]
            found = bool(np.any(window == cls))
            where = {"x": x, "y": y, "search_radius": r}
        else:
            x0, x1 = sorted((_nx(entry["u0"], w), _nx(entry["u1"], w)))
            y0, y1 = sorted((_ny(entry["v0"], h), _ny(entry["v1"], h)))
            x1, y1 = min(w - 1, x1 + 1), min(h - 1, y1 + 1)
            window = auto[y0:y1, x0:x1]
            found = bool(np.any(window == cls))
            where = {"x0": x0, "y0": y0, "x1": x1 - 1, "y1": y1 - 1}

        if found:
            seen[cls]["hit"] += 1
        else:
            misses.append({"expected": cls, "got": "no matching pixel", **where})
    return seen, misses
