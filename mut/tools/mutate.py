#!/usr/bin/env python3
"""Apply one source mutation to a sandbox copy of a module.

Usage:
  mutate.py --list                     # print the catalogue
  mutate.py --check                    # verify every find-string matches exactly once
  mutate.py <id> <src> <out>           # write mutated copy of <src> to <out>

No third-party modules; stdlib only.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CATALOGUE = os.path.join(HERE, os.pardir, "mutants.json")


def load():
    with open(CATALOGUE, "r", encoding="utf-8") as f:
        return json.load(f)


def apply_one(mut, text):
    find = mut["find"]
    repl = mut["replace"]
    n = text.count(find)
    if n != 1:
        raise SystemExit(
            "mutant %s: find-string occurs %d times (need exactly 1)\n"
            "--- find ---\n%s" % (mut["id"], n, find))
    return text.replace(find, repl, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("id", nargs="?")
    ap.add_argument("src", nargs="?")
    ap.add_argument("out", nargs="?")
    a = ap.parse_args()
    cat = load()
    if a.list:
        for m in cat["mutants"]:
            print("%-26s %-8s %-10s %s" % (m["id"], m["module"], m["operator"], m["title"]))
        return 0
    if a.check:
        bad = 0
        for m in cat["mutants"]:
            src = os.path.join(cat["root"], m["module"], "src", m["module"] + ".c")
            with open(src, "r", encoding="utf-8") as f:
                text = f.read()
            n = text.count(m["find"])
            if n != 1:
                bad += 1
                print("BAD  %-26s matches=%d" % (m["id"], n))
            else:
                line = text[:text.index(m["find"])].count("\n") + 1
                print("ok   %-26s %s:%d" % (m["id"], m["module"] + ".c", line))
        return 1 if bad else 0
    if not (a.id and a.src and a.out):
        ap.error("need <id> <src> <out>, or --check/--list")
    m = next((x for x in cat["mutants"] if x["id"] == a.id), None)
    if not m:
        raise SystemExit("no such mutant: " + a.id)
    with open(a.src, "r", encoding="utf-8") as f:
        text = f.read()
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(apply_one(m, text))
    return 0


if __name__ == "__main__":
    sys.exit(main())
