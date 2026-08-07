#!/usr/bin/env python3
"""Reading and writing of SAND-64 differential characteristics.

A characteristic file is a plain text file of the form

    # any comment line
    # characteristic 1
      0 0012 0000 0019 0000
      1 0002 0000 0012 0000
      ...
     12 008A 0000 0090 0000

    # characteristic 2
      ...

Each data line is  r  dX^r(hi) dX^r(lo)  dY^r(hi) dY^r(lo)  in hexadecimal,
so that a characteristic of R rounds occupies R+1 lines.
"""


def fmt(word):
    """32-bit word -> 'HHHH HHHH'."""
    return f"{word >> 16:04X} {word & 0xFFFF:04X}"


def read_trails(path, num_rounds=None):
    """Return the list of characteristics stored in `path`.

    Each characteristic is a list of (dX^r, dY^r) pairs of 32-bit integers,
    indexed by r = 0, 1, ..., R.  Incomplete records (which can only arise if a
    file was truncated) are dropped when `num_rounds` is given.
    """
    trails, cur = [], None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# characteristic"):
                if cur:
                    trails.append(cur)
                cur = []
            elif line and not line.startswith("#"):
                p = line.split()
                if len(p) >= 5 and p[0].isdigit():
                    cur.append((int(p[1] + p[2], 16), int(p[3] + p[4], 16)))
    if cur:
        trails.append(cur)
    if num_rounds is not None:
        trails = [t for t in trails if len(t) == num_rounds + 1]
    return trails


def write_trails(path, trails, header=()):
    with open(path, "w") as f:
        for line in header:
            f.write(f"# {line}\n")
        for idx, t in enumerate(trails, 1):
            f.write(f"\n# characteristic {idx}\n")
            for r, (dx, dy) in enumerate(t):
                f.write(f"{r:>3} {fmt(dx)} {fmt(dy)}\n")
