import numpy as np

P8 = [7, 4, 1, 6, 3, 0, 5, 2]
P16 = [14, 15, 8, 9, 2, 3, 12, 13, 6, 7, 0, 1, 10, 11, 4, 5]


class SAND:
    def __init__(self, n):
        self.n = n
        self.w = n // 4
        self.wm = np.uint64((1 << self.w) - 1)
        self.full = np.uint64((1 << n) - 1) if n < 64 else np.uint64(0xFFFFFFFFFFFFFFFF)
        self.p = P8 if n == 32 else P16
        self.colmask = [sum(1 << (self.w * i + j) for i in range(4)) for j in range(self.w)]

    def rows(self, x):
        w = self.w
        return [(x >> np.uint64(w * i)) & self.wm for i in range(4)]

    def join(self, r):
        w = self.w
        out = np.zeros_like(r[0])
        for i in range(4):
            out |= (r[i] & self.wm) << np.uint64(w * i)
        return out

    def rot(self, x):
        w = self.w
        r = self.rows(x)
        return self.join([((v << np.uint64(1)) | (v >> np.uint64(w - 1))) & self.wm for v in r])

    def G0(self, x):
        x0, x1, x2, x3 = self.rows(x)
        y0 = (x3 & x2) ^ x0
        y3 = (y0 & x1) ^ x3
        return self.join([y0, x1, x2, y3])

    def G1(self, x):
        x0, x1, x2, x3 = self.rows(x)
        y2 = (x3 & x1) ^ x2
        y1 = (y2 & x0) ^ x1
        return self.join([x0, y1, y2, x3])

    def P(self, x):
        out = np.zeros_like(x)
        for j, pj in enumerate(self.p):
            m = x & np.uint64(self.colmask[j])
            s = pj - j
            out |= (m << np.uint64(s)) if s >= 0 else (m >> np.uint64(-s))
        return out

    def f(self, x):
        return self.P(self.G0(x) ^ self.G1(self.rot(x)))


def rnd(rng, n, size):
    v = rng.integers(0, 2**63, size=size, dtype=np.uint64) ^ (rng.integers(0, 2, size=size, dtype=np.uint64) << np.uint64(63))
    if n < 64:
        v &= np.uint64((1 << n) - 1)
    return v
