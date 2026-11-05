import numpy as np, math
from sand_np import SAND, rnd
rng=np.random.default_rng(1)
T64=[(0x00900000,0x008C0000),(0x00010000,0x00900000),(0,0x00010000),(0x00010000,0),(0x00900000,0x00010000),(0x008C0100,0x00900000),(0x000E1000,0x008C0100),(0x00950000,0x000E1000),(0x00510000,0x00950000),(0x00280000,0x00510000),(0x00380000,0x00280000),(0x00480000,0x00380000),(0x00540000,0x00480000)]
L64=[0,-4,-6,-6,-8,-12,-21,-30,-38,-44,-48,-54,-58]
T128x=[0x1115,0x3C0,0x114,0xC,0x10,0,0x10,0xC,0x114,0x3C0,0x1115,0xCC0C,0x5400,0xC03F,0x4410,0xC000,0x400,0xC003,0x510,0xC0CC]
T128y=[0xCC0C]+T128x[:-1]
L128=[0,-10,-18,-24,-28,-30,-30,-32,-36,-42,-50,-60,-72,-78,-94,-100,-104,-106,-114,-120]
def chk(S,X,Y,L,N=1<<22):
    for r in range(len(X)-1):
        w=L[r]-L[r+1]
        x=rnd(rng,S.n,N); d=np.uint64(X[r])
        out=S.f(x)^S.f(x^d)^np.uint64(Y[r])
        cnt=int(np.sum(out==np.uint64(X[r+1])))
        est = -math.log2(cnt/N) if cnt else float('inf')
        print(r,w,round(est,2),cnt)
chk(SAND(32),[t[0] for t in T64],[t[1] for t in T64],L64)
chk(SAND(64),T128x,T128y,L128)
