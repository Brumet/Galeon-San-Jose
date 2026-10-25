import numpy as np
from scipy.spatial import cKDTree
s = np.load("scan_aligned.npy").astype(np.float64); ax = np.load("axis.npy"); xb, pr = np.load("profile_adapt.npy")
x = s[:,0]; cyx = np.polyval(ax[0], x); dy = s[:,1]-cyx; dz = s[:,2]-np.polyval(ax[1], x)
r = np.hypot(dy, dz) - np.interp(x, xb, pr)
m = (x > 100) & (x < 132) & (r > 0.4) & (dz > 0) & (np.abs(np.degrees(np.arctan2(dy, dz))) < 60)
A = s[m & (dy > 0)]; B = s[m & (dy < 0)]
print("A(+y)", len(A), "B(-y)", len(B), "A x", A[:,0].min().round(2), A[:,0].max().round(2), "B x", B[:,0].min().round(2), B[:,0].max().round(2))
Bm = B.copy(); Bm[:,1] = 2 * np.polyval(ax[0], Bm[:,0]) - Bm[:,1]      # espejo en el plano vertical del eje
def icp(P, Q, it=40):
    tq = cKDTree(Q); Rt = np.eye(3); t = Q.mean(0) - P.mean(0); cur = P + t
    for k in range(it):
        d, i = tq.query(cur, workers=-1); th = np.percentile(d, 85); g = d < th
        a, b = cur[g], Q[i[g]]; ca, cb = a.mean(0), b.mean(0)
        U, _, Vt = np.linalg.svd((a - ca).T @ (b - cb)); Rk = (U @ Vt).T
        if np.linalg.det(Rk) < 0: Vt[-1] *= -1; Rk = (U @ Vt).T
        cur = (cur - ca) @ Rk.T + cb
    d, i = tq.query(cur, workers=-1); return cur, d
Bal, d = icp(Bm, A)
print("dist B'->A pct", np.percentile(d, [50, 75, 90, 95, 99]).round(3))
np.save("asaA.npy", A.astype(np.float32)); np.save("asaB_on_A.npy", Bal.astype(np.float32)); np.save("asaB_d.npy", d.astype(np.float32))
