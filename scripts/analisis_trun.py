import numpy as np
s = np.load("scan_aligned.npy").astype(np.float64); c = np.load("scan_cyl.npy").astype(np.float64)
x, rho, th = c.T
P = np.load("prof_sectors.npy"); xb = P[0]; pr = P[3]
ok = ~np.isnan(pr); pi = np.interp(x, xb[ok], pr[ok])
ax = np.load("axis.npy"); cz = np.polyval(ax[1], s[:,0]); cy = np.polyval(ax[0], s[:,0])
def fit_circle(u, v):
    A = np.c_[2*u, 2*v, np.ones_like(u)]; b = u*u + v*v
    q, *_ = np.linalg.lstsq(A, b, rcond=None); return q[0], q[1], np.sqrt(q[2] + q[0]**2 + q[1]**2)
for side in (1, -1):
    m = (x > 108) & (x < 138) & (rho > pi + 0.8) & (np.sign(s[:,1]-cy) == side) & (np.abs(np.abs(th) - 90) < 45)
    p = s[m]; yy = np.abs(p[:,1] - cy[m])
    print("side", side, "n", m.sum(), "y range", yy.min().round(2), yy.max().round(2))
    # cuerpo cilíndrico: yy entre ρ+1 y extremo-0.4
    yend = np.percentile(yy, 99.5)
    k = (yy > yy.min() + 1.0) & (yy < yend - 0.5)
    u, v = p[k,0], p[k,2] - cz[m][k]
    xc, zc, r = fit_circle(u, v)
    for _ in range(3):
        d = np.hypot(u-xc, v-zc); kk = np.abs(d - r) < 0.3; xc, zc, r = fit_circle(u[kk], v[kk])
    ang = np.degrees(np.arctan2(u - xc, v - zc)); h, _ = np.histogram(ang, bins=24, range=(-180, 180))
    print(f"  xc {xc:.3f} zc_rel {zc:.3f} r {r:.3f} yend {yend:.3f} rms {np.std(np.hypot(u-xc,v-zc)[np.abs(np.hypot(u-xc,v-zc)-r)<0.3]-r):.3f}")
    print("  coverage by 15deg (0=top):", h)
    # cara extrema
    f = yy > yend - 0.3
    print("  end face pts", f.sum(), "radial max", np.hypot(p[f,0]-xc, p[f,2]-cz[m][f]-zc).max().round(2))
