import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
c = np.load("scan_cyl.npy").astype(np.float64); x, rho, th = c.T; ath = np.abs(th)
step = 0.02; xb = np.arange(0, 220 + step, step)
trun = (x > 106) & (x < 136)
sel = ((ath > 35) & (ath < 88) & ~trun) | (trun & (ath > 33) & (ath < 62))
idx = np.floor(x[sel]/step).astype(int); o = np.argsort(idx); ii = idx[o]; rr = rho[sel][o]
u, st, n = np.unique(ii, return_index=True, return_counts=True)
p = np.full(len(xb), np.nan)
for k, a, m in zip(u, st, n):
    if m >= 10: p[k] = np.median(rr[a:a+m])
ok = ~np.isnan(p); p = np.interp(xb, xb[ok], p[ok])
ps = savgol_filter(p, 25, 3)          # ventana 0.5 mm
ps = gaussian_filter1d(ps, 2.5)       # sigma 0.05 mm
old = np.load("prof_sectors.npy")
np.save("profile_final.npy", np.vstack([xb, ps]))
fig, axs = plt.subplots(3,1, figsize=(26,12))
for ax,(lo,hi) in zip(axs, [(0,20),(105,135),(185,210)]):
    m=(xb>=lo)&(xb<=hi); ax.plot(xb[m], p[m], lw=.6, label="mediana"); ax.plot(xb[m], ps[m], lw=1.2, label="final")
    mo=(old[0]>=lo)&(old[0]<=hi); ax.plot(old[0][mo], old[3][mo], lw=.6, label="v1"); ax.grid(lw=.3); ax.set_xticks(np.arange(lo,hi,0.5))
axs[0].legend(); plt.tight_layout(); plt.savefig("perfil2.png", dpi=65)
d = np.diff(ps); print("noise rms raw", np.std(p - ps).round(4))
