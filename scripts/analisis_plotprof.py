import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
c = np.load("scan_cyl.npy").astype(np.float64); x, rho, th = c.T; ath = np.abs(th)
step = 0.05; xb = np.arange(0, 220 + step, step)
def prof(sel):
    idx = np.floor(x[sel]/step).astype(int); o = np.argsort(idx); ii = idx[o]; rr = rho[sel][o]
    u, st, n = np.unique(ii, return_index=True, return_counts=True)
    p = np.full(len(xb), np.nan)
    for k, a, m in zip(u, st, n):
        if m >= 8: p[k] = np.median(rr[a:a+m])
    return p
fig, axs = plt.subplots(3, 1, figsize=(26, 14))
sects = [(0,15),(15,35),(35,60),(60,85),(95,120),(120,140)]
P = {}
for a, b in sects:
    P[(a,b)] = prof((ath >= a) & (ath < b))
for ax, (lo, hi) in zip(axs, [(0, 75), (70, 145), (140, 222)]):
    for k, p in P.items(): ax.plot(xb, p, lw=0.7, label=f"{k}")
    ax.set_xlim(lo, hi); m = (xb>=lo)&(xb<=hi); ax.set_ylim(np.nanmin(P[(60,85)][m])-1, np.nanmax(P[(60,85)][m])+2); ax.grid(True, lw=0.3)
    ax.set_xticks(np.arange(lo, hi, 2))
axs[0].legend()
plt.tight_layout(); plt.savefig("perfil.png", dpi=70)
np.save("prof_sectors.npy", np.vstack([xb] + [P[k] for k in sects]))
