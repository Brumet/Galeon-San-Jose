import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d, maximum_filter1d
xb, ps = np.load("profile_final.npy")
step = xb[1] - xb[0]
light = ps
heavy = gaussian_filter1d(ps, 0.6 / step)
mid = gaussian_filter1d(ps, 0.25 / step)
curv = np.abs(np.gradient(np.gradient(mid, step), step))
slope = np.abs(np.gradient(mid, step))
ss = lambda a, b, v: np.clip((v - a) / (b - a), 0, 1)
feat = maximum_filter1d(np.maximum(ss(0.25, 0.8, curv), ss(0.1, 0.2, slope)), int(1.6 / step))   # cerca de anillos
w = np.clip(feat, 0, 1); w = gaussian_filter1d(w, 0.2 / step)
out = w * light + (1 - w) * heavy
print("diff light-heavy flat rms", np.std((light - heavy)[w < 0.1]).round(4), "flat frac", (w < 0.1).mean().round(2))
np.save("profile_adapt.npy", np.vstack([xb, out]))
fig, axs = plt.subplots(3,1, figsize=(26,12))
for ax,(lo,hi) in zip(axs, [(0,25),(120,150),(180,210)]):
    m=(xb>=lo)&(xb<=hi); ax.plot(xb[m], light[m], lw=.6); ax.plot(xb[m], out[m], lw=1); ax2=ax.twinx(); ax2.plot(xb[m], w[m], 'r', lw=.4); ax.grid(lw=.3)
plt.tight_layout(); plt.savefig("perfil3.png", dpi=65)
