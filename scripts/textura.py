"""Mapa de textura real del escaneo en coordenadas (x, θ): H = detalle fino (< ~1 mm), L = bultos (> ~1 mm)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
c = np.load("scan_cyl.npy").astype(np.float64); x, rho, th = c.T
s = np.load("scan_aligned.npy").astype(np.float64); ax = np.load("axis.npy")
xb, prof = np.load("profile_adapt.npy"); R = lambda v: np.interp(v, xb, prof)
dy = s[:,1] - np.polyval(ax[0], s[:,0]); dz = s[:,2] - np.polyval(ax[1], s[:,0])
TRUN = {+1: (117.765, -3.399, 4.523, 21.233), -1: (117.048, -5.076, 4.391, 21.098)}
trun = np.zeros(len(x), bool)
for sd, (xc, zc, rt, ye) in TRUN.items():
    trun |= (np.sign(dy) == sd) & (np.abs(dy) > 5) & (np.hypot(x - xc, dz - zc) < rt + 1.2)
res = rho - R(x)
ok = ~trun & (np.abs(th) <= 113) & (x > 0.4) & (x < 208.3) & (np.abs(res) < 1.2)
DX, DT = 0.05, 0.2
nx, nt = int(round(220 / DX)) + 1, int(round(360 / DT)) + 1
ix = np.clip(np.round(x[ok] / DX).astype(int), 0, nx - 1); it = np.clip(np.round((th[ok] + 180) / DT).astype(int), 0, nt - 1)
S = np.zeros((nx, nt), np.float64); W = np.zeros((nx, nt), np.float64)
np.add.at(S, (ix, it), res[ok]); np.add.at(W, (ix, it), 1)
def nconv(sx, st):
    a = gaussian_filter(S, (sx, st), mode=("nearest", "wrap")); b = gaussian_filter(W, (sx, st), mode=("nearest", "wrap"))
    return a / np.maximum(b, 1e-9), b
fine, wf = nconv(0.8, 0.8)
coarse, wc = nconv(1.0 / DX, np.degrees(1.0 / 12) / DT)
H = fine - coarse
H[wf < 0.15] = np.nan
H = np.clip(H, -0.25, 0.25)
# huecos (muñones, bordes): copiar textura desplazada 15 mm en x
for sh in (300, -300, 600, -600):
    m = np.isnan(H); H[m] = np.roll(H, sh, axis=0)[m]
H[np.isnan(H)] = 0
L = np.where(wc > 0.05, coarse, 0)
from scipy.ndimage import median_filter
np.save("tex_H.npy", H.astype(np.float32)); np.save("tex_Hm.npy", median_filter(H, size=(3, 3)).astype(np.float32))  # sin moteado del borde del escaneo; np.save("tex_L.npy", L.astype(np.float32))
print("H std (zona valida)", np.nanstd(H[wf >= 0.15]).round(4), "p1/p99", np.percentile(H[wf >= 0.15], [1, 99]).round(3))
cols = slice(int((180 + 30) / DT), int((180 + 113) / DT))
plt.figure(figsize=(22, 6)); plt.imshow(H[:, cols].T, cmap="gray", vmin=-0.12, vmax=0.12, aspect="auto", origin="lower", extent=[0, 220, 30, 113])
plt.xlabel("x mm"); plt.ylabel("θ"); plt.tight_layout(); plt.savefig("tex_H.png", dpi=70)
