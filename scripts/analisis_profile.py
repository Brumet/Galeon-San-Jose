import numpy as np
s = np.load("scan_aligned.npy").astype(np.float64)
sl = np.load("slices.npy")
# eje curvo suave: polinomio grado 3 robusto sobre centros de rodajas válidas (fuera muñones/extremos)
ok = (sl[:,0] > 8) & (sl[:,0] < 200) & ~((sl[:,0] > 105) & (sl[:,0] < 140))
cy = np.polyfit(sl[ok,0], sl[ok,1], 3); cz = np.polyfit(sl[ok,0], sl[ok,2], 3)
for j, c in ((1, cy), (2, cz)):
    print("res", np.round(sl[ok,j] - np.polyval(c, sl[ok,0]), 3)[::8])
np.save("axis.npy", np.vstack([cy, cz]))
x = s[:,0]; dy = s[:,1] - np.polyval(cy, x); dz = s[:,2] - np.polyval(cz, x)
rho = np.hypot(dy, dz); th = np.degrees(np.arctan2(dy, dz))  # 0 arriba
ath = np.abs(th)
np.save("scan_cyl.npy", np.c_[x, rho, th].astype(np.float32))
# perfil: mediana de rho por rodaja de 0.05 mm, ángulos laterales 35..100 grados
step = 0.05; xb = np.arange(0, 220 + step, step)
idx = np.floor(x / step).astype(int)
sel = (ath > 35) & (ath < 100)
trun = (x > 106) & (x < 136)
sel2 = sel & ~(trun & (ath > 45) & (ath < 135))
sel2 |= trun & (ath > 20) & (ath < 45)
prof = np.full(len(xb), np.nan); cnt = np.zeros(len(xb), int)
order = np.argsort(idx[sel2]); ii = idx[sel2][order]; rr = rho[sel2][order]
u, start, c = np.unique(ii, return_index=True, return_counts=True)
for k, a, n in zip(u, start, c):
    if n >= 15: prof[k] = np.median(rr[a:a+n]); cnt[k] = n
np.save("profile_raw.npy", np.c_[xb, prof, cnt])
print("profile nan", np.isnan(prof).sum(), "of", len(prof))
for X in range(0, 221, 10): print(X, np.round(prof[int(X/step)-1 if X else 0], 3))
# residuo por ángulo
pi = np.interp(x, xb[~np.isnan(prof)], prof[~np.isnan(prof)])
res = rho - pi
m = (x > 5) & (x < 200) & ~trun
for a in range(0, 180, 10):
    k = m & (ath >= a) & (ath < a + 10)
    if k.sum(): print(f"th {a:3d}-{a+10:3d} n {k.sum():7d} med {np.median(res[k]):+.3f} p10 {np.percentile(res[k],10):+.3f} p90 {np.percentile(res[k],90):+.3f}")
