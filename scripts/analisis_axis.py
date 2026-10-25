import numpy as np
s = np.load("scan_aligned.npy").astype(np.float64)
x = s[:,0]
# primer eje aproximado: centro de bbox YZ por rodaja, usando puntos con |y|<17 (fuera muñones)
def fit_circle(y, z):
    A = np.c_[2*y, 2*z, np.ones_like(y)]; b = y*y + z*z
    c, *_ = np.linalg.lstsq(A, b, rcond=None); r = np.sqrt(c[2] + c[0]**2 + c[1]**2); return c[0], c[1], r
bins = np.arange(0, 220.01, 1.0)
rows = []
for a, b in zip(bins[:-1], bins[1:]):
    m = (x >= a) & (x < b)
    p = s[m]
    if len(p) < 500: continue
    cy, cz, r = fit_circle(p[:,1], p[:,2])
    for _ in range(4):  # robusto: descartar puntos lejos del círculo (muñones, delfines)
        rr = np.hypot(p[:,1]-cy, p[:,2]-cz); k = np.abs(rr - np.median(rr)) < 0.6
        cy, cz, r = fit_circle(p[k,1], p[k,2])
    th = np.degrees(np.arctan2(p[:,1]-cy, p[:,2]-cz))  # 0 = arriba, ±180 abajo
    h, _ = np.histogram(np.abs(th), bins=36, range=(0,180))
    empty = np.where(h == 0)[0]
    gap = 180 - (empty.min()*5) if len(empty) else 0
    rows.append((a+0.5, cy, cz, r, gap))
rows = np.array(rows)
np.save("slices.npy", rows)
for r in rows[::5]: print(np.round(r, 3))
# eje = recta robusta en y(x), z(x)
for j, n in ((1,'y'),(2,'z')):
    p = np.polyfit(rows[:,0], rows[:,j], 1); res = rows[:,j] - np.polyval(p, rows[:,0])
    k = np.abs(res) < 0.15; p = np.polyfit(rows[k,0], rows[k,j], 1); print(n, p, "inliers", k.mean())
