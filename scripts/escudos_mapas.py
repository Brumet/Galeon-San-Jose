"""Despliega la banda de escudos (parte superior de la caña) y genera mapas de relieve para leer inscripciones.
Uso: python escudos_mapas.py <trabajo> <salida>   (lee scan_cyl.npy y profile_adapt.npy)"""
import numpy as np, sys, os
from scipy.ndimage import gaussian_filter
from PIL import Image
T, OUT = sys.argv[1:3]
c = np.load(os.path.join(T, "scan_cyl.npy")).astype(np.float64); x, rho, th = c.T
xb, pr = np.load(os.path.join(T, "profile_adapt.npy")); R = lambda v: np.interp(v, xb, pr)
X0, X1, TH = 16.0, 104.0, 42.0
m = (x > X0) & (x < X1) & (np.abs(th) < TH)
x, rho, th = x[m], rho[m], th[m]
s = np.radians(th) * R(x)                       # arco en mm (desplegado)
DX = 0.02
ix = ((x - X0) / DX).astype(int); S0 = s.min(); js = ((s - S0) / DX).astype(int)
H = np.zeros((ix.max() + 1, js.max() + 1)); W = np.zeros_like(H)
res = rho - R(x)
np.add.at(H, (ix, js), res); np.add.at(W, (ix, js), 1)
nc = lambda sg: gaussian_filter(H, sg) / np.maximum(gaussian_filter(W, sg), 1e-9)
fino, base = nc(1.6), nc(1.5 / DX)          # el escaneo tiene ~0.06 mm entre puntos: σ 1.6 celdas rellena la rejilla de 0.02
hp = fino - base                                # relieve sin la forma del cañón
valid = gaussian_filter(W, 3.0) > 0.02
hp[~valid] = 0
def norm(a, lo, hi):
    return (np.clip((a - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)
def guarda(a, nombre):                         # filas = arco (arriba del cañón al centro), columnas = x
    Image.fromarray(np.flipud(a.T)).save(os.path.join(OUT, nombre))
# contraste local: cada zona se normaliza con su propio relieve (σ ≈ 2 mm)
def local(a, sg=2.0 / DX):
    mu = gaussian_filter(a, sg); sd = np.sqrt(np.maximum(gaussian_filter((a - mu) ** 2, sg), 1e-12))
    return np.where(valid, (a - mu) / sd, 0)
hpl = local(hp)
guarda(norm(hpl, -2.5, 2.5), "relieve_altura.png")
gy, gx = np.gradient(gaussian_filter(hp, 0.8), DX)
for nombre, az in (("luz_NO", 315), ("luz_NE", 45), ("luz_SO", 225), ("luz_SE", 135)):
    a = np.radians(az); alt = np.radians(25)
    sl = np.arctan(np.hypot(gx, gy) * 1.0); asp = np.arctan2(-gy, gx)
    sh = np.sin(alt) * np.cos(sl) + np.cos(alt) * np.sin(sl) * np.cos(a - asp)
    guarda(norm(local(sh), -2.5, 2.5), f"relieve_{nombre}.png")
lap = gaussian_filter(hp, 1.2) - gaussian_filter(hp, 4)
guarda(norm(local(-lap), -2.5, 2.5), "relieve_bordes.png")
np.save(os.path.join(T, "escudos_hp.npy"), hp.astype(np.float32))
print("rejilla", hp.shape, "mm x:", X0, X1, "arco:", round(S0, 2), round(S0 + hp.shape[1] * DX, 2))
