"""Rehace la parte de abajo del cañón (franja sin datos del escaneo) girando el perfil real.

Entradas en `work` (generadas por fondo_analisis.py / perfil): scan_aligned.npy, axis.npy, prof_sectors.npy.
1. Carga el escaneo original con normales, lo lleva al sistema de la reconstrucción.
2. Quita los puntos de abajo (|θ| >= TH_CUT) y funde los de la banda TH_BLEND..TH_CUT hacia el perfil.
3. Añade puntos sintéticos: superficie de revolución del perfil real + parte baja del muñón +Y.
4. Poisson apantallado (mismos parámetros que la reconstrucción original).
"""
import pymeshlab as ml, numpy as np, sys, os, time

def log(*a): print(*a, flush=True)
src, work = sys.argv[1], sys.argv[2]
TH_BLEND, TH_CUT = 95.0, 112.0          # grados desde arriba
X0, X1 = 0.4, 206.6                      # tramo de revolución (boca .. inicio del cascabel)
DS = 0.04                                # separación de puntos sintéticos (mm)
# muñones, valores iniciales (trun.py): xc, zc relativo al eje, radio, extremo en |y|; se reajustan con la mitad de arriba
TRUN = {+1: (117.793, -3.636, 4.704, 21.682), -1: (117.105, -5.175, 4.489, 21.366)}
TA_BLEND, TA_CUT = 80.0, 100.0           # ángulo en el muñón (0 = arriba) a partir del cual se recrea

t0 = time.time()
ms = ml.MeshSet(); ms.load_new_mesh(src)
ms.meshing_remove_duplicate_vertices(); ms.meshing_remove_duplicate_faces()
ms.meshing_remove_connected_component_by_face_number(mincomponentsize=500)
ms.compute_normal_per_vertex()
V = ms.current_mesh().vertex_matrix(); N = ms.current_mesh().vertex_normal_matrix()
log("scan", len(V), round(time.time() - t0, 1))
M = np.array([[0, 1, 0], [0, 0, -1], [-1, 0, 0]], float)
V = V @ M.T; N = N @ M.T
V += np.load(os.path.join(work, "scan_aligned.npy")).astype(np.float64).min(0) - V.min(0)

ax = np.load(os.path.join(work, "axis.npy"))
cy = lambda x: np.polyval(ax[0], x); cz = lambda x: np.polyval(ax[1], x)
xb, prof = np.load(os.path.join(work, "profile_adapt.npy"))   # perfil_final.py + perfil_adaptativo.py
# molduras: se afila el perfil medido (máscara de enfoque 1D, σ 0.12 mm) para que los anillos tengan canto marcado
from scipy.ndimage import gaussian_filter1d
K_MOLD = 1.0
prof0 = prof.copy()
prof = prof0 + K_MOLD * (prof0 - gaussian_filter1d(prof0, 0.12 / (xb[1] - xb[0])))
DELTA = prof - prof0                                       # corrección que se aplica también a la piel real
R = lambda x: np.interp(x, xb, prof)
# canto del plato trasero: el escaneo deja láminas sueltas; se recrea completo (todo el contorno) con el perfil
RIM = (204.3, 206.6)
# plato trasero (cara casi vertical): ahí se corrige en x en vez de en radio
XS0 = 204.6   # (con X1 = 205.8 el plato trasero queda con datos reales del escaneo)
_m = (xb >= XS0) & (xb <= X1)
_pm = np.minimum.accumulate(prof[_m]) - np.arange(_m.sum()) * 1e-7   # estrictamente decreciente
Xinv = lambda r: np.interp(r, _pm[::-1], xb[_m][::-1])

# textura real (textura.py): H detalle fino, L bultos; rejilla x 0.05 mm, θ 0.2° desde -180
TH = np.load(os.path.join(work, "tex_Hm.npy")); TL = np.load(os.path.join(work, "tex_L.npy"))
from scipy.ndimage import map_coordinates
def tex(G, xx, tt):
    return map_coordinates(G, [xx / 0.05, (tt + 180) / 0.2], order=1, mode="nearest")
def pingpong(v, a, b):
    t = np.mod(v - a, 2 * (b - a)); return a + np.where(t > b - a, 2 * (b - a) - t, t)
# zonas donantes de corrosión uniforme (sin escudos, letras ni concreciones)
DON = {+1: ((35, 82), (52, 108)), -1: ((36, 82), (-105, -55))}
def donor(xx, tt, side):
    (xa, xb_), (ta, tb_) = DON[side]
    u = np.abs(tt) - TH_CUT                               # 0..68° bajo la unión
    xs_ = pingpong(xx + (0 if side > 0 else 11.0), xa, xb_)
    ts_ = pingpong(ta + u * 0.92, ta, tb_)
    return tex(TH, xs_, ts_)
TILE_X, TILE_U, FADE = 30.0, 26.0, 0.25     # parches (mm, grados) y fracción de fundido entre parches
def bottom_tex(xx, tt):
    """Textura de la parte de abajo: parches aleatorios de las zonas donantes, sin repetición visible.
    Mezcla entre parches con pesos normalizados en energía (no baja la intensidad de la textura)."""
    rng = np.random.default_rng(1654)
    u = np.abs(tt) - TH_BLEND + 4 + 120                   # grados desde el inicio del fundido (+120: positivo también arriba, en el canto)
    kx0 = np.floor(xx / TILE_X).astype(int); ku0 = np.floor(u / TILE_U).astype(int)
    nkx, nku = int(xx.max() / TILE_X) + 3, int(200 / TILE_U) + 3
    side = rng.integers(0, 2, (nkx + 1, nku + 1)) * 2 - 1
    offx = rng.random((nkx + 1, nku + 1)); offu = rng.random((nkx + 1, nku + 1))
    acc = np.zeros(len(xx)); wsum = np.zeros(len(xx))
    for dx in (-1, 0, 1):
        for du in (-1, 0, 1):
            kx = np.clip(kx0 + dx, 0, nkx); ku = np.clip(ku0 + du, 0, nku)
            # peso: meseta en el parche con bordes suaves que se solapan con los vecinos
            fx = (xx - kx * TILE_X) / TILE_X; fu = (u - ku * TILE_U) / TILE_U
            wx = smooth(-FADE, FADE, fx) * (1 - smooth(1 - FADE, 1 + FADE, fx))
            wu = smooth(-FADE, FADE, fu) * (1 - smooth(1 - FADE, 1 + FADE, fu))
            w = wx * wu
            if not np.any(w > 0): continue
            sd = side[kx, ku]
            out = np.zeros(len(xx))
            for sgn in (+1, -1):
                m = (w > 0) & (sd == sgn)
                if not m.any(): continue
                (xa, xb_), (ta, tb_) = DON[sgn]
                span_x = TILE_X * (1 + 2 * FADE); span_u = TILE_U * (1 + 2 * FADE)
                x0 = xa + offx[kx[m], ku[m]] * max(xb_ - xa - span_x, 0)
                t0 = ta + offu[kx[m], ku[m]] * max(tb_ - ta - span_u, 0)
                out[m] = tex(TH, x0 + (fx[m] + FADE) * TILE_X, t0 + (fu[m] + FADE) * TILE_U)
            acc += w * out; wsum += w * w
    return acc / np.sqrt(np.maximum(wsum, 1e-9))

def cyl(p):
    dy = p[:, 1] - cy(p[:, 0]); dz = p[:, 2] - cz(p[:, 0])
    return np.hypot(dy, dz), np.degrees(np.arctan2(dy, dz)), dy, dz

def smooth(a, b, v):
    t = np.clip((v - a) / (b - a), 0, 1); return t * t * (3 - 2 * t)

def in_trunnion(x, dy, dz, side, pad):
    xc, zc, r, ye = TRUN[side]
    return (np.sign(dy) == side) & (np.abs(dy) > 5) & (np.abs(dy) < ye + pad) & (np.hypot(x - xc, dz - zc) < r + pad)

rho, th, dy, dz = cyl(V); x = V[:, 0]; ath = np.abs(th)
trun_any = in_trunnion(x, dy, dz, +1, 0.6) | in_trunnion(x, dy, dz, -1, 0.6)
body = (x > X0) & (x < X1) & ~trun_any
keep = ~(body & (ath >= TH_CUT))
keep &= ~((x > RIM[0]) & (x < RIM[1]) & ~trun_any)          # canto del plato: fuera todo lo real
# banda de fusión: llevar el radio real hacia el perfil
bl = body & (ath >= TH_BLEND) & (ath < TH_CUT)
w = smooth(TH_BLEND, TH_CUT, ath[bl])
sx = smooth(XS0 - 0.3, XS0 + 0.3, x[bl])                  # 0 caña, 1 plato trasero
rn = rho[bl] - w * (1 - sx) * tex(TL, x[bl], th[bl])      # quita bultos, deja la textura
xn = x[bl] + w * sx * (Xinv(rho[bl]) - x[bl])
V[bl, 0] = xn
V[bl, 1] = cy(xn) + rn * np.sin(np.radians(th[bl])); V[bl, 2] = cz(xn) + rn * np.cos(np.radians(th[bl]))
# muñones: reajuste con la mitad de arriba, se recrea la mitad de abajo y la cara extrema
def fit_circle(u, v):
    A = np.c_[2 * u, 2 * v, np.ones_like(u)]; q, *_ = np.linalg.lstsq(A, u * u + v * v, rcond=None)
    return q[0], q[1], np.sqrt(q[2] + q[0] ** 2 + q[1] ** 2)
nbt = 0
for sd in (+1, -1):
    xc, zc, rt, ye = TRUN[sd]
    tb = in_trunnion(x, dy, dz, sd, 0.6)
    ta = np.degrees(np.arctan2(x - xc, dz - zc)); ata = np.abs(ta)
    fz = tb & (np.abs(dy) > ye - 1.0) & (np.hypot(x - xc, dz - zc) < rt - 0.6) & (ata < 90)
    ye = float(np.median(np.abs(dy[fz])))
    k = tb & (np.abs(dy) > R(x) + 1.0) & (np.abs(dy) < ye - 0.6) & (ata < 95)
    for _ in range(4):
        xc, zc, rt = fit_circle(x[k], dz[k])
        k &= np.abs(np.hypot(x - xc, dz - zc) - rt) < 0.3
    TRUN[sd] = (xc, zc, rt, ye)
    log("muñón", sd, "xc", round(xc, 3), "zc", round(zc, 3), "r", round(rt, 3), "cara", round(ye, 3))
    tb = in_trunnion(x, dy, dz, sd, 0.6)
    ta = np.degrees(np.arctan2(x - xc, dz - zc)); ata = np.abs(ta)
    face = tb & (np.abs(dy) > ye - 0.35)
    side_s = tb & ~face
    keep &= ~(side_s & (ata >= TA_CUT)) & ~face
    bt = side_s & (ata >= TA_BLEND) & (ata < TA_CUT)
    w = smooth(TA_BLEND, TA_CUT, ata[bt]); rr = np.hypot(x[bt] - xc, dz[bt] - zc); rn = rr + w * (rt - rr)
    V[bt, 0] = xc + rn * np.sin(np.radians(ta[bt]))
    V[bt, 2] = cz(V[bt, 0]) + zc + rn * np.cos(np.radians(ta[bt]))
    nbt += bt.sum()
log("real kept", keep.sum(), "removed", (~keep).sum(), "blended", bl.sum() + nbt)
# delfines de las asas: consenso de las dos gemelas (delfines_consenso.py)
DELFINES = "--sin-delfines" not in sys.argv
if DELFINES:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from delfines_consenso import consenso
    kd, DP, DN, diag = consenso(V, N, cy, cz, R, log)
    keep &= kd
    np.savez(os.path.join(work, "delfines_diag.npz"), **{k: v for k, v in diag.items()})
else:
    DP, DN = np.zeros((0, 3)), np.zeros((0, 3))
V, N = V[keep], N[keep]

# --- puntos sintéticos: revolución del perfil ---
xs = np.arange(X0, X1, 0.005); rs = R(xs)
seg = np.hypot(np.diff(xs), np.diff(rs)); s = np.r_[0, np.cumsum(seg)]
sq = np.arange(0, s[-1], DS); xq = np.interp(sq, s, xs); rq = np.interp(sq, s, rs)
tx = np.gradient(xq); tr = np.gradient(rq); tn = np.hypot(tx, tr); nx, nr = -tr / tn, tx / tn
SP, SN = [], []
for xi, ri, nxi, nri in zip(xq, rq, nx, nr):
    full = RIM[0] < xi < RIM[1]
    t0 = 0.0 if full else TH_CUT
    n = int(np.ceil(2 * np.radians(180 - t0) * ri / DS))
    a = np.radians(np.linspace(t0, 360 - t0, n, endpoint=not full))
    sa, ca = np.sin(a), np.cos(a)
    SP.append(np.c_[np.full(n, xi), cy(xi) + ri * sa, cz(xi) + ri * ca]); SN.append(np.c_[np.full(n, nxi), nri * sa, nri * ca])
SP = np.vstack(SP); SN = np.vstack(SN)
_, _, sdy, sdz = cyl(SP)
inside = np.zeros(len(SP), bool)
for sd in (+1, -1):
    xc_, zc_, r_, ye_ = TRUN[sd]
    inside |= (np.sign(sdy) == sd) & (np.abs(sdy) < ye_) & (np.hypot(SP[:, 0] - xc_, sdz - zc_) < r_)
SP, SN = SP[~inside], SN[~inside]
log("synthetic body", len(SP))

# --- muñones: cilindro (mitad de abajo) y cara extrema completa ---
TP, TN = [], []
for sd in (+1, -1):
    xc, zc, rt, ye = TRUN[sd]
    na = int(np.ceil(np.radians(360 - 2 * TA_CUT) * rt / DS))
    a = np.radians(np.linspace(TA_CUT, 360 - TA_CUT, na))
    for yy in np.arange(6.0, ye, DS):
        px = xc + rt * np.sin(a); pzr = zc + rt * np.cos(a)
        p = np.c_[px, cy(px) + sd * yy, cz(px) + pzr]
        m = np.hypot(yy, pzr) > R(px)                     # solo lo que queda fuera del cuerpo
        TP.append(p[m]); TN.append(np.c_[np.sin(a), np.zeros(na), np.cos(a)][m])
    for t in np.r_[0, np.arange(DS, rt + DS / 2, DS)]:
        k = max(1, int(np.ceil(2 * np.pi * t / DS)))
        a2 = np.linspace(0, 2 * np.pi, k, endpoint=False)
        px = xc + t * np.sin(a2)
        TP.append(np.c_[px, cy(px) + sd * ye, cz(px) + zc + t * np.cos(a2)]); TN.append(np.tile([0, sd, 0], (k, 1)))
TP = np.vstack(TP); TN = np.vstack(TN).astype(float)
log("synthetic trunnion", len(TP))

allP = np.vstack([V, SP, TP, DP]); allN = np.vstack([N, SN, TN, DN])
allN /= np.linalg.norm(allN, axis=1, keepdims=True)
np.save(os.path.join(work, "synth_only.npy"), np.vstack([SP, TP]).astype(np.float32))
ms2 = ml.MeshSet(); ms2.add_mesh(ml.Mesh(vertex_matrix=allP, v_normals_matrix=allN))
log("poisson start", len(allP), round(time.time() - t0, 1))
# pointweight alto y pocas muestras por nodo: la superficie sigue más de cerca los puntos (bordes menos redondeados)
PW = 10 if "--fiel" in sys.argv else 6; SPN = 1.0 if "--fiel" in sys.argv else 1.5
log("poisson pointweight", PW, "samplespernode", SPN)
ms2.generate_surface_reconstruction_screened_poisson(depth=12, samplespernode=SPN, pointweight=PW, preclean=True)
ms2.meshing_remove_connected_component_by_diameter(mincomponentdiag=ml.PercentageValue(10))
m2 = ms2.current_mesh(); PV = m2.vertex_matrix().copy(); PF = m2.face_matrix()
rho, th, dy, dz = cyl(PV); x = PV[:, 0]; ath = np.abs(th)
nt = np.zeros(len(PV), bool)
for sd in (+1, -1): nt |= in_trunnion(x, dy, dz, sd, 0.9)
w = smooth(TH_BLEND - 3, TH_BLEND + 9, ath) * smooth(X0, X0 + 0.3, x) * (1 - smooth(X1 - 0.3, X1, x)) * ~nt
w = np.maximum(w, smooth(RIM[0], RIM[0] + 0.3, x) * (1 - smooth(RIM[1] - 0.3, RIM[1], x)) * ~nt)
sx = smooth(XS0 - 0.3, XS0 + 0.3, x)
w *= np.where(sx < 0.5, np.abs(rho - R(x)) < 0.35, np.abs(x - Xinv(rho)) < 0.35)   # solo la piel de revolución
hb = bottom_tex(x, th)
rn = rho + w * (1 - sx) * (R(x) + hb - rho)
# piel real (fuera de escudos, asas, muñones y canto): mismo afilado de molduras
piel = (np.abs(rho - R(x)) < 0.35) & ~nt & ~((x > 14) & (x < 106) & (ath < 55)) & (x < RIM[0])
rn += (1 - w) * (1 - sx) * np.interp(x, xb, DELTA) * piel
xn = x + w * sx * (Xinv(rho) + hb - x)
PV[:, 0] = xn
PV[:, 1] = cy(xn) + rn * np.sin(np.radians(th)); PV[:, 2] = cz(xn) + rn * np.cos(np.radians(th))
npj = int((w > 0).sum())
for sd in (+1, -1):
    xc, zc, rt, ye = TRUN[sd]
    rho, th, dy, dz = cyl(PV); x = PV[:, 0]
    tb = in_trunnion(x, dy, dz, sd, 0.4)
    ta = np.degrees(np.arctan2(x - xc, dz - zc)); rr = np.hypot(x - xc, dz - zc)
    ws = smooth(TA_CUT, TA_CUT + 10, np.abs(ta)) * tb * (np.abs(dy) > R(x) + 0.4) * (np.abs(dy) < ye - 0.35) * (np.abs(rr - rt) < 0.25)
    ht = donor(60 + np.abs(dy) * 1.0, np.sign(ta) * (TH_CUT + (np.abs(ta) - TA_CUT) * rt / 13.0), +1 if sd > 0 else -1)
    rn = rr + ws * (rt + ht - rr)
    PV[:, 0] = np.where(ws > 0, xc + rn * np.sin(np.radians(ta)), PV[:, 0])
    PV[:, 2] = np.where(ws > 0, cz(PV[:, 0]) + zc + rn * np.cos(np.radians(ta)), PV[:, 2])
    wf = tb * (np.abs(np.abs(dy) - ye) < 0.15) * (rr < rt - 0.35)
    hf = donor(70 + (x - xc) * 1.0, TH_CUT + 30 + np.degrees((dz - zc) / 13.0), -1 if sd > 0 else +1)
    PV[wf > 0, 1] = cy(PV[wf > 0, 0]) + sd * (ye + hf[wf > 0])
    npj += int((ws > 0).sum() + (wf > 0).sum())
log("proyectados", npj)
ms2 = ml.MeshSet(); ms2.add_mesh(ml.Mesh(vertex_matrix=PV, face_matrix=PF))
tm = ms2.get_topological_measures(); g = ms2.get_geometric_measures()
b = ms2.current_mesh().bounding_box()
log("poisson faces", ms2.current_mesh().face_number(), "boundary", tm['boundary_edges'], "nm", tm['non_two_manifold_edges'],
    "parts", tm['connected_components_number'], "vol", round(g.get('mesh_volume', 0)), "bbox", np.round(b.min(), 3), np.round(b.max(), 3),
    round(time.time() - t0, 1))
ms2.save_current_mesh(os.path.join(work, "canon_fondo_full.ply"), binary=True)
log("saved")
