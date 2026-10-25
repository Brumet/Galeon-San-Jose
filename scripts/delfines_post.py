"""Delfines de las asas "como salieron del molde", aplicado sobre la malla final (sin mezclar nubes).

Cada asa conserva su malla. Sus vértices se desplazan sobre la normal hacia el consenso con su gemela
(reflejada y alineada por ICP):
  - coinciden (|dif| < UMBRAL)           -> mitad de camino (promedio de las dos mediciones)
  - la gemela queda por dentro (-UMBRAL)  -> este vértice es un grumo: se hunde hasta la gemela
  - la gemela sobresale (+UMBRAL)         -> el grumo es de la gemela: no se mueve
Uso: python delfines_post.py entrada.ply salida.ply work_dir
"""
import pymeshlab as ml, numpy as np, sys, os
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delfines_consenso import _icp, X_ASAS

UMBRAL = 0.18          # mm
R_SUAVE = 0.18         # mm, suavizado de la diferencia (quita ruido de medición, no detalle del molde)
PIE = (0.5, 1.2)       # mm sobre la caña: transición del pie (0 = cuerpo intacto, 1 = asa)
R_CAMPO = 0.35         # mm, suavizado del campo de desplazamiento (evita pliegues en las transiciones)
MAX_D = 0.4            # mm, desplazamiento máximo

src, dst, work = sys.argv[1:4]
ms = ml.MeshSet(); ms.load_new_mesh(src); ms.compute_normal_per_vertex()
V = ms.current_mesh().vertex_matrix().copy(); N = ms.current_mesh().vertex_normal_matrix(); F = ms.current_mesh().face_matrix()
ax = np.load(os.path.join(work, "axis.npy")); xb, pr = np.load(os.path.join(work, "profile_adapt.npy"))
cy = lambda x: np.polyval(ax[0], x); cz = lambda x: np.polyval(ax[1], x); R = lambda x: np.interp(x, xb, pr)
x = V[:, 0]; dy = V[:, 1] - cy(x); dz = V[:, 2] - cz(x)
r = np.hypot(dy, dz) - R(x); th = np.degrees(np.arctan2(dy, dz))
zona = (x > X_ASAS[0]) & (x < X_ASAS[1]) & (r > PIE[0]) & (dz > 0) & (np.abs(th) < 60)
iA = np.where(zona & (dy > 0))[0]; iB = np.where(zona & (dy < 0))[0]
mir = lambda P: np.c_[P[:, 0], 2 * cy(P[:, 0]) - P[:, 1], P[:, 2]]
mirn = lambda Nn: Nn * [1, -1, 1]
Rt, t = _icp(mir(V[iB]), V[iA])
print(f"asas A {len(iA)} B {len(iB)} vértices")

def desplaza(P, NP, Q):
    """Desplazamiento de P (normales NP) hacia el consenso con Q (gemela, ya en el marco de P)."""
    d, i = cKDTree(Q).query(P, workers=-1)
    sd = ((Q[i] - P) * NP).sum(1); sd[d > 0.8] = 0
    tp = cKDTree(P); nb = tp.query_ball_point(P, R_SUAVE, workers=-1)
    s = np.array([sd[n].mean() for n in nb])
    t_in = np.clip((-s - UMBRAL * 0.7) / (UMBRAL * 0.6), 0, 1)       # grumo propio
    t_out = np.clip((s - UMBRAL * 0.7) / (UMBRAL * 0.6), 0, 1)       # grumo de la gemela
    k = 0.5 * (1 - t_in) * (1 - t_out) + 1.0 * t_in                  # fracción de s a recorrer
    print(f"  |dif| mediana {np.median(np.abs(s)):.3f} mm, grumo propio {(t_in > 0.5).mean()*100:.1f} %, grumo gemela {(t_out > 0.5).mean()*100:.1f} %")
    conf = np.clip(1 - (d - 0.25) / 0.25, 0, 1)                      # sin gemela fiable cerca -> no se toca
    D = NP * np.clip(k * s * conf, -MAX_D, MAX_D)[:, None]
    nb2 = tp.query_ball_point(P, R_CAMPO, workers=-1)                  # campo suave: sin pliegues
    return np.array([D[n].mean(0) for n in nb2])

QA = mir(V[iB]) @ Rt.T + t                       # gemela B en el marco de A
DA = desplaza(V[iA], N[iA], QA)
QB = mir((V[iA] - t) @ Rt)                       # gemela A en el marco de B
DB = desplaza(V[iB], N[iB], QB)
w = np.clip((r - PIE[0]) / (PIE[1] - PIE[0]), 0, 1)
# cascabel: simétrico izquierda-derecha (girado ~7° sobre el eje): consenso con su propio reflejo
X_CASC = (208.6, 209.4)                          # arranque desde el plato trasero (rampa)
iC = np.where(x > X_CASC[0])[0]
PC = V[iC]; ref = PC[PC[:, 0] > X_CASC[1]]
Rc, tc = _icp(mir(ref), ref)
DC = desplaza(PC, N[iC], mir(PC) @ Rc.T + tc)
wc = np.clip((x - X_CASC[0]) / (X_CASC[1] - X_CASC[0]), 0, 1)
V[iA] += DA * w[iA, None]; V[iB] += DB * w[iB, None]; V[iC] += DC * wc[iC, None]

# realce suave del relieve medido (máscara de enfoque multibanda): amplifica surcos y bordes existentes,
# sin la banda más fina (grano del escáner) para no amplificar ruido
KF, KM = 0.6, 0.5      # banda de bordes (suavizado 1->4 pasadas) y de escamas/volutas (4->20)
def realce(sel, peso):
    idx = np.where(sel)[0]; rm = -np.ones(len(V), int); rm[idx] = np.arange(len(idx))
    f = rm[F[sel[F].all(1)]]
    def lisa(it):
        m = ml.MeshSet(); m.add_mesh(ml.Mesh(vertex_matrix=V[idx], face_matrix=f))
        m.apply_coord_laplacian_smoothing(stepsmoothnum=it, boundary=False, cotangentweight=True)
        return m.current_mesh().vertex_matrix()
    S1, S4, S20 = lisa(1), lisa(4), lisa(20)
    V[idx] += peso[idx, None] * (KF * (S1 - S4) + KM * (S4 - S20))
asas_m = (x > X_ASAS[0] - 2) & (x < X_ASAS[1] + 2) & (dz > 0) & (np.abs(th) < 70) & (r > -0.5)
def metal(sel, peso, fn):
    """Aplica una variante 'metal' (metal2.py: C o D) al submodelo sel, mezclada con peso."""
    idx = np.where(sel)[0]; rm = -np.ones(len(V), int); rm[idx] = np.arange(len(idx))
    f = rm[F[sel[F].all(1)]]
    m = ml.MeshSet(); m.add_mesh(ml.Mesh(vertex_matrix=V[idx], face_matrix=f))
    m.meshing_remove_unreferenced_vertices()               # el filtro de dos pasos exige malla compacta
    keep = sel.copy(); keep[idx] = np.isin(np.arange(len(idx)), np.unique(f))
    idx2 = np.where(keep)[0]; nv = fn(m.current_mesh().vertex_matrix(), m.current_mesh().face_matrix())
    V[idx2] += peso[idx2, None] * (nv - V[idx2])
MODO = next((a.split("=")[1] for a in sys.argv if a.startswith("--metal=")), None)
if MODO:
    from metal2 import C as _C, D as _D
    fn = {"C": _C, "D": _D}[MODO]
    metal(asas_m, w * zona, fn); metal(x > X_CASC[0] - 0.6, wc, fn)
    print("variante metal", MODO, "aplicada a asas y cascabel")
elif "--sin-realce" not in sys.argv:
    realce(asas_m, w * zona)
    realce(x > X_CASC[0] - 0.6, wc)
    print("realce aplicado a asas y cascabel")
out = ml.MeshSet(); out.add_mesh(ml.Mesh(vertex_matrix=V, face_matrix=F))
tm = out.get_topological_measures()
print("boundary", tm['boundary_edges'], "nm", tm['non_two_manifold_edges'], "parts", tm['connected_components_number'])
out.save_current_mesh(dst, binary=True)
# control: diferencia entre las dos asas después
A2, B2 = V[iA], mir(V[iB]) @ Rt.T + t
d, _ = cKDTree(B2).query(A2, workers=-1)
print("dif entre asas después: mediana", np.median(d).round(3), "p95", np.percentile(d, 95).round(3))
