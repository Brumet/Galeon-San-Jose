"""Delfín "como salió del molde": consenso de las dos asas gemelas, solo con datos medidos del escaneo.

Las asas son del mismo molde y están espejadas respecto al plano vertical del eje.
1. Se refleja el asa −Y y se alinea (ICP rígido) sobre la +Y.
2. Diferencia firmada entre ambas, suavizada (~0.4 mm) para detectar grumos localizados, no matices de detalle.
3. Donde coinciden se conservan las dos mediciones (Poisson las promedia). Donde una sobresale más de
   UMBRAL respecto a la otra (concreción), se descartan sus puntos en esa zona.
4. El consenso se coloca en los dos lados (espejado y con la transformación inversa).
"""
import numpy as np
from scipy.spatial import cKDTree

UMBRAL = 0.18        # mm de grumo a partir del cual manda la otra asa
R_SUAVE = 0.4        # mm, radio para suavizar la diferencia
R_CORTE = 0.6        # mm sobre la caña: por encima es "asa"; por debajo queda el cuerpo real
X_ASAS = (100.0, 135.0)


def _icp(P, Q, it=50):
    tq = cKDTree(Q); cur = P + (Q.mean(0) - P.mean(0))
    for _ in range(it):
        d, i = tq.query(cur, workers=-1); g = d < np.percentile(d, 85)
        a, b = cur[g], Q[i[g]]; ca, cb = a.mean(0), b.mean(0)
        U, _, Vt = np.linalg.svd((a - ca).T @ (b - cb))
        D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
        Rk = Vt.T @ D @ U.T
        cur = (cur - ca) @ Rk.T + cb
    # transformación total P -> cur (Procrustes exacto)
    cp, cc = P.mean(0), cur.mean(0)
    U, _, Vt = np.linalg.svd((P - cp).T @ (cur - cc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    Rt = Vt.T @ D @ U.T
    return Rt, cc - cp @ Rt.T


def _sd(P, NP, Q):
    """Distancia firmada de Q respecto a P, sobre la normal de P (+ = Q más afuera)."""
    d, i = cKDTree(Q).query(P, workers=-1)
    return ((Q[i] - P) * NP).sum(1), d


def _suave(P, v, r):
    t = cKDTree(P); nb = t.query_ball_point(P, r, workers=-1)
    return np.array([np.mean(v[n]) for n in nb])


def consenso(V, N, cy, cz, R, log=print):
    """Devuelve (keep_mask, P_nuevos, N_nuevos): quita las asas reales y añade el consenso en los dos lados."""
    x = V[:, 0]; dy = V[:, 1] - cy(x); dz = V[:, 2] - cz(x)
    r = np.hypot(dy, dz) - R(x); th = np.degrees(np.arctan2(dy, dz))
    asa = (x > X_ASAS[0]) & (x < X_ASAS[1]) & (r > R_CORTE) & (dz > 0) & (np.abs(th) < 60)
    iA = np.where(asa & (dy > 0))[0]; iB = np.where(asa & (dy < 0))[0]
    A, NA = V[iA], N[iA]
    Bm = V[iB].copy(); Bm[:, 1] = 2 * cy(Bm[:, 0]) - Bm[:, 1]
    NBm = N[iB].copy(); NBm[:, 1] *= -1
    Rt, t = _icp(Bm, A)
    B, NB = Bm @ Rt.T + t, NBm @ Rt.T
    sdA, dA = _sd(A, NA, B)          # + : B sobresale respecto a A
    sdB, dB = _sd(B, NB, A)          # + : A sobresale respecto a B
    sA = _suave(A, np.where(dA < 0.6, sdA, 0), R_SUAVE)
    sB = _suave(B, np.where(dB < 0.6, sdB, 0), R_SUAVE)
    keepA = ~(sA < -UMBRAL)          # A tiene grumo (B queda por dentro)
    keepB = ~(sB < -UMBRAL)
    log(f"asas: A {len(A)} B {len(B)} | mediana |dif| {np.median(np.abs(sdA[dA < 0.6])):.3f} mm | "
        f"descartados por grumo A {(~keepA).sum()} B {(~keepB).sum()}")
    C = np.vstack([A[keepA], B[keepB]]); NC = np.vstack([NA[keepA], NB[keepB]])
    # lado −Y: inversa de la alineación y espejo
    Cb = (C - t) @ Rt; NCb = NC @ Rt
    Cb[:, 1] = 2 * cy(Cb[:, 0]) - Cb[:, 1]; NCb[:, 1] *= -1
    keep = ~asa
    return keep, np.vstack([C, Cb]), np.vstack([NC, NCb]), dict(A=A, B=B, sA=sA, sB=sB, keepA=keepA, keepB=keepB, Rt=Rt, t=t)
