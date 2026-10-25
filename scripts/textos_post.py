"""Textos legibles con cantos vivos, a partir del relieve medido (sin tipografía ni letras inventadas).

En cada zona de texto se despliega la malla en (x, arco), se separa el relieve de las letras del fondo y se
umbraliza con contraste local. Se descartan las manchas pequeñas (corrosión) y cada trazo medido se convierte
en una meseta de altura uniforme con borde nítido, en su posición exacta.
Zona: inscripción de la culata (IVAN TER HORST ME FE[CIT] / SEVILLA ANNO 1654), el texto legible.
Uso: python textos_post.py entrada.ply salida.ply work_dir
"""
import pymeshlab as ml, numpy as np, sys, os
from scipy.ndimage import gaussian_filter, label, map_coordinates, find_objects

DX = 0.02
MAX_TRAZO = 1.4      # mm: un trazo de letra no es más largo (las letras miden ~0.9 mm)
ZONAS = [  # nombre, x mín/máx (mm), |θ| mín/máx (°), lado (+1, −1 o 0 = ambos), altura de letra final (mm), área mínima de trazo (mm²)
    ("culata", 201.0, 204.2, 12, 110, 0, 0.15, 0.06),
    # la placa se probó y se descartó: su relieve es tan débil que las letras salían partidas (se queda con R2)
]
src, dst, work = sys.argv[1:4]
ms = ml.MeshSet(); ms.load_new_mesh(src)
V = ms.current_mesh().vertex_matrix().copy(); F = ms.current_mesh().face_matrix()
ax = np.load(os.path.join(work, "axis.npy")); xb, pr = np.load(os.path.join(work, "profile_adapt.npy"))
cy = lambda x: np.polyval(ax[0], x); cz = lambda x: np.polyval(ax[1], x); R0 = lambda x: np.interp(x, xb, pr)
ss = lambda a, b, v: np.clip((v - a) / (b - a), 0, 1) ** 2 * (3 - 2 * np.clip((v - a) / (b - a), 0, 1))

for nombre, xa, xb_, ta, tb, lado, prof, amin in ZONAS:
    x = V[:, 0]; dy = V[:, 1] - cy(x); dz = V[:, 2] - cz(x); rho = np.hypot(dy, dz)
    th = np.degrees(np.arctan2(dy, dz)); res = rho - R0(x)
    sel = (x > xa) & (x < xb_) & (np.abs(th) >= ta) & (np.abs(th) <= tb) & (np.abs(res) < 0.6)
    if lado: sel &= np.sign(th) == lado
    i = np.where(sel)[0]; s = np.radians(th[i]) * R0(x[i])
    gx = ((x[i] - xa) / DX); s0 = s.min(); gs = (s - s0) / DX
    H = np.zeros((int(gx.max()) + 2, int(gs.max()) + 2)); W = np.zeros_like(H)
    np.add.at(H, (gx.astype(int), gs.astype(int)), res[i]); np.add.at(W, (gx.astype(int), gs.astype(int)), 1)
    nc = lambda sg: gaussian_filter(H, sg) / np.maximum(gaussian_filter(W, sg), 1e-9)
    fino, base = nc(1.0), nc(0.25 / DX)
    rel = fino - base
    sd = np.sqrt(np.maximum(gaussian_filter(rel ** 2, 0.5 / DX), 1e-12))
    z = rel / sd
    M = ss(0.4, 1.1, z) * (rel > 0.02)
    lab, n = label(M > 0.5)
    area = np.bincount(lab.ravel()) * DX * DX
    # solo trazos con tamaño de letra: fuera manchas pequeñas y crestas largas (bordes de molduras)
    ext = np.zeros(n + 1)
    for k, sl in enumerate(find_objects(lab), start=1):
        if sl is not None: ext[k] = max(sl[0].stop - sl[0].start, sl[1].stop - sl[1].start) * DX
    buenos = np.where((area >= amin) & (ext <= MAX_TRAZO))[0]; buenos = buenos[buenos > 0]
    M = M * np.isin(lab, buenos)
    M = gaussian_filter(M, 0.7)                                  # borde nítido pero sin escalones de rejilla
    destino = base + prof * M + 0.25 * rel * (1 - M)               # meseta de letra + algo de textura en el fondo
    D = destino - fino
    d = map_coordinates(D, [gx, gs], order=1, mode="nearest")
    borde = ss(0, 0.2, x[i] - xa) * ss(0, 0.2, xb_ - x[i]) * ss(0, 3, np.abs(th[i]) - ta + (3 if ta == 0 else 0)) * ss(0, 3, tb - np.abs(th[i]))
    rn = rho[i] + d * borde
    V[i, 1] = cy(x[i]) + rn * np.sin(np.radians(th[i])); V[i, 2] = cz(x[i]) + rn * np.cos(np.radians(th[i]))
    print(f"{nombre}: {len(i)} vértices, {len(buenos)} trazos, cobertura de letra {np.mean(M > 0.5) * 100:.1f} %", flush=True)

out = ml.MeshSet(); out.add_mesh(ml.Mesh(vertex_matrix=V, face_matrix=F))
tm = out.get_topological_measures()
print("boundary", tm['boundary_edges'], "nm", tm['non_two_manifold_edges'], "parts", tm['connected_components_number'])
out.save_current_mesh(dst, binary=True)
