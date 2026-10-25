"""Escudos y texto de la culata: profundiza el relieve medido (sin suavizar, para no perder legibilidad).

Variante R2 elegida en el paso 4: V + 0.7*(L1 - L4) + 0.9*(L4 - L40), con Lk = laplaciano de k pasadas.
Así se marcan las letras y las molduras de los escudos y el relieve gana profundidad para el dorado.
Uso: python escudos_post.py entrada.ply salida.ply work_dir
"""
import pymeshlab as ml, numpy as np, sys, os

KF, KM = 0.7, 0.9
src, dst, work = sys.argv[1:4]
ms = ml.MeshSet(); ms.load_new_mesh(src)
V = ms.current_mesh().vertex_matrix().copy(); F = ms.current_mesh().face_matrix()
ax = np.load(os.path.join(work, "axis.npy")); xb, pr = np.load(os.path.join(work, "profile_adapt.npy"))
x = V[:, 0]; dy = V[:, 1] - np.polyval(ax[0], x); dz = V[:, 2] - np.polyval(ax[1], x)
r = np.hypot(dy, dz) - np.interp(x, xb, pr); ath = np.abs(np.degrees(np.arctan2(dy, dz)))
ramp = lambda v, a, b: np.clip((v - a) / (b - a), 0, 1)
# banda de los cinco escudos (parte superior de la caña) y faja de la culata con la inscripción
w_esc = ramp(x, 16, 18) * (1 - ramp(x, 102, 104)) * (1 - ramp(ath, 42, 50))
w_cul = ramp(x, 195.5, 196.5) * (1 - ramp(x, 204.0, 204.6)) * (1 - ramp(ath, 112, 118))
peso = np.maximum(w_esc, w_cul) * (np.abs(r) < 1.0)        # solo la piel, no asas ni muñones
sel = peso > 0
for zona in (sel & (x < 150), sel & (x > 150)):
    idx = np.where(zona)[0]; rm = -np.ones(len(V), int); rm[idx] = np.arange(len(idx))
    m = ml.MeshSet(); m.add_mesh(ml.Mesh(vertex_matrix=V[idx], face_matrix=rm[F[zona[F].all(1)]]))
    m.meshing_remove_unreferenced_vertices()
    f = rm[F[zona[F].all(1)]]; used = np.unique(f); idx = idx[used]
    V0 = m.current_mesh().vertex_matrix(); F0 = m.current_mesh().face_matrix()
    def lap(it):
        q = ml.MeshSet(); q.add_mesh(ml.Mesh(vertex_matrix=V0, face_matrix=F0))
        q.apply_coord_laplacian_smoothing(stepsmoothnum=it, boundary=False, cotangentweight=True)
        return q.current_mesh().vertex_matrix()
    L1, L4, L40 = lap(1), lap(4), lap(40)
    V[idx] += peso[idx, None] * (KF * (L1 - L4) + KM * (L4 - L40))
    print("zona", "escudos" if x[idx].mean() < 150 else "culata", len(idx), "vértices", flush=True)
out = ml.MeshSet(); out.add_mesh(ml.Mesh(vertex_matrix=V, face_matrix=F))
tm = out.get_topological_measures()
print("boundary", tm['boundary_edges'], "nm", tm['non_two_manifold_edges'], "parts", tm['connected_components_number'])
out.save_current_mesh(dst, binary=True)
