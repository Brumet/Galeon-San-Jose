"""Alisa las láminas del borde del escaneo en la cúpula trasera, junto al arranque inferior del cascabel.
Solo actúa en la mitad de abajo de la cúpula (x > 207.4) y no toca el volumen del cascabel.
Uso: python cupula_post.py entrada.ply salida.ply work_dir"""
import pymeshlab as ml, numpy as np, sys, os
src, dst, work = sys.argv[1:4]
ms = ml.MeshSet(); ms.load_new_mesh(src); V = ms.current_mesh().vertex_matrix().copy(); F = ms.current_mesh().face_matrix()
ax = np.load(os.path.join(work, "axis.npy"))
x = V[:, 0]; dy = V[:, 1] - np.polyval(ax[0], x); dz = V[:, 2] - np.polyval(ax[1], x)
ath = np.abs(np.degrees(np.arctan2(dy, dz))); rho = np.hypot(dy, dz)
casc = (np.abs(dy) < 2.4) & (rho < 9.5) & (x > 208.0)          # tubo del cascabel
zona = (x > 207.2) & (x < 209.2) & (ath > 45) & (rho > 4.0) & ~casc
peso = np.clip((x - 207.2) / 0.4, 0, 1) * np.clip((ath - 45) / 25, 0, 1)
idx = np.where(zona)[0]; rm = -np.ones(len(V), int); rm[idx] = np.arange(len(idx))
m = ml.MeshSet(); m.add_mesh(ml.Mesh(vertex_matrix=V[idx], face_matrix=rm[F[zona[F].all(1)]]))
m.apply_coord_laplacian_smoothing(stepsmoothnum=60, boundary=False, cotangentweight=True)
L = m.current_mesh().vertex_matrix()
V[idx] += peso[idx, None] * (L - V[idx])
o = ml.MeshSet(); o.add_mesh(ml.Mesh(vertex_matrix=V, face_matrix=F))
print("cúpula: vértices alisados", len(idx)); o.save_current_mesh(dst, binary=True)
