import pymeshlab as ml, numpy as np, sys
for tag, f in (("antes", "canon_base_v12.ply"), ("despues", "canon_fondo_full.ply")):
    ms = ml.MeshSet(); ms.load_new_mesh(f); V = ms.current_mesh().vertex_matrix(); F = ms.current_mesh().face_matrix()
    for side, sg in (("masY", 1), ("menosY", -1)):
        m = (V[:,0] > 112) & (V[:,0] < 133) & (V[:,2] > 29) & (np.sign(V[:,1]) == sg) & (np.abs(V[:,1]) < 9)
        kf = m[F].all(1); idx = np.where(m)[0]; rm = -np.ones(len(V), int); rm[idx] = np.arange(len(idx))
        o = ml.MeshSet(); o.add_mesh(ml.Mesh(vertex_matrix=V[m], face_matrix=rm[F[kf]]))
        o.save_current_mesh(f"cab_{side}_{tag}.ply", binary=True)
