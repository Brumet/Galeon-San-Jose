import pymeshlab as ml, numpy as np, struct, sys
ax = np.load("axis.npy"); D, C, OUT = sys.argv[1:4]
def zona(src, name, cond, target):
    ms = ml.MeshSet(); ms.load_new_mesh(src); V = ms.current_mesh().vertex_matrix(); F = ms.current_mesh().face_matrix()
    m = cond(V); kf = m[F].all(1); idx = np.where(m)[0]; rm = -np.ones(len(V), int); rm[idx] = np.arange(len(idx))
    o = ml.MeshSet(); o.add_mesh(ml.Mesh(vertex_matrix=V[m], face_matrix=rm[F[kf]])); o.meshing_remove_unreferenced_vertices()
    if o.current_mesh().face_number() > target:
        o.meshing_decimation_quadric_edge_collapse(targetfacenum=target, qualitythr=0.5, preserveboundary=True, preservenormal=True, optimalplacement=True, planarquadric=True)
    V2 = o.current_mesh().vertex_matrix(); F2 = o.current_mesh().face_matrix().astype(np.uint32)
    mn = V2.min(0); sc = (V2.max(0) - mn) / 65535.0; q = np.round((V2 - mn) / sc).astype(np.uint16)
    with open(f"{OUT}/{name}.wasm", "wb") as f:
        f.write(struct.pack("<II6f", len(q), len(F2), *mn, *sc)); f.write(q.tobytes()); f.write(F2.tobytes())
    print(name, len(F2), flush=True)
th = lambda V: np.degrees(np.arctan2(V[:,1] - np.polyval(ax[0], V[:,0]), V[:,2] - np.polyval(ax[1], V[:,0])))
zona(D, "completo", lambda V: np.ones(len(V), bool), 450000)
zona(D, "asas_D", lambda V: (V[:,0] > 110) & (V[:,0] < 136) & (V[:,2] > 26), 450000)
zona(C, "asas_C", lambda V: (V[:,0] > 110) & (V[:,0] < 136) & (V[:,2] > 26), 450000)
zona(D, "cascabel", lambda V: V[:,0] > 203.5, 450000)
zona(D, "escudos", lambda V: (V[:,0] > 14) & (V[:,0] < 106) & (np.abs(th(V)) < 60), 650000)
zona(D, "culata", lambda V: (V[:,0] > 194) & (V[:,0] < 206), 550000)
