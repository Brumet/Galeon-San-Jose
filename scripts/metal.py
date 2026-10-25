"""Delfines 'metal fundido': suavizado que conserva aristas + afilado de banda media + pizca de textura fina."""
import pymeshlab as ml, numpy as np
VARIANTES = {"A": dict(ts=3, thr=30, km=0.8, kf=0.25), "B": dict(ts=8, thr=35, km=1.2, kf=0.3)}
def lap(V, F, it):
    m = ml.MeshSet(); m.add_mesh(ml.Mesh(vertex_matrix=V, face_matrix=F))
    m.apply_coord_laplacian_smoothing(stepsmoothnum=it, boundary=False, cotangentweight=True); return m.current_mesh().vertex_matrix()
def metal(V0, F, ts, thr, km, kf):
    m = ml.MeshSet(); m.add_mesh(ml.Mesh(vertex_matrix=V0, face_matrix=F))
    m.apply_coord_two_steps_smoothing(stepsmoothnum=ts, normalthr=thr, stepnormalnum=20, stepfitnum=20)
    TS = m.current_mesh().vertex_matrix()
    S4, S20 = lap(TS, F, 4), lap(TS, F, 20)
    return TS + km * (S4 - S20) + kf * (V0 - lap(V0, F, 1))
if __name__ == "__main__":
    for crop in ("m_casc",):
        ms = ml.MeshSet(); ms.load_new_mesh(crop + ".ply"); ms.meshing_remove_unreferenced_vertices(); ms.save_current_mesh(crop + ".ply", binary=True); V0 = ms.current_mesh().vertex_matrix(); F = ms.current_mesh().face_matrix()
        for k, p in VARIANTES.items():
            o = ml.MeshSet(); o.add_mesh(ml.Mesh(vertex_matrix=metal(V0, F, **p), face_matrix=F)); o.save_current_mesh(f"{crop}_{k}.ply", binary=True)
        print(crop, "ok", flush=True)
