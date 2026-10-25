import pymeshlab as ml, numpy as np
from metal import lap
def C(V0, F):   # aristas: filtro de normales muy suave + afilado medio
    m = ml.MeshSet(); m.add_mesh(ml.Mesh(vertex_matrix=V0, face_matrix=F))
    m.apply_coord_two_steps_smoothing(stepsmoothnum=1, normalthr=20, stepnormalnum=6, stepfitnum=6)
    TS = m.current_mesh().vertex_matrix()
    return TS + 0.8 * (lap(TS, F, 4) - lap(TS, F, 20)) + 0.4 * (V0 - lap(V0, F, 1))
def D(V0, F):   # HC laplaciano (conserva volumen) + afilado medio
    m = ml.MeshSet(); m.add_mesh(ml.Mesh(vertex_matrix=V0, face_matrix=F))
    for _ in range(2): m.apply_coord_hc_laplacian_smoothing()
    H = m.current_mesh().vertex_matrix()
    return H + 1.0 * (lap(H, F, 4) - lap(H, F, 20)) + 0.3 * (V0 - lap(V0, F, 1))
if __name__ == "__main__":
    for crop in ("m_asaN", "m_asaP", "m_casc"):
        ms = ml.MeshSet(); ms.load_new_mesh(crop + ".ply"); V0 = ms.current_mesh().vertex_matrix(); F = ms.current_mesh().face_matrix()
        for k, f in (("C", C), ("D", D)):
            o = ml.MeshSet(); o.add_mesh(ml.Mesh(vertex_matrix=f(V0, F), face_matrix=F)); o.save_current_mesh(f"{crop}_{k}.ply", binary=True)
        print(crop, "ok", flush=True)
