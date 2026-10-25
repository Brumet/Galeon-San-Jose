import pymeshlab as ml, numpy as np, os, time
from scipy.spatial import cKDTree
def log(*a): print(*a, flush=True)
ms = ml.MeshSet(); ms.load_new_mesh("canon_fondo_full.ply")
V = ms.current_mesh().vertex_matrix()
s = np.load("scan_aligned.npy").astype(np.float64); ax = np.load("axis.npy")
dy = s[:,1]-np.polyval(ax[0], s[:,0]); dz = s[:,2]-np.polyval(ax[1], s[:,0]); top = np.abs(np.degrees(np.arctan2(dy, dz))) < 90
d, _ = cKDTree(V).query(s[top][::5], workers=-1)
log("escaneo (mitad de arriba) -> malla nueva: mediana", np.median(d).round(4), "p95", np.percentile(d,95).round(4), "p99.9", np.percentile(d,99.9).round(3))
zmin = V[:,2].min(); log("zmin", zmin)
ms.compute_matrix_from_translation(traslmethod='XYZ translation', axisx=0, axisy=0, axisz=-zmin)
def topo(n):
    tm = ms.get_topological_measures(); b = ms.current_mesh().bounding_box()
    log(f"[{n}] faces {ms.current_mesh().face_number()} boundary {tm['boundary_edges']} nm {tm['non_two_manifold_edges']} parts {tm['connected_components_number']} bbox {np.round(b.min(),2)} {np.round(b.max(),2)}")
os.makedirs("out", exist_ok=True)
topo("v2_MAXIMO"); ms.save_current_mesh("out/canon_v2_resina_MAXIMO.stl", binary=True)
base = ms.current_mesh_id()
for target, name in ((4_000_000, "canon_v2_resina_4M"), (300_000, "canon_v2_fusion_300k"), (100_000, "canon_v2_fusion_100k")):
    ms.set_current_mesh(base); ms.generate_copy_of_current_mesh()
    ms.meshing_decimation_quadric_edge_collapse(targetfacenum=target, qualitythr=0.5, preserveboundary=True, preservenormal=True, preservetopology=True, optimalplacement=True, planarquadric=True, autoclean=True)
    ms.meshing_repair_non_manifold_edges(); ms.meshing_repair_non_manifold_vertices(); ms.meshing_close_holes(maxholesize=100); ms.meshing_re_orient_faces_coherently()
    topo(name); ms.save_current_mesh(f"out/{name}.stl", binary=True); ms.save_current_mesh(f"out/{name}.3mf")
