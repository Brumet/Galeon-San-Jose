import pymeshlab as ml, numpy as np
from scipy.spatial import cKDTree
ms = ml.MeshSet(); ms.load_new_mesh("canon_fondo_full.ply")
V = ms.current_mesh().vertex_matrix(); F = ms.current_mesh().face_matrix()
ax = np.load("axis.npy"); xb, pr = np.load("profile_adapt.npy")
x = V[:,0]; dy = V[:,1]-np.polyval(ax[0],x); dz = V[:,2]-np.polyval(ax[1],x)
r = np.hypot(dy,dz) - np.interp(x, xb, pr); th = np.degrees(np.arctan2(dy,dz))
s = np.load("scan_aligned.npy").astype(np.float64)
regions = {"asas": (x > 98) & (x < 132) & (np.abs(th) < 70) & (r > -0.5) & (dz > 0),
           "cascabel": (x > 203)}
for name, m in regions.items():
    sub = s[(s[:,0] > x[m].min()-1) & (s[:,0] < x[m].max()+1)]
    d, _ = cKDTree(sub).query(V[m], workers=-1)
    print(name, "verts", m.sum(), "dist pct", np.percentile(d,[50,90,95,99]).round(3), "frac>0.15", (d>0.15).mean().round(3))
    keepf = m[F].all(1); idx = np.where(m)[0]; remap = -np.ones(len(V), int); remap[idx] = np.arange(len(idx))
    f2 = remap[F[keepf]]
    t = np.clip((d - 0.05) / 0.25, 0, 1)
    col = np.c_[0.85 + 0.15*t, 0.68*(1-t) + 0.1*t, 0.32*(1-t) + 0.1*t, np.ones_like(t)]
    out = ml.MeshSet(); out.add_mesh(ml.Mesh(vertex_matrix=V[m], face_matrix=f2, v_color_matrix=col))
    out.save_current_mesh(f"{name}_diag.ply", save_vertex_color=True, binary=True)
    np.save(f"{name}_bbox.npy", np.r_[V[m].min(0), V[m].max(0)])
    print(name, "bbox", V[m].min(0).round(2), V[m].max(0).round(2))
