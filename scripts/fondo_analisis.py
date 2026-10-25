"""Paso 1: alinear escaneo original con la reconstrucción y marcar qué vértices son 'inventados' (sin datos)."""
import numpy as np, sys, os, time
from scipy.spatial import cKDTree

def log(*a): print(*a, flush=True)

def read_stl(p):
    with open(p, "rb") as f:
        f.seek(80); n = np.frombuffer(f.read(4), np.uint32)[0]
        d = np.frombuffer(f.read(n * 50), dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", 9), ("a", "<u2")]))
    return d["v"].reshape(-1, 3, 3)

root = sys.argv[1]; work = sys.argv[2]
t = time.time()
scan = read_stl(os.path.join(root, "Sin título.stl")).reshape(-1, 3).astype(np.float64)
scan = np.unique(scan, axis=0)
log("scan verts", len(scan), round(time.time() - t, 1))
M = np.array([[0, 1, 0], [0, 0, -1], [-1, 0, 0]], float)
scan = scan @ M.T
rec_tri = read_stl(os.path.join(root, "STL_preparado", "canon_resina_4M.stl")).astype(np.float64)
rec = rec_tri.reshape(-1, 3)
log("scan bbox", scan.min(0).round(3), scan.max(0).round(3))
log("rec  bbox", rec.min(0).round(3), rec.max(0).round(3))
# traslación inicial por bbox, luego refinamiento por vecinos cercanos (solo traslación)
off = rec.min(0) - scan.min(0)
tree_src = scan + off
for it in range(3):
    tr = cKDTree(tree_src)
    sub = rec[::37]
    d, i = tr.query(sub)
    m = d < 0.05
    delta = np.median(sub[m] - tree_src[i[m]], axis=0)
    tree_src = tree_src + delta
    log("iter", it, "delta", delta.round(4), "inliers", m.mean().round(3))
tr = cKDTree(tree_src)
uv, inv = np.unique(rec, axis=0, return_inverse=True)
d, _ = tr.query(uv, workers=-1)
log("dist pct", np.percentile(d, [50, 80, 90, 95, 99]).round(4))
np.save(os.path.join(work, "scan_aligned.npy"), tree_src.astype(np.float32))
np.save(os.path.join(work, "rec_uv.npy"), uv.astype(np.float32))
np.save(os.path.join(work, "rec_dist.npy"), d.astype(np.float32))
np.save(os.path.join(work, "rec_faces.npy"), inv.reshape(-1, 3).astype(np.int32))
log("done", round(time.time() - t, 1))
