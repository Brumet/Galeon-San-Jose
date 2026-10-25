import numpy as np
from scipy.spatial import cKDTree
s = np.load("scan_aligned.npy").astype(np.float64); uv = np.load("rec_uv.npy").astype(np.float64)
# deshacer offset previo: volver a bbox del escaneo y aplicar offset correcto
s = s - s.min(0)  # x 0..220, y 0..43.5, z 0..36.8
off = np.array([uv[:,0].min(), -(s[:,1].max())/2, uv[:,2].max() - s[:,2].max()])
s = s + off
sub = uv[::31]
for th in (1.0, 0.3, 0.1, 0.05, 0.03):
    tr = cKDTree(s); d, i = tr.query(sub, workers=-1); m = d < th
    delta = np.median(sub[m] - s[i[m]], axis=0); s += delta
    print(th, delta.round(4), m.mean().round(3), flush=True)
tr = cKDTree(s); d, _ = tr.query(uv, workers=-1)
print("pct", np.percentile(d, [50, 60, 70, 80, 90, 95, 99]).round(4))
np.save("scan_aligned.npy", s.astype(np.float32)); np.save("rec_dist.npy", d.astype(np.float32))
