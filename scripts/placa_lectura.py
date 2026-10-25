"""Puntúa textos candidatos por renglón de la placa contra el relieve medido (correlación normalizada)."""
import numpy as np, sys
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter, zoom
hp = np.load("escudos_hp.npy"); DX = 0.02; X0 = 16.0; S0 = -8.82
FONT = "C:/Windows/Fonts/palab.ttf"
def banda(xc, h=0.55, s_a=-6.0, s_b=6.0):
    a, b = int((xc - h - X0) / DX), int((xc + h - X0) / DX); c, d = int((s_a - S0) / DX), int((s_b - S0) / DX)
    P = hp[a:b, c:d]; P = P - gaussian_filter(P, 8)
    return (P - P.mean()) / (P.std() + 1e-9)
def render(txt, cap_px):
    f = ImageFont.truetype(FONT, int(cap_px * 1.42))
    w = int(f.getlength(txt)) + 4; im = Image.new("L", (w, int(cap_px * 2)), 0)
    ImageDraw.Draw(im).text((2, int(cap_px * 0.5)), txt, fill=255, font=f, anchor="lt")
    a = np.asarray(im, float); rows = np.where(a.max(1) > 0)[0]; cols = np.where(a.max(0) > 0)[0]
    return a[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1] / 255.0
def score(M, txt, cap_mm=0.78):
    best = (-1, None)
    for cap in (cap_mm * 0.9, cap_mm, cap_mm * 1.1):
        T0 = render(txt, cap / DX)
        for wfrac in np.linspace(0.75, 1.0, 6):              # ancho del renglón respecto a la placa útil
            W = int(M.shape[1] * wfrac); T = zoom(T0, (1, W / T0.shape[1]), order=1)
            T = gaussian_filter(T, 1.0); T = (T - T.mean()) / (T.std() + 1e-9)
            for dx in range(0, M.shape[0] - T.shape[0] + 1, 2):
                for ds in range(0, M.shape[1] - T.shape[1] + 1, 6):
                    c = (M[dx:dx + T.shape[0], ds:ds + T.shape[1]] * T).mean()
                    if c > best[0]: best = (c, (cap, wfrac, dx, ds))
    return best
if __name__ == "__main__":
    cand = {72.6: ["DON DIEGO PHELIPE", "DON DIEGO FELIPE", "MARQVES DE LEGANES", "CAPITAN GENERAL"],
            73.8: ["DE GVZMAN MARQS", "PE DE GVZMAN MARQS", "DE GVZMAN MARQVES", "DE GVZMAN MARQ", "MARQVES DE LEGANES", "DE LA ARTILLERIA"],
            74.9: ["DE LEGANES CAP", "S DE LEGANES CAP", "DE LEGANES CAPN", "DE LEGANES CAPITAN", "DE LA ARTILLERIA", "GENERAL DE LA ART"],
            75.95: ["GENERAL DE LA ART", "GL DE LA ARTILLA", "DE LA ARTILLERIA", "GENERAL DE ESPAÑA", "DE ESPAÑA", "ARTILLERIA DE ESP"]}
    for xc, lst in cand.items():
        M = banda(xc)
        res = sorted(((score(M, t)[0], t) for t in lst), reverse=True)
        print(f"renglón x={xc}: " + " | ".join(f"{t} {c:.3f}" for c, t in res), flush=True)
