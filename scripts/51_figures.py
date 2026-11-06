#!/usr/bin/env python3
"""Step 51: draw the data figures from the panel tables written by 50_figure_data.py (matplotlib only).

Every figure is drawn on a canvas in pixel units (760 px = 180 mm; 1 px = 0.237 mm = 0.671 pt;
11 px text = 7.4 pt).  No figure title line is drawn (titles belong in the legends);
panel letters and short panel headings are kept.  Supplementary Figs S1, S2 and S8 (one panel each) are 88 mm wide,
all others 180 mm.  Every plotted number is read from results/figures/data/*.tsv.

Not generated here (schematics): Fig. 1a (approach) and Fig. 7 (model).  For Fig. 1a, an image
placed at results/figures/schematics/Figure_1a.png is embedded if present, otherwise a labelled frame is left.
Fig. 6 and Supp. Fig. S10 are genome-browser figures made by scripts 32 and 43.

Output: results/figures/files/<Figure_1..5, Supplementary_Figure_S1..S9>.pdf (vector, TrueType fonts: pdf.fonttype 42)
and .png (300 dpi).
Usage: python scripts/51_figures.py [--only Figure_2 Supplementary_Figure_S5] [--data DIR] [--out DIR]
"""
import argparse
import math
import os
import sys
import traceback

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle  # noqa: E402

MM_PER_PX = 180.0 / 760.0
PT = MM_PER_PX * 72 / 25.4  # points per drawing pixel
K, RED, NAVY, TEAL, CYAN, LAV, SAL, GRY = "#000000", "#D55E00", "#0072B2", "#009E73", "#56B4E9", "#CC79A7", "#E69F00", "#8C8C8C"
LGR, PALE = "#D9D9D9", "#EDEDED"
DATA, OUT = "results/figures/data", "results/figures/files"
WARN = []


def pick_font():
    for name in ("Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"):
        try:
            p = fm.findfont(fm.FontProperties(family=name), fallback_to_default=False)
            if p:
                return name, p
        except Exception:
            continue
    return "DejaVu Sans", fm.findfont("DejaVu Sans")


FONT, FONT_PATH = pick_font()
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": [FONT, "DejaVu Sans"], "font.size": 12 * PT,
    "mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": f"{FONT}:italic", "mathtext.bf": f"{FONT}:bold",
    "mathtext.sf": FONT, "mathtext.cal": FONT, "mathtext.tt": FONT, "mathtext.default": "regular", "lines.scale_dashes": False, "axes.unicode_minus": True,
    "lines.solid_capstyle": "butt", "lines.dash_capstyle": "butt", "pdf.compression": 6})


# --------------------------------------------------------------------------------------------- canvas
class Canvas:
    """A figure whose data coordinates are the drawing's pixels (y downwards), with an offset stack for the
    drawings' transform='translate(dx,dy)' groups."""

    def __init__(self, w_px, h_px, dy=0.0):
        self.w, self.h = w_px, h_px
        self.fig = plt.figure(figsize=(w_px * MM_PER_PX / 25.4, h_px * MM_PER_PX / 25.4))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, w_px)
        self.ax.set_ylim(h_px, 0)
        self.ax.axis("off")
        self.off = [(0.0, dy)]
        self.z = 2
        self.drawn, self.skipped = 0, []  # panels drawn / panels left as a frame (missing table)

    def push(self, dx=0.0, dy=0.0):
        x, y = self.off[-1]
        self.off.append((x + dx, y + dy))

    def pop(self):
        self.off.pop()

    def xy(self, x, y):
        ox, oy = self.off[-1]
        return x + ox, y + oy

    def T(self, x, y, s, fs=12, ha="left", fw="normal", st="normal", c=K, rot=0, va="baseline", alpha=1.0):
        x, y = self.xy(x, y)
        ha = {"start": "left", "end": "right", "middle": "center"}.get(ha, ha)
        self.ax.text(x, y, s, fontsize=fs * PT, ha=ha, va=va, fontweight=fw, fontstyle=st, color=c, rotation=rot,
                     rotation_mode="anchor", alpha=alpha, zorder=10)

    def L(self, x1, y1, x2, y2, c=K, lw=1.0, dash=None, alpha=1.0, z=3):
        (x1, y1), (x2, y2) = self.xy(x1, y1), self.xy(x2, y2)
        ln, = self.ax.plot([x1, x2], [y1, y2], color=c, lw=lw * PT, alpha=alpha, zorder=z, solid_capstyle="butt")
        if dash:
            ln.set_dashes([d * PT for d in dash])

    def PL(self, pts, c=K, lw=1.0, dash=None, z=3, alpha=1.0):
        pts = [self.xy(x, y) for x, y in pts]
        if not pts:
            return
        ln, = self.ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, lw=lw * PT, zorder=z, alpha=alpha,
                           solid_joinstyle="miter")
        if dash:
            ln.set_dashes([d * PT for d in dash])

    def R(self, x, y, w, h, fc="none", ec="none", lw=1.0, rx=0.0, alpha=1.0, dash=None, z=2):
        x, y = self.xy(x, y)
        if rx:
            p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={rx}", fc=fc, ec=ec, lw=lw * PT,
                               alpha=alpha, zorder=z)
        else:
            p = Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw * PT, alpha=alpha, zorder=z)
        if dash:
            p.set_linestyle((0, [d * PT for d in dash]))
        self.ax.add_patch(p)

    def C(self, x, y, r, fc=K, ec="none", lw=0.0, alpha=1.0, ealpha=None, z=5):
        x, y = self.xy(x, y)
        fcol = matplotlib.colors.to_rgba(fc, alpha) if fc != "none" else "none"
        ecol = matplotlib.colors.to_rgba(ec, alpha if ealpha is None else ealpha) if ec != "none" else "none"
        self.ax.add_patch(Circle((x, y), r, fc=fcol, ec=ecol, lw=lw * PT, zorder=z))

    def P(self, pts, fc="none", ec=K, lw=1.0, z=5, closed=True):
        pts = [self.xy(x, y) for x, y in pts]
        self.ax.add_patch(Polygon(pts, closed=closed, fc=fc, ec=ec, lw=lw * PT, zorder=z, joinstyle="miter"))

    def tri(self, x, y, up, half, down, ec, lw=1.6, fc="#FFFFFF", z=6):
        """Open triangle as in the drawings: apex (x, y-up), base corners (x +- half, y+down)."""
        self.P([(x, y - up), (x + half, y + down), (x - half, y + down)], fc=fc, ec=ec, lw=lw, z=z)

    def image(self, path, x0, y0, x1, y1):
        img = plt.imread(path)
        (x0, y0), (x1, y1) = self.xy(x0, y0), self.xy(x1, y1)
        h_w = img.shape[0] / img.shape[1]
        w = x1 - x0
        h = min(y1 - y0, w * h_w)
        self.ax.imshow(img, extent=(x0, x0 + h / h_w, y0 + h, y0), zorder=1, interpolation="lanczos")
        self.ax.set_xlim(0, self.w)
        self.ax.set_ylim(self.h, 0)

    def save(self, name, outdir):
        os.makedirs(outdir, exist_ok=True)
        pdf, png = os.path.join(outdir, name + ".pdf"), os.path.join(outdir, name + ".png")
        self.fig.savefig(pdf)
        self.fig.savefig(png, dpi=300)
        plt.close(self.fig)
        return pdf, png


# ---------------------------------------------------------------------------------------------- data
class MissingTable(FileNotFoundError):
    pass


def D(name):
    p = os.path.join(ARGS.data, name)
    if not os.path.exists(p):
        raise MissingTable(p)
    return pd.read_csv(p, sep="\t", low_memory=False)


class panel:
    """`with panel(cv, "Figure 2c", box, letter):` around the drawing of one panel of a multi-panel figure.  The block
    reads its tables with D() before drawing anything; if one is missing (its panel failed or was skipped in
    50_figure_data.py), the panel is left as a labelled dashed frame, a SKIP line is printed and the figure goes on."""

    def __init__(self, cv, name, box, let=None):
        self.cv, self.name, self.box, self.let = cv, name, box, let

    def __enter__(self):
        self.depth = len(self.cv.off)
        return self

    def __exit__(self, et, e, tb):
        cv = self.cv
        if et is None:
            cv.drawn += 1
            return False
        if not issubclass(et, MissingTable):
            return False
        x, y, w, h = self.box
        if self.let:
            letter(cv, *self.let)
        cv.R(x, y, w, h, ec=LGR, lw=1, dash=(4, 3))
        cv.T(x + w / 2, y + h / 2, f"{self.name} not drawn: {os.path.basename(str(e))} missing", 10, ha="middle", c=GRY)
        del cv.off[self.depth:]
        cv.skipped.append(self.name)
        print(f"    SKIP {self.name}: {e} not found (that panel is missing or failed in 50_figure_data.py; "
              f"the rest of the figure is drawn)")
        return True


def row(df, **key):
    m = pd.Series(True, index=df.index)
    for k, v in key.items():
        m &= df[k].astype(str) == str(v)
    x = df[m]
    if len(x) != 1:
        raise KeyError(f"{len(x)} rows for {key}")
    return x.iloc[0]


def val(df, col, **key):
    return float(row(df, **key)[col])


def fint(v):
    return f"{int(round(float(v))):,}"


def minus(s):
    return s.replace("-", "−")


def fmt_p(p, style="std"):
    """p as written in the figures: 3 × 10⁻⁵ below 0.001; 0.002 (1 significant digit) below 0.01;
    else 2 decimals (style 'std') or 1 significant digit (style '1sig', Fig. 4a)."""
    if p is None or not np.isfinite(p):
        return ""
    if p < 0.001:
        m, e = f"{p:.0e}".split("e")
        return f"{int(m)} × 10$^{{{int(e)}}}$"
    if p < 0.01:
        return f"{p:.1g}" if f"{p:.1g}" != "0.01" else "0.010"
    if style == "1sig":
        return f"{p:.1g}"
    return "1.0" if p >= 0.995 else f"{p:.2f}"


def check_range(fig, what, v, lo, hi):
    """Values outside the axis range of the figure are clipped to it and reported."""
    if v is not None and np.isfinite(v) and (v < lo or v > hi):
        WARN.append((fig, what, float(v), lo, hi))
        return min(max(v, lo), hi)
    return v


def warn_lines(ws):
    out = {}
    for w in ws:
        if isinstance(w, str):
            out.setdefault((w,), []).append(None)
        else:
            out.setdefault((w[0], w[3], w[4]), []).append(w)
    lines = []
    for k, v in out.items():
        if len(k) == 1:
            lines.append(k[0])
        else:
            ex = ", ".join(f"{w[1]} {w[2]:.3g}" for w in v[:4]) + (", ..." if len(v) > 4 else "")
            lines.append(f"{k[0]}: {len(v)} value(s) outside the drawn axis range [{k[1]}, {k[2]}], clipped ({ex})")
    return lines


def wrap(s, n):
    out = []
    for w in s.split(" "):
        if out and len(out[-1] + " " + w) <= n:
            out[-1] += " " + w
        else:
            out.append(w)
    return out


def letter(cv, x, y, s):
    cv.T(x, y, s, 14, fw="bold")


def ax_left(cv, x0, y0, y1, ticks, Y, fmt=str, fs=11, tick=4):
    cv.L(x0, y0, x0, y1)
    for t in ticks:
        cv.L(x0 - tick, Y(t), x0, Y(t))
        cv.T(x0 - 7, Y(t) + 4, fmt(t), fs, ha="end")


def ax_bottom(cv, x0, x1, y, ticks, X, fmt=str, fs=11, dy=16):
    cv.L(x0, y, x1, y)
    for t in ticks:
        cv.L(X(t), y, X(t), y + 4)
        cv.T(X(t), y + dy, fmt(t), fs, ha="middle")


def gnum(v):
    """tick label like the drawings (0.5, 1, 1.5)"""
    return f"{v:g}"


# ======================================================================================== Figure 1
def figure1():
    cv = Canvas(760, 904, -32)
    # ---- a: schematic (not data) ----
    letter(cv, 24, 58, "a")
    cv.T(44, 58, "Approach: capped start sites inside introns, on the host gene’s strand, away from exons")
    img = "results/figures/schematics/Figure_1a.png"
    if os.path.exists(img):
        cv.image(img, 40, 64, 746, 301)
    else:
        cv.R(44, 70, 696, 229, fc="none", ec=LGR, lw=1, dash=(4, 3))
        cv.T(392, 188, "panel a: schematic, not generated from data (insert the drawing, or place it at "
                       "results/figures/schematics/Figure_1a.png)", 10, ha="middle", c=GRY)
    # ---- b (the drawing shifts panels b and c down to leave room below panel a) ----
    cv.push(0, 44)
    # ---- b: waffle ----
    with panel(cv, "Figure 1b (left)", (44, 300, 290, 250), (24, 290, "b")):
        W = D("fig1_b_waffle.tsv")
        col = {"itss": RED, "tss": NAVY, "inter": TEAL, "anti": LAV}
        total = float(W.total.iloc[0])
        letter(cv, 24, 290, "b")
        cv.T(44, 290, f"{fint(total)} K562 PRO-cap start sites (1 square = 0.5%)")
        order = ["tss", "anti", "inter", "other", "itss"]
        labels = {"tss": "at annotated transcript starts", "anti": "antisense to genes", "inter": "between genes",
                  "other": "elsewhere in genes, on their strand", "itss": "in introns, host strand, ≥1 kb from exons"}
        i = 0
        wx0, wy0, cs, gap = 44, 316, 11, 2
        for c in order:
            for _ in range(int(val(W, "squares", category=c))):
                cv.R(wx0 + (i % 20) * (cs + gap), wy0 + (i // 20) * (cs + gap), cs, cs, fc=col.get(c, LGR), rx=1.5)
                i += 1
        if i != 200:
            WARN.append(f"Figure 1b: {i} waffle squares instead of 200")
        for j, c in enumerate(order):
            cv.R(wx0, 456 + j * 17, 10, 10, fc=col.get(c, LGR), rx=1.5)
            cv.T(wx0 + 16, 465 + j * 17, f"{labels[c]} ({val(W, 'pct', category=c):.1f}%)", 11, fw="bold" if c == "itss" else "normal")
    # ---- b: flow ----
    with panel(cv, "Figure 1b (right)", (348, 300, 400, 290)):
        F = D("fig1_b_flow.tsv")
        tx, cw, ch = 348, 254, 46
        cy = lambda i: 384 + i * 72
        c1, c2 = tx + cw - 86, tx + cw - 30
        n = lambda line, step: val(F, "n", line=line, step=step)
        cv.T(tx, 290, "Intronic start sites through the analysis")
        cv.T(c1, 302, "K562", 10, ha="middle", fw="bold")
        cv.T(c2, 302, "HCT116", 10, ha="middle", fw="bold")
        cv.R(tx, 308, cw, 30, fc="#F4F4F4", ec=LGR, rx=6)
        cv.T(tx + 12, 327, "All PRO-cap start sites", 11)
        allk = n("K562", "all_starts")
        cv.T(c1, 328, fint(allk) if np.isfinite(allk) else "–", 13, ha="middle", fw="bold")
        cv.T(c2, 328, "–", 11, ha="middle")
        cv.T(tx + cw - 8, 350, "genome-wide breakdown for K562 only", 9.5, ha="end")
        cv.L(tx + 16, 338, tx + 16, cy(0) - 3, c=RED, lw=1.5)
        main = [("intronic", "Intronic start sites", "in introns, on the host strand, ≥1 kb from exons"),
                ("capped", "Passed the cap test", "capped-start test passed"),
                ("capped_active", "Confirmed start, active host gene", "host gene transcribed")]
        side = [("near_annot", "Near an annotated intronic lncRNA", "of all intronic sites, ≤500 nt from an annotated intronic lncRNA", 0),
                ("capped_silent", "Confirmed start, silent host gene", "host gene silent", 1)]
        fillc = lambda i: "#F7DCCB" if i == 0 else "#FBEDE4" if i == 1 else "#FEF8F4"
        for i, (step, lab, crit) in enumerate(main):
            y = cy(i)
            cv.T(tx + 28, y - 10, crit, 10, st="italic")
            cv.PL([(tx + 12, y - 7), (tx + 16, y - 1), (tx + 20, y - 7)], c=RED, lw=1.5)
            cv.R(tx, y, cw, ch, fc=fillc(i), ec=LGR, rx=6)
            cv.R(tx, y, 5, ch, fc=RED, rx=2, alpha=1 - i * 0.3, z=3)
            lns = wrap(lab, 22)
            for j, l in enumerate(lns):
                cv.T(tx + 14, y + 27 - (len(lns) - 1) * 6.5 + j * 13, l, 11, fw="bold" if i == 0 else "normal")
            cv.T(c1, y + 29, fint(n("K562", step)), 15, ha="middle", fw="bold")
            cv.T(c2, y + 29, fint(n("HCT116", step)), 15, ha="middle", fw="bold")
            if i < len(main) - 1:
                cv.L(tx + 16, y + ch, tx + 16, cy(i + 1) - 3, c=RED, lw=1.5)
        sx, sw, sh = tx + cw + 20, 122, 74
        prev = 0
        for step, lab, crit, frm in side:
            ym = cy(frm) + ch / 2
            cl = wrap(crit, 22)
            yt = max(ym - sh + 12, prev + 12 + len(cl) * 11)
            prev = yt + sh
            for j, l in enumerate(cl):
                cv.T(sx + 8, yt - 8 - (len(cl) - 1 - j) * 11, l, 10, st="italic")
            cv.L(tx + cw, ym, sx - 3, ym, c=RED, lw=1.5)
            cv.PL([(sx - 9, ym - 4), (sx - 3, ym), (sx - 9, ym + 4)], c=RED, lw=1.5)
            cv.R(sx, yt, sw, sh, fc=fillc(2), ec=LGR, rx=6)
            cv.R(sx, yt, 5, sh, fc=RED, rx=2, alpha=0.4, z=3)
            for j, l in enumerate(wrap(lab, 18)):
                cv.T(sx + 12, yt + 16 + j * 13, l, 11)
            cv.T(sx + 38, yt + 48, "K562", 9.5, ha="middle", fw="bold")
            cv.T(sx + 92, yt + 48, "HCT116", 9.5, ha="middle", fw="bold")
            cv.T(sx + 38, yt + 66, fint(n("K562", step)), 15, ha="middle", fw="bold")
            cv.T(sx + 92, yt + 66, fint(n("HCT116", step)), 15, ha="middle", fw="bold")
    cv.pop()
    # ---- c ----
    cv.push(0, 56)
    with panel(cv, "Figure 1c", (44, 616, 700, 260), (24, 606, "c")):
        B, S = D("fig1_c.tsv"), D("fig1_c_summary.tsv")
        letter(cv, 24, 606, "c")
        cv.T(44, 606, f"Validation: start sites detected for {fint(val(S, 'n_nested_pc', line='K562'))} protein-coding genes nested in introns")
        pw, py0, py1 = 250, 662, 800
        for line, x0 in (("K562", 82), ("HCT116", 450)):
            PX = lambda k: x0 + k / 5 * pw
            PY = lambda v: py1 - v / 100 * (py1 - py0)
            cv.T(x0, 638, line, fw="bold")
            cv.T(x0 + 52, 638, f"{fint(val(S, 'n_clear_detected', line=line))} of {fint(val(S, 'n_clear', line=line))} "
                               f"clearly expressed detected", 11)
            ax_left(cv, x0, py0, py1, [0, 50, 100], PY)
            ax_bottom(cv, x0, x0 + pw, py1, range(6), PX)
            cv.T(x0 - 34, (py0 + py1) / 2, "% detected", 11, ha="middle", rot=90)
            cv.T(x0 + pw / 2, py1 + 32, "expression bin (0 = none, 5 = highest)", 11, ha="middle")
            for kind in ("alone", "nested"):
                pts = B[(B.line == line) & (B.kind == kind)].sort_values("bin")
                xy = [(PX(r.bin), PY(r.pct)) for r in pts.itertuples() if np.isfinite(r.pct)]
                if kind == "nested":
                    cv.PL(xy, c=CYAN, lw=2.2)
                    for x, y in xy:
                        cv.C(x, y, 4, fc=CYAN, ec=K, lw=0.8)
                else:
                    cv.PL(xy, c=GRY, lw=1.4, dash=(5, 3))
                    for x, y in xy:
                        cv.R(x - 3.5, y - 3.5, 7, 7, fc="#FFFFFF", ec=GRY, lw=1.2, z=5)
        cv.L(82, 862, 106, 862, c=CYAN, lw=2.2)
        cv.C(94, 862, 4, fc=CYAN, ec=K, lw=0.8)
        cv.T(112, 866, "nested protein-coding genes (either strand)", 11)
        cv.L(420, 862, 444, 862, c=GRY, lw=1.4, dash=(5, 3))
        cv.R(428.5, 858.5, 7, 7, fc="#FFFFFF", ec=GRY, lw=1.2, z=5)
        cv.T(450, 866, "stand-alone genes, same expression", 11)
    cv.pop()
    return cv


# ======================================================================================== Figure 2
def figure2():
    cv = Canvas(760, 1010, -32)
    # ---- a ----
    with panel(cv, "Figure 2a", (44, 66, 340, 266), (24, 58, "a")):
        A = D("fig2_a.tsv")
        letter(cv, 24, 58, "a")
        cv.T(44, 58, "Histone marks within 1 kb (K562, HCT116)")
        groups = [("host", "host gene promoters", K), ("int", "internal promoters with a capped start", NAVY),
                  ("itss", "intronic start sites, confirmed, active host", RED), ("rand", "random intronic positions", GRY)]
        marks = ["H3K27ac", "H3K4me3", "H3K4me1"]
        aw, ay1 = 110, 230
        ayb = lambda i: 92 + i * 48
        for line, ax in (("K562", 112), ("HCT116", 262)):
            AX = lambda v: ax + v / 100 * aw
            cv.T(ax + aw / 2, 80, line, 11, ha="middle", fw="bold")
            ax_bottom(cv, ax, ax + aw, ay1, [0, 50, 100], AX)
            for i, mk in enumerate(marks):
                for j, (g, _, cc) in enumerate(groups):
                    y, v = ayb(i) + j * 11, val(A, "pct", line=line, mark=mk, group=g)
                    cv.L(AX(0), y, AX(v), y, c=LGR, lw=1, z=2)
                    if line == "K562":
                        cv.C(AX(v), y, 3.6, fc=cc, ec=K, lw=0.5)
                    else:
                        cv.tri(AX(v), y, 4.2, 4, 3, cc, lw=1.4)
                    cv.T(AX(v) + 7, y + 3.3, f"{v:.0f}", 9.5)
        for i, mk in enumerate(marks):
            cv.T(100, ayb(i) + 20, mk, 12, ha="end")
        cv.T((112 + 262 + aw) / 2, ay1 + 31, "% of sites with a peak within 1 kb", 11, ha="middle")
        for j, (g, lab, cc) in enumerate(groups):
            y = ay1 + 48 + j * 12
            cv.C(64, y - 4, 3.6, fc=cc, ec=K, lw=0.5)
            cv.P([(78, y - 8.2), (82, y - 1), (74, y - 1)], fc="#FFFFFF", ec=cc, lw=1.4)
            cv.T(90, y, f"{lab} (n = {fint(val(A, 'n', line='K562', mark='H3K27ac', group=g))} K562, "
                        f"{fint(val(A, 'n', line='HCT116', mark='H3K27ac', group=g))} HCT116)", 11)
        cv.C(64, 322, 3.6, fc="#7F7F7F", ec=K, lw=0.5)
        cv.P([(78, 317.8), (82, 325), (74, 325)], fc="#FFFFFF", ec="#7F7F7F", lw=1.4)
        cv.T(90, 326, "circles: K562; triangles: HCT116", 11)
    # ---- b ----
    with panel(cv, "Figure 2b", (420, 66, 320, 200), (400, 58, "b")):
        TE = D("fig2_b.tsv")
        letter(cv, 400, 58, "b")
        cv.T(420, 58, "Transposable elements at the start site")
        bx0, bw = 560, 170
        BX = lambda e: bx0 + (e + 3) / 5 * bw
        ax_bottom(cv, bx0, bx0 + bw, 206, [-3, -2, -1, 0, 1, 2], BX)
        cv.L(BX(0), 74, BX(0), 206, lw=0.8, dash=(2, 3))
        cv.T(bx0 + bw / 2, 238, "log2 (intronic sites / random positions)", 11, ha="middle")
        for i, (c, lab) in enumerate((("none", "no transposable element"), ("alu", "Alu"), ("l1", "L1"),
                                      ("dna", "DNA transposons"), ("ltr", "LTR / ERV"))):
            y = 86 + i * 26
            cv.T(bx0 - 10, y + 4, lab, 11, ha="end")
            for line in ("K562", "HCT116"):
                e = check_range("Figure 2b", f"{line} {c}", val(TE, "log2_ratio", line=line, cls=c), -3, 2)
                cc = RED if e > 0 else NAVY
                if line == "K562":
                    cv.C(BX(e), y - 4, 5, fc=cc, ec=K, lw=0.7)
                else:
                    cv.tri(BX(e), y, 1.5, 5, 8, cc, lw=1.6)
        cv.C(566, 254, 4.5, fc="#7F7F7F", ec=K, lw=0.7)
        cv.T(575, 258, "K562", 11)
        cv.P([(626, 249), (630.5, 257.5), (621.5, 257.5)], fc="#FFFFFF", ec="#7F7F7F", lw=1.5)
        cv.T(636, 258, "HCT116", 11)
    # ---- c ----
    with panel(cv, "Figure 2c", (44, 290, 700, 345), (24, 282, "c")):
        cv.push(0, 68)
        M = D("fig2_c.tsv")
        letter(cv, 24, 282, "c")
        cv.T(44, 282, "Known transcription factor motifs within 150 nt of intronic start sites (HOMER; top 25 per cell line)")
        fc = {"GATA": RED, "AP-1": NAVY, "ETS": TEAL, "MARE": SAL, "KLF": LAV}
        vw, vy0, vy1, RK = 260, 316, 520, 1.35
        for line, x0, vmax, ticks in (("K562", 80, 90, [0, 30, 60, 90]), ("HCT116", 440, 800, [0, 200, 400, 600, 800])):
            rs = M[M.line == line].copy()
            top = rs.top_in_family.astype(str) == "True"
            VX = lambda e: x0 + e / 4.5 * vw
            VY = lambda lp: vy1 - lp / vmax * (vy1 - vy0)
            cv.T(x0 + vw / 2, vy0 - 14, line, 11, ha="middle", fw="bold")
            ax_left(cv, x0, vy0 + 6, vy1, ticks, VY)
            cv.L(x0, vy0 + 6, x0, vy0)  # the drawing's y axis starts 6 px below the top tick
            ax_bottom(cv, x0, x0 + vw, vy1, [0, 1, 2, 3, 4], VX)
            cv.T(x0 + vw / 2, vy1 + 32, "log2 enrichment over random intronic positions", 11, ha="middle")
            cv.T(x0 - 40, (vy0 + vy1) / 2, "−log10 p", 11, ha="middle", rot=90)
            rs["tflag"] = top.astype(int)
            for r in rs.sort_values("tflag", kind="stable").itertuples():
                pale = r.family == "AP-1" and not r.tflag
                e = check_range("Figure 2c", f"{line} {r.name} enrichment", r.log2_enrichment, 0, 4.5)
                lp = check_range("Figure 2c", f"{line} {r.name} -log10 p", r.log10p, 0, vmax)
                cv.C(VX(e), VY(lp), RK * math.sqrt(r.pct_target), fc=fc.get(r.family, "#BDBDBD"), alpha=0.3 if pale else 0.75,
                     ec=K, ealpha=0.45 if pale else 1.0, lw=0.6, z=5 + r.tflag)
            offs = ({"GATA": (-14, -10, "end"), "AP-1": (22, -6, "start"), "ETS": (-16, -8, "end")} if line == "K562"
                    else {"AP-1": (-30, -8, "end"), "MARE": (12, 3, "start"), "ETS": (12, -26, "start")})
            dflt = (12, 4, "start") if line == "K562" else (-6, -30, "middle")
            for r in rs[top].itertuples():
                cx, cy_, rr = VX(r.log2_enrichment), VY(r.log10p), RK * math.sqrt(r.pct_target)
                dx, dy, an = offs.get(r.family, dflt)
                dl = math.hypot(dx, dy)
                cv.L(cx + dx / dl * rr, cy_ + dy / dl * rr, cx + dx, cy_ + dy, lw=0.6, z=8)
                tx_ = cx + dx if an == "middle" else (cx + dx - 2 if an == "end" else cx + dx + 2)
                ty_ = cy_ + dy - 3 if an == "middle" else cy_ + dy + 3.5
                cv.T(tx_, ty_, r.name, 10, ha=an, fw="bold")
        cv.T(80, 586, "motif family:", 11)
        for x, c, lab in ((172, RED, "GATA"), (230, NAVY, "AP-1"), (286, TEAL, "ETS"), (336, SAL, "NRF2 / NF-E2 / BACH (MARE)"),
                          (520, LAV, "KLF"), (572, "#BDBDBD", "other")):
            cv.C(x, 582, 5, fc=c, alpha=0.75, ec=K, ealpha=1, lw=0.6)
            cv.T(x + 9, 586, lab, 11)
        cv.T(80, 612, "circle area ∝ share of sites with the motif:", 11)
        for i, v in enumerate((10, 25, 50)):
            r = RK * math.sqrt(v)
            cv.C(340 + i * 52, 607, r, fc="#BDBDBD", alpha=0.75, ec=K, ealpha=1, lw=0.6)
            cv.T(340 + i * 52 + r + 4, 611, f"{v}%", 11)
        cv.T(80, 632, "pale circles: other AP-1 motifs (variants of the labelled ones)", 11)
        cv.pop()
    # ---- d ----
    with panel(cv, "Figure 2d", (44, 636, 700, 295), (24, 628, "d")):
        cv.push(0, 104)
        V, O = D("fig2_d_venn.tsv"), D("fig2_d_or.tsv")
        letter(cv, 24, 628, "d")
        cv.T(44, 628, "Start sites used in both cell lines versus one (overlap: all sites passing the cap test; "
                      "odds ratios: sites in active hosts)")
        nk, nh = val(V, "n_confirmed", line="K562"), val(V, "n_confirmed", line="HCT116")
        sk, shh = val(V, "n_shared", line="K562"), val(V, "n_shared", line="HCT116")
        kr = 41.4 / math.sqrt(1346)  # drawing scale: radius = kr * sqrt(number of sites)
        r1, r2 = kr * math.sqrt(nk), kr * math.sqrt(nh)
        target = math.pi * kr ** 2 * (sk + shh) / 2

        def lens(d):
            if d >= r1 + r2:
                return 0.0
            if d <= abs(r1 - r2):
                return math.pi * min(r1, r2) ** 2
            a = r1 ** 2 * math.acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
            b = r2 ** 2 * math.acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2))
            c = 0.5 * math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
            return a + b - c

        lo_, hi_ = abs(r1 - r2), r1 + r2
        for _ in range(60):
            mid = (lo_ + hi_) / 2
            lo_, hi_ = (mid, hi_) if lens(mid) > target else (lo_, mid)
        dd = (lo_ + hi_) / 2
        ex_, ey = 118, 760
        cv.C(ex_, ey, r1, fc=CYAN, alpha=0.35, ec=K, ealpha=1, lw=0.8, z=3)
        cv.C(ex_ + dd, ey, r2, fc=SAL, alpha=0.35, ec=K, ealpha=1, lw=0.8, z=3)
        cv.T(ex_ - r1 + 4, ey - r1 - 8, f"K562 ({fint(nk)})", 11, fw="bold")
        cv.T(ex_ + dd - 20, ey - r2 - 8, f"HCT116 ({fint(nh)})", 11, fw="bold")
        cv.T(ex_ - 18, ey + 4, fint(val(V, "n_specific", line="K562")), 11, ha="middle")
        cv.L(ex_ + 28, ey + 12, ex_ + 2, ey + 52, lw=0.6)
        cv.T(ex_ - 4, ey + 64, f"$\\bf{{{fint(sk)}\\ /\\ {fint(shh)}}}$ shared".replace(",", "{,}"), 11, ha="middle")
        cv.T(ex_ + dd + 26, ey + 4, fint(val(V, "n_specific", line="HCT116")), 11, ha="middle")
        cv.T(ex_ + dd / 2, ey + r2 + 22, "area ∝ number of sites", 11, ha="middle")
        fx0, fw = 425, 250
        FX = lambda v: fx0 + (math.log2(v) + 3) / 5 * fw
        cv.T(fx0 + fw / 2, 656, "odds ratio, shared vs specific, at equal strength (CMH)", 11, ha="middle")
        ax_bottom(cv, fx0, fx0 + fw, 846, [0.125, 0.25, 0.5, 1, 2, 4], FX, fmt=gnum)
        cv.T(FX(1) - 5, 878, "← more common in cell-type-specific sites", 10, ha="end")
        cv.T(FX(1) + 5, 878, "more common in shared sites →", 10)
        cv.T(fx0 + fw / 2, 900, "odds ratio (log scale): shared vs cell-type-specific sites,", 10.5, ha="middle")
        cv.T(fx0 + fw / 2, 914, "Cochran–Mantel–Haenszel across thirds of PRO-cap reads", 10.5, ha="middle")
        cv.T(fx0 + fw / 2, 928, "whiskers: 95% CI (Robins–Breslow–Greenland)", 10.5, ha="middle")
        cv.L(FX(1), 668, FX(1), 846, lw=0.8, dash=(2, 3))
        rows_ = [("conserved", "K562", "conserved (phyloP > 1)"), ("conserved", "HCT116", ""),
                 ("transposon", "K562", "in a transposon"), ("transposon", "HCT116", "")]
        for i, (f, line, lab) in enumerate(rows_):
            y = 700 + i * 34 + (16 if i >= 2 else 0)
            r = row(O, line=line, feature=f)
            o_, lo, hi = (check_range("Figure 2d", f"{line} {f}", v, 0.125, 4) for v in (r["or"], r.ci_low, r.ci_high))
            if lab:
                cv.T(fx0 - 10, y + 10, lab, 11, ha="end")
            cv.L(FX(lo), y, FX(hi), y, lw=1.2)
            cv.L(FX(lo), y - 3, FX(lo), y + 3, lw=1.2)
            cv.L(FX(hi), y - 3, FX(hi), y + 3, lw=1.2)
            if line == "K562":
                cv.C(FX(o_), y, 5, fc=NAVY, ec=K, lw=0.7)
            else:
                cv.tri(FX(o_), y, 5.5, 5, 4, NAVY, lw=1.6)
            cv.T(fx0 + fw + 6, y + 4, f"p = {fmt_p(r.p)}", 10)
        cv.pop()
    return cv


# ======================================================================================== Figure 3
def figure3():
    cv = Canvas(760, 668, -32)
    gc = {"itss": RED, "int": NAVY, "rand": GRY}
    # ---- a ----
    with panel(cv, "Figure 3a", (44, 66, 330, 280), (24, 58, "a")):
        A = D("fig3_a.tsv")
        ax0, ax1, ay0, ay1 = 80, 330, 86, 270
        AX = lambda t: ax0 + 30 + t * ((ax1 - ax0 - 60) / 2)
        AY = lambda v: ay1 - v / 80 * (ay1 - ay0)
        letter(cv, 24, 58, "a")
        cv.T(44, 58, "Stable RNA at equal initiation strength")
        cv.T(44, 74, "sites at least 3 kb upstream of the next host exon", 11)
        rp = val(A, "pct", line="K562", group="rand", third="all")
        cv.L(ax0, AY(rp), ax1, AY(rp), c=GRY, lw=1, dash=(3, 3))
        cv.T(ax1 + 4, AY(rp) + 4, f"random {rp:.0f}%", 10, c=GRY)
        ax_left(cv, ax0, ay0, ay1, [0, 20, 40, 60, 80], AY)
        cv.L(ax0, ay1, ax1, ay1)
        cv.T(ax0 - 34, (ay0 + ay1) / 2, "% sites with stable RNA", 11, ha="middle", rot=90)
        for t, th in enumerate(("weak", "medium", "strong")):
            cv.L(AX(t), ay1, AX(t), ay1 + 4)
            cv.T(AX(t), ay1 + 16, th, 11, ha="middle")
        cv.T((ax0 + ax1) / 2, ay1 + 32, "initiation strength (third of PRO-cap reads)", 11, ha="middle")
        for g in ("itss", "int"):
            for line in ("K562", "HCT116"):
                pts = [(AX(t), AY(check_range("Figure 3a", f"{line} {g} {th}", val(A, "pct", line=line, group=g, third=th), 0, 80)))
                       for t, th in enumerate(("weak", "medium", "strong"))]
                cv.PL(pts, c=gc[g], lw=1.8, dash=None if line == "K562" else (5, 3))
                for x, y in pts:
                    if line == "K562":
                        cv.C(x, y, 4.2, fc=gc[g], ec=K, lw=0.7)
                    else:
                        cv.tri(x, y, 5, 4.5, 3.5, gc[g], lw=1.6)
        cv.L(80, 316, 100, 316, c=RED, lw=1.8)
        cv.T(106, 320, "intronic start sites", 11)
        cv.L(210, 316, 230, 316, c=NAVY, lw=1.8)
        cv.T(236, 320, "internal promoters", 11)
        cv.C(90, 334, 4, fc="#7F7F7F", ec=K, lw=0.7)
        cv.T(106, 338, "K562 (solid)", 11)
        cv.P([(220, 329), (224.5, 337.5), (215.5, 337.5)], fc="#FFFFFF", ec="#7F7F7F", lw=1.6)
        cv.T(236, 338, "HCT116 (dashed)", 11)
    # ---- b ----
    with panel(cv, "Figure 3b", (414, 66, 330, 230), (394, 58, "b")):
        Bc = D("fig3_b.tsv")
        PROD, LG = "#CC79A7", "#E3E3E3"
        letter(cv, 394, 58, "b")
        cv.T(414, 58, "Single sites after DRB release (HCT116, both replicates)")
        cv.T(414, 74, "each square = 1% of sites", 11)
        cs, cg, wy0 = 8.6, 1.6, 112

        def show(v):
            return f"{v:.1f}" if v < 3 else f"{int(math.floor(v + 0.5))}"

        for g, x0, name in (("rand", 404, "Random positions"), ("int", 520, "Internal promoters"), ("itss", 636, "Intronic start sites")):
            r = row(Bc, group=g)
            prod, short = int(math.floor(r.pct_productive + 0.5)), int(math.floor(r.pct_short + 0.5))
            cv.T(x0 + 50, wy0 - 10, name, 11, ha="middle", fw="bold" if g == "itss" else "normal")
            for i in range(100):
                c = PROD if i < prod else SAL if i < prod + short else LG
                cv.R(x0 + (i % 10) * (cs + cg), wy0 + (i // 10) * (cs + cg), cs, cs, fc=c)
            cv.T(x0 + 50, wy0 + 114, f"productive {show(r.pct_productive)}%", 10.5, ha="middle")
            cv.T(x0 + 50, wy0 + 127, f"short burst {show(r.pct_short)}%", 10.5, ha="middle")
        cv.R(418, 262, 10, 10, fc=PROD)
        cv.T(432, 271, "productive (signal 2–8 kb)", 11)
        cv.R(592, 262, 10, 10, fc=SAL)
        cv.T(606, 271, "short burst (0–2 kb)", 11)
        cv.R(418, 280, 10, 10, fc=LG)
        cv.T(432, 289, "no detectable response", 11)
    # ---- c ----
    with panel(cv, "Figure 3c", (44, 380, 700, 300), (24, 372, "c")):
        P = D("fig3_c.tsv")
        letter(cv, 24, 372, "c")
        cv.T(44, 372, "Nascent RNA after DRB release relative to 40 min, around sites beyond the host’s wave front (HCT116)")
        pw, py0, py1 = 280, 400, 600
        PY = lambda v: py1 - v / 1.6 * (py1 - py0)
        for m, x0 in ((10, 80), (20, 440)):
            PXf = lambda kb: x0 + (kb + 10) / 20 * pw
            cv.T(x0 + pw / 2, py0 - 6, f"{m} min after release", 11, ha="middle", fw="bold")
            ax_left(cv, x0, py0, py1, [0, 0.4, 0.8, 1.2, 1.6], PY, fmt=lambda t: f"{t:.1f}")
            ax_bottom(cv, x0, x0 + pw, py1, [-10, -5, 0, 5, 10], PXf, fmt=lambda t: f"+{t}" if t > 0 else f"{t}")
            cv.L(PXf(0), py0, PXf(0), py1, lw=0.8, dash=(2, 3))
            cv.T(x0 + pw / 2, py1 + 32, "kb from the site (direction of transcription)", 11, ha="middle")
            if m == 10:
                cv.T(x0 - 34, (py0 + py1) / 2, "nascent RNA, 10 or 20 min / 40 min", 11, ha="middle", rot=90)
            for s in ("rand", "int", "itss"):
                q = P[(P.minutes == m) & (P.group == s)].sort_values("kb_mid")
                pts = [(PXf(r.kb_mid), PY(check_range("Figure 3c", f"{m} min {s}", r.ratio, 0, 1.6))) for r in q.itertuples()]
                cv.PL(pts, c=gc[s], lw=1.3 if s == "rand" else 2, dash=(4, 3) if s == "rand" else None)
                for x, y in pts:
                    cv.C(x, y, 1.6 if s == "rand" else 2.4, fc=gc[s], z=6)
        nn = lambda m, g: fint(P[(P.minutes == m) & (P.group == g)].n_sites.iloc[0])
        cv.L(80, 672, 100, 672, c=RED, lw=2)
        cv.T(106, 676, f"intronic start sites (n = {nn(10, 'itss')} / {nn(20, 'itss')})", 11)
        cv.L(330, 672, 350, 672, c=NAVY, lw=2)
        cv.T(356, 676, f"internal promoters (n = {nn(10, 'int')} / {nn(20, 'int')})", 11)
        cv.L(580, 672, 600, 672, c=GRY, lw=1.3, dash=(4, 3))
        cv.T(606, 676, "random positions", 11)
    return cv


# ======================================================================================== Figure 4
def figure4():
    cv = Canvas(760, 658, -32)
    gcol = {"rand": GRY, "host": K, "int": NAVY, "inter": TEAL, "silent": SAL}
    # ---- a ----
    with panel(cv, "Figure 4a", (44, 66, 410, 285), (24, 58, "a")):
        A = D("fig4_a.tsv")
        letter(cv, 24, 58, "a")
        cv.T(44, 58, "RNA after the site rises without MTREX (K562)")
        ax0, aw = 230, 160
        AX = lambda v: ax0 + v / 50 * aw
        ax_bottom(cv, ax0, ax0 + aw, 282, [0, 25, 50], AX, dy=16)
        cv.T(ax0 + aw / 2, 314, "% sites above control threshold", 11, ha="middle")
        cv.L(AX(5), 76, AX(5), 282, lw=0.8, dash=(2, 3))
        cv.T(AX(5) + 3, 74, "random positions, 5% by construction", 10)
        names = [("rand", "Random intronic positions"), ("int", "Internal promoters"), ("itss", "Intronic start sites, active host"),
                 ("prom", "… promoter-like"), ("enh", "… enhancer-like"), ("silent", "Intronic start sites, silent host")]
        RK = 14 / math.sqrt(A.n.max())
        for i, (g, lab) in enumerate(names):
            r = row(A, group=g)
            y = 92 + i * 34
            rad = max(3, RK * math.sqrt(r.n))
            v = check_range("Figure 4a", g, r.pct, 0, 50)
            cv.T(ax0 - 10, y + 4, f"{lab} ({fint(r.n)})", 11, ha="end")
            cv.C(AX(v), y, rad, fc=gcol.get(g, RED), alpha=0.55 if g in ("prom", "enh") else 1, ec=K, ealpha=1, lw=0.7)
            cv.T(AX(v) + rad + 4, y + 4, f"{r.pct:.0f}%", 11)
            if g != "rand":
                cv.T(ax0 + aw + 14, y + 4, f"p = {fmt_p(r.p, '1sig')}", 10)
        cv.T(44, 334, "Knockout check: upstream antisense RNAs (PROMPTs) rose only without MTREX, and only partly.", 10)
        cv.T(44, 347, "Dot area ∝ number of sites; p: one-sided Mann–Whitney test.", 10)
    # ---- b ----
    with panel(cv, "Figure 4b", (490, 66, 255, 255), (470, 58, "b")):
        Bb = D("fig4_b.tsv")
        letter(cv, 470, 58, "b")
        cv.T(490, 58, "Splice-donor site before the first PAS")
        cv.T(490, 74, "share of sites, first kb of the transcribed strand", 11)
        dx0, dx1 = 490, 740
        DX = lambda v: dx0 + (v - 30) / 50 * (dx1 - dx0)
        dc = {"host": K, "int": NAVY, "intall": LAV, "rand": GRY, "itss": RED}
        for i, (g, lab) in enumerate((("host", "Host gene promoters"), ("int", "Internal promoters in use"),
                                      ("intall", "All annotated internal promoters"), ("itss", "Intronic start sites"),
                                      ("rand", "Random intronic positions"))):
            y = 98 + i * 38
            k_ = check_range("Figure 4b", f"K562 {g}", val(Bb, "pct", line="K562", group=g), 30, 80)
            h_ = check_range("Figure 4b", f"HCT116 {g}", val(Bb, "pct", line="HCT116", group=g), 30, 80)
            cv.T(dx0, y, lab, 11, fw="bold" if g == "itss" else "normal")
            cv.T(dx1, y, f"{int(math.floor(k_ + 0.5))}% / {int(math.floor(h_ + 0.5))}%", 10, ha="end")
            cv.L(DX(k_), y + 9.5, DX(h_), y + 19, c="#BDBDBD", lw=2)
            cv.C(DX(k_), y + 9.5, 5, fc=dc[g], ec=K, lw=0.7)
            cv.tri(DX(h_), y + 19, 5.5, 5, 4, dc[g], lw=1.6)
        ax_bottom(cv, dx0, dx1, 284, [30, 40, 50, 60, 70, 80], DX)
        cv.T(dx0, 316, "% of sites", 11)
        cv.C(620, 312, 4.5, fc="#7F7F7F", ec=K, lw=0.7)
        cv.T(629, 316, "K562", 11)
        cv.P([(677, 307), (681.5, 316), (672.5, 316)], fc="#FFFFFF", ec="#7F7F7F", lw=1.5)
        cv.T(688, 316, "HCT116", 11)
    # ---- c ----
    with panel(cv, "Figure 4c", (44, 386, 700, 290), (24, 378, "c")):
        Cc = D("fig4_c.tsv")
        letter(cv, 24, 378, "c")
        cv.T(44, 378, "Polyadenylation signals (PAS) in the first kilobase: transcribed strand (filled) versus opposite strand (open)")
        bw = 190
        for line, x0 in (("K562", 180), ("HCT116", 500)):
            BX = lambda v: x0 + v / 2 * bw
            cv.T(x0 + bw / 2, 402, line, 11, ha="middle", fw="bold")
            cv.T(x0 + bw + 30, 402, "ratio", 10, ha="middle")
            ax_bottom(cv, x0, x0 + bw, 606, [0, 0.5, 1, 1.5, 2], BX, fmt=gnum)
            cv.T(x0 + bw / 2, 638, "PAS per kb", 11, ha="middle")
            for i, (g, lab) in enumerate((("host", "Host gene promoters"), ("int", "Internal promoters"), ("itss", "Intronic start sites"),
                                          ("rand", "Random intronic positions"), ("inter", "Start sites between genes"))):
                r = row(Cc, line=line, group=g)
                s, o = (check_range("Figure 4c", f"{line} {g}", v, 0, 2) for v in (r.pas_transcribed, r.pas_opposite))
                y = 424 + i * 38
                ov = abs(BX(s) - BX(o)) < 13
                ys, yo = (y - 5, y + 5) if ov else (y, y)
                c = gcol.get(g, RED)
                if line == "K562":
                    cv.T(x0 - 12, y + 4, lab, 11, ha="end", fw="bold" if g == "itss" else "normal")
                cv.L(BX(s), ys, BX(o), yo, c="#BDBDBD", lw=3)
                cv.C(BX(o), yo, 5, fc="#FFFFFF", ec=c, lw=1.8)
                cv.C(BX(s), ys, 5, fc=c, ec=K, lw=0.7, z=6)
                cv.T(x0 + bw + 30, y + 4, f"{r.ratio:.2f}", 11, ha="middle", fw="bold" if g == "itss" else "normal")
        cv.T(44, 660, "Ratio = transcribed / opposite strand; ratio below 1: fewer polyadenylation signals in the direction of transcription than in the opposite direction.", 10)
        cv.T(44, 673, "Host promoters: opposite strand = PROMPT direction.", 10)
    return cv


# ======================================================================================== Figure 5
def figure5():
    cv = Canvas(760, 618, -32)
    gc = {"act": RED, "sil": SAL, "inter": TEAL, "rand": GRY, "prom": NAVY}
    glab = {"act": "Intronic start sites, active host", "sil": "Intronic start sites, silent host",
            "inter": "Start sites between genes", "rand": "Random intronic positions, active host", "prom": "Internal promoters"}
    # ---- a ----
    with panel(cv, "Figure 5a", (44, 66, 700, 300), (24, 58, "a")):
        A = D("fig5_a.tsv")
        letter(cv, 24, 58, "a")
        cv.T(44, 58, "Stable RNA per initiation event (unannotated sites passing the same cap test;")
        cv.T(44, 73, "internal promoters and random intronic positions for reference)")
        cv.T(44, 89, "share of sites with RNA after the start above the level before it, 95% CI; active hosts: sites ≥ 3 kb "
                     "before the next exon", 11)
        cv.push(0, 15)
        cv.T(640, 96, "median ×10$^{6}$", 10, ha="middle")
        cv.T(640, 108, "K562 / HCT116", 10, ha="middle")
        ax0, aw = 250, 300
        AX = lambda v: ax0 + v / 100 * aw
        ax_bottom(cv, ax0, ax0 + aw, 316, [0, 25, 50, 75, 100], AX)
        cv.L(AX(50), 112, AX(50), 316, lw=0.8, dash=(2, 3))
        cv.T(AX(50), 106, "half of sites", 10, ha="middle")
        cv.T(ax0 + aw / 2, 348, "% of sites with a positive excess", 11, ha="middle")

        mcol = "median_x1e6"

        def med(v):
            if not np.isfinite(v):
                return ""
            s = f"{v:.0f}" if abs(v) >= 100 else f"{v:.1f}"
            return "0.0" if s == "-0.0" else minus(s)

        for i, g in enumerate(("act", "sil", "inter", "prom", "rand")):
            y = 128 + i * 40
            k, h = row(A, line="K562", group=g), row(A, line="HCT116", group=g)
            cv.T(ax0 - 12, y + 8, glab[g], 11, ha="end")
            cv.L(AX(k.pct_lo), y, AX(k.pct_hi), y, lw=1.2)
            cv.C(AX(k.pct_positive), y, 5, fc=gc[g], ec=K, lw=0.7)
            cv.L(AX(h.pct_lo), y + 15, AX(h.pct_hi), y + 15, lw=1.2)
            cv.tri(AX(h.pct_positive), y + 15, 5.5, 5, 4, gc[g], lw=1.6)
            if g != "rand":
                cv.T(640, y + 12, f"{med(k[mcol])} / {med(h[mcol])}", 11, ha="middle", fw="bold" if g == "prom" else "normal")
        cv.push(0, 40)
        cv.C(256, 330, 4.5, fc="#7F7F7F", ec=K, lw=0.7)
        cv.T(265, 334, "K562", 11)
        cv.P([(316, 325), (320.5, 333.5), (311.5, 333.5)], fc="#FFFFFF", ec="#7F7F7F", lw=1.5)
        cv.T(326, 334, "HCT116", 11)
        cv.pop()
        cv.pop()
    # ---- b ----
    with panel(cv, "Figure 5b", (44, 382, 700, 190), (24, 374, "b")):
        cv.push(0, 55)
        Bn = D("fig5_b.tsv")
        cv.T(400, 566, "median per PRO-cap read ×1000 (log scale)", 11, ha="middle")
        cv.T(400, 581, "A confidence limit of 0 cannot be shown on a log scale and is drawn at the left edge of the axis.", 10, ha="middle")
        letter(cv, 24, 374, "b")
        cv.T(44, 374, "Nascent RNA 2–8 kb downstream per initiation event, 10 min after DRB release (HCT116)")
        cv.T(44, 390, "median excess per PRO-cap read ×1000, 95% bootstrap CI", 11)
        bx0, bw, lo_, hi_ = 250, 300, 0.03, 20
        BX = lambda v: bx0 + (math.log10(max(v, lo_)) - math.log10(lo_)) / (math.log10(hi_) - math.log10(lo_)) * bw
        ax_bottom(cv, bx0, bx0 + bw, 528, [0.03, 0.1, 0.3, 1, 3, 10], BX, fmt=gnum)
        for i, g in enumerate(("act", "sil", "inter", "prom")):
            r = row(Bn, group=g)
            y = 414 + i * 30
            m, lo, hi = (check_range("Figure 5b", g, v, -np.inf, hi_) for v in (r.median_x1e3, r.lo_x1e3, r.hi_x1e3))
            cv.T(bx0 - 12, y + 4, f"{glab[g]} ({fint(r.n)})", 11, ha="end")
            cv.L(BX(lo), y, BX(hi), y, lw=1.2)
            cv.L(BX(lo), y - 4, BX(lo), y + 4, lw=1.2)
            cv.L(BX(hi), y - 4, BX(hi), y + 4, lw=1.2)
            cv.C(BX(m), y, 5, fc=gc[g], ec=K, lw=0.7)
            cv.T(BX(hi) + 10, y + 4, f"{m:.2f}", 11)
        cv.pop()
    return cv


# =================================================================== Supplementary Figs S1, S2, S8 (88 mm)
W88 = 88.0 / MM_PER_PX  # 371.6 px


def legend_kh(cv, x, y):
    cv.C(x, y - 4, 4.5, fc="#7F7F7F", ec=K, lw=0.7)
    cv.T(x + 9, y, "K562", 10.5)
    cv.P([(x + 52, y - 9), (x + 56.5, y - 0.5), (x + 47.5, y - 0.5)], fc="#FFFFFF", ec="#7F7F7F", lw=1.5)
    cv.T(x + 62, y, "HCT116", 10.5)


def figS1():
    cv = Canvas(W88, 318)
    T = D("figS1.tsv")
    col = {"rand": GRY, "annot": LAV, "int": NAVY, "itss": RED}
    x0, x1 = 158, 298
    X = lambda v: x0 + (math.log10(max(v, 0.1)) + 1) / 3 * (x1 - x0)
    cv.T(336, 30, "K562 / HCT116", 9.5, ha="middle", fw="bold")
    labs = {"rand": ["Random intronic", "positions"], "annot": ["Annotated", "intronic RNAs"],
            "int": ["Internal promoters", "(annotated)"], "itss": ["Intronic PRO-cap", "peak summits"]}
    for i, g in enumerate(("rand", "annot", "int", "itss")):
        y = 52 + i * 54
        k, h = row(T, line="K562", group=g), row(T, line="HCT116", group=g)
        for j, l in enumerate(labs[g]):
            cv.T(x0 - 10, y + j * 12, l, 10.5, ha="end")
        cv.T(x0 - 10, y + 25, f"n = {fint(k.n)} / {fint(h.n)}", 9.5, ha="end")
        cv.L(X(k.lo), y + 2, X(k.hi), y + 2, lw=1.2)
        cv.C(X(k.pct), y + 2, 4.5, fc=col[g], ec=K, lw=0.7)
        cv.L(X(h.lo), y + 16, X(h.hi), y + 16, lw=1.2)
        cv.tri(X(h.pct), y + 16, 5, 4.5, 4, col[g], lw=1.5)
        f = lambda v: f"{v:.2f}" if v < 1 else f"{v:.1f}"
        cv.T(336, y + 13, f"{f(k.pct)} / {f(h.pct)}%", 10.5, ha="middle", fw="bold" if g == "itss" else "normal")
    ax_bottom(cv, x0, x1, 262, [0.1, 1, 10, 100], X, fmt=lambda t: f"{t:g}%", fs=10.5)
    cv.T((x0 + x1) / 2, 292, "% of sites passing (log scale)", 10.5, ha="middle")
    legend_kh(cv, 160, 312)
    return cv


def figS2():
    cv = Canvas(W88, 300)
    T, Cr = D("figS2.tsv"), D("figS2_cor.tsv")
    x0, x1 = 162, 300
    X = lambda v: x0 + (v - 50) / 50 * (x1 - x0)
    cv.T(338, 26, "passed", 9.5, ha="middle", fw="bold")
    names = {"weak": "Weakest third", "medium": "Middle third", "strong": "Strongest third", "all": "All sites"}
    for i, th in enumerate(("weak", "medium", "strong", "all")):
        r = row(T, third=th)
        y = 48 + i * 50 + (14 if th == "all" else 0)
        if th == "all":
            cv.L(16, y - 24, W88 - 16, y - 24, c=LGR)
        cv.T(x0 - 10, y, names[th], 10.5, ha="end", fw="bold" if th == "all" else "normal")
        if th == "all":
            cv.T(x0 - 10, y + 13, f"n = {fint(r.n)}", 9.5, ha="end")
        else:
            cv.T(x0 - 10, y + 13, f"median {fint(r.median_reads)} PRO-cap reads", 9.5, ha="end")
            cv.T(x0 - 10, y + 25, f"n = {fint(r.n)}", 9.5, ha="end")
        p, lo, hi = (check_range("Supplementary Figure S2", th, v, 50, 100) for v in (r.pct, r.lo, r.hi))
        cv.L(X(lo), y + 4, X(hi), y + 4, lw=1.2)
        cv.C(X(p), y + 4, 5, fc=RED, alpha=1 if th == "all" else 0.85, ec=K, ealpha=1, lw=0.7)
        cv.T(338, y + 8, f"{r.pct:.1f}%", 10.5, ha="middle", fw="bold" if th == "all" else "normal")
    ax_bottom(cv, x0, x1, 234, [50, 60, 70, 80, 90, 100], X, fs=10.5)
    cv.T((x0 + x1) / 2, 264, "% passing in the RAD21 line (95% CI)", 10.5, ha="middle")
    c = Cr.iloc[0]
    cv.T(16, 284, "PRO-cap reads at sites passing in both lines:", 9.5)
    cv.T(16, 296, f"Spearman r = {c.spearman_r:.2f} (n = {fint(c.n)}).", 9.5)
    return cv


def figS8():
    cv = Canvas(W88, 346)
    T = D("figS8.tsv")
    x0, x1 = 150, 352
    X = lambda v: x0 + (v - 0.4) / 0.9 * (x1 - x0)
    legend_kh(cv, 196, 22)
    cv.T(x1, 38, "whiskers: 95% bootstrap CI", 9.5, ha="end")
    labs = {"inter": ["Pairs between genes", "(orientation arbitrary)"], "silent": ["Pairs in introns", "of silent genes"],
            "active": ["Pairs in introns", "of active genes"], "k36": ["… with H3K36me3", "within 1 kb"],
            "nok36": ["… without H3K36me3"]}
    for i, g in enumerate(("inter", "silent", "active", "k36", "nok36")):
        y = 60 + i * 46
        k, h = row(T, line="K562", group=g), row(T, line="HCT116", group=g)
        cc = TEAL if g == "inter" else SAL if g == "silent" else RED
        for j, l in enumerate(labs[g]):
            cv.T(x0 - 10, y + 4 + j * 12, l, 10.5, ha="end")
        for r, yy, kind in ((k, y, "k"), (h, y + 13, "h")):
            m, lo, hi = (check_range("Supplementary Figure S8", f"{g}", v, 0.4, 1.3) for v in (r.median_ratio, r.ci_low, r.ci_high))
            cv.L(X(lo), yy, X(hi), yy, lw=1.2)
            cv.L(X(lo), yy - 3, X(lo), yy + 3, lw=1.2)
            cv.L(X(hi), yy - 3, X(hi), yy + 3, lw=1.2)
            if kind == "k":
                cv.C(X(m), yy, 5, fc=cc, ec=K, lw=0.7, z=6)
            else:
                cv.tri(X(m), yy, 5.5, 5.5, 4.5, cc, lw=1.6)
    ax_bottom(cv, X(0.4), X(1.3), 290, [0.4, 0.6, 0.8, 1.0, 1.2], X, fmt=lambda t: f"{t:.1f}", fs=10.5)
    cv.L(X(1), 48, X(1), 290, lw=0.8, dash=(2, 3))
    cv.T((x0 + x1) / 2, 324, "ratio of the two arms", 10.5, ha="middle")
    cv.T((x0 + x1) / 2, 337, "(gene strand / opposite strand)", 10.5, ha="middle")
    return cv


# ============================================================================== Supplementary S3-S7, S9
def figS3():
    cv = Canvas(760, 357, -43)
    T = D("figS3.tsv")
    T = T[T.shown.astype(str) == "True"]
    col = {"10min": "#0072B2", "20min": "#56B4E9", "30min": "#999999", "DRB10": "#0072B2"}
    pw, y0, y1 = 270, 92, 300
    for line, x0, ymax, ticks, lab, series in (
            ("HCT116", 80, 3.2, [0, 1, 2, 3], "HCT116, DRB/TTchem-seq / 40 min", ["10min", "20min", "30min"]),
            ("K562", 450, 2.4, [0, 0.5, 1, 1.5, 2], "K562, BruDRB-seq 10 min / Bru-seq", ["DRB10"])):
        X = lambda kb: x0 + kb / 100 * pw
        Y = lambda v: y1 - v / ymax * (y1 - y0)
        letter(cv, x0 - 56, 76, "a" if line == "HCT116" else "b")
        cv.T(x0 - 36, 76, lab, fw="bold")
        ax_left(cv, x0, y0, y1, ticks, Y, fmt=gnum)
        ax_bottom(cv, x0, x0 + pw, y1, [0, 20, 40, 60, 80, 100], X)
        cv.T(x0 + pw / 2, y1 + 32, "kb from the host gene’s start", 11, ha="middle")
        cv.L(x0, Y(1), x0 + pw, Y(1), lw=0.8, dash=(2, 3))
        cv.T(x0 - 40, (y0 + y1) / 2, "nascent RNA / reference", 11, ha="middle", rot=90)
        for s in series:
            q = T[(T.line == line) & (T.label == s)].sort_values("kb_mid")
            pts = [(X(r.kb_mid), Y(check_range("Supplementary Figure S3", f"{line} {s}", r.ratio, 0, ymax))) for r in q.itertuples()]
            cv.PL(pts, c=col[s], lw=2)
            for x, y in pts:
                cv.C(x, y, 2.2, fc=col[s], z=6)
    for x, c, lab, dash, lw in ((80, "#0072B2", "10 min after release", None, 2), (240, "#56B4E9", "20 min", None, 2),
                                (330, "#999999", "30 min", None, 2), (420, K, "equal to the reference", (2, 3), 0.8)):
        cv.L(x, 372, x + (24 if dash else 20), 372, c=c, lw=lw, dash=dash)
        cv.T(x + (30 if dash else 26), 376, lab, 11)
    return cv


def figS4():
    cv = Canvas(760, 344, -36)
    P, Bt = D("figS4_a.tsv"), D("figS4_b.tsv")
    sc = {"itss": RED, "int": NAVY, "rand": GRY, "sil": SAL}
    ax0, aw, ay0, ay1 = 70, 280, 92, 280
    AX = lambda kb: ax0 + (kb + 10) / 20 * aw
    AY = lambda v: ay1 - (v - 0.6) / 0.6 * (ay1 - ay0)
    letter(cv, 24, 62, "a")
    cv.T(44, 62, "Nascent RNA around sites beyond the wave front")
    ax_left(cv, ax0, ay0, ay1, [0.6, 0.8, 1.0, 1.2], AY, fmt=lambda t: f"{t:.1f}")
    ax_bottom(cv, ax0, ax0 + aw, ay1, [-10, -5, 0, 5, 10], AX, fmt=lambda t: f"+{t}" if t > 0 else f"{t}")
    cv.L(AX(0), ay0, AX(0), ay1, lw=0.8, dash=(2, 3))
    cv.T(ax0 + aw / 2, ay1 + 32, "kb from the site (direction of transcription)", 11, ha="middle")
    cv.T(ax0 - 40, (ay0 + ay1) / 2, "DRB 10 min / Bru-seq", 11, ha="middle", rot=90)
    nm = {"itss": "Intronic start sites, active host", "int": "Internal promoters", "rand": "Random intronic positions"}
    for s in ("rand", "int", "itss"):
        q = P[P.group == s].sort_values("kb_mid")
        pts = [(AX(r.kb_mid), AY(check_range("Supplementary Figure S4a", s, r.ratio, 0.6, 1.2))) for r in q.itertuples()]
        cv.PL(pts, c=sc[s], lw=1.3 if s == "rand" else 2, dash=(4, 3) if s == "rand" else None)
        for x, y in pts:
            cv.C(x, y, 1.6 if s == "rand" else 2.4, fc=sc[s], z=6)
    for i, s in enumerate(("itss", "int", "rand")):
        n = fint(P[P.group == s].n_sites.iloc[0])
        cv.L(70, 336 + i * 16, 90, 336 + i * 16, c=sc[s], lw=1.3 if s == "rand" else 2, dash=(4, 3) if s == "rand" else None)
        cv.T(96, 340 + i * 16, f"{nm[s]} (n = {n})", 11)
    letter(cv, 384, 62, "b")
    cv.T(404, 62, "Sites above the random-position threshold")
    bx0, bw = 560, 170
    BX = lambda v: bx0 + v / 40 * bw
    ax_bottom(cv, BX(0), BX(40), 280, [0, 20, 40], BX, fmt=lambda t: f"{t}%")
    cv.L(BX(5), 92, BX(5), 280, lw=0.8, dash=(2, 3))
    cv.T(BX(5) + 3, 90, "chance", 10)
    for i, (g, lab) in enumerate((("rand", "Random intronic positions"), ("int", "Internal promoters"),
                                  ("itss", "Intronic start sites, active host"), ("sil", "Intronic start sites, silent host"))):
        r = row(Bt, group=g)
        y = 112 + i * 44
        cv.T(bx0 - 10, y + 2, lab, 11, ha="end")
        cv.T(bx0 - 10, y + 15, f"n = {fint(r.n)}" + ("" if g == "rand" else f"; p = {fmt_p(r.p)}"), 10, ha="end")
        p, lo, hi = (check_range("Supplementary Figure S4b", g, v, 0, 40) for v in (r.pct, r.lo, r.hi))
        cv.L(BX(lo), y + 6, BX(hi), y + 6, lw=1.2)
        cv.C(BX(p), y + 6, 5, fc=sc[g], ec=K, lw=0.7)
    return cv


def figS5():
    cv = Canvas(760, 536, -24)
    A, Bt = D("figS5_a.tsv"), D("figS5_b.tsv")
    verdict = {"MTREX": "rise: knockout worked (partly)", "EXOSC10": "no rise; gene starts fall", "EXOSC5": "no rise in PROMPTs"}
    letter(cv, 24, 60, "a")
    cv.T(44, 60, f"Knockout check at {fint(val(A, 'n_promoters', ko='MTREX'))} active host promoters: RNA change "
                 f"(log2, knockout / control), median")
    ax0, aw = 160, 260
    AX = lambda v: ax0 + (v + 0.2) / 0.6 * aw
    ax_bottom(cv, AX(-0.2), AX(0.4), 200, [-0.2, 0, 0.2, 0.4], AX, fmt=lambda t: f"+{t:.1f}" if t > 0 else f"{t:.1f}")
    cv.L(AX(0), 78, AX(0), 200, lw=0.8, dash=(2, 3))
    cv.T(ax0 + aw / 2, 232, "log2 fold change", 11, ha="middle")
    cv.T(440, 84, "PROMPT higher than gene start", 10, fw="bold")
    for i, ko in enumerate(("MTREX", "EXOSC10", "EXOSC5")):
        r = row(A, ko=ko)
        y = 102 + i * 36
        pr, ge = (check_range("Supplementary Figure S5a", ko, v, -0.2, 0.4) for v in (r.prompt_median, r.gene_median))
        cv.T(ax0 - 12, y + 4, ko, 12, ha="end", fw="bold")
        cv.L(AX(ge), y, AX(pr), y, c="#BDBDBD", lw=3)
        cv.C(AX(ge), y, 5, fc="#FFFFFF", ec=GRY, lw=1.6)
        cv.C(AX(pr), y, 5, fc=NAVY, ec=K, lw=0.7)
        cv.T(440, y + 4, f"{r.pct_prompt_higher:.0f}% of promoters; {verdict[ko]}", 11)
    cv.C(170, 252, 4.5, fc=NAVY, ec=K, lw=0.7)
    cv.T(179, 256, "upstream antisense RNA (PROMPT, 200–2,000 nt upstream)", 11)
    cv.C(520, 252, 4.5, fc="#FFFFFF", ec=GRY, lw=1.5)
    cv.T(529, 256, "first 3 kb of the gene", 11)
    letter(cv, 24, 296, "b")
    cv.T(44, 296, "Rise in RNA after each kind of site: % of sites above the 95th percentile of random intronic positions")
    gcl = {"int": NAVY, "sil": SAL}
    gn = {"int": "Internal promoters", "itss": "Intronic start sites, active host", "prom": "… promoter-like",
          "enh": "… enhancer-like", "sil": "Intronic start sites, silent host"}
    bw = 120
    for ko, x0 in (("MTREX", 250), ("EXOSC10", 435), ("EXOSC5", 620)):
        BX = lambda v: x0 + v / 50 * bw
        cv.T(x0 + bw / 2, 322, ko, fw="bold", ha="middle")
        ax_bottom(cv, x0, x0 + bw, 520, [0, 25, 50], BX, fmt=lambda t: f"{t}%")
        cv.L(BX(5), 332, BX(5), 520, lw=0.8, dash=(2, 3))
        for i, g in enumerate(("int", "itss", "prom", "enh", "sil")):
            r = row(Bt, ko=ko, group=g)
            y = 346 + i * 36
            v = check_range("Supplementary Figure S5b", f"{ko} {g}", r.pct, 0, 50)
            cv.C(BX(v), y, 5, fc=gcl.get(g, RED), alpha=0.6 if g in ("prom", "enh") else 1, ec=K, ealpha=1, lw=0.7)
            if v < 25:
                cv.T(BX(v) + 8, y + 4, f"p = {fmt_p(r.p)}", 10)
            else:
                cv.T(BX(v) - 8, y + 4, f"p = {fmt_p(r.p)}", 10, ha="end")
    for i, g in enumerate(("int", "itss", "prom", "enh", "sil")):
        cv.T(236, 346 + i * 36 + 4, gn[g], 11, ha="end")
    cv.T(24, 552, "Dotted lines: 5% expected by chance. p: one-sided Mann–Whitney tests against random intronic positions.", 10)
    return cv


def figS6():
    cv = Canvas(760, 352, -48)
    T = D("figS6.tsv")
    gc = {"act": RED, "sil": SAL, "inter": TEAL}
    dg = {"act": -6, "sil": 0, "inter": 6}
    pw, y0, y1 = 230, 96, 290
    PY = lambda v: y1 - v / 40 * (y1 - y0)
    for line, x0 in (("K562", 90), ("HCT116", 450)):
        PX = lambda t: x0 + 30 + t * (pw - 60) / 2
        cv.T(x0, 80, line, fw="bold")
        ax_left(cv, x0, y0, y1, [0, 10, 20, 30, 40], PY)
        cv.L(x0, y1, x0 + pw, y1)
        for t, th in enumerate(("weak", "medium", "strong")):
            cv.L(PX(t), y1, PX(t), y1 + 4)
            cv.T(PX(t), y1 + 16, th, 11, ha="middle")
        cv.T(x0 + pw / 2, y1 + 32, "third of PRO-cap strength", 11, ha="middle")
        if line == "K562":
            cv.T(x0 - 36, (y0 + y1) / 2, "% of sites with stable RNA", 11, ha="middle", rot=90)
        ref = val(T, "pct", line=line, group="prom", third="all")
        cv.L(x0, PY(ref), x0 + pw, PY(ref), c=NAVY, lw=1.2, dash=(2, 3))
        cv.T(x0 + 6, PY(ref) - 5, f"internal promoters, {ref:.1f}%", 10)
        for g in ("inter", "sil", "act"):
            pts = [(PX(t) + dg[g], PY(check_range("Supplementary Figure S6", f"{line} {g}", val(T, "pct", line=line, group=g, third=th), 0, 40)))
                   for t, th in enumerate(("weak", "medium", "strong"))]
            cv.PL(pts, c=gc[g], lw=2, dash=(5, 3) if g == "sil" else None)
            for x, y in pts:
                cv.C(x, y, 4, fc=gc[g], ec=K, lw=0.6)
    for x, g, lab in ((90, "act", "intronic, active host"), (270, "sil", "intronic, silent host"), (450, "inter", "between genes")):
        cv.L(x, 364, x + 22, 364, c=gc[g], lw=2, dash=(5, 3) if g == "sil" else None)
        cv.C(x + 11, 364, 3.5, fc=gc[g], ec=K, lw=0.6)
        cv.T(x + 28, 368, lab, 11)
    return cv


def figS7():
    cv = Canvas(760, 615, -30)
    T, O = D("figS7.tsv"), D("figS7_oe.tsv")
    gc = {"same": RED, "anti": NAVY, "stand": GRY}
    pw, ph = 250, 170
    for ty, y0, lab in (("lncRNA", 92, "long non-coding genes"), ("protein_coding", 362, "protein-coding genes")):
        for line, x0 in (("K562", 90), ("HCT116", 470)):
            X = lambda b: x0 + b / 5 * pw
            Y = lambda v: y0 + ph - v / 100 * ph
            cv.T(x0, y0 - 30, f"{line}, {lab}", 12, fw="bold")
            ax_left(cv, x0, y0, y0 + ph, [0, 50, 100], Y)
            ax_bottom(cv, x0, x0 + pw, y0 + ph, range(6), X)
            cv.T(x0 + pw / 2, y0 + ph + 32, "expression bin (0 = none, 5 = highest)", 11, ha="middle")
            cv.T(x0 - 36, y0 + ph / 2, "% with own start", 11, ha="middle", rot=90)
            for g in ("stand", "anti", "same"):
                q = T[(T.line == line) & (T.type == ty) & (T.group == g)].sort_values("bin")
                q = q[np.isfinite(q.pct)]
                pts = [(X(r.bin), Y(r.pct)) for r in q.itertuples()]
                cv.PL(pts, c=gc[g], lw=1.4 if g == "stand" else 2, dash=(5, 3) if g == "stand" else None)
                for (x, y), r in zip(pts, q.itertuples()):
                    small = r.n < 10
                    cv.C(x, y, 3 if small else 4, fc="#FFFFFF" if small else gc[g], ec=gc[g], lw=1.4)
            for i, g in enumerate(("same", "anti")):
                oe = val(O, "oe", line=line, type=ty, group=g)
                cv.T(x0 + i * 120, y0 - 12, f"{'host strand' if g == 'same' else 'opposite strand'} o/e {oe:.2f}", 11, c=gc[g])
    cv.L(90, 606, 114, 606, c=GRY, lw=1.4, dash=(5, 3))
    cv.T(120, 610, "stand-alone genes", 11)
    cv.L(280, 606, 304, 606, c=RED, lw=2)
    cv.T(310, 610, "nested, host strand", 11)
    cv.L(470, 606, 494, 606, c=NAVY, lw=2)
    cv.T(500, 610, "nested, opposite strand", 11)
    cv.C(96, 630, 3, fc="#FFFFFF", ec=K, lw=1.2)
    cv.T(106, 634, "open symbols: fewer than 10 genes in the bin", 11)
    return cv


def figS9():
    cv = Canvas(760, 325, -20)
    T = D("figS9.tsv")
    BL = ["1–2", "2–3", "3–5", "5–10", "10–20", ">20"]
    pw, y0, y1 = 180, 92, 252
    gcol = {"itss": RED, "rand": GRY}
    for x0, let, title, sub, meas, mx, ticks, shade in (
            (70, "a", "Apparent stable RNA, all distances", "% of sites above the original threshold", "stable", 40, [0, 10, 20, 30, 40], True),
            (325, "b", "Initiation strength", "median capped reads, intronic start sites", "capped_median", 60, [0, 20, 40, 60], False),
            (575, "c", "Response after DRB release", "% of sites above random 95th percentile", "drb", 20, [0, 5, 10, 15, 20], False)):
        X = lambda b: x0 + 12 + b * (pw - 24) / 5
        Y = lambda v: y1 - v / mx * (y1 - y0)
        letter(cv, x0 - 46, 52, let)
        cv.T(x0 - 30, 52, title)
        cv.T(x0 - 30, 67, sub, 10.5)
        if shade:
            cv.R(X(0) - 10, y0, X(1) - X(0) + 20, y1 - y0, fc="#EFEFEF", z=1)
            cv.T((X(0) + X(1)) / 2, y1 - 7, "excluded", 10, ha="middle")
        ax_left(cv, x0, y0, y1, ticks, Y, fs=10.5)
        cv.L(x0, y1, x0 + pw, y1)
        for b_, l in enumerate(BL):
            cv.L(X(b_), y1, X(b_), y1 + 4)
            cv.T(X(b_), y1 + 15, l, 10, ha="middle")
        cv.T(x0 + pw / 2, y1 + 30, "kb to the next host exon", 10.5, ha="middle")
        for g in ("rand", "itss"):
            for line in ("K562", "HCT116"):
                q = T[(T.measure == meas) & (T.group == g) & (T.line == line)].sort_values("bin")
                q = q[np.isfinite(q.value)]
                if not len(q):
                    continue
                pts = [(X(r.bin), Y(check_range("Supplementary Figure S9", f"{meas} {line} {g}", r.value, 0, mx))) for r in q.itertuples()]
                cv.PL(pts, c=gcol[g], lw=1.6, dash=None if line == "K562" else (4, 3))
                for x, y in pts:
                    if line == "K562":
                        cv.C(x, y, 3.8, fc=gcol[g], ec=K, lw=0.6)
                    else:
                        cv.tri(x, y, 4.6, 4.2, 3.2, gcol[g], lw=1.5)
    cv.L(70, 322, 92, 322, c=RED, lw=1.6)
    cv.T(98, 326, "intronic start sites", 11)
    cv.L(230, 322, 252, 322, c=GRY, lw=1.6)
    cv.T(258, 326, "random intronic positions", 11)
    cv.C(430, 322, 3.8, fc="#7F7F7F", ec=K, lw=0.6)
    cv.T(440, 326, "K562 (solid)", 11)
    cv.P([(550, 317.4), (554.2, 325.2), (545.8, 325.2)], fc="#FFFFFF", ec="#7F7F7F", lw=1.5)
    cv.T(561, 326, "HCT116 (dashed)", 11)
    return cv


FIGURES = [("Figure_1", figure1), ("Figure_2", figure2), ("Figure_3", figure3), ("Figure_4", figure4), ("Figure_5", figure5),
           ("Supplementary_Figure_S1", figS1), ("Supplementary_Figure_S2", figS2), ("Supplementary_Figure_S3", figS3),
           ("Supplementary_Figure_S4", figS4), ("Supplementary_Figure_S5", figS5), ("Supplementary_Figure_S6", figS6),
           ("Supplementary_Figure_S7", figS7), ("Supplementary_Figure_S8", figS8), ("Supplementary_Figure_S9", figS9)]


def main():
    global ARGS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--data", default=DATA)
    ap.add_argument("--out", default=OUT)
    ARGS = ap.parse_args()
    print(f"51_figures.py  matplotlib {matplotlib.__version__}  font {FONT} ({FONT_PATH})")
    if FONT not in ("Arial", "Liberation Sans", "Helvetica"):
        print("  WARNING: neither Arial nor Liberation Sans found; install one or the text will use DejaVu Sans")
    ok = 0
    for name, f in FIGURES:
        if ARGS.only and name not in ARGS.only:
            continue
        n0 = len(WARN)
        try:
            cv = f()
            if cv.skipped and not cv.drawn:
                raise MissingTable(", ".join(cv.skipped) + ": no panel table")
            w_mm, h_mm = cv.w * MM_PER_PX, cv.h * MM_PER_PX
            pdf, png = cv.save(name, ARGS.out)
            tag = "[ok]" if not cv.skipped else "[partial]"
            print(f"{tag} {name}: {w_mm:.0f} x {h_mm:.0f} mm -> {pdf}, {png}"
                  + (f"  (not drawn: {', '.join(cv.skipped)})" if cv.skipped else ""))
            ok += 1
        except FileNotFoundError as e:
            plt.close("all")
            print(f"[SKIP] {name}: {e} not found (its panel is missing or failed in 50_figure_data.py; "
                  f"see results/figures/data/_status.tsv)")
        except Exception as e:  # one figure must not stop the others
            plt.close("all")
            tb = traceback.format_exc().strip().splitlines()
            print(f"[error] {name}: {type(e).__name__}: {e}\n    " + "\n    ".join(tb[-5:]))
        for w in warn_lines(WARN[n0:]):
            print(f"    warning: {w}")
    print(f"{ok} figures written to {ARGS.out}/ (Fig. 1a and Fig. 7 are schematics, not generated here)")


if __name__ == "__main__":
    main()
