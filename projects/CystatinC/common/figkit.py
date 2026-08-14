"""shared figure helpers -- collision-free layer labelling."""
import numpy as np
COL = {"prism": "#eef2f7", "Al2O3": "#cfe3f7", "Cu": "#e08a4a",
       "Ni": "#9aa0a6", "ZnS": "#f4e08a", "sense": "#dff0e4"}
PRETTY = {"Al2O3": r"Al$_2$O$_3$", "sense": "analyte"}


def _declash(mids, lo, hi, minsep):
    """push overlapping label positions apart, keeping them inside [lo, hi]."""
    order = np.argsort(mids)
    pos = np.array(mids, dtype=float)
    s = pos[order]
    for k in range(1, len(s)):
        if s[k] - s[k - 1] < minsep:
            s[k] = s[k - 1] + minsep
    over = s[-1] - (hi - 0.02 * (hi - lo))
    if over > 0:
        s -= over
        for k in range(len(s) - 2, -1, -1):
            if s[k + 1] - s[k] < minsep:
                s[k] = s[k + 1] - minsep
    pos[order] = s
    return pos


def stack_strip(ax, bands, lo, hi, fs=7, minfrac=0.11, show_t=True):
    """
    narrow panel: coloured layer bands on the left, labels on the right,
    joined by leader lines. label positions are de-collided, so thin films
    can never write on top of one another.
    """
    keep = [(nm, a, b) for nm, a, b in bands if b > lo and a < hi]
    for nm, a, b in keep:
        ax.axhspan(max(a, lo), min(b, hi), xmin=0.0, xmax=0.42,
                   color=COL[nm], ec="k", lw=0.4)
    mids = [0.5 * (max(a, lo) + min(b, hi)) for nm, a, b in keep]
    pos = _declash(mids, lo, hi, minfrac * (hi - lo))
    for (nm, a, b), m, p in zip(keep, mids, pos):
        txt = PRETTY.get(nm, nm)
        if show_t and nm not in ("prism", "sense"):
            txt += "\n%g nm" % (b - a)
        ax.plot([0.42, 0.52], [m, p], color="0.45", lw=0.5,
                clip_on=False, solid_capstyle="butt")
        ax.text(0.56, p, txt, fontsize=fs, va="center", ha="left",
                linespacing=1.3, color="0.15")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(hi, lo)
    ax.set_xticks([])
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(left=False, labelleft=False)
