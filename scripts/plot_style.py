"""
plot_style.py — Shared academic-paper style for all analytical plots.

Import once at the top of each plot script:

    from plot_style import apply_academic_style, COLORS, format_pvalue, add_source_note
    apply_academic_style()
"""
import matplotlib.pyplot as plt


COLORS = {
    "treated":     "#C0392B",  # red — treatment / post-period / actual
    "control":     "#1F4E79",  # deep blue — control / pre-period
    "neutral":     "#7F7F7F",  # gray — placebos, null, secondary
    "accent":      "#D4A017",  # gold — highlights
    "covid_shade": "#F2D7D5",  # very light pink — COVID period
    "pre_shade":   "#EAF2F8",  # very light blue — pre-treatment period
    "grid":        "#E5E5E5",
    "spine":       "#666666",
    "text_muted":  "#555555",

    # Legacy aliases used by some scripts (so imports don't break mid-refactor)
    "blue":        "#1F4E79",
    "red":         "#C0392B",
    "gray":        "#7F7F7F",
    "green":       "#4DAF4A",
    "orange":      "#F4A582",
    "purple":      "#8856A7",
}


def apply_academic_style():
    """Set rcParams for a clean academic-paper look. Call once per script."""
    plt.rcParams.update({
        "font.family":       "sans-serif",
        "font.sans-serif":   ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size":         11,
        "axes.titlesize":    13,
        "axes.titleweight":  "bold",
        "axes.labelsize":    11,
        "xtick.labelsize":   10,
        "ytick.labelsize":   10,
        "legend.fontsize":   10,
        "figure.dpi":        200,
        "savefig.dpi":       200,
        "savefig.bbox":      "tight",
        "savefig.facecolor": "white",
        "figure.facecolor":  "white",
        "axes.facecolor":    "white",
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.edgecolor":    COLORS["spine"],
        "axes.linewidth":    0.8,
        "axes.axisbelow":    True,
        "grid.color":        COLORS["grid"],
        "grid.linewidth":    0.5,
        "grid.alpha":        1.0,
        "axes.grid":         False,
    })


def format_pvalue(p):
    """Format a p-value the academic way: 'p < 0.001' for small p, else 'p = 0.123'."""
    if p is None:
        return "p = n/a"
    if p < 0.001:
        return "p < 0.001"
    if p < 0.01:
        return f"p = {p:.3f}"
    return f"p = {p:.3f}"


def sig_stars(p):
    """Return significance stars for a p-value."""
    if p is None:
        return ""
    if p < 0.01:
        return "***"
    if p < 0.05:
        return "**"
    if p < 0.1:
        return "*"
    return ""


def add_source_note(fig, text, *, sample_n=None, y=0.005, x=0.01):
    """Write a small italic gray source note at the bottom-left of the figure."""
    full = text if not sample_n else f"{text}  |  {sample_n}"
    fig.text(x, y, full, fontsize=9, color=COLORS["text_muted"],
             style="italic", ha="left", va="bottom")


def add_n_label(ax, text, loc="upper left"):
    """Place a small sample-size box on the axes."""
    coords = {
        "upper left":  (0.02, 0.97, "left",  "top"),
        "upper right": (0.98, 0.97, "right", "top"),
        "lower left":  (0.02, 0.03, "left",  "bottom"),
        "lower right": (0.98, 0.03, "right", "bottom"),
    }[loc]
    x, y, ha, va = coords
    ax.text(x, y, text, transform=ax.transAxes, fontsize=9,
            color=COLORS["text_muted"], ha=ha, va=va,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor=COLORS["grid"], linewidth=0.5, alpha=0.9))


def clean_spines(ax):
    """Hide top/right spines and lighten left/bottom — call after plotting."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(COLORS["spine"])
        ax.spines[side].set_linewidth(0.8)


SOURCE_DEFAULT = ("Source: City of Chicago Data Portal — Violence Reduction "
                  "homicide & ShotSpotter alerts (2009–2023).")
