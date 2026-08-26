"""Plotting helpers (Agg backend, consistent style)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def savefig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def newfig(w=7.0, h=4.5):
    return plt.subplots(figsize=(w, h))
