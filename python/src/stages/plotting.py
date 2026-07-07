"""Diagram plotting for the stages package.

Matplotlib is an optional dependency (install the ``stages-thermo[plot]``
extra); everything else in ``stages`` works without it. The import is
deferred into the functions so ``import stages.plotting`` alone doesn't pull
matplotlib in.

The M1 surface is the McCabe–Thiele staircase diagram. Later milestones add
the H–x–y (Ponchon–Savarit) diagram, Gilliland charts, and column profile /
convergence-history plots.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover - typing only
    from matplotlib.axes import Axes

    from stages import EquilibriumCurve, McCabeThieleResult, TotalRefluxResult

__all__ = ["plot_mccabe_thiele", "plot_total_reflux"]


def _require_pyplot() -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - exercised without extra
        raise ImportError(
            "matplotlib is required for plotting — install the "
            "'stages-thermo[plot]' extra"
        ) from exc
    return plt


def _draw_frame(ax: "Axes", curve: "EquilibriumCurve") -> None:
    """The parts every x-y diagram shares: equilibrium curve + diagonal."""
    ax.plot(curve.x, curve.y, lw=2, label="equilibrium $y^*(x)$", zorder=3)
    ax.plot([0, 1], [0, 1], color="0.4", lw=1, label="$y = x$")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("$x$ — light component in liquid")
    ax.set_ylabel("$y$ — light component in vapor")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)


def plot_mccabe_thiele(
    result: "McCabeThieleResult",
    curve: "EquilibriumCurve",
    ax: "Axes | None" = None,
    show_rmin: bool = False,
) -> "Axes":
    """Draw the full McCabe–Thiele diagram for a construction result.

    Parameters
    ----------
    result:
        The construction returned by :func:`stages.mccabe_thiele`.
    curve:
        The equilibrium curve the construction was run on.
    ax:
        Existing matplotlib axes to draw into; a new figure is created when
        omitted.
    show_rmin:
        Also draw the limiting (minimum-reflux) operating line through the
        controlling pinch point, dashed.

    Returns
    -------
    matplotlib.axes.Axes
        The axes drawn into — call ``ax.figure.savefig(...)`` or let the
        notebook render it.
    """
    plt = _require_pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))

    _draw_frame(ax, curve)

    spec = result.spec
    xd, xb, zf = spec.x_distillate, spec.x_bottoms, spec.z_feed
    xi, yi = result.intersection

    # Operating lines, drawn as the segments actually used by the staircase.
    ax.plot([xd, xi], [xd, yi], color="tab:red", lw=1.5, label="rectifying line")
    ax.plot([xb, xi], [xb, yi], color="tab:green", lw=1.5, label="stripping line")
    # q-line from the feed point on the diagonal to the intersection.
    ax.plot([zf, xi], [zf, yi], color="tab:purple", lw=1.5, ls="-.", label="$q$-line")

    # The staircase itself.
    sx = [p[0] for p in result.staircase]
    sy = [p[1] for p in result.staircase]
    ax.plot(sx, sy, color="tab:blue", lw=1.0, drawstyle="default", label="stages")

    # Stage numbers at the equilibrium-curve corners.
    for stage in result.stages:
        ax.annotate(
            str(stage.index),
            (stage.x, stage.y),
            textcoords="offset points",
            xytext=(-8, 4),
            fontsize=8,
            color="tab:blue",
        )

    # Anchor compositions on the diagonal.
    for label, xv in (("$x_B$", xb), ("$z_F$", zf), ("$x_D$", xd)):
        ax.plot([xv], [xv], marker="o", ms=4, color="k")
        ax.annotate(label, (xv, xv), textcoords="offset points", xytext=(6, -12))

    if show_rmin:
        px, py = result.rmin.pinch
        ax.plot(
            [xd, px],
            [xd, py],
            color="tab:red",
            lw=1.2,
            ls="--",
            label=f"limiting line ($R_{{min}}$ = {result.rmin.r_min:.3f})",
        )
        ax.plot([px], [py], marker="x", ms=8, color="tab:red")

    ax.set_title(
        f"McCabe–Thiele: N = {result.n_stages:.2f} stages, "
        f"feed stage {result.feed_stage}, R = {spec.reflux:.3f}"
    )
    ax.legend(loc="upper left", fontsize=9)
    return ax


def plot_total_reflux(
    result: "TotalRefluxResult",
    curve: "EquilibriumCurve",
    ax: "Axes | None" = None,
) -> "Axes":
    """Draw the total-reflux (minimum stages) staircase against the diagonal.

    Parameters mirror :func:`plot_mccabe_thiele`; the operating lines
    collapse onto ``y = x`` so only the staircase and frame are drawn.
    """
    plt = _require_pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))

    _draw_frame(ax, curve)

    sx = [p[0] for p in result.staircase]
    sy = [p[1] for p in result.staircase]
    ax.plot(sx, sy, color="tab:blue", lw=1.0, label="stages at total reflux")
    for stage in result.stages:
        ax.annotate(
            str(stage.index),
            (stage.x, stage.y),
            textcoords="offset points",
            xytext=(-8, 4),
            fontsize=8,
            color="tab:blue",
        )
    ax.set_title(f"Total reflux: $N_{{min}}$ = {result.n_min:.2f} stages")
    ax.legend(loc="upper left", fontsize=9)
    return ax
