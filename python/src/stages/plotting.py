"""Diagram plotting for the stages package.

Matplotlib is an optional dependency (install the ``stages-thermo[plot]``
extra); everything else in ``stages`` works without it. The import is
deferred into the functions so ``import stages.plotting`` alone doesn't pull
matplotlib in.

The M1 surface is the McCabe–Thiele staircase diagram; M2 adds the
enthalpy–composition (H–x–y) diagram and the Ponchon–Savarit construction on
it. Later milestones add Gilliland charts and column profile /
convergence-history plots.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover - typing only
    from matplotlib.axes import Axes

    from stages import (
        EnthalpyCurve,
        EquilibriumCurve,
        McCabeThieleResult,
        PonchonSavaritResult,
        TotalRefluxResult,
    )

__all__ = [
    "plot_hxy",
    "plot_mccabe_thiele",
    "plot_ponchon_savarit",
    "plot_total_reflux",
]


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


# --------------------------------------------------------------------------- #
# Ponchon–Savarit — the enthalpy–composition (H–x–y) diagram                  #
# --------------------------------------------------------------------------- #
#
# Unlike the McCabe–Thiele frame, this puts **molar enthalpy on the y-axis** and
# composition on the x-axis, and is NOT a unit square — so no
# ``set_aspect("equal")`` (enthalpy spans tens of MJ/kmol while composition
# spans [0, 1]).


def _draw_hxy_frame(ax: "Axes", curve: "EnthalpyCurve") -> None:
    """The parts every H–x–y diagram shares: saturated-liquid and
    saturated-vapor enthalpy curves."""
    ax.plot(
        curve.x,
        curve.h_liq,
        color="tab:blue",
        lw=2,
        label="saturated liquid $h_L(x)$",
        zorder=3,
    )
    ax.plot(
        curve.y,
        curve.h_vap,
        color="tab:red",
        lw=2,
        label="saturated vapor $H_V(y)$",
        zorder=3,
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("composition — light component (mole fraction)")
    ax.set_ylabel(r"molar enthalpy $H$ [kJ/kmol]")
    ax.grid(True, alpha=0.3)


def plot_hxy(
    curve: "EnthalpyCurve",
    ax: "Axes | None" = None,
    n_tie_lines: int = 0,
) -> "Axes":
    """Draw the bare enthalpy–composition (H–x–y) diagram for a curve.

    Parameters
    ----------
    curve:
        The :class:`stages.EnthalpyCurve` to draw.
    ax:
        Existing axes to draw into; a new figure is created when omitted.
    n_tie_lines:
        If > 0, overlay this many equilibrium tie lines (connecting each
        liquid point to its equilibrium vapor point), evenly spaced in ``x`` —
        a purely illustrative aid to reading the diagram.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _require_pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))
    _draw_hxy_frame(ax, curve)

    if n_tie_lines > 0:
        xs = curve.x
        step = max(1, len(xs) // (n_tie_lines + 1))
        drawn_label = False
        for i in range(step, len(xs) - 1, step):
            xi = xs[i]
            yi = curve.y_of_x(xi)
            ax.plot(
                [xi, yi],
                [curve.h_liquid_of_x(xi), curve.h_vapor_of_y(yi)],
                color="0.6",
                lw=0.8,
                zorder=2,
                label="tie lines" if not drawn_label else None,
            )
            drawn_label = True

    ax.legend(loc="upper left", fontsize=9)
    return ax


def plot_ponchon_savarit(
    result: "PonchonSavaritResult",
    curve: "EnthalpyCurve",
    ax: "Axes | None" = None,
    show_poles: bool = True,
) -> "Axes":
    """Draw the full Ponchon–Savarit construction on the H–x–y diagram.

    Draws the saturated-liquid/vapor curves, the stepped tie lines (equilibrium
    stages) and pole (operating) lines, the feed point, and — when
    ``show_poles`` is set — the two difference points ``Δ_D`` / ``Δ_B`` with the
    pole lines extended to them (which lets matplotlib autoscale the enthalpy
    axis to include the poles above/below the diagram).

    Parameters
    ----------
    result:
        The construction returned by :func:`stages.ponchon_savarit`.
    curve:
        The enthalpy curve the construction was run on.
    ax:
        Existing axes to draw into; a new figure is created when omitted.
    show_poles:
        Also plot ``Δ_D`` and ``Δ_B`` and the pole lines reaching them.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _require_pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))
    _draw_hxy_frame(ax, curve)

    tie = result.tie_lines  # list of ((x_n, h_L), (y_n, H_V))

    # Equilibrium tie lines (one per stage).
    for k, (liq, vap) in enumerate(tie):
        ax.plot(
            [liq[0], vap[0]],
            [liq[1], vap[1]],
            color="tab:green",
            lw=1.0,
            zorder=2,
            label="tie lines (stages)" if k == 0 else None,
        )
        ax.annotate(
            str(k + 1),
            (liq[0], liq[1]),
            textcoords="offset points",
            xytext=(-9, -10),
            fontsize=8,
            color="tab:green",
        )

    # Pole (operating) lines: L_n → V_{n+1} for successive stages.
    for n in range(len(tie) - 1):
        liq_n = tie[n][0]
        vap_next = tie[n + 1][1]
        ax.plot(
            [liq_n[0], vap_next[0]],
            [liq_n[1], vap_next[1]],
            color="tab:purple",
            lw=0.9,
            ls="-.",
            zorder=2,
            label="pole (operating) lines" if n == 0 else None,
        )

    # Feed point on the diagram.
    fx, fh = result.feed_point
    ax.plot([fx], [fh], marker="s", ms=6, color="k", zorder=4)
    ax.annotate(
        "$F$", (fx, fh), textcoords="offset points", xytext=(6, -4), fontsize=10
    )

    if show_poles:
        dd_x, dd_h = result.delta_d
        db_x, db_h = result.delta_b
        # The poles.
        ax.plot([dd_x], [dd_h], marker="^", ms=9, color="tab:red", zorder=5)
        ax.annotate(
            r"$\Delta_D$",
            (dd_x, dd_h),
            textcoords="offset points",
            xytext=(6, 0),
            fontsize=11,
            color="tab:red",
        )
        ax.plot([db_x], [db_h], marker="v", ms=9, color="tab:blue", zorder=5)
        ax.annotate(
            r"$\Delta_B$",
            (db_x, db_h),
            textcoords="offset points",
            xytext=(6, 0),
            fontsize=11,
            color="tab:blue",
        )
        # The Δ_D–F–Δ_B collinear balance line.
        ax.plot(
            [dd_x, db_x],
            [dd_h, db_h],
            color="0.5",
            lw=0.8,
            ls=":",
            zorder=1,
            label=r"$\Delta_D$–$F$–$\Delta_B$ line",
        )
        # Extend the first/last pole lines up to the poles.
        v1 = tie[0][1]
        ax.plot([dd_x, v1[0]], [dd_h, v1[1]], color="tab:red", lw=0.7, ls="--", zorder=1)
        l_last = tie[-1][0]
        ax.plot(
            [db_x, l_last[0]],
            [db_h, l_last[1]],
            color="tab:blue",
            lw=0.7,
            ls="--",
            zorder=1,
        )

    ax.set_title(
        f"Ponchon–Savarit: N = {result.n_stages:.2f} stages, "
        f"feed stage {result.feed_stage}, R = {result.spec.reflux:.3f}\n"
        f"$Q_C/F$ = {result.q_condenser:,.0f}, "
        f"$Q_R/F$ = {result.q_reboiler:,.0f} kJ/kmol feed"
    )
    ax.legend(loc="center left", fontsize=8)
    return ax
