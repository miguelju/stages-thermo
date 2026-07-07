"""stages — a staged-separation (distillation) column solver.

High-level Python interface to the stages-thermo engine (Rust core via PyO3).
Walks the pedagogical ladder of column methods — McCabe–Thiele,
Ponchon–Savarit, the Fenske–Underwood–Gilliland shortcut, and the rigorous
MESH solvers — and exposes a granular, batch-capable API. All thermodynamics
come from ``vle-thermo``; this package adds none of its own.

Status
------
Milestone 1 (column model + McCabe–Thiele). The binary layer is live:

- :class:`ThermoSystem` — Peng–Robinson (φ-φ) and van Laar (γ-φ) systems
  over vle-thermo's component database.
- :class:`EquilibriumCurve` — the y*(x) curve from thermodynamics, from a
  constant relative volatility, or from raw data points.
- :func:`mccabe_thiele`, :func:`rmin`, :func:`total_reflux`, :func:`n_vs_r`
  — the full McCabe–Thiele construction with pinch/R_min detection.
- :class:`BinaryColumn`, :class:`Feed` — the binary column model.
- :mod:`stages.plotting` — the staircase diagram (needs matplotlib; install
  the ``stages-thermo[plot]`` extra).

Units: temperature K, pressure kPa (absolute), compositions are mole
fractions of the *light* (more volatile) component, which is listed first.

>>> import stages
>>> sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
>>> curve = stages.EquilibriumCurve.from_thermo(sys, pressure=101.325)
>>> result = stages.mccabe_thiele(
...     curve, x_distillate=0.95, x_bottoms=0.05, z_feed=0.5, reflux=2.0
... )
>>> result.n_stages          # doctest: +SKIP
11.5...

Column methods land milestone by milestone — see the repo's ``ROADMAP.md``.
"""

from __future__ import annotations

# The native extension. Importing it exercises the Rust shared object; if the
# wheel was built without the `python` feature or the abi3 target mismatched,
# this fails loudly at import time.
from stages._engine import (
    BinaryColumn,
    CondenserKind,
    EquilibriumCurve,
    Feed,
    Line,
    McCabeThieleResult,
    McCabeThieleSpec,
    RminResult,
    StagePoint,
    ThermoSystem,
    TotalRefluxResult,
    mccabe_thiele,
    n_vs_r,
    rmin,
    smoke_bubble_temperature,
    total_reflux,
    version,
)

try:  # Prefer the installed distribution's version; fall back for source trees.
    from importlib.metadata import PackageNotFoundError
    from importlib.metadata import version as _dist_version

    __version__ = _dist_version("stages-thermo")
except (ImportError, PackageNotFoundError):  # pragma: no cover
    __version__ = version()

__all__ = [
    "BinaryColumn",
    "CondenserKind",
    "EquilibriumCurve",
    "Feed",
    "Line",
    "McCabeThieleResult",
    "McCabeThieleSpec",
    "RminResult",
    "StagePoint",
    "ThermoSystem",
    "TotalRefluxResult",
    "__version__",
    "mccabe_thiele",
    "n_vs_r",
    "rmin",
    "smoke_bubble_temperature",
    "total_reflux",
    "version",
]
