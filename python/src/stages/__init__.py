"""stages — a staged-separation (distillation) column solver.

High-level Python interface to the stages-thermo engine (Rust core via PyO3).
Walks the pedagogical ladder of column methods — McCabe–Thiele,
Ponchon–Savarit, the Fenske–Underwood–Gilliland shortcut, and the rigorous
MESH solvers — and exposes a granular, batch-capable API. All thermodynamics
come from ``vle-thermo``; this package adds none of its own.

Status
------
Milestone 0 (repo bootstrap). The only surface today is the vle-thermo smoke
path, which proves the PyO3 boundary end-to-end. Column methods land milestone
by milestone — see the repo's ``ROADMAP.md``.

>>> import stages
>>> stages.__version__          # doctest: +SKIP
'0.0.1'
>>> t = stages.smoke_bubble_temperature()   # methanol/water bubble T [K] via vle-thermo
>>> 280.0 < t < 400.0           # doctest: +SKIP
True
"""

from __future__ import annotations

# The native extension. Importing it exercises the Rust shared object; if the
# wheel was built without the `python` feature or the abi3 target mismatched,
# this fails loudly at import time.
from stages._engine import smoke_bubble_temperature, version

try:  # Prefer the installed distribution's version; fall back for source trees.
    from importlib.metadata import PackageNotFoundError
    from importlib.metadata import version as _dist_version

    __version__ = _dist_version("stages-thermo")
except (ImportError, PackageNotFoundError):  # pragma: no cover
    __version__ = version()

__all__ = ["__version__", "smoke_bubble_temperature", "version"]
