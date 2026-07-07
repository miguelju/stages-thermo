# stages-thermo

**A staged-separation (distillation) learning library and fast steady-state
column solver, built on [`vle-thermo`](https://pypi.org/project/vle-thermo/).**

`import stages` — walk the full pedagogical ladder of column methods
(McCabe–Thiele → Ponchon–Savarit → Fenske–Underwood–Gilliland shortcut →
rigorous MESH), each implemented from scratch and anchored to its textbook
equations, with a granular, batch-capable API ("numpy for distillation
columns"). All thermodynamics come from `vle-thermo`; this package adds none of
its own.

```sh
pip install stages-thermo           # import name is `stages`
pip install "stages-thermo[plot]"   # + matplotlib for the staircase diagrams
```

> **Status:** `0.1.x` ships the first rung of the ladder — the binary
> McCabe–Thiele layer (Milestone 1): equilibrium curves from real
> thermodynamics, minimum reflux by geometric pinch detection (tangent
> pinches included), stage stepping with Murphree efficiency, total reflux,
> N(R), staircase plotting, and the binary column material balances.
> Multicomponent and rigorous MESH solvers land milestone by milestone — see
> the repo's `ROADMAP.md`. The API may still move before `1.0`.

```python
import stages

# Real thermodynamics from vle-thermo: Peng–Robinson benzene–toluene at 1 atm
# (light component first; units K / kPa absolute / mole fractions).
sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
curve = stages.EquilibriumCurve.from_thermo(sys, 101.325)

# Minimum reflux by pinch detection — tangent pinches included …
r = stages.rmin(curve, x_distillate=0.95, x_bottoms=0.05, z_feed=0.50, q=1.0)

# … then the full construction. Rich result objects, never bare numbers:
design = stages.mccabe_thiele(curve, 0.95, 0.05, 0.50, reflux=1.5 * r.r_min)
print(f"N = {design.n_stages:.2f} stages, feed stage {design.feed_stage}, "
      f"R_min = {r.r_min:.3f}")

design.stages          # every (x, y) stage corner
design.staircase       # the full polyline, ready to plot
design.rectifying      # operating lines as slope/intercept
r.pinch, r.tangent     # where the pinch sits, and whether it's a tangent pinch

# The staircase diagram (requires the [plot] extra):
# from stages import plotting
# plotting.plot_mccabe_thiele(design, curve, show_rmin=True)
```

Strongly non-ideal systems go through the γ-φ route — the same construction,
different thermodynamics:

```python
# van Laar methanol–water (the system validated in vle's Chapter IV notebooks).
mw = stages.ThermoSystem.van_laar(["methanol", "water"], 0.5853, 0.3458)
curve_mw = stages.EquilibriumCurve.from_thermo(mw, 101.325)
```

The executable learning path lives in the repo's `notebooks/` —
`01-mccabe-thiele.ipynb` designs benzene–toluene and methanol–water columns
end-to-end, with exercises.

The native core is a Rust crate (`stages-thermo` on crates.io) with PyO3
bindings; wheels are abi3 (`cp310-abi3-*`), so one wheel per (OS, arch) covers
CPython 3.10+.

## License

MIT © Miguel Roberto Jackson Ugueto
