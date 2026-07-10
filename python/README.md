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

> **Status:** `0.2.x` ships the first two rungs of the ladder — the binary
> **McCabe–Thiele** (Milestone 1) and **Ponchon–Savarit** (Milestone 2) layers:
> equilibrium and enthalpy–composition (H–x–y) curves from real thermodynamics,
> minimum reflux by geometric pinch detection (tangent pinches included), stage
> stepping with Murphree efficiency, total reflux, N(R), the energy-exact
> difference-point construction (with condenser/reboiler duties), the NRTL γ-φ
> model, per-phase enthalpies, and the diagram plots. Multicomponent and
> rigorous MESH solvers land milestone by milestone — see the repo's
> `ROADMAP.md`. The API may still move before `1.0`.

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

Rung 2 — **Ponchon–Savarit** — closes the energy balance on the
enthalpy–composition (H–x–y) diagram, so it also returns the condenser and
reboiler duties (which McCabe–Thiele cannot):

```python
# H–x–y curve: saturated-liquid and -vapor enthalpies alongside y*(x).
ec = stages.EnthalpyCurve.from_thermo(sys, 101.325)
ps = stages.ponchon_savarit(ec, x_distillate=0.95, x_bottoms=0.05,
                            z_feed=0.50, reflux=1.5)
print(f"N = {ps.n_stages:.2f} stages, feed stage {ps.feed_stage}")
print(f"Q_C/F = {ps.q_condenser:,.0f}, Q_R/F = {ps.q_reboiler:,.0f} kJ/kmol feed")
ps.delta_d, ps.delta_b   # the two difference points (poles), (x, H) in kJ/kmol

# NRTL for strongly non-ideal aqueous-organic systems (a12/a21 in kJ/kmol):
aw = stages.ThermoSystem.nrtl(["ammonia", "water"], a12=-1800.0, a21=-1200.0, alpha=0.2)
# from stages import plotting; plotting.plot_ponchon_savarit(ps, ec)
```

The executable learning path lives in the repo's `notebooks/` —
`01-mccabe-thiele.ipynb` and `02-ponchon-savarit.ipynb` design benzene–toluene,
methanol–water, and ammonia–water columns end-to-end, with exercises.

The native core is a Rust crate (`stages-thermo` on crates.io) with PyO3
bindings; wheels are abi3 (`cp310-abi3-*`), so one wheel per (OS, arch) covers
CPython 3.10+.

## License

MIT © Miguel Roberto Jackson Ugueto
