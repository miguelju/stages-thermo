# stages-thermo — Staged-Separation Learning Library & Column Solver

A staged-separation (distillation) learning library and fast steady-state column
solver, built on [`vle-thermo`](https://pypi.org/project/vle-thermo/). Rust core,
Python bindings, executable notebooks.

> **Early-stage, first two rungs live.** v0.2.0 ships the McCabe–Thiele
> (Milestone 1) and Ponchon–Savarit (Milestone 2) binary layers: real-thermo
> equilibrium and enthalpy–composition (H–x–y) curves, the full graphical
> constructions (tangent-pinch R_min; difference-point energy balances), the
> NRTL γ-φ model, and the executable `01-mccabe-thiele.ipynb` /
> `02-ponchon-savarit.ipynb`. Later rungs are still planned — see
> **[Status](#status)** and [ROADMAP.md](ROADMAP.md).

## About This Project

Two deliberately intertwined goals:

1. **A learning repo that opens the black box.** Chemical-engineering school
   taught *that* columns are solved by "inside-out seeded by Fenske–Underwood–
   Gilliland" but never *how*; professional life delegated it to a process
   simulator. This repo walks the full pedagogical ladder — McCabe–Thiele →
   Ponchon–Savarit → FUG shortcut → rigorous MESH — with **every method
   implemented from scratch**, documented against the textbook equations, and
   demonstrated in executable notebooks. The rigorous-solver notebooks show the
   machinery (residuals, Jacobian structure, convergence history), not just the
   answer.
2. **A fast, granular, production-quality column engine** — "numpy for
   distillation columns." Rust core, Python bindings, a batch API, and every
   intermediate quantity queryable (so an MCP server can later answer "what's the
   temperature on stage 7?" or "how many stages at 1.3·Rmin?").

All thermodynamics (K-values, enthalpies, derivatives) come **entirely from
`vle-thermo`**. This project adds no thermo of its own — it is the first
downstream consumer of that crate, which pressure-tests its public API.

## The Pedagogical Ladder

The ladder **is** the implementation order — successive rungs share kernels (the
tridiagonal solve, bubble-point-per-stage), so it is also a sane build sequence.

| Rung | Method | System | What it teaches |
|---|---|---|---|
| 1 | **McCabe–Thiele** | binary | The equilibrium stage, operating lines, q-line, reflux, pinch, stage-stepping |
| 2 | **Ponchon–Savarit** | binary | Energy balances enter; H–x–y diagram; why CMO is an approximation |
| 3 | **FUG shortcut** (Fenske, Underwood, Gilliland, Kirkbride, Winn) | multicomponent | Nmin, Rmin, N(R), feed stage — and how simulators seed rigorous runs |
| 4 | **Wang–Henke bubble-point** | multicomponent | The MESH equations, tridiagonal structure, tearing, why it fails wide-boiling |
| 5 | **Burningham–Otto sum-rates** | absorbers | The flipped tearing for wide-boiling systems |
| 6 | **Naphtali–Sandholm simultaneous correction** | multicomponent, nonideal | Full Newton on MESH — **the flagship solver** |
| 7 | Robustness layer (homotopy, pseudo-transient) | hard columns | What simulators do when Newton fails |
| 8 | **Inside-out** (Boston/Russell) *(stretch)* | large columns | The method school mentioned; why it was the commercial default |

Rungs 1–3 are the "learning repo" ask; rungs 4–8 are the "solution engine" ask.

## Install

> **Pre-1.0.** v0.2.0 ships the McCabe–Thiele and Ponchon–Savarit binary
> layers. The API may still move before 1.0 — track progress in
> [ROADMAP.md](ROADMAP.md).

### Python (PyPI)

```sh
pip install stages-thermo
```

Distribution name is `stages-thermo`, import name is `stages` (like `vle-thermo`
→ `vle`):

```python
import stages
```

### Rust (crates.io)

```sh
cargo add stages-thermo
```

Both track the same version and build from the same source tree. `stages-thermo`
depends on `vle-thermo` (crates.io / PyPI, ≥ 0.11) for all thermodynamics.

## Quickstart — rungs 1–2, shipping today

**Every method returns a rich, inspectable result object**, never bare numbers:

```python
import stages

# Real thermodynamics from vle-thermo: Peng–Robinson benzene–toluene at 1 atm.
sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
curve = stages.EquilibriumCurve.from_thermo(sys, pressure=101.325)

# Minimum reflux by geometric pinch detection (tangent pinches included) …
r = stages.rmin(curve, x_distillate=0.95, x_bottoms=0.05, z_feed=0.50, q=1.0)

# … and the full construction: stages, feed stage, operating lines, staircase.
design = stages.mccabe_thiele(curve, 0.95, 0.05, 0.50, reflux=1.5 * r.r_min)
print(f"N = {design.n_stages:.1f}, feed stage {design.feed_stage}, R_min = {r.r_min:.3f}")

from stages import plotting
plotting.plot_mccabe_thiele(design, curve, show_rmin=True)   # the classic diagram
```

The same construction runs unchanged on a γ-φ activity-model curve
(`stages.ThermoSystem.van_laar(["methanol", "water"], 0.5853, 0.3458)`), a
constant-α idealization (`EquilibriumCurve.constant_alpha(2.5)`), or digitized
literature data (`EquilibriumCurve.from_points(...)`) — the solvers never know
which thermodynamics produced the curve. Later rungs (rigorous multicomponent
columns, the batch `solve_batch` layer) follow the same rich-result design; see
[PLAN.md](PLAN.md) §8.

## Status

**Milestones 0, 1, 2, and 4 are complete** (v0.2.0): the McCabe–Thiele and
Ponchon–Savarit binary layers are live on top of vle-thermo 0.11.x, with the
executable [`notebooks/01-mccabe-thiele.ipynb`](notebooks/01-mccabe-thiele.ipynb)
and [`notebooks/02-ponchon-savarit.ipynb`](notebooks/02-ponchon-savarit.ipynb).

| | State |
|---|---|
| M0 — Repo bootstrap | **Complete** (scaffold, CI, `0.0.1` stubs, docs split) |
| M1 — Column model + McCabe–Thiele | **Complete — v0.1.0** (equilibrium curves, R_min/pinch, stage stepping, Murphree, N(R), staircase plots, 📓 01) |
| M2 — Ponchon–Savarit | **Complete — v0.2.0** (H–x–y enthalpy curves, NRTL γ-φ + per-phase enthalpy adapter, difference-point construction, duties, 📓 02) |
| M4 — Upstream vle-thermo derivative release | **Complete** (vle-thermo v0.11: NRTL + ammonia, `k_values_with_derivs`, γ-φ enthalpy, 25-compound Rust DB) |
| M3, M5–M10 — The remaining ladder → 1.0 | Pending |
| M11–M12 — Inside-out, MCP server *(stretch)* | Pending |

Per-milestone detail and hour estimates: [ROADMAP.md](ROADMAP.md) and
[TODO.md](TODO.md). The full technical spec (method selection, thermo adapter,
granular API, validation strategy, references) is [PLAN.md](PLAN.md).

## Relationship to vle-thermo

stages-thermo is the **first downstream consumer** of
[`vle-thermo`](https://pypi.org/project/vle-thermo/) — a modern Rust + Python VLE
engine (22+ EOS, 5 activity models, 11 mixing rules, the modern flash suite, all
with exact composition derivatives). stages-thermo:

- gets **all** K-values, enthalpies, and derivatives from vle-thermo, funneled
  through a single adapter module (see [CLAUDE.md](CLAUDE.md) §"The thermo
  adapter");
- reuses vle's entire proven toolchain (workspace layout, maturin build,
  cibuildwheel CI, `_batch` + rayon + GIL-release pattern, criterion benches);
- **pressure-tests vle-thermo's public API** — the gaps it surfaced (∂K/∂T,
  ∂H/∂T, packaged γ-φ enthalpy, an expanded Rust-side component DB) fed back
  as vle-thermo's derivative release (v0.9.x, stages-thermo's Milestone 4), and
  the NRTL activity model + ammonia component (v0.11, needed for M2's
  ammonia–water showcase); the engine now builds on `vle-thermo = "0.11"` with
  the `component-db` feature.

There is **no Rust staged-column / MESH library on crates.io** — this is
greenfield.

## Documentation

- [PLAN.md](PLAN.md) — the source-of-truth technical plan (vision, method
  selection, repo structure, thermo adapter, granular API, validation,
  references)
- [ROADMAP.md](ROADMAP.md) — milestones M0–M12 with hour estimates
- [TODO.md](TODO.md) — actionable task lists per milestone
- [CLAUDE.md](CLAUDE.md) — conventions (units, citations, GPL caution, PyO3 rule,
  notebook conventions, signing flow)
- `docs/theory/` + `docs/references.md` — the equations as implemented, cited
  (ACS style)

## License

MIT — see [LICENSE](LICENSE). Built with [Claude Code](https://claude.ai/code) as
a development partner; each milestone records the exact Claude model that executed
it in [ROADMAP.md](ROADMAP.md) and the commit trailers.
