# stages-thermo — Staged-Separation Learning Library & Column Solver

A staged-separation (distillation) learning library and fast steady-state column
solver, built on [`vle-thermo`](https://pypi.org/project/vle-thermo/). Rust core,
Python bindings, executable notebooks.

> **Early-stage.** Only the Milestone 0 bootstrap scaffolding exists today. The
> methods below describe the planned library — see **[Status](#status)** and
> [ROADMAP.md](ROADMAP.md).

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

> **Name-holding stubs only (pre-1.0).** The `0.0.1` releases below just reserve
> the names on both registries while the library is built. Nothing is usable yet
> — track progress in [ROADMAP.md](ROADMAP.md).

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
depends on `vle-thermo` (crates.io / PyPI, v0.8.1) for all thermodynamics.

## Quickstart (sketch — the planned API)

The API is designed so **every method returns a rich, inspectable result
object**, never bare numbers:

```python
import stages

# Rung 1 — McCabe–Thiele on a binary, thermo from vle-thermo.
mt = stages.mccabe_thiele(["benzene", "toluene"], x_feed=0.5, reflux=1.5, eos="PR")
print(mt.n_stages, mt.r_min, mt.pinch)          # result carries the curve, op lines, q-line, every (x,y) corner

# Rung 6 — a rigorous multicomponent column (flagship Naphtali–Sandholm solver).
col = stages.Column.debutanizer()                # multi-feed, side draws, duties, pressure profile
sol = col.solve(spec={"reflux": 2.5, "distillate": 40.0})
print(sol.stage(7).temperature)                  # everything about stage 7
print(sol.mass_balance_closure())                # audit residual
print(sol.report.iterations, sol.report.residual_history)   # open the black box

# The batch layer — "numpy for columns": one FFI crossing, GIL released, parallel.
import numpy as np
sweep = stages.solve_batch(col, reflux=np.linspace(1.2, 4.0, 500))   # N-vs-R curve in one cell
```

*(Illustrative — no code ships yet. See [ROADMAP.md](ROADMAP.md) for what exists.)*

## Status

**Milestone 0 (repo bootstrap) is in progress; everything else is pending.** No
solver code exists yet — only the scaffolding, CI, and these planning documents.

| | State |
|---|---|
| M0 — Repo bootstrap | **In progress** (scaffold, CI, name-holding `0.0.1` stubs, docs split) |
| M1–M10 — The ladder (McCabe–Thiele → 1.0) | Pending |
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
- **pressure-tests vle-thermo's public API** — gaps it surfaces (∂K/∂T, ∂H/∂T,
  packaged γ-φ enthalpy, an expanded Rust-side component DB) feed back as a
  dedicated vle-thermo derivative release (0.9.x), which stages-thermo's
  Milestone 4 tracks.

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
