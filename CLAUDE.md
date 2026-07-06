# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository.

**stages-thermo** is a staged-separation (distillation) learning library and
fast steady-state column solver, built on **`vle-thermo`** (crates.io / PyPI,
MIT, currently v0.8.1). It adds no thermodynamics of its own — it is the first
downstream consumer of vle-thermo. Conventions here are carried over from the
vle repo with names adjusted (`vle` → `stages`, import name `stages`, crate/PyPI
name `stages-thermo`).

## Project Overview

Two intertwined goals (PLAN.md §1):

1. **A learning repo that opens the black box** — walk the full pedagogical
   ladder (McCabe–Thiele → Ponchon–Savarit → FUG shortcut → rigorous MESH), every
   method implemented from scratch, documented against textbook equations, and
   demonstrated in executable notebooks. The rigorous-solver notebooks must show
   the machinery (residuals, Jacobian structure, convergence history), not just
   the answer.
2. **A fast, granular, production-quality column engine** — "numpy for
   distillation columns." Rust core, Python bindings, batch API, every
   intermediate quantity queryable.

The milestone order **is** the pedagogical ladder — see [ROADMAP.md](ROADMAP.md),
[TODO.md](TODO.md), and the source-of-truth [PLAN.md](PLAN.md).

## Target Architecture

```
engine/     — Rust crate `stages-thermo` (lib `stages_thermo`), PyO3 bindings via maturin
python/     — Python package `stages-thermo`, import name `stages` (high-level API, plotting)
notebooks/  — Jupyter notebooks, one per milestone (the learning path)
docs/       — theory/ (equations as implemented, cited) + references.md
```

**Build chain:** Rust (`engine/`) → PyO3/maturin → Python native module → Python
wrapper (`python/`) → Jupyter notebooks. The engine depends on `vle-thermo` from
crates.io (pinned minor version); the Python package depends on `vle-thermo` from
PyPI. Full repo layout: [PLAN.md](PLAN.md) §6.

## The thermo adapter — exactly one module talks to vle-thermo

**Design rule (PLAN.md §7): exactly one module — `engine/src/thermo.rs` — imports
`vle-thermo`. Everything else consumes a `ThermoProvider`-style interface.** Use a
concrete struct first; make it a trait **only** if/when a surrogate model (the
inside-out milestone) or a mock needs to slot in — do not abstract prematurely.

- The adapter is where every K-value, enthalpy, and derivative request to
  vle-thermo is funneled. If a needed quantity is missing upstream (e.g. ∂K/∂T,
  ∂H/∂T, γ-φ packaged enthalpy in vle-thermo 0.8.1), the adapter fills the gap
  with a **central finite-difference interim** until the upstream analytic
  version lands (M4 / vle-thermo 0.9.x). The FD path doubles as a correctness
  oracle for the analytic Jacobian later.
- **Reference-state discipline**: pick one convention for the whole column
  (vle-thermo's `t_ref`/`p_ref`/`h_ref` machinery), set it once in the adapter,
  never per-stage.

## Canonical Units

**Identical to vle** — non-negotiable, since vle-thermo's public API speaks these:

- Temperature: **K** (absolute)
- Pressure: **kPa** (absolute — **never** gauge)
- Energy (molar): **kJ/kmol**
- Entropy (molar): **kJ/(kmol·K)**
- Amount: **kmol**
- Gas constant R: **8.31451 kJ/(kmol·K)**

**Every Rust function and Python wrapper function** that accepts or returns a
physical quantity MUST state the units in its doc comment.

```rust
/// Per-stage bubble temperature at the given stage pressure.
///
/// # Arguments
/// * `pressure` — Stage pressure in **kPa** (absolute)
///
/// # Returns
/// Bubble temperature in **K**
pub fn stage_bubble_temperature(pressure: f64, /* ... */) -> f64 { ... }
```

User-facing Python accepts unit strings via `vle.units` (e.g. `Q_(300, "K")`,
`"1 atm"`); engine-internal functions always take absolute kPa / K.

## Educational Comments + Textbook Anchoring

Two layers, both required (PLAN.md §12):

1. **Beginner-friendly Rust/Python idiom explanations** — Miguel's standing
   preference. Explain the language mechanics where they aren't obvious.
2. **Textbook anchoring** — every solver module's doc comment states the method,
   its reference in ACS form, and **maps code symbols to the book's symbols**,
   e.g.:

   ```rust
   // Wang–Henke, S&H eq. 10-25: A_j = L_{j-1}  (sub-diagonal of the tridiagonal M-balance)
   // B_j = -(L_j + V_j K_{ij}) ; C_j = V_{j+1} K_{i,j+1}
   ```

Never paraphrase equations in docs/notebooks — render them with LaTeX and cite
the source.

## Reference Citation Requirements

Code implementing a published method **cites it** in the module-level doc comment
in ACS form (the analogue of vle's Pascal-source rule). The full list lives in
[`docs/references.md`](docs/references.md) (ACS style, PLAN.md §14) with a
reference-to-code mapping. Method → paper:

- **Wang–Henke bubble-point** (`rigorous/bubble_point.rs`) — Wang & Henke, *Hydrocarbon Process.* **1966**, 45(8), 155
- **Burningham–Otto sum-rates** (`rigorous/sum_rates.rs`) — Burningham & Otto, *Hydrocarbon Process.* **1967**, 46(10), 163
- **BP↔narrow / SR↔wide pairing** (docs) — Friday & Smith, *AIChE J.* **1964**, 10, 698
- **Naphtali–Sandholm simultaneous correction** (`rigorous/naphtali_sandholm.rs`, `jacobian.rs`) — Naphtali & Sandholm, *AIChE J.* **1971**, 17, 148
- **Tridiagonal / tearing** — Tomich, *AIChE J.* **1970**, 16, 229
- **Inside-out + tridiagonal stabilization** (`rigorous/inside_out.rs`, `numerics/thomas.rs`) — Boston & Sullivan, *Can. J. Chem. Eng.* **1974**, 52, 52; Boston & Britt, *Comput. Chem. Eng.* **1978**, 2, 109; Boston, ACS Symp. Ser. 124, **1980**; Russell, *Chem. Eng.* **1983** (Oct 17), 53
- **Nonsmooth inside-out inner loop** — Watson, Vikse, Gundersen & Barton, *IECR* **2017**, 56, 960
- **Thermodynamic homotopy** (`rigorous/homotopy.rs`) — Vickery & Taylor, *AIChE J.* **1986**, 32, 547; Wayburn & Seader, *Comput. Chem. Eng.* **1987**, 11, 7
- **Pseudo-transient relaxation + Newton switchover** (`rigorous/relaxation.rs`) — Pattison & Baldea, *AIChE J.* **2014**, 60, 4104; Ketchum, *Chem. Eng. Sci.* **1979**, 34, 387
- **Rate-based reference (non-goal, cite in docs)** — Krishnamurthy & Taylor, *AIChE J.* **1985**, 31, 449/456
- **Gilliland closed form** (`shortcut/gilliland.rs`) — Molokanov et al. **1972**
- **Textbook set** — Seader, Henley & Roper, *Separation Process Principles* (Ch. 10 = the implementation-level bible); Holland, *Fundamentals of Multicomponent Distillation*, 1981; Kister, *Distillation Design*, 1992; Doherty & Malone, *Conceptual Design of Distillation Systems*, 2001; King, *Separation Processes*; Górak & Sorensen (eds.), *Distillation: Fundamentals and Principles*, 2014; Taylor & Kooijman, *The ChemSep Book* (chemsep.org)

### GPL caution (load-bearing)

**DWSIM is GPLv3.** Its `RigorousColumnSolvers/` are textbook-annotated
Wang–Henke / sum-rates / Naphtali–Sandholm — **consult them as literature, never
translate or port their code** into this MIT project. **BioSTEAM** (`stage.py`)
is permissively licensed (NCSA) — its inside-out and simultaneous-correction
implementations may be **mimicked freely**. ChemSep-LITE (free binary) and DWSIM
are the cross-validation oracles (run the same specs, compare numbers), which is
license-clean regardless.

## Method / Algorithm Choices

The rationale behind the solver design (PLAN.md §5) — apply these when
implementing:

- **Primary rigorous solver = Naphtali–Sandholm simultaneous correction**:
  damped full Newton on the stage-grouped MESH system, unknowns `(v_1j..v_Cj,
  T_j, l_1j..l_Cj)` per stage, **exact analytic block-tridiagonal Jacobian**,
  block-Thomas elimination (LU-factorize the pivot blocks, **never invert**),
  Armijo line search + physical safeguards (clamp ΔT to ~10 K/step, work in
  ln-component-flows for positivity). Chosen over inside-out because vle-thermo's
  property calls are cheap Rust and its derivatives are analytic and nearly free
  — an exact cheap Jacobian makes full Newton the simplest correct code path.
- **Seeding is the whole game** (Newton's weakness is convergence radius): FUG
  shortcut → N, feed stage, R → linear T profile → constant-molal-overflow flows
  → compositions from a feed flash or 2–3 Wang–Henke sweeps. The seeding pipeline
  is itself pedagogy.
- **Fallback ladder when Newton stalls** (M9): thermodynamic homotopy
  (ideal-Raoult → rigorous λ-continuation) first, then pseudo-transient
  relaxation with a switch back to Newton once residuals shrink.
- **Tridiagonal kernels**: scalar Thomas (`numerics/thomas.rs`) with
  Boston–Sullivan-style near-cancellation stabilization; block-tridiagonal LU
  (`numerics/block_thomas.rs`) using nalgebra blocks, no inversion.
- **Underwood roots**: interval-bracketed (Brent per adjacent-volatility
  interval); budget test time for distributed components (companion-matrix
  eigenvalue formulation is the robust fallback).
- **Derivatives**: exact/analytic once M4 lands; FD only as the interim adapter
  path **and** as a permanent test oracle for the analytic Jacobian (the same
  pattern vle uses).
- **Hot-path rules** (from M6/M10 on): no heap allocation inside iteration loops;
  criterion benches guard regressions.

## PyO3 Bindings Rule

**Every milestone that adds Rust functionality must expose the new public
functions or types as PyO3 bindings in the same commit series.** Pure-Rust-
without-Python is not acceptable — Python is a first-class consumer.

- New `#[pyfunction]`s go in `engine/src/py_bindings.rs` (+ `py_column.rs`),
  feature-gated behind `python` (same shape as vle).
- New public types get `#[cfg_attr(feature = "python", pyo3::pyclass(...))]`.
- Add at least one `python/tests/` test per binding, exercising it through the
  built wheel (CI runs pytest against the cibuildwheel wheel, so a missing
  binding is a hard failure).

## Notebook Conventions

Every milestone-level notebook (`notebooks/NN-<name>.ipynb`) MUST follow a
professional structure so the collection reads as a coherent learning path.

**Required sections (top to bottom):**

1. **Title + one-sentence motivation** (H1 + lead paragraph).
2. **Setup (optional)** — a markdown cell followed by a code cell containing
   exactly a commented `# %pip install --upgrade stages-thermo` (and
   `vle-thermo`), left **commented out** so the notebook executes top-to-bottom
   without it.
3. **Context** — quote/paraphrase the relevant textbook method or `docs/theory/`
   page, with a relative link back. Blockquote direct quotes; render equations
   with LaTeX (never paraphrase them).
4. **What was built in this milestone** — short prose pointing at the modules /
   structs the reader will call.
5. **Worked example** — one fully-executed example end-to-end, matching a
   literature table where possible.
6. **User exercises — at least 2** — each with a problem statement, a `# TODO:`
   template cell, and a hidden/collapsed solution (`<details>` block or a
   bottom "Solutions" section).
7. **References** — cross-links to `docs/theory/`, `docs/references.md`, and the
   PLAN section describing the algorithm.

**Other requirements:**

- All cells execute top-to-bottom in a fresh kernel. Verify with
  `~/miniconda3/envs/stages/bin/jupyter nbconvert --to notebook --execute`
  before committing.
- `import matplotlib.pyplot as plt` + inline `%matplotlib inline`.
- Import units as `from vle.units import ureg, Q_`; express inputs with explicit
  units, e.g. `T = Q_(300, "K")`.
- **Pin numeric expectations** (literature table values, mass-balance closure) in
  assertion cells so regressions surface as a failing notebook, not silent drift.
- The rigorous-solver notebooks must show the machinery (residuals, Jacobian
  sparsity, convergence history) — that's the "open the black box" requirement.

## Python Environment (conda `stages` env — mandatory)

**All Python work in this repo goes through the dedicated `stages` conda env.
Never invoke a bare `python`/`python3`/`pip`/`pytest`, and never create a `.venv`
in the repo.** Use the env's binaries directly by absolute path (more reliable
from non-interactive shells than `conda activate`):

- `~/miniconda3/envs/stages/bin/python` — running any Python script or one-liner
- `~/miniconda3/envs/stages/bin/pytest` — running `python/tests/`
- `~/miniconda3/envs/stages/bin/maturin` — building/installing the PyO3 wheel (`maturin develop` from `python/`)
- `~/miniconda3/envs/stages/bin/jupyter` — executing notebooks (`nbconvert --execute`)

If the env is missing, recreate it with conda (`conda create -n stages
python=3.12` + `pip install maturin pytest jupyter matplotlib numpy vle-thermo`)
rather than falling back to the system Python. The env installs **`vle-thermo`
from PyPI** — the notebooks use both (`vle.System` for thermo exploration,
`stages` for columns). Note: `conda env list` from a non-interactive shell
sometimes misses envs — check `ls ~/miniconda3/envs/` directly.

## Milestone Reiteration Before Execution

**Before executing a milestone, re-read its ROADMAP/TODO/PLAN entries and restate
the plan** (deliverables, module list, validation cases, the notebook) to Miguel
for review before writing code. Each milestone follows a plan-then-execute cycle:
plan → review → execute (code + tests + docs) → validate against literature/
cross-simulator cases → commit. Do not start coding a milestone from memory of
the plan; reiterate it first.

## Milestone Tracking Rules (model-attribution lines)

**When completing a milestone**, you MUST record the exact LLM model that
executed it in **three** places:

1. **ROADMAP.md** — add `*Executed by Claude Code using <model name and
   version>*` under the milestone header.
2. **TODO.md** — the same line under the milestone section header.
3. **Git commit message** — a `Co-Authored-By` trailer naming the model.

The model name must be the exact model powering the session (e.g. `Claude Opus
4.8 (1M context)`, `Claude Fable 5`). M0 was executed by **Claude Opus 4.8 (1M
context)**; future milestones are attributed as they land.

## Release & Push Rules

### Doc-sync-before-push

**Before every `git push` or release**, review and update all documentation to
reflect current state:

1. **README.md** — feature/ladder status, structure, status line
2. **ROADMAP.md** — check off completed milestones, model-attribution line
3. **TODO.md** — check off tasks, update the summary table + estimates
4. **PLAN.md** — update only if architecture / method selection / milestone scope
   actually changed (it is the source-of-truth spec; keep it authoritative)
5. **CLAUDE.md** — update if new conventions, paths, or tools were introduced
6. **`python/README.md`** (the **PyPI** long-description for `stages-thermo`) —
   update whenever the change affects the Python-facing story (new/removed public
   API, status/version language, install steps). Every snippet must run verbatim
   against the current wheel — execute it (`~/miniconda3/envs/stages/bin/python`)
   before committing.
7. **`engine/README.md`** (the **crates.io** page for `stages-thermo`, via
   `readme = "README.md"` in `engine/Cargo.toml`) — update whenever the change
   affects the Rust-crate story. Any `rust` code block must compile — verify it
   before committing.
8. **`docs/theory/` + `docs/references.md`** — update when a method's
   documentation or citation changed.

**Package-page docs are immutable per published version.** PyPI and crates.io
render the README **bundled with each release**; editing `python/README.md` /
`engine/README.md` does **not** refresh the live page — only publishing a new
version does (bump `[workspace.package] version` in the root `Cargo.toml` +
`version` in `python/pyproject.toml`, tag `v<x.y.z>` → `release.yml`). Batch
doc-only README fixes into the next release rather than tagging solely for a typo,
unless the staleness is materially misleading. The two package READMEs, the root
`README.md`, and both `description` fields must tell a mutually consistent
version/status story.

Do NOT push until all documentation accurately reflects the current state.

### fmt-before-push (`hooks/pre-push`)

The repo ships a versioned pre-push hook that runs `cargo fmt --check` and
**blocks the push on any diff** — mirroring the first step of CI's `lint` job.
Activate it once per clone with `git config core.hooksPath hooks` (a local
setting, not committed). **Always run `cargo fmt --check` before pushing** even if
the hook is active — never `--no-verify`. (clippy runs in CI only, matching vle,
if the local pyo3/rustc mismatch reproduces here.)

This repo carries **no private infrastructure** (distribution = crates.io + PyPI +
notebooks only), so the private-data gate is a placeholder — keep the rule but it
should never fire. Miguel's public professional addresses on `migueljackson.dev`
are safe to commit.

### YubiKey signing flow (Miguel taps; Claude never commits)

Claude cannot sign — every commit needs Miguel at the keyboard to tap the
YubiKey. The pattern:

1. **Claude writes the commit message to `/tmp/<name>-commit.txt`** via `Write`
   (inline heredocs hit the zsh `EOF`-must-be-flush-left trap; the file sidesteps
   it). Single-line trivial messages can use `git commit -m`.
2. **Miguel runs `git commit -F /tmp/<name>-commit.txt`** and taps the YubiKey
   (single tap, ~3s window).
3. **Claude runs `git push`** — needs no touch (the signature is on the commit
   object).
4. **Claude verifies:** `git log -1 --format='%h %s%nVerified: %G?  Signer:
   %GS'` → expect `Verified: G`.

**Don't background `git push`** waiting for a tap — the Bash tool has no presence
channel and it'll hang. **Never `--no-verify` / `--no-gpg-sign`** to bypass a
stuck signing step (the server-side `require-signed-commits` rejects unsigned
anyway); fix the underlying issue. All commit `Co-Authored-By` / author fields use
Miguel's public identity.

## Domain Context

This is a chemical-engineering / separations codebase. Standard notation: stages
numbered top-down (1 = condenser, N = reboiler), K_ij (K-value of component i on
stage j), x/y (liquid/vapor mole fractions), L/V (liquid/vapor molar flows),
R (reflux ratio), Nmin/Rmin (minimum stages/reflux), α (relative volatility),
Murphree efficiency, MESH (Material / Equilibrium / Summation / entHalpy
equations). Out of scope for v1 (documented non-goals, PLAN.md §4): rate-based/
nonequilibrium stages, reactive distillation, dividing-wall columns, tray
hydraulics — the stage-equation layout should merely not preclude them.

## Validation Strategy (PLAN.md §9)

Four layers: (1) kernel unit tests (`thomas`/`block_thomas` vs dense solves;
Underwood bracketing on constructed volatility sets; Fenske/Gilliland vs hand
calcs); (2) **literature tables** as the primary benchmark, pinned in assertion
cells (Seader/Henley/Roper Ch. 10 worked examples — Wang–Henke per-iteration
tables, Burningham–Otto absorber, Naphtali–Sandholm example, debutanizer;
Holland stage-by-stage tables; nonideal methanol–water / 2-propanol–water from
vle Ch. IV); (3) cross-simulator checks vs DWSIM + ChemSep-LITE (~1% target);
(4) invariant/property tests on every converged solution (mass-balance closure
≤1e-9 relative, ΣxΣy=1, energy closure, y=Kx consistency).
