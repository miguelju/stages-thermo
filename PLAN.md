# stages-thermo — Project Plan

**A staged-separation (distillation) learning library and fast steady-state column solver, built on `vle-thermo`.**

*Planning document, 2026-07-05. Prepared by Claude Fable 5 for execution by Claude Opus 4.8. This file is the handoff spec: it becomes the new repo's `PLAN.md`, and Milestone 0 splits its milestone/task sections into `ROADMAP.md` and `TODO.md` following the vle conventions.*

---

## 1. Vision

Two goals, deliberately intertwined:

1. **A learning repo that opens the black box.** Miguel's stated motivation: chemical engineering school (25+ years ago) taught *that* columns are solved by "inside-out seeded by Fenske–Underwood–Gilliland" but never *how*; professional life delegated it to a process simulator. This repo walks the full pedagogical ladder — McCabe–Thiele → Ponchon–Savarit → FUG shortcut → rigorous MESH — with every method implemented from scratch, documented against the textbook equations, and demonstrated in executable notebooks. The rigorous-solver notebooks must show the machinery (residuals, Jacobian structure, convergence history), not just the answer.

2. **A fast, granular, production-quality column engine** — "numpy for distillation columns." Rust core, Python bindings, batch API, every intermediate quantity queryable (so an MCP server can later answer questions like "what's the temperature on stage 7?" or "how many stages at 1.3·Rmin?").

The thermodynamics (K-values, enthalpies, derivatives) come entirely from **`vle-thermo`** (crates.io / PyPI, MIT, currently v0.8.1). This project adds no thermo of its own — it is the first downstream consumer of vle-thermo, which will pressure-test that crate's public API.

---

## 2. Name

**Recommendation: `stages-thermo`** (Rust crate `stages-thermo`, PyPI `stages-thermo`, Python import `stages`).

- **Free on both PyPI and crates.io** (verified 2026-07-05 via registry APIs; register 0.0.1 stubs in Milestone 0 to hold both, same as vle-thermo does).
- Sibling branding with `vle-thermo` — the `*-thermo` family does the disambiguation work, and "vle-thermo + stages-thermo" reads as one project family.
- Names the core concept: the theoretical equilibrium stage.
- Import-name note: PyPI package `stages` exists but is a dead `1.0.dev0` placeholder from years ago; the import collision risk for `import stages` is negligible, and matching vle's pattern (`vle-thermo` → `import vle`) is worth it. If Opus prefers zero risk, use `import stages_thermo`.

Alternatives checked (all free on **both** registries unless noted): `etapas` (Spanish "stages" — distinctive, ties to the Spanish-language research heritage; the strongest alternative), `mesh-thermo` (names the mathematics; instantly meaningful to ChemEs), `cascade-thermo`, `destila`, `platos`, `equistage`, `trays`. Rejected: `mccabe` (PyPI = the flake8 complexity plugin on every dev machine), `distill`/`distil` (taken + ML-culture collision), `stages` alone (PyPI squatted), `cascada`/`rectify` (split availability).

---

## 3. Language & stack

**Rust + PyO3/maturin + Python wrapper + Jupyter notebooks — the identical stack to vle.** No real decision to make here: the requirement is "fast, interfaces with Python, numpy for distillation columns," which is exactly what the vle stack already delivers, and reusing it means reusing the entire proven toolchain (workspace layout, maturin build, cibuildwheel CI, `_batch` + rayon + GIL-release pattern, criterion benches, the conda-env discipline).

Notable confirmation from the ecosystem survey: **there is no Rust staged-column/MESH library anywhere** — crates.io has nothing. This is greenfield; stages-thermo would be the first.

- Engine crate depends on `vle-thermo` from **crates.io** (pinned minor version). For local co-development, a `[patch.crates-io] vle-thermo = { path = "../vle/engine" }` entry in a git-ignored `.cargo/config.toml` or a documented patch workflow — but CI always builds against the published crate, which keeps vle-thermo honest about publishing what downstream needs.
- Python package depends on `vle-thermo` from PyPI (the notebooks use both: `vle.System` for thermo exploration, `stages` for columns).
- `nalgebra` for the dense (2C+1)×(2C+1) blocks; `num-dual` (already a vle dependency) for any derivative not yet analytic; `thiserror` for errors; same edition/rust-version as vle.

---

## 4. Scope — what the library computes

### The pedagogical ladder (also the execution order)

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

Rungs 1–3 are exactly the "learning repo" ask; rungs 4–8 are the "solution engine" ask. They share kernels (the Thomas tridiagonal solve appears in rungs 4, 5, and 8; bubble-point-per-stage in rung 4 reuses vle-thermo's `bubble_temperature` machinery), so the ladder is also a sane implementation order.

### Column model (the domain object everything shares)

`Column`: N stages numbered top-down (1 = condenser, N = reboiler), each stage with optional feed `F_j, z_j, T/P/q`, optional vapor/liquid side draws `W_j / U_j`, heat duty `Q_j`, stage pressure `P_j` (profile, not single value). Condenser: total or partial. Reboiler: partial (kettle). Specs: two of {reflux ratio R, distillate rate D, boilup ratio, bottoms rate B, condenser/reboiler duty, component purity or recovery in a product} — a proper spec system with residual equations, not hardcoded R+D.

Out of scope for v1 (documented as such): rate-based/nonequilibrium stages (Taylor–Krishna — model size ~3×, needs hardware correlations orthogonal to thermo; the stage-equation trait should be designed so an NEQ stage could slot in later), reactive distillation, dividing-wall/interlinked columns, tray hydraulics/efficiency correlations (a constant Murphree efficiency knob is cheap and worth including in the MESH model, though).

---

## 5. Rigorous method selection (the research answer)

**Primary solver: Naphtali–Sandholm-type simultaneous correction** — damped full Newton on the stage-grouped MESH system, unknowns per stage `(v_1j..v_Cj, T_j, l_1j..l_Cj)` (component molar flows + temperature, 2C+1 per stage), **exact analytic block-tridiagonal Jacobian**, block-Thomas elimination (LU-factorize the pivot blocks, never invert), Armijo line search + physical safeguards (clamp ΔT to ~10 K/step, work in ln-component-flows for positivity).

**Why this over inside-out** (the answer to Miguel's implicit question): inside-out won commercially (Aspen RadFrac and HYSYS default to it) because in 1974–1990 a rigorous K/H evaluation was the dominant cost and derivatives had to be finite-differenced — the inner loop's cheap surrogate thermo (Kb model + constant relative volatilities + linear enthalpies) minimized rigorous calls. **vle-thermo inverts that premise**: property calls are cheap Rust and the derivatives (∂lnφ̂ᵢ/∂nⱼ already public; ∂K/∂T to be added) are analytic and nearly free. With an exact cheap Jacobian, full Newton gives quadratic terminal convergence, one code path for every spec set, and dramatically simpler code than two nested loops with a surrogate model between them. This is also what ChemSep (the academic gold standard, Taylor/Kooijman) uses throughout, and what Aspen itself falls back to for hard nonideal columns. The 1971 criticism of Naphtali–Sandholm ("large number of partial derivatives, large storage") is precisely what modern hardware + analytic derivatives nullify: one Newton iteration costs O(N·(2C+1)³) — microseconds for N≤200, C≤20.

**Seeding** (the FUG connection from school, made real): FUG shortcut → N, feed stage, R estimate → linear T profile between feed-flash bubble/dew estimates → constant-molal-overflow flows → compositions from a feed flash or 2–3 Wang–Henke sweeps (nearly free once the tridiagonal kernel exists). Newton's weakness is convergence radius; seeding is the whole game, and the seeding pipeline is itself the pedagogy.

**Fallback ladder when Newton stalls** (Milestone: robustness layer):
1. **Thermodynamic homotopy** (Vickery & Taylor 1986): solve the column with ideal-Raoult K and ideal enthalpies first (easy), then continue λ: 0→1 blending in the rigorous model. Trivially implementable since vle-thermo exposes both ideal and rigorous paths.
2. **Pseudo-transient relaxation** (modern citable form: Pattison & Baldea 2014): backward-Euler pseudo-time on holdup-augmented MESH; huge convergence basin, slow endgame, so switch back to Newton once residuals shrink (the Ketchum 1979 combination).

**Optional later**: inside-out (Boston 1980 / Russell 1983) as a stretch milestone — worth having for the historical/pedagogical payoff ("the method school never explained") and for very large columns; if built, use the nonsmooth inner loop of Watson/Vikse/Gundersen/Barton (IECR 2017) so stages that dry out don't kill it. Not on the critical path.

**Reference implementations to study** (survey result): **BioSTEAM** `biosteam/units/stage.py` (Python, permissive NCSA license — the only complete open inside-out, plus a modern simultaneous-correction with analytic block-tridiagonal Jacobian; mimic freely) and **DWSIM** `DWSIM.UnitOperations/.../RigorousColumnSolvers/` (VB.NET, textbook-annotated Wang–Henke/sum-rates/Naphtali–Sandholm — but **GPLv3: read as literature, never port code** into this MIT project). ChemSep-LITE (free binary) and DWSIM are the cross-validation oracles. IDAES shows the equation-oriented/IPOPT alternative — good to cite in docs as "the other modern approach," not what we build.

---

## 6. Repo structure

Mirrors vle so every convention transfers:

```
stages-thermo/
├── Cargo.toml                  # workspace: [engine]; workspace.package version/license (MIT)
├── PLAN.md                     # this document
├── ROADMAP.md                  # milestones (split from §10, vle format + model-attribution lines)
├── TODO.md                     # tasks per milestone (split from §11)
├── CLAUDE.md                   # conventions (seeded from §12)
├── README.md                   # project front page
├── PUBLISHING.md               # release flow (copy vle's, adjust names)
├── LICENSE                     # MIT
├── hooks/pre-push              # cargo fmt --check gate (copy vle's; activate via core.hooksPath)
├── .github/workflows/          # ci.yml, release.yml — copy vle's shape (SHA-pinned, paths-filtered,
│                               #   idempotent publish probes, cibuildwheel)
├── engine/                     # Rust crate `stages-thermo`, lib `stages_thermo`
│   ├── Cargo.toml              # deps: vle-thermo (crates.io), nalgebra, num-dual, thiserror,
│   │                           #   smallvec; feature "python" → pyo3/numpy/rayon (abi3-py310)
│   ├── README.md               # crates.io page (same immutable-per-version rules as vle)
│   ├── benches/                # criterion: block_thomas, ns_iteration, full column solves
│   └── src/
│       ├── lib.rs
│       ├── thermo.rs           # ADAPTER over vle-thermo (see §7) — the only module that imports it
│       ├── column/
│       │   ├── model.rs        # Column, Stage, Feed, SideDraw, CondenserKind, PressureProfile
│       │   ├── specs.rs        # Spec enum + residual forms (R, D, B, boilup, Q, purity, recovery)
│       │   ├── profiles.rs     # Profiles { t, p, l, v, x, y, q } — THE result object (§8)
│       │   └── init.rs         # seeding: FUG → T-profile → CMO flows → feed-flash compositions
│       ├── binary/
│       │   ├── equilibrium.rs  # y*(x) curve generation via vle-thermo bubble_temperature/pressure
│       │   ├── mccabe_thiele.rs# op lines, q-line, stage stepping, pinch, Rmin, N(R), efficiency
│       │   └── ponchon_savarit.rs # H–x–y construction, tie lines, difference points, stage stepping
│       ├── shortcut/
│       │   ├── fenske.rs       # Nmin + non-key distribution
│       │   ├── underwood.rs    # bracketed root(s) between adjacent volatilities (Brent/interval)
│       │   ├── gilliland.rs    # Molokanov closed form
│       │   ├── kirkbride.rs    # feed-stage split
│       │   ├── winn.rs         # K_LK = β·K_HK^δ variant (β,δ regressed from vle-thermo — showcase)
│       │   └── fug.rs          # orchestrated FUG(K) design: specs → (N, feed stage, R, splits)
│       ├── rigorous/
│       │   ├── mesh.rs         # MESH residual assembly (shared by ALL rigorous solvers)
│       │   ├── bubble_point.rs # Wang–Henke
│       │   ├── sum_rates.rs    # Burningham–Otto
│       │   ├── naphtali_sandholm.rs  # flagship Newton SC
│       │   ├── jacobian.rs     # analytic block-tridiagonal Jacobian assembly
│       │   ├── homotopy.rs     # Vickery–Taylor thermodynamic continuation
│       │   ├── relaxation.rs   # pseudo-transient fallback
│       │   └── inside_out.rs   # (stretch, late milestone)
│       ├── numerics/
│       │   ├── thomas.rs       # scalar tridiagonal (Thomas), with the Boston–Sullivan-style
│       │   │                   #   stabilization / clamping for near-cancellation
│       │   └── block_thomas.rs # block-tridiagonal LU elimination (nalgebra blocks, no inversion)
│       ├── types.rs            # errors (thiserror), SolveReport (iterations, residual history)
│       └── py_bindings.rs      # + py_column.rs — feature-gated, same shape as vle
├── python/                     # Python package `stages-thermo`, import `stages`
│   ├── pyproject.toml          # maturin backend; deps: vle-thermo, numpy; extras: matplotlib
│   ├── README.md               # PyPI page
│   ├── src/stages/
│   │   ├── __init__.py         # Column, solve results, shortcut fns, McCabeThiele, PonchonSavarit
│   │   ├── plotting.py         # the diagrams: staircase, H-x-y, Gilliland, profiles, convergence
│   │   └── data/               # example column definitions (JSON): debutanizer, absorber, ...
│   └── tests/                  # pytest vs the wheel (mass-balance closure, literature tables)
├── notebooks/                  # the learning path, one per milestone (§10) — vle notebook rules
└── docs/
    ├── theory/                 # markdown per rung: the equations as implemented, cited
    └── references.md           # full ACS-style reference list (§14) + reference-to-code mapping
```

Not carried over from vle: `deploy/` (no hub deployment; distribution = crates.io + PyPI + notebooks only), `legacy/` (no legacy code — though see §13 on checking the thesis/VB6 for any staged-column remnants), `units/` crate (reuse `vle-units` via vle-thermo; user-facing Python takes `Q_`/unit strings through `vle.units`).

---

## 7. The thermo adapter (`engine/src/thermo.rs`) and upstream vle-thermo gaps

**Design rule: exactly one module talks to vle-thermo.** Everything else consumes a `ThermoProvider`-style interface (concrete struct first; make it a trait only if/when a surrogate model — inside-out — or a mock needs to slot in. The inside-out milestone would make it a trait; until then don't abstract prematurely).

What the adapter needs per stage evaluation, vs. what vle-thermo 0.8.1 provides (API-mapping result):

| Need | vle-thermo today | Gap handling |
|---|---|---|
| K(T,P,x,y) | `flash::system::k_values` — public, γ-φ AND φ-φ, 22 EOS/5 activity/11 mixing rules | ✅ use directly |
| Bubble/dew T,P | `flash::bubble::*`, `flash::dew::*` | ✅ (drives McCabe–Thiele curve + WH temperature loop) |
| PT / PH flash | `flash_isothermal_warm` (warm-start K!), `flash_adiabatic` | ✅ (feed flashes; PH is cubic-only — fine for v1) |
| ∂lnφ̂ᵢ/∂nⱼ | `mixture::d_ln_phi_d_n` — public, analytic/dual | ✅ → ∂K/∂x columns of the Jacobian |
| H per phase (φ-φ) | `energy::phase_enthalpy_entropy` — public | ✅ |
| H per phase (γ-φ) | **not packaged** — building blocks public (`ideal_enthalpy_mix` + `excess_h_s`) | **adapter assembles it** (day 1); upstream later |
| ∂K/∂T, ∂lnφ/∂T, ∂lnφ/∂P | **missing** | interim: central FD in the adapter (cheap, one extra K call); **upstream analytic** in the vle-thermo derivative release |
| ∂H/∂T (real Cp), ∂H/∂nⱼ | **missing** (only ideal-gas `ideal_cp`) | same: FD interim → upstream analytic/dual (num-dual over `mixture_params::<D>` which is already generic over dual scalars) |
| Rust component DB | **missing** — JSON DB is Python-side only (15 compounds), lacks `cp_coeffs`; Rust consumers hand-build `Component`s | stages-thermo ships its own loader over a vendored JSON initially; **upstream**: move DB into the crate behind a feature, add `cp_coeffs` |
| Components for the classic examples | toluene, ethanol, acetone, chloroform, iC4, iC5, nC8+ absent | **upstream DB expansion** (needed by M1's benzene–toluene!) |

**Upstream strategy**: these gaps become a dedicated milestone (M4, "vle-thermo derivative release", target vle-thermo 0.9.x) executed in the vle repo under its own rules (PyO3 bindings same commit, docs sync, notebook if user-facing). The FD-interim adapter means stages-thermo milestones M1–M3 and even M5–M6 **do not block** on upstream; M8's exact-Jacobian Naphtali–Sandholm is where analytic derivatives genuinely pay, so M4 must land before M8 (FD Jacobian works meanwhile — it's a correctness oracle for the analytic one anyway, the same pattern vle uses).

Reference-state discipline: pick one convention for the whole column (vle-thermo's `t_ref/p_ref/h_ref` machinery), set it once in the adapter, never per-stage.

---

## 8. Granular API — designed for the MCP server from day 1

The MCP server itself is a stretch milestone, but the API shape it needs costs nothing if designed in now:

- **Every method returns a rich result object, never bare numbers.** `McCabeThieleResult` carries the equilibrium curve samples, both operating lines, q-line, every stage's (x,y) corner, Rmin + pinch location, N. `FugResult` carries Nmin, all Underwood roots, Rmin, the Gilliland point, feed stage, per-component D/B splits. `ColumnSolution` carries per-stage `Profiles` (T, P, L, V, x, y, K, H_L, H_V), duties (Qc, Qr), and a `SolveReport` (method, iterations, residual-norm history, damping steps, homotopy path if used).
- **Solvers are inspectable**: an iteration-callback / trace mode records intermediate profiles, so a notebook can animate how the temperature profile relaxes to solution — this is the "open the black box" requirement, and it's also the MCP answer to "why didn't it converge".
- **Stage-level queries** on `ColumnSolution`: `solution.stage(7)` → everything about stage 7; `solution.mass_balance_closure()`, `solution.energy_balance_closure()` → audit residuals.
- **Batch API** (the "numpy for columns"): `solve_batch` over parameter sweeps (R, feed composition, P) with rayon + GIL release, NaN-row-on-failure — identical contract to vle-thermo's `_batch` methods. This is what makes N-vs-R curves and operating-envelope maps one-liner notebook cells.
- Python surface stays thin: dataclass-like wrappers + `plotting.py`; all computation in Rust.

---

## 9. Validation & testing strategy

Four layers, same philosophy as vle's Chapter-IV benchmark:

1. **Kernel unit tests**: `thomas`/`block_thomas` vs dense nalgebra solves on random diagonally-dominant systems; Underwood root-bracketing on constructed volatility sets (including k distributed components → k−1 roots); Fenske/Gilliland against hand calculations.
2. **Literature tables** (the primary correctness benchmark, pinned in assertion cells like vle):
   - Seader/Henley/Roper *Separation Process Principles* Ch. 10: the Wang–Henke worked example (per-iteration T/V tables — unit-tests the BP loop itself, not just the endpoint), the Burningham–Otto absorber example, the Naphtali–Sandholm example, the debutanizer/depropanizer exercises.
   - Holland *Fundamentals of Multicomponent Distillation* (full text on archive.org): stage-by-stage answer tables.
   - Nonideal set: methanol–water (van Laar — already validated thermo in vle Ch. IV), 2-propanol–water (Wilson — same), then acetone–methanol–water extractive and ethanol–water–benzene azeotropic as the hard robustness cases.
3. **Cross-simulator checks**: DWSIM (open, scriptable) and ChemSep-LITE (free) on identical specs — agreement targets ~1% on product compositions/duties, with thermo-package differences documented.
4. **Invariant/property tests** on every converged solution regardless of source: component mass balance closure ≤1e-9 relative, ΣxΣy=1 per stage, energy closure, y=Kx consistency — cheap and catch everything structural.

Performance: criterion benches from M6 on (block-Thomas, one NS iteration, full debutanizer solve); regression-guarded like vle. Target headline (for the README, once measured): full rigorous solve of a 40-stage, 10-component column in ~milliseconds.

---

## 10. Milestones (→ ROADMAP.md)

Estimates assume the vle working style (focused sessions, tests + notebook + docs per milestone). Total ≈ **130–190 h** to M10; stretch milestones extra.

- **M0 — Repo bootstrap** *(4–6 h)*: scaffold per §6; split this plan into ROADMAP/TODO/CLAUDE.md/README; CI (`ci.yml` lint+test+wheel, `release.yml` with idempotent publish probes — copy vle's, SHA-pinned); register `stages-thermo` 0.0.1 stubs on PyPI + crates.io; apply the four-ruleset GitHub hardening pattern (protect-main / require-signed-commits / protect-version-tags / require-signed-tags); `hooks/pre-push` fmt gate; conda env decision (dedicated `stages` env: python 3.12, maturin, pytest, jupyter, matplotlib, numpy, vle-thermo).
- **M1 — Column model + McCabe–Thiele** *(14–20 h)*: `column/model.rs` (binary-sufficient subset), `binary/equilibrium.rs` + `mccabe_thiele.rs`; Rmin via pinch detection (tangent pinches included), N(R), Murphree efficiency, total-reflux; plotting (staircase diagram); **needs upstream toluene+ethanol in the component DB — do that tiny vle-thermo PR first**. 📓 `01-mccabe-thiele.ipynb`: benzene–toluene (PR EOS) and methanol–water (van Laar, ties to vle Ch. IV).
- **M2 — Ponchon–Savarit** *(10–14 h)*: adapter's γ-φ enthalpy assembly (ideal + excess) lands here; H–x–y diagram construction from vle-thermo enthalpies, difference points, tie-line stepping; comparison cell: P–S vs M–T stage counts vs CMO error. 📓 `02-ponchon-savarit.ipynb`.
- **M3 — FUG shortcut** *(12–16 h)*: Fenske (+ non-key distribution), Underwood with interval-bracketed roots (Brent per (α_HK, α_LK) interval; handle distributed components), Gilliland–Molokanov, Kirkbride, Winn (β,δ regressed from vle-thermo K's — nice engine showcase); orchestrated `fug.rs` design function. 📓 `03-shortcut-design.ipynb`: depropanizer design, FUG vs (later) rigorous comparison teaser.
- **M4 — Upstream: vle-thermo derivative release (0.9.x, in the vle repo)** *(12–18 h)*: analytic/dual ∂lnφ/∂T, ∂lnφ/∂P → `k_values_with_derivs`; real-mixture Cp (∂H/∂T) and partial-molar ∂H/∂n via `num-dual` over the already-generic `mixture_params`; packaged γ-φ `phase_enthalpy_entropy`; Rust-side component DB with `cp_coeffs` + the new components; PyO3 bindings + tests per vle's M5+ rule; publish. *(Can start any time; must precede M8.)*
- **M5 — MESH infrastructure** *(10–14 h)*: full `column/model.rs` (multi-feed, side draws, duties, pressure profile, condenser kinds), spec system with residual forms, `rigorous/mesh.rs` residual assembly, `numerics/thomas.rs`, `column/init.rs` seeding pipeline (FUG → linear T → CMO → feed-flash x). No solver yet — tests assert residuals vanish on hand-constructed consistent states.
- **M6 — Wang–Henke bubble-point** *(12–16 h)*: tridiagonal composition step (with clamping/renormalization), per-stage bubble-T via vle-thermo, energy-balance flow update; convergence on ΣΔT²; trace mode. Validate against the Seader/Henley per-iteration tables. Criterion benches start. 📓 `04-mesh-and-bubble-point.ipynb` — *the "how a column is actually solved" notebook, part 1*.
- **M7 — Sum-rates** *(8–12 h)*: Burningham–Otto for absorbers/strippers (no condenser/reboiler specs); the Friday–Smith narrow/wide-boiling pairing story in docs; absorber validation case. 📓 `05-absorbers-sum-rates.ipynb`.
- **M8 — Naphtali–Sandholm (flagship)** *(20–28 h)*: stage-grouped variables (ln component flows + T), analytic block-tridiagonal Jacobian from M4 derivatives (FD Jacobian kept as test oracle), block-Thomas, Armijo line search + step clamps, spec equations in the Newton system, seeding from M5 pipeline (option: 2–3 WH sweeps). Validate: S/H/R NS example, debutanizer vs DWSIM/ChemSep, methanol–water + 2-propanol–water nonideal columns. 📓 `06-naphtali-sandholm.ipynb` — *part 2: Jacobian sparsity plot, quadratic-convergence plot, "what the simulator was doing all those years".*
- **M9 — Robustness layer** *(12–18 h)*: Vickery–Taylor thermodynamic homotopy (ideal→rigorous λ-continuation) + pseudo-transient fallback with Newton switchover; automatic escalation policy in a top-level `solve()`; hard-case suite (extractive acetone–methanol–water, azeotropic ethanol–water–benzene if thermo supports it — else document the boundary). 📓 `07-when-newton-fails.ipynb`.
- **M10 — Batch API + performance + 1.0 polish** *(10–14 h)*: `solve_batch` (rayon, GIL-release, NaN-on-fail), parameter sweeps; bench-guarded optimization pass (no allocation in Newton loop — vle hot-path rules); README headline numbers; API docs pass; version 1.0.0. 📓 `08-numpy-for-columns.ipynb`: N-vs-R maps, feed-sensitivity heatmaps.
- **M11 *(stretch)* — Inside-out**: Boston/Russell Kb-surrogate method, ideally with the Watson–Barton nonsmooth inner loop; head-to-head benchmark vs NS across the validation suite — the empirical answer to "was inside-out ever going to win here?". 📓 `09-inside-out.ipynb`.
- **M12 *(stretch)* — MCP server**: thin server (likely Python package `stages-mcp`) exposing: solve column from JSON spec, query stage/profile/duty, shortcut design, "explain convergence" from SolveReport. The §8 API means this is mostly plumbing.

Every milestone: notebook per vle conventions (setup cell, context, worked example, ≥2 exercises with hidden solutions, pinned assertions), docs sync, model-attribution lines in ROADMAP/TODO, YubiKey-signed commits by Miguel.

---

## 11. TODO seed (→ TODO.md)

M0 checklist detail (the rest expand the same way at execution time):

- [ ] `gh repo create miguelju/stages-thermo` (private initially?), MIT LICENSE, .gitignore (Rust+Python+Jupyter)
- [ ] Cargo workspace + `engine/` skeleton compiling against crates.io `vle-thermo = "0.8"` with a smoke test (bubble T of methanol–water matches vle's known value)
- [ ] `python/` maturin skeleton; `import stages; stages.__version__` wheel test
- [ ] Split PLAN → ROADMAP.md / TODO.md / CLAUDE.md / README.md
- [ ] CI: `ci.yml` (fmt+clippy, cargo test, wheel build+pytest; paths-filtered, `**/*.md` ignore), `release.yml` (tag-triggered, registry probes → gated idempotent publish); all actions SHA-pinned; workflow self-test paths
- [ ] Publish 0.0.1 name-holding stubs to PyPI + crates.io
- [ ] Rulesets (4-pattern) + `sha_pinning_required` + repo pre-push hooks (`core.hooksPath hooks`; fmt gate; private-data gate probably unneeded — no infra in this repo — but keep the placeholder-only rule)
- [ ] Create `stages` conda env; document absolute-path binary rule
- [ ] Project memory dir for Claude sessions in the new repo slug

(M1–M12 task lists derive directly from §10; expand each into checkbox items in TODO.md at M0 time, with the estimates from §10.)

---

## 12. Conventions to carry into the new repo's CLAUDE.md

- **Units**: identical canonical set to vle (K, kPa abs, kJ/kmol, kmol, R = 8.31451) — non-negotiable since vle-thermo's API speaks these; same units-in-docstrings rule for every function; user-facing Python accepts unit strings via `vle.units`.
- **Educational comments**: beginner-friendly Rust/Python idiom explanations (Miguel's standing preference) **plus** textbook anchoring: every solver module's doc comment states the method, its reference in ACS form, and maps code symbols to the book's symbols (e.g. `// S&H eq. 10-25: A_j = L_{j-1}`).
- **Citation rule** (adapted from vle's Pascal rule): code implementing a published method cites it — Wang & Henke 1966, Naphtali & Sandholm 1971, Boston & Sullivan 1974, Boston 1980, Russell 1983, Burningham & Otto 1967, Vickery & Taylor 1986, Pattison & Baldea 2014, Watson et al. 2017, Molokanov 1972, plus the textbook set (Seader/Henley/Roper; Holland; Kister; Doherty & Malone; King). Full list → `docs/references.md`. **GPL caution in CLAUDE.md**: DWSIM is GPLv3 — consult as literature, never translate its code.
- **PyO3 bindings same-commit rule**, notebook conventions, milestone-reiteration-before-execution, model-attribution lines, phase/milestone sync invariants, fmt-before-push, YubiKey signing flow (Miguel taps; Claude never commits), doc-sync-before-push list (README/ROADMAP/TODO/PLAN + both package READMEs with the immutable-per-published-version rule): all copied from vle verbatim with names adjusted.
- **Conda env**: dedicated `stages` env (`~/miniconda3/envs/stages/bin/{python,maturin,pytest,jupyter}`), same absolute-path discipline; installs `vle-thermo` from PyPI.

---

## 13. Risks & open questions

1. **Underwood for distributed components** is the classic implementation trap (multiple roots, poles at every α) — the interval-bracketing design in §6 handles it, but budget test time; a companion-matrix eigenvalue formulation (CES 2012) is the robust fallback.
2. **Wide-boiling + γ-φ interactions**: the validation ladder deliberately climbs ideal-HC → mildly nonideal → azeotropic; expect M9 to reveal thermo-side issues (e.g., activity-model T-extrapolation) that feed back as vle-thermo issues — that's a feature (downstream pressure-testing), but it's schedule risk.
3. **Component DB coverage** gates which literature cases run verbatim; the M1/M4 upstream additions must include ideal-gas Cp coefficients (needed for every enthalpy balance), which the current JSON lacks entirely — sourcing those (e.g., from the thesis appendix or Poling/Prausnitz) is real work, not typing.
4. **Aij/kij parameters for the nonideal validation columns**: vle Ch. IV gives methanol–water and 2-propanol–water; extractive/azeotropic cases need parameter hunting. Document sources per case in `docs/references.md`.
5. **Naming check at M0**: re-verify registry availability at execution time (a 404 today isn't a reservation) and register immediately.
6. **Legacy tie-in worth checking** (nice-to-have): whether the 1999 thesis / VB6 code contains any staged-column or bubble/dew-cascade remnants worth citing in the teaching docs — it would close the personal-history loop from thesis → vle → stages.
7. **Scope discipline**: rate-based, reactive, dividing-wall are *documented* non-goals for 1.0. The stage-equation layout should merely not preclude them.

## 14. Core references (→ docs/references.md, ACS style)

Wang & Henke, *Hydrocarbon Process.* **1966**, 45(8), 155 · Friday & Smith, *AIChE J.* **1964**, 10, 698 (why BP↔narrow, SR↔wide) · Burningham & Otto, *Hydrocarbon Process.* **1967**, 46(10), 163 · Tomich, *AIChE J.* **1970**, 16, 229 · Naphtali & Sandholm, *AIChE J.* **1971**, 17, 148 · Boston & Sullivan, *Can. J. Chem. Eng.* **1974**, 52, 52 · Boston & Britt, *Comput. Chem. Eng.* **1978**, 2, 109 · Boston, ACS Symp. Ser. 124, **1980** · Russell, *Chem. Eng.* **1983** (Oct 17), 53 · Vickery & Taylor, *AIChE J.* **1986**, 32, 547 · Wayburn & Seader, *Comput. Chem. Eng.* **1987**, 11, 7 · Ketchum, *Chem. Eng. Sci.* **1979**, 34, 387 · Pattison & Baldea, *AIChE J.* **2014**, 60, 4104 · Watson, Vikse, Gundersen & Barton, *IECR* **2017**, 56, 960 · Krishnamurthy & Taylor, *AIChE J.* **1985**, 31, 449/456 · Molokanov et al. **1972** (Gilliland closed form) · Seader, Henley & Roper, *Separation Process Principles*, 3rd/4th ed. (Ch. 10 = the implementation-level bible) · Holland, *Fundamentals of Multicomponent Distillation*, 1981 · Kister, *Distillation Design*, 1992 · Doherty & Malone, *Conceptual Design of Distillation Systems*, 2001 · Górak & Sorensen (eds.), *Distillation: Fundamentals and Principles*, 2014 · Taylor & Kooijman, *The ChemSep Book* (chemsep.org).
