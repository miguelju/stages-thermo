# Project Roadmap

High-level milestones for **stages-thermo** — a staged-separation (distillation)
learning library and fast steady-state column solver, built on `vle-thermo`. For
actionable tasks with time estimates, see [TODO.md](TODO.md). For the full
technical rationale (method selection, thermo adapter, references), see
[PLAN.md](PLAN.md).

The milestone order **is** the pedagogical ladder: each rung is a distillation
method implemented from scratch, documented against its textbook equations, and
demonstrated in an executable notebook. Rungs 1–3 are the "learning repo" ask;
rungs 4–8 are the "fast column engine" ask. They share kernels (the tridiagonal
solve, bubble-point-per-stage), so the ladder is also a sane implementation
order.

## The pedagogical ladder (§4 of PLAN.md)

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

Total ≈ **130–190 h** to M10; the two stretch milestones (M11, M12) are extra.

---

## Milestone 0: Repo Bootstrap ~ *(4–6 h)* — **complete**
**Goal**: Scaffold the repo per PLAN.md §6, split this plan into ROADMAP / TODO /
CLAUDE.md / README.md, stand up CI, and hold both registry names.
*Executed by Claude Code using Claude Opus 4.8 (1M context)*

- [✓] Scaffold the workspace per §6 (Cargo workspace + `engine/` skeleton, `python/` maturin skeleton)
- [✓] Split PLAN → ROADMAP.md / TODO.md / CLAUDE.md / README.md
- [✓] CI: `ci.yml` (fmt + clippy, cargo test, wheel build + pytest; paths-filtered, `**/*.md` ignore) and `release.yml` (tag-triggered, idempotent publish probes) — copy vle's shape, SHA-pinned
- [✓] Register `stages-thermo` 0.0.1 name-holding stubs on PyPI + crates.io
- [✓] Apply the four-ruleset GitHub hardening pattern (protect-main / require-signed-commits / protect-version-tags / require-signed-tags) + `sha_pinning_required`
- [✓] `hooks/pre-push` fmt gate (`core.hooksPath hooks`)
- [✓] Create the dedicated `stages` conda env; document the absolute-path binary rule
- [✓] Create the Claude project-memory dir for the new repo slug

## Milestone 1: Column Model + McCabe–Thiele ~ *(14–20 h)* — **complete**
**Goal**: The binary-sufficient `Column` model, the equilibrium curve from
`vle-thermo`, and the full McCabe–Thiele method with staircase plotting.
*Executed by Claude Code using Claude Fable 5*

> The upstream toluene + ethanol component-DB gap was closed by vle-thermo's
> 0.9.x derivative release (M4, executed in the vle repo) — the engine now
> builds on `vle-thermo = "0.9"` with the `component-db` feature.
> 📓 `01-mccabe-thiele.ipynb`: benzene–toluene (PR EOS) and methanol–water
> (van Laar, ties to vle Ch. IV via the Table 4.6 reproduction).

- [✓] `column/model.rs` — binary-sufficient subset (`BinaryColumn`, `Feed`, `CondenserKind`, material balances)
- [✓] `binary/equilibrium.rs` — y*(x) curve via vle-thermo bubble T/P (+ `constant_alpha` / `from_points` constructors as test oracle and literature-data path)
- [✓] `binary/mccabe_thiele.rs` — operating lines, q-line, stage stepping, pinch, Rmin, N(R), Murphree efficiency, total-reflux
- [✓] Rmin via pinch detection (tangent pinches included, both sections; validated on an ethanol–water tangent pinch)
- [✓] Plotting — staircase diagram (`stages.plotting`, McCabe–Thiele + total-reflux)
- [✓] 📓 `01-mccabe-thiele.ipynb` (benzene–toluene, methanol–water; 3 exercises with hidden solutions; pinned assertions)

## Milestone 2: Ponchon–Savarit ~ *(10–14 h)* — **complete**
**Goal**: Energy balances enter — the γ-φ enthalpy assembly in the adapter, the
H–x–y diagram, and stage stepping by difference points.
*Executed by Claude Code using Claude Opus 4.8 (1M context)*

> Consumed vle-thermo **0.11** (NRTL activity model + ammonia component, shipped
> upstream in the vle repo); the engine pin moved `0.9` → `0.11`.
> 📓 `02-ponchon-savarit.ipynb`: benzene–toluene consistency check
> (P–S ≈ M–T ≈ 12.2 stages), methanol–water CMO-error gap, and the
> ammonia–water two-route showcase (NRTL-computed chart vs reference-data chart).

- [✓] Adapter per-phase enthalpy (`ThermoSystem::phase_enthalpy`, wrapping vle-thermo's γ-φ / φ-φ enthalpy) + NRTL constructor + one reference state (`t_ref`/`p_ref`)
- [✓] H–x–y diagram construction (`EnthalpyCurve`: `from_thermo` computes it, `from_points` feeds reference data)
- [✓] Difference points (poles Δ_D/Δ_B), tie-line + pole-line stepping, per-mole-feed duties (`binary/ponchon_savarit.rs`)
- [✓] Comparison: P–S vs M–T stage counts vs CMO error (near-ideal agree; nonideal diverge)
- [✓] PyO3 bindings (`nrtl`, `phase_enthalpy`, `EnthalpyCurve`, `ponchon_savarit`) + `plot_hxy`/`plot_ponchon_savarit`
- [✓] 📓 `02-ponchon-savarit.ipynb` (3 worked examples + 2 exercises, executes top-to-bottom)

## Milestone 3: FUG Shortcut ~ *(12–16 h)*
**Goal**: The Fenske–Underwood–Gilliland shortcut suite — the "how simulators
seed rigorous runs" connection, made real.

- [ ] `shortcut/fenske.rs` — Nmin + non-key distribution
- [ ] `shortcut/underwood.rs` — interval-bracketed roots (Brent per (α_HK, α_LK) interval; handle distributed components)
- [ ] `shortcut/gilliland.rs` — Gilliland–Molokanov closed form
- [ ] `shortcut/kirkbride.rs` — feed-stage split
- [ ] `shortcut/winn.rs` — K_LK = β·K_HK^δ (β,δ regressed from vle-thermo K's — engine showcase)
- [ ] `shortcut/fug.rs` — orchestrated FUG(K) design function
- [ ] 📓 `03-shortcut-design.ipynb` (depropanizer design; FUG vs rigorous teaser)

## Milestone 4: Upstream — vle-thermo Derivative Release (0.9.x, in the vle repo) ~ *(12–18 h)* — **complete**
**Goal**: Close the derivative/database gaps that stages-thermo needs. Executed
in the **vle repo** under its own rules (PyO3 bindings same commit, docs sync,
notebook if user-facing). *Landed as vle-thermo v0.9.0/v0.9.1 (published on
crates.io + PyPI) before stages-thermo M1, so the FD-interim adapter path was
never needed for the component DB.*

> This is vle's own Milestone 12 — see the vle repo's tracking docs for the
> execution record and model attribution.

- [✓] Analytic/dual ∂lnφ/∂T, ∂lnφ/∂P → `k_values_with_derivs`
- [✓] Real-mixture Cp (∂H/∂T) and partial-molar ∂H/∂n via `num-dual` over the generic `mixture_params`
- [✓] Packaged γ-φ `phase_enthalpy_entropy`
- [✓] Rust-side component DB (`component-db` feature, 24 compounds) with `cp_coeffs` + the new components (toluene, ethanol, acetone, chloroform, iC4, iC5, nC8–nC10)
- [✓] PyO3 bindings + tests per vle's M5+ rule; published (v0.9.1 current — includes the Wong–Sandler departure-enthalpy patch)

## Milestone 5: MESH Infrastructure ~ *(10–14 h)*
**Goal**: The full column model, the spec system, MESH residual assembly, the
tridiagonal kernel, and the seeding pipeline — everything the rigorous solvers
share, with no solver yet.

- [ ] Full `column/model.rs` — multi-feed, side draws, duties, pressure profile, condenser kinds
- [ ] `column/specs.rs` — spec system with residual forms (R, D, B, boilup, Q, purity, recovery)
- [ ] `rigorous/mesh.rs` — MESH residual assembly (shared by all rigorous solvers)
- [ ] `numerics/thomas.rs` — scalar tridiagonal solve
- [ ] `column/init.rs` — seeding pipeline (FUG → linear T → CMO flows → feed-flash x)
- [ ] Tests assert residuals vanish on hand-constructed consistent states (no solver)

## Milestone 6: Wang–Henke Bubble-Point ~ *(12–16 h)*
**Goal**: The first rigorous solver — the bubble-point tearing method — plus the
first criterion benches. *"How a column is actually solved," part 1.*

- [ ] Tridiagonal composition step (with clamping / renormalization)
- [ ] Per-stage bubble-T via vle-thermo
- [ ] Energy-balance flow update; convergence on ΣΔT²
- [ ] Trace mode (records intermediate profiles)
- [ ] Validate against the Seader/Henley per-iteration tables
- [ ] Criterion benches start
- [ ] 📓 `04-mesh-and-bubble-point.ipynb`

## Milestone 7: Sum-Rates ~ *(8–12 h)*
**Goal**: The flipped tearing for wide-boiling systems — Burningham–Otto for
absorbers and strippers.

- [ ] `rigorous/sum_rates.rs` — Burningham–Otto (no condenser/reboiler specs)
- [ ] The Friday–Smith narrow/wide-boiling pairing story in docs
- [ ] Absorber validation case
- [ ] 📓 `05-absorbers-sum-rates.ipynb`

## Milestone 8: Naphtali–Sandholm (flagship) ~ *(20–28 h)*
**Goal**: Full damped Newton on the stage-grouped MESH system with an exact
analytic block-tridiagonal Jacobian — **the flagship solver**. *"What the
simulator was doing all those years," part 2.*

> Depends on M4 (analytic derivatives). The FD Jacobian is kept as a correctness
> oracle for the analytic one.

- [ ] Stage-grouped variables (ln component flows + T)
- [ ] `rigorous/jacobian.rs` — analytic block-tridiagonal Jacobian from M4 derivatives (FD Jacobian kept as oracle)
- [ ] `numerics/block_thomas.rs` — block-tridiagonal LU elimination (no inversion)
- [ ] `rigorous/naphtali_sandholm.rs` — Armijo line search + step clamps; spec equations in the Newton system
- [ ] Seeding from the M5 pipeline (option: 2–3 WH sweeps)
- [ ] Validate: S/H/R NS example, debutanizer vs DWSIM/ChemSep, methanol–water + 2-propanol–water nonideal columns
- [ ] 📓 `06-naphtali-sandholm.ipynb` (Jacobian sparsity plot, quadratic-convergence plot)

## Milestone 9: Robustness Layer ~ *(12–18 h)*
**Goal**: What simulators do when Newton fails — thermodynamic homotopy and
pseudo-transient relaxation, behind an automatic escalation policy.

- [ ] `rigorous/homotopy.rs` — Vickery–Taylor thermodynamic continuation (ideal → rigorous λ)
- [ ] `rigorous/relaxation.rs` — pseudo-transient fallback with Newton switchover
- [ ] Automatic escalation policy in a top-level `solve()`
- [ ] Hard-case suite (extractive acetone–methanol–water; azeotropic ethanol–water–benzene if thermo supports it, else document the boundary)
- [ ] 📓 `07-when-newton-fails.ipynb`

## Milestone 10: Batch API + Performance + 1.0 Polish ~ *(10–14 h)*
**Goal**: The "numpy for columns" batch layer, a bench-guarded optimization pass,
and version **1.0.0**.

- [ ] `solve_batch` — rayon, GIL-release, NaN-on-fail; parameter sweeps
- [ ] Bench-guarded optimization pass (no allocation in the Newton loop — vle hot-path rules)
- [ ] README headline numbers (once measured)
- [ ] API docs pass
- [ ] Version **1.0.0**
- [ ] 📓 `08-numpy-for-columns.ipynb` (N-vs-R maps, feed-sensitivity heatmaps)

## Milestone 11 *(stretch)*: Inside-Out
**Goal**: The Boston/Russell Kb-surrogate method — "the method school never
explained" — with a head-to-head benchmark vs Naphtali–Sandholm.

- [ ] `rigorous/inside_out.rs` — Boston/Russell Kb-surrogate method, ideally with the Watson–Barton nonsmooth inner loop
- [ ] Head-to-head benchmark vs NS across the validation suite
- [ ] 📓 `09-inside-out.ipynb`

## Milestone 12 *(stretch)*: MCP Server
**Goal**: A thin server (likely a Python package `stages-mcp`) over the §8
granular API — mostly plumbing, because the result objects were designed for it
from day 1.

- [ ] Solve column from JSON spec
- [ ] Query stage / profile / duty
- [ ] Shortcut design endpoint
- [ ] "Explain convergence" from `SolveReport`

---

**Every milestone**: a notebook per vle conventions (setup cell, research/context
snippets, worked example, ≥2 exercises with hidden solutions, pinned assertions),
docs sync before push, a model-attribution line added here and in `TODO.md`, and
YubiKey-signed commits by Miguel.

**Status key**: `[✓]` complete · `[ ]` not started · `[~]` in progress
