# TODO

Actionable tasks with rough time estimates. Grouped by [ROADMAP.md](ROADMAP.md)
milestone. Check off items as they're completed. Time estimates assume working
with Claude Code in the vle working style (focused sessions, tests + notebook +
docs per milestone). The M0 checklist is expanded verbatim from PLAN.md §11; the
M1–M12 lists derive from the deliverables in PLAN.md §10.

---

## Milestone 0: Repo Bootstrap *(4–6 h)* — **complete**
*Executed by Claude Code using Claude Opus 4.8 (1M context)*

- [✓] **Create the repo** — `gh repo create miguelju/stages-thermo` (private initially?), MIT LICENSE, `.gitignore` (Rust + Python + Jupyter)
- [✓] **Cargo workspace + `engine/` skeleton** — compiling against crates.io `vle-thermo = "0.8"` with a smoke test (bubble T of methanol–water matches vle's known value)
- [✓] **`python/` maturin skeleton** — `import stages; stages.__version__` wheel test
- [✓] **Split PLAN → ROADMAP.md / TODO.md / CLAUDE.md / README.md**
- [✓] **CI** — `ci.yml` (fmt + clippy, cargo test, wheel build + pytest; paths-filtered, `**/*.md` ignore), `release.yml` (tag-triggered, registry probes → gated idempotent publish); all actions SHA-pinned; workflow self-test paths
- [✓] **Publish 0.0.1 name-holding stubs** to PyPI + crates.io
- [✓] **Rulesets** (4-pattern) + `sha_pinning_required` + repo pre-push hooks (`core.hooksPath hooks`; fmt gate; private-data gate probably unneeded — no infra in this repo — but keep the placeholder-only rule)
- [✓] **Create `stages` conda env** — document the absolute-path binary rule
- [✓] **Project memory dir** for Claude sessions in the new repo slug

## Milestone 1: Column Model + McCabe–Thiele *(14–20 h)* — **complete**
*Executed by Claude Code using Claude Fable 5*

Binary-sufficient column model plus the full McCabe–Thiele method. The
component-DB prerequisite (toluene + ethanol) landed upstream via M4
(vle-thermo 0.9.x, `component-db` feature) before this milestone started, so
the separate tiny PR was unnecessary.

- [✓] **`column/model.rs` (binary subset)** (~3–4h) — `BinaryColumn`, `Feed`, `CondenserKind`, single pressure, material-balance product rates (the full multicomponent `Column`/`Stage` arrives at M5)
- [✓] **`binary/equilibrium.rs`** (~2–3h) — y*(x) curve generation via vle-thermo `bubble_temperature`/`bubble_pressure`; `constant_alpha` + `from_points` constructors (analytic test oracle / literature data)
- [✓] **`binary/mccabe_thiele.rs`** (~4–6h) — operating lines, q-line, stage stepping, pinch, Rmin, N(R), Murphree efficiency, total-reflux
- [✓] **Rmin via pinch detection** (~2h) — tangent pinches included (both sections; rectifying max-slope + stripping min-slope converted through the feed-section balances)
- [✓] **Staircase plotting** (~2h) — `plotting.py` staircase + total-reflux diagrams
- [✓] **Upstream component DB** (~1h) — superseded by M4: vle-thermo 0.9.x ships the 24-compound Rust DB (toluene, ethanol, benzene, … with `cp_coeffs`); engine dep bumped to `vle-thermo = "0.9"` + `component-db`
- [✓] **📓 `01-mccabe-thiele.ipynb`** (~2–3h) — benzene–toluene (PR EOS), methanol–water (van Laar, ties to vle Ch. IV via Table 4.6); 3 exercises with hidden solutions; pinned assertions; executes top-to-bottom

## Milestone 2: Ponchon–Savarit *(10–14 h)*

Energy balances enter. The adapter's γ-φ enthalpy assembly lands here.

- [ ] **Adapter γ-φ enthalpy assembly** (~3–4h) — ideal + excess, assembled in `thermo.rs` from vle-thermo building blocks
- [ ] **H–x–y diagram construction** (~2–3h) — from vle-thermo enthalpies
- [ ] **`binary/ponchon_savarit.rs`** (~3–4h) — difference points, tie-line stepping, stage stepping
- [ ] **Comparison cell** (~1h) — P–S vs M–T stage counts vs CMO error
- [ ] **📓 `02-ponchon-savarit.ipynb`** (~2–3h)

## Milestone 3: FUG Shortcut *(12–16 h)*

Fenske–Underwood–Gilliland–Kirkbride–Winn, orchestrated into a design function.

- [ ] **`shortcut/fenske.rs`** (~1–2h) — Nmin + non-key distribution
- [ ] **`shortcut/underwood.rs`** (~3–4h) — interval-bracketed roots (Brent per (α_HK, α_LK) interval; handle distributed components)
- [ ] **`shortcut/gilliland.rs`** (~1h) — Molokanov closed form
- [ ] **`shortcut/kirkbride.rs`** (~1h) — feed-stage split
- [ ] **`shortcut/winn.rs`** (~1–2h) — K_LK = β·K_HK^δ (β,δ regressed from vle-thermo K's)
- [ ] **`shortcut/fug.rs`** (~2–3h) — orchestrated FUG(K) design: specs → (N, feed stage, R, splits)
- [ ] **📓 `03-shortcut-design.ipynb`** (~2–3h) — depropanizer design, FUG vs (later) rigorous comparison teaser

## Milestone 4: Upstream — vle-thermo Derivative Release (0.9.x) *(12–18 h)* — **complete**

Executed in the **vle repo** under its own rules (this is vle's Milestone 12;
see that repo's tracking docs for the execution record). Landed as vle-thermo
**v0.9.0/v0.9.1, published on crates.io + PyPI** — before stages-thermo M1, so
M1 built directly on the new DB. Full task detail lives in vle's
`DERIVATIVE_RELEASE_PLAN.md` / `TODO.md`; the stages-thermo-facing deliverables:

- [✓] **Analytic/dual ∂lnφ/∂T, ∂lnφ/∂P** → `k_values_with_derivs`
- [✓] **Real-mixture Cp (∂H/∂T) + partial-molar ∂H/∂n** via `num-dual` over the generic `mixture_params`
- [✓] **Packaged γ-φ `phase_enthalpy_entropy`**
- [✓] **Rust-side component DB** (`component-db` feature, 24 compounds) with `cp_coeffs` + the new components (toluene, ethanol, acetone, chloroform, isobutane, isopentane, n-octane, n-nonane, n-decane)
- [✓] **PyO3 bindings + tests** per vle's M5+ rule; published (v0.9.1 current, incl. the Wong–Sandler departure-enthalpy patch)

## Milestone 5: MESH Infrastructure *(10–14 h)*

The full column model, the spec system, MESH assembly, the tridiagonal kernel,
and seeding — no solver yet.

- [ ] **Full `column/model.rs`** (~3–4h) — multi-feed, side draws, duties, pressure profile, condenser kinds
- [ ] **`column/specs.rs`** (~2–3h) — spec enum + residual forms (R, D, B, boilup, Q, purity, recovery)
- [ ] **`rigorous/mesh.rs`** (~2–3h) — MESH residual assembly (shared by all rigorous solvers)
- [ ] **`numerics/thomas.rs`** (~1–2h) — scalar tridiagonal (Thomas) with clamping for near-cancellation
- [ ] **`column/init.rs`** (~2h) — seeding pipeline (FUG → linear T → CMO flows → feed-flash x)
- [ ] **Residual-vanishing tests** (~1h) — assert residuals vanish on hand-constructed consistent states

## Milestone 6: Wang–Henke Bubble-Point *(12–16 h)*

The first rigorous solver — bubble-point tearing. Criterion benches start.

- [ ] **Tridiagonal composition step** (~3–4h) — with clamping / renormalization
- [ ] **Per-stage bubble-T** (~2h) — via vle-thermo
- [ ] **Energy-balance flow update** (~2–3h) — convergence on ΣΔT²
- [ ] **Trace mode** (~1–2h) — records intermediate profiles for animation / debugging
- [ ] **Validate vs Seader/Henley per-iteration tables** (~2h)
- [ ] **Criterion benches** (~1h) — block-Thomas, one WH iteration, full column solve
- [ ] **📓 `04-mesh-and-bubble-point.ipynb`** (~2–3h) — "how a column is actually solved," part 1

## Milestone 7: Sum-Rates *(8–12 h)*

Burningham–Otto for absorbers/strippers — the flipped tearing for wide-boiling.

- [ ] **`rigorous/sum_rates.rs`** (~3–4h) — Burningham–Otto (no condenser/reboiler specs)
- [ ] **Friday–Smith pairing story in docs** (~1h) — why BP↔narrow, SR↔wide
- [ ] **Absorber validation case** (~2h)
- [ ] **📓 `05-absorbers-sum-rates.ipynb`** (~2–3h)

## Milestone 8: Naphtali–Sandholm (flagship) *(20–28 h)*

Full damped Newton with an exact analytic block-tridiagonal Jacobian. Depends on
M4; FD Jacobian kept as oracle.

- [ ] **Stage-grouped variables** (~2–3h) — ln component flows + T, positivity by construction
- [ ] **`rigorous/jacobian.rs`** (~5–7h) — analytic block-tridiagonal Jacobian from M4 derivatives (FD Jacobian kept as test oracle)
- [ ] **`numerics/block_thomas.rs`** (~3–4h) — block-tridiagonal LU elimination (nalgebra blocks, no inversion)
- [ ] **`rigorous/naphtali_sandholm.rs`** (~4–6h) — Armijo line search + step clamps; spec equations in the Newton system
- [ ] **Seeding from M5 pipeline** (~1–2h) — option: 2–3 WH sweeps
- [ ] **Validation** (~3–4h) — S/H/R NS example, debutanizer vs DWSIM/ChemSep, methanol–water + 2-propanol–water nonideal columns
- [ ] **📓 `06-naphtali-sandholm.ipynb`** (~2–3h) — Jacobian sparsity plot, quadratic-convergence plot

## Milestone 9: Robustness Layer *(12–18 h)*

What simulators do when Newton fails.

- [ ] **`rigorous/homotopy.rs`** (~3–4h) — Vickery–Taylor thermodynamic continuation (ideal → rigorous λ)
- [ ] **`rigorous/relaxation.rs`** (~3–4h) — pseudo-transient fallback with Newton switchover
- [ ] **Automatic escalation policy** (~2–3h) — top-level `solve()` that escalates Newton → homotopy → pseudo-transient
- [ ] **Hard-case suite** (~2–3h) — extractive acetone–methanol–water; azeotropic ethanol–water–benzene if thermo supports it, else document the boundary
- [ ] **📓 `07-when-newton-fails.ipynb`** (~2–3h)

## Milestone 10: Batch API + Performance + 1.0 Polish *(10–14 h)*

The "numpy for columns" layer, a bench-guarded optimization pass, and v1.0.0.

- [ ] **`solve_batch`** (~3–4h) — rayon, GIL-release, NaN-row-on-failure; parameter sweeps (R, feed composition, P)
- [ ] **Optimization pass** (~3–4h) — no allocation in the Newton loop (vle hot-path rules), bench-guarded
- [ ] **README headline numbers** (~1h) — once measured (target: 40-stage, 10-component solve in ~ms)
- [ ] **API docs pass** (~1–2h)
- [ ] **Version 1.0.0** (~0.5h) — bump + tag
- [ ] **📓 `08-numpy-for-columns.ipynb`** (~2–3h) — N-vs-R maps, feed-sensitivity heatmaps

## Milestone 11 *(stretch)*: Inside-Out

Boston/Russell Kb-surrogate method + a head-to-head benchmark vs NS.

- [ ] **`rigorous/inside_out.rs`** — Boston/Russell Kb-surrogate method, ideally with the Watson–Barton nonsmooth inner loop
- [ ] **Head-to-head benchmark** vs Naphtali–Sandholm across the validation suite
- [ ] **📓 `09-inside-out.ipynb`**

## Milestone 12 *(stretch)*: MCP Server

A thin server (likely Python package `stages-mcp`) over the §8 granular API.

- [ ] **Solve column from JSON spec**
- [ ] **Query stage / profile / duty**
- [ ] **Shortcut design endpoint**
- [ ] **"Explain convergence" from `SolveReport`**

---

## Summary

| Milestone | Est. | Status |
|-----------|------|--------|
| 0. Repo Bootstrap | ~4–6h | **Complete** |
| 1. Column Model + McCabe–Thiele | ~14–20h | **Complete** (v0.1.0) |
| 2. Ponchon–Savarit | ~10–14h | Not started |
| 3. FUG Shortcut | ~12–16h | Not started |
| 4. Upstream vle-thermo Derivative Release (0.9.x) | ~12–18h | **Complete** *(in the vle repo; vle-thermo v0.9.1 published)* |
| 5. MESH Infrastructure | ~10–14h | Not started |
| 6. Wang–Henke Bubble-Point | ~12–16h | Not started |
| 7. Sum-Rates | ~8–12h | Not started |
| 8. Naphtali–Sandholm (flagship) | ~20–28h | Not started |
| 9. Robustness Layer | ~12–18h | Not started |
| 10. Batch API + Performance + 1.0 | ~10–14h | Not started |
| **Total to M10** | **~130–190h** | |
| 11. Inside-Out | *stretch* | Not started |
| 12. MCP Server | *stretch* | Not started |

Every active milestone's total includes its milestone notebook (~2–3h) and the
docs-sync pass (per CLAUDE.md's doc-sync-before-push list). Estimates assume the
vle working style.
