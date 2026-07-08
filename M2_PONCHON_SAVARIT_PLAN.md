# Plan: Milestone 2 — Ponchon–Savarit (with the ammonia–water side-by-side)

> ⛔ **BLOCKED — do not start until the upstream vle-thermo work ships.**
> This milestone consumes a new **NRTL** activity model and an **ammonia** component that do
> not yet exist in vle-thermo. Execute the vle-repo plan first
> (`/Users/migueljackson/dev/vle/NRTL_AMMONIA_PLAN.md`), publish **vle-thermo 0.11.0**, then
> bump `engine/Cargo.toml` here to `vle-thermo = "0.11"` and proceed. Everything below assumes
> that release is live.

## Context

Milestone 2 is Ponchon–Savarit — the enthalpy–composition (H–x–y) method, where energy balances
enter and constant-molal-overflow (CMO) stops being assumed. The canonical teaching system is
**ammonia–water**, because its heat of mixing is so large that CMO (hence McCabe–Thiele) is badly
wrong — you *need* the energy-exact method.

Because a computed γ-φ chart can't pixel-match the historical ammonia–water chart, the notebook does
**both routes side by side** — the educational payload:
- **(a) compute** the H–x–y chart from NRTL and step off stages — "how a chart is built from thermo";
- **(b) reference data** — feed Ibrahim–Klein / Tillner-Roth enthalpy+VLE points through an
  enthalpy-augmented `from_points` path and run the pure graphical construction — "the textbook, faithfully".

The full design rationale (NRTL over UNIQUAC/Helmholtz) lives in the vle plan; the pedagogical write-up
is the lesson section at the bottom of this doc (destined for `docs/theory/ponchon-savarit.md`).

## B1. Adapter enthalpy (`engine/src/thermo.rs`)

- Add const reference-state fields `t_ref` (298.15 K), `p_ref` (101.325 kPa) to `ThermoSystem`, set once in
  the `peng_robinson` / `van_laar` constructors (PLAN §7 reference-state discipline — never per-stage).
- Add `phase_enthalpy(&self, t, p, comp, phase: PhaseId) -> Result<f64>` wrapping
  `vle_thermo::flash::system::phase_enthalpy_entropy(spec, t, p, comp, phase, t_ref, p_ref, &[], &[])`
  via the existing `with_spec` closure (~L196); `&[]` h_ref/s_ref = zero datum (P–S uses differences).
  Import `vle_thermo::eos::PhaseId` (not re-exported at crate root). Mirror the `check_composition` +
  `map_err(StagesError::Thermo)` idiom. Integration-test the γ-φ path (empty `sat_models` slice).

## B2. Enthalpy-augmented curve (`engine/src/binary/equilibrium.rs`)

- Add saturated-liquid `h_liq: Vec<f64>` (at T_bubble(x)) and saturated-vapor `h_vap: Vec<f64>`
  (equilibrium vapor y*(x) at the same T) arrays, with `h_liquid_of_x` / `h_vapor_of_y` built on the
  existing free `interp` fn (~L280). Prefer a **separate `EnthalpyCurve`** (wraps an `EquilibriumCurve` +
  the two enthalpy arrays) so synthetic constant-α curves aren't forced to carry enthalpy.
- **Route (b) path:** a `from_points`-with-enthalpy constructor (mirror `from_points` ~L170 validation)
  taking literature `x, y, T, h_liq, h_vap` — how the notebook feeds Ibrahim–Klein / Tillner-Roth data.

## B3. Ponchon–Savarit (`engine/src/binary/ponchon_savarit.rs`, new)

- `PonchonSavaritSpec` (x_D, x_B, z_F, R or duties, condenser kind), `PonchonSavaritResult`
  (stages, difference points Δ_D/Δ_B, n_stages, feed_stage, Q_C, Q_R); reuse `StagePoint`.
- Difference points: Δ_D = (x_D, h_D + Q_C/D), Δ_B = (x_B, h_B − Q_R/B); R fixes Δ_D height.
- **Stepping mirrors `mccabe_thiele.rs` (~L481)** but replaces the CMO operating line with a **pole-line**
  through Δ, and the horizontal step with a **tie line** on the H–x–y diagram. Reuse the sign-change-scan +
  linear-root idiom (`q_line_curve_intersection` ~L256), the fractional-last-stage + `PINCH_PROGRESS`
  guards, and `MAX_STAGES`.
- Expose enough for the notebook to put P–S N against M–T N on the same spec (CMO-error comparison).

## B4. Bindings, plotting, wiring

- `engine/src/binary/mod.rs` — `pub mod ponchon_savarit;` + re-exports.
- `engine/src/py_bindings.rs` — `phase_enthalpy` method (add a `parse_phase` helper like `parse_condenser`
  ~L53), `ponchon_savarit` `#[pyfunction]`, `PonchonSavaritResult` + `EnthalpyCurve` pyclasses, enthalpy
  accessors, `from_points`-with-enthalpy static method; register via `add_class` / `add_function`.
- `python/src/stages/__init__.py` — import + `__all__` (~L43, L71).
- `python/src/stages/plotting.py` — `plot_hxy` + `plot_ponchon_savarit` with a **new `_draw_hxy_frame`**
  (H on y-axis, composition on x-axis, saturated-liquid/vapor curves, tie lines, poles, stage staircase;
  **not** the unit-square frame, no `set_aspect("equal")`). Add to `__all__`.
- Tests: `engine/tests/m2_ponchon_savarit.rs` (benzene–toluene + methanol–water invariants);
  `python/tests/test_ponchon_savarit.py` (every binding exercised, per the PyO3 rule).

## B5. Docs + notebook

- `docs/theory/ponchon-savarit.md` — equations-as-implemented (difference points, tie-line stepping, energy
  balances mapped to S&H symbols) + a header diagram (`docs/theory/img/ponchon-savarit.png`, generated from
  the library) + **the technical-lesson section below**. Must follow the **GitHub-math rendering rules** in
  CLAUDE.md (no `&` in math, no `\,\;\!`, no `$…$` inside emphasis, no `($…$)` double-hug; use `\ast`).
  Verify with the codified greps + a MathJax round-trip.
- `docs/references.md` — add Ponchon/Savarit + the ammonia–water refs (below) + a code-mapping bullet.
- `notebooks/02-ponchon-savarit.ipynb` — mirror `01`'s structure. Worked examples:
  (1) **benzene–toluene** near-ideal consistency check (P–S N ≈ M–T N ≈ 12.22, pinned);
  (2) **methanol–water** CMO-error demonstration (P–S vs M–T diverge, pin the gap);
  (3) **ammonia–water side by side** — route (a) NRTL-computed H–x–y + stepping, route (b)
  Ibrahim–Klein/Tillner-Roth reference points via `from_points`-with-enthalpy + same construction; a
  comparison cell where the gap between (a) and (b) *is* the lesson. ≥2 exercises w/ hidden solutions,
  pinned assertions; execute top-to-bottom via nbconvert before commit.

## Reuse (verbatim) from M1

The free `interp` fn (`equilibrium.rs:~280`), the `with_spec` closure (`thermo.rs:~196`), the
sign-change-scan + linear-root idiom (`mccabe_thiele.rs:~277`), the stepping-loop skeleton with
fractional-last-stage + pinch guards (`mccabe_thiele.rs:~481`), the `to_py_err` / `parse_condenser`
binding idioms (`py_bindings.rs:~41,53`), and the deferred-matplotlib `_draw_frame` plotting idiom.

## Decisions baked in (flag to change)

1. NRTL α storage = **option B** (parallel `alpha` matrix, ternary-capable) — set in the vle plan.
2. NRTL **τ = energy/(RT)** first; `τ = a + b/T` only if the NH₃–H₂O fit is inadequate.
3. NH₃–H₂O params = **published (Aspen/DECHEMA) preferred**; accuracy bar = few-% at moderate P.
4. Ammonia–water **stays in M2** via the two-route notebook (needs the vle release first).

## Verification (end-to-end)

`cargo test --workspace`; `maturin develop` (from `python/`) then
`~/miniconda3/envs/stages/bin/pytest python/tests/test_ponchon_savarit.py`; execute the notebook via
`~/miniconda3/envs/stages/bin/jupyter nbconvert --to notebook --execute notebooks/02-ponchon-savarit.ipynb`
with pinned assertions (benzene–toluene consistency; methanol–water CMO-error gap; ammonia–water route-a-vs-b).
Doc render: run the CLAUDE.md GitHub-math greps + MathJax round-trip on `docs/theory/ponchon-savarit.md`.
Cross-check: P–S invariants (mass + energy closure on every converged construction). Milestone close:
model-attribution lines in ROADMAP.md + TODO.md + commit trailer; doc-sync pass; YubiKey-signed commits.

---

## The technical lesson — drop-in prose for `docs/theory/ponchon-savarit.md`

> Concise, references inline, GitHub-math-safe. Title/level to fit the doc when created.

### Ammonia–water: how the textbook did it, and the limits of a computed chart

*Motivating the two ammonia–water constructions in the notebook — one computed from a model, one from reference data.*

**The textbook case was done on charts, not on-the-fly thermo.** Ponchon (1921) and Savarit (1922)
devised the enthalpy–composition method as a *graphical* procedure run on a pre-drawn H–x–y diagram for
one specific binary: nothing is computed during the construction — you locate the difference points
(poles) and step off stages with a straightedge on a printed chart. For ammonia–water that chart is the
Merkel–Bošnjaković enthalpy–concentration diagram (Bošnjaković, *Technische Thermodynamik*, 1935), built
from **experimental data** — vapor–liquid equilibrium plus **calorimetry** (heats of solution/mixing,
heat capacities, latent heats) — reduced onto one diagram. Its modern refinements are the correlations
still used today: Scatchard et al. (1947), Macriss et al. (1964), the ASHRAE/EES formulation of
Ibrahim & Klein (1993), and the reference Helmholtz-energy EOS of Tillner-Roth & Friend (1998). So the
"textbook thermo" for ammonia–water *is* a data artifact — a chart/correlation — and Ponchon–Savarit is
the graphical construction performed on it.

**Why ammonia–water was the showcase.** Its heat of mixing is enormous. Constant molal overflow — the
assumption behind McCabe–Thiele — requires equal molar latent heats and negligible enthalpy of mixing so
that liquid and vapor flows stay constant along each section. Ammonia–water violates that badly, so
McCabe–Thiele returns the wrong stage count *and* the wrong internal flows. Ponchon–Savarit closes the
energy balance exactly through the chart, so it is the method you genuinely need here — which is precisely
why textbooks introduce it on ammonia–water.

**Why a computed NRTL γ-φ chart frays.** Regenerating the H–x–y diagram from a model splits it: activity
coefficients (NRTL) for the liquid, a cubic EOS for the vapor, enthalpies from ideal-gas terms plus a
liquid excess-enthalpy and a vapor departure. The weak links, in order: (1) **the vapor at elevated
pressure** — the classic charts live at several bar, where the cubic-EOS vapor enthalpy departure is only
approximate; this dominates the error and is not the activity model's fault; (2) **the reference state for
a dissolving light gas** — γ-φ references the liquid to pure liquid ammonia, an awkward datum at generator
conditions; (3) **the curvature of the excess enthalpy** — real ammonia–water heat of mixing is large,
asymmetric, and strongly temperature-dependent across the full range, and a few-parameter NRTL cannot match
it at the water-rich and ammonia-rich ends at once. Net: NRTL gets a few percent at moderate pressure — a
real improvement over Wilson or van Laar — but not a pixel-match to the Bošnjaković chart.

**The "UNIQUAC is in the good models" trap.** It is tempting to reach for UNIQUAC because the celebrated
ammonia–water model uses it — but that model (Thomsen & Rasmussen, 1999) is *extended* UNIQUAC: plain
UNIQUAC plus a Debye–Hückel term plus speciation (NH3 + H2O to NH4+ and OH-). All the accuracy comes from
those additions, not from the local-composition kernel. On a single binary fit, plain UNIQUAC has only two
adjustable energy parameters (its r and q are fixed structural constants) against NRTL's three, and its
size-asymmetry advantage is wasted on two small molecules — so plain UNIQUAC is, if anything, slightly
*less* flexible than NRTL here and no closer to the chart.

**What actually would get closer,** in increasing order of effort: (1) feed the reference correlation's
data (Ibrahim–Klein / Tillner-Roth) directly and run the pure graphical construction on it — matches the
textbook chart, needs no new model; (2) implement *extended* UNIQUAC or a Helmholtz-energy EOS — the models
that genuinely reproduce the chart, but each a large undertaking.

**Why this project does not build the specialized models.** Their distinguishing capability serves nothing
else on the roadmap. Extended UNIQUAC's electrolyte/speciation machinery is used by exactly one planned
system (ammonia–water); every other case — hydrocarbons, alcohol–water, acetone–water, the extractive and
azeotropic ternaries — is neutral and needs none of it. A Helmholtz EOS is either single-system
(Tillner-Roth) or a from-scratch reference backbone nothing planned requires (the cubic EOS already meets
the roughly 1 percent cross-simulator target, and the validation oracles use cubic EOS themselves). NRTL,
by contrast, is general infrastructure: it improves every aqueous-organic column on the ladder and is the
standard model for the later extractive and azeotropic cases. So the honest engineering choice is to add
NRTL once (broad benefit) and reproduce the ammonia–water textbook chart from reference data rather than
build single-use thermodynamics.

**The two constructions in this notebook.** Route (a) computes the H–x–y diagram from NRTL and steps off
stages — teaching how a chart is built from thermo, and showing the fray. Route (b) feeds the
Ibrahim–Klein / Tillner-Roth reference points through the enthalpy-augmented `from_points` path and runs the
same construction — reproducing the textbook faithfully and cleanly separating the method from the thermo
model. Shown side by side, the gap between them is the lesson.

### References to add to `docs/references.md` (ACS style — verify page numbers at execution)

- Ponchon, M. Étude graphique de la distillation fractionnée. *La Technique Moderne* **1921**, 13, 20 and 55.
- Savarit, R. *Arts et Métiers* **1922** (enthalpy–composition construction). *(verify volume/pages)*
- Bošnjaković, F. *Technische Thermodynamik*; Theodor Steinkopff: Dresden, **1935**.
- Scatchard, G.; Epstein, L. F.; Warburton, J.; Cody, P. J. Thermodynamic properties of saturated liquid and
  vapor of ammonia–water mixtures. *Refrig. Eng.* **1947**, 53, 413. *(verify)*
- Macriss, R. A.; Eakin, B. E.; Ellington, R. T.; Huebler, J. *Physical and Thermodynamic Properties of
  Ammonia–Water Mixtures*; Institute of Gas Technology Research Bulletin 34, **1964**.
- Ibrahim, O. M.; Klein, S. A. Thermodynamic properties of ammonia–water mixtures. *ASHRAE Trans.* **1993**,
  99 (1), 1495.
- Tillner-Roth, R.; Friend, D. G. A Helmholtz free energy formulation of the thermodynamic properties of the
  mixture {water + ammonia}. *J. Phys. Chem. Ref. Data* **1998**, 27, 63.
- Renon, H.; Prausnitz, J. M. Local compositions in thermodynamic excess functions for liquid mixtures (NRTL).
  *AIChE J.* **1968**, 14, 135.
- Abrams, D. S.; Prausnitz, J. M. Statistical thermodynamics of liquid mixtures (UNIQUAC). *AIChE J.* **1975**,
  21, 116.
- Thomsen, K.; Rasmussen, P. Modeling of vapor–liquid–solid equilibrium in gas–aqueous electrolyte systems
  (extended UNIQUAC). *Chem. Eng. Sci.* **1999**, 54, 1787.
- Wilson, G. M. Vapor–liquid equilibrium. XI. A new expression for the excess free energy of mixing.
  *J. Am. Chem. Soc.* **1964**, 86, 127.
- Code mapping: `Ponchon–Savarit construction → engine/src/binary/ponchon_savarit.rs` (Ponchon 1921 /
  Savarit 1922; S&H Ch. 7 energy-balance treatment); `NRTL → vle-thermo activity.rs` (Renon & Prausnitz 1968).
