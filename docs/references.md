# References

ACS-style reference list for the methods implemented in stages-thermo, plus the
reference→code mapping (filled in as each method module lands). Sourced from
`PLAN.md` §14. Every solver module cites the paper it implements in its doc
comment (see `CLAUDE.md`, "Reference Citation Requirements").

## Core method references

- **McCabe & Thiele**, *Ind. Eng. Chem.* **1925**, 17(6), 605–611. — the graphical binary method (M1).
- **Ponchon**, M. *La Technique Moderne* **1921**, 13, 20 and 55. — the enthalpy–composition (H–x–y) construction (M2).
- **Savarit**, R. *Arts et Métiers* **1922**. — the enthalpy–composition construction (M2). *(verify volume/pages)*
- **Renon & Prausnitz**, *AIChE J.* **1968**, 14, 135. — NRTL local-composition activity model (M2 liquid model, in vle-thermo).
- **Wang & Henke**, *Hydrocarbon Process.* **1966**, 45(8), 155. — bubble-point MESH method (M6).
- **Friday & Smith**, *AIChE J.* **1964**, 10, 698. — why bubble-point ↔ narrow-boiling, sum-rates ↔ wide-boiling.
- **Burningham & Otto**, *Hydrocarbon Process.* **1967**, 46(10), 163. — sum-rates method for absorbers (M7).
- **Tomich**, *AIChE J.* **1970**, 16, 229. — simultaneous θ-method.
- **Naphtali & Sandholm**, *AIChE J.* **1971**, 17, 148. — simultaneous-correction Newton on grouped MESH (M8, flagship).
- **Boston & Sullivan**, *Can. J. Chem. Eng.* **1974**, 52, 52. — inside-out class of methods; the Thomas-stabilization ideas (M11).
- **Boston & Britt**, *Comput. Chem. Eng.* **1978**, 2, 109.
- **Boston**, ACS Symp. Ser. 124, **1980**. — inside-out (M11).
- **Russell**, *Chem. Eng.* **1983** (Oct 17), 53. — inside-out, engineering write-up (M11).
- **Vickery & Taylor**, *AIChE J.* **1986**, 32, 547. — thermodynamic homotopy continuation (M9).
- **Wayburn & Seader**, *Comput. Chem. Eng.* **1987**, 11, 7. — homotopy-continuation for difficult columns.
- **Ketchum**, *Chem. Eng. Sci.* **1979**, 34, 387. — relaxation + Newton combination (M9).
- **Pattison & Baldea**, *AIChE J.* **2014**, 60, 4104. — pseudo-transient continuation (M9).
- **Watson, Vikse, Gundersen & Barton**, *IECR* **2017**, 56, 960. — nonsmooth inside-out inner loop (M11).
- **Krishnamurthy & Taylor**, *AIChE J.* **1985**, 31, 449/456. — rate-based (nonequilibrium) model; documented non-goal for 1.0.
- **Molokanov et al.**, **1972**. — Gilliland closed-form correlation (M3).

## Textbook set

- **Seader, Henley & Roper**, *Separation Process Principles*, 3rd/4th ed. — Ch. 10 is the implementation-level bible; the primary validation-table source (M6–M8).
- **Holland**, *Fundamentals of Multicomponent Distillation*, 1981. — stage-by-stage answer tables (full text on archive.org).
- **Kister**, *Distillation Design*, 1992.
- **Doherty & Malone**, *Conceptual Design of Distillation Systems*, 2001.
- **Górak & Sørensen (eds.)**, *Distillation: Fundamentals and Principles*, 2014.
- **Taylor & Kooijman**, *The ChemSep Book* (chemsep.org). — ChemSep is the academic gold-standard cross-validation oracle.

## Ammonia–water enthalpy–composition data (M2 showcase)

The classic Ponchon–Savarit teaching system. Its enthalpy–concentration chart is
a *data artifact* (VLE + calorimetry reduced onto one diagram), not an on-the-fly
model — see the lesson in [`theory/ponchon-savarit.md`](theory/ponchon-savarit.md).
Page numbers to be firmed up at citation time.

- **Merkel, F.; Bošnjaković, F.** *Diagramme und Tabellen zur Berechnung der
  Absorptions-Kältemaschinen*; Springer: Berlin, **1929**. — the original
  NH₃–H₂O enthalpy–concentration diagrams.
- **Bošnjaković, F.** *Technische Thermodynamik*; Theodor Steinkopff: Dresden, **1935**.
- **Scatchard, G.; Epstein, L. F.; Warburton, J.; Cody, P. J.** Thermodynamic
  properties of saturated liquid and vapor of ammonia–water mixtures. *Refrig.
  Eng.* **1947**, 53, 413. *(verify)*
- **Macriss, R. A.; Eakin, B. E.; Ellington, R. T.; Huebler, J.** *Physical and
  Thermodynamic Properties of Ammonia–Water Mixtures*; Institute of Gas
  Technology Research Bulletin 34, **1964**.
- **Ibrahim, O. M.; Klein, S. A.** Thermodynamic properties of ammonia–water
  mixtures. *ASHRAE Trans.* **1993**, 99 (1), 1495. — the ASHRAE/EES formulation.
- **Tillner-Roth, R.; Friend, D. G.** A Helmholtz free energy formulation of the
  thermodynamic properties of the mixture {water + ammonia}. *J. Phys. Chem. Ref.
  Data* **1998**, 27, 63.
- **Abrams & Prausnitz**, *AIChE J.* **1975**, 21, 116. — UNIQUAC (cited in the
  lesson for the "UNIQUAC trap"; not used here).
- **Thomsen, K.; Rasmussen, P.** *Chem. Eng. Sci.* **1999**, 54, 1787. —
  *extended* UNIQUAC (the electrolyte model whose accuracy the lesson attributes
  to its additions, not the local-composition kernel).

## Reference implementations studied (not ported)

- **BioSTEAM** — `biosteam/units/stage.py` (Python, NCSA permissive): the only complete open inside-out plus a modern simultaneous-correction with analytic block-tridiagonal Jacobian. **Mimic freely** (license-compatible).
- **DWSIM** — `DWSIM.UnitOperations/.../RigorousColumnSolvers/` (VB.NET, textbook-annotated Wang–Henke / sum-rates / Naphtali–Sandholm). **GPLv3 — read as literature, never port code** into this MIT project. Used as a cross-validation oracle only.
- **ChemSep-LITE** (free binary) — cross-validation oracle.
- **IDAES** — the equation-oriented / IPOPT alternative; cited in docs as "the other modern approach," not built here.

## Parameter sources (validation cases)

- **Orbey & Sandler** (via vle Chapter IV, Tables 4.5/4.6) — methanol(1)–water(2)
  van Laar Λ₁₂ = 0.5853, Λ₂₁ = 0.3458; the Table 4.6 bubble pressures at 298 K
  are reproduced in `notebooks/01-mccabe-thiele.ipynb`.
- **Perry's classic van Laar set** — ethanol(1)–water(2) A₁₂ = 1.6798,
  A₂₁ = 0.9227 (notebook exercise 1: azeotrope + tangent pinch).
- **Ammonia–water NRTL (illustrative)** — NH₃(1)–H₂O(2) with `aij[0][1] = −1800`,
  `aij[1][0] = −1200` kJ/kmol (energy convention `gᵢⱼ − gⱼⱼ`) and α = 0.2
  (carried over from vle Milestone 14). Signs are physically correct (negative
  deviation, exothermic mixing); magnitudes are illustrative, not a certified
  regression — used for the route-(a) NRTL chart in
  `notebooks/02-ponchon-savarit.ipynb`. The route-(b) reference chart uses
  representative Bošnjaković/Ibrahim–Klein-style H–x–y data (see the notebook).

## Reference → code mapping

Populated as modules land. Format: `<method>` → `engine/src/<path>` (implements `<citation>`).

- McCabe–Thiele construction → `engine/src/binary/mccabe_thiele.rs` (implements
  McCabe & Thiele 1925; equations per S&H Ch. 7: 7-9, 7-12, 7-24, 7-26).
- Equilibrium curve y*(x) → `engine/src/binary/equilibrium.rs` (bubble-point
  sweep via vle-thermo; constant-α form S&H eq. 7-13 as test oracle).
- Binary column balances → `engine/src/column/model.rs` (S&H eqs. 7-2/7-3).
- Ponchon–Savarit construction → `engine/src/binary/ponchon_savarit.rs`
  (implements Ponchon 1921 / Savarit 1922; S&H Ch. 7 energy-balance treatment).
- Enthalpy–composition curve (H–x–y) → `engine/src/binary/equilibrium.rs`
  (`EnthalpyCurve`); per-phase molar enthalpy → `engine/src/thermo.rs`
  (`ThermoSystem::phase_enthalpy`, wrapping vle-thermo's γ-φ / φ-φ enthalpy).
- NRTL activity model → vle-thermo `activity.rs` (Renon & Prausnitz 1968);
  exposed to stages via `ThermoSystem::nrtl`.
- Theory write-up: [`theory/mccabe-thiele.md`](theory/mccabe-thiele.md).
