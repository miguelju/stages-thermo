# References

ACS-style reference list for the methods implemented in stages-thermo, plus the
reference→code mapping (filled in as each method module lands). Sourced from
`PLAN.md` §14. Every solver module cites the paper it implements in its doc
comment (see `CLAUDE.md`, "Reference Citation Requirements").

## Core method references

- **McCabe & Thiele**, *Ind. Eng. Chem.* **1925**, 17(6), 605–611. — the graphical binary method (M1).
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

## Reference → code mapping

Populated as modules land. Format: `<method>` → `engine/src/<path>` (implements `<citation>`).

- McCabe–Thiele construction → `engine/src/binary/mccabe_thiele.rs` (implements
  McCabe & Thiele 1925; equations per S&H Ch. 7: 7-9, 7-12, 7-24, 7-26).
- Equilibrium curve y*(x) → `engine/src/binary/equilibrium.rs` (bubble-point
  sweep via vle-thermo; constant-α form S&H eq. 7-13 as test oracle).
- Binary column balances → `engine/src/column/model.rs` (S&H eqs. 7-2/7-3).
- Theory write-up: [`theory/mccabe-thiele.md`](theory/mccabe-thiele.md).
