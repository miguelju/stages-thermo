//! Binary (two-component) staged-separation methods — rungs 1–2 of the
//! pedagogical ladder (PLAN §4).
//!
//! Everything in this module works on the **light component's** mole
//! fraction: `x` in the liquid, `y` in the vapor, with the more volatile
//! component listed first. This is the classical McCabe–Thiele convention
//! (Seader, Henley & Roper, *Separation Process Principles*, Ch. 7).

pub mod equilibrium;
pub mod mccabe_thiele;
pub mod ponchon_savarit;

pub use equilibrium::{EnthalpyCurve, EquilibriumCurve};
pub use mccabe_thiele::{
    McCabeThieleResult, McCabeThieleSpec, RminResult, StagePoint, TotalRefluxResult, mccabe_thiele,
    n_vs_r, rmin, total_reflux,
};
pub use ponchon_savarit::{PonchonSavaritResult, PonchonSavaritSpec, ponchon_savarit};
