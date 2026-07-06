//! The thermodynamics adapter — **the only module that imports `vle_thermo`**
//! (PLAN §7).
//!
//! Every other module in this crate reaches thermodynamics through this
//! boundary, never through `vle_thermo` directly. Keeping the dependency
//! surface in one file means:
//!
//! - a surrogate thermo model (the inside-out inner loop, a stretch milestone)
//!   or a test mock can replace it without touching any solver, and
//! - the exact set of vle-thermo entry points stages-thermo relies on is
//!   auditable in one place — which is what turns this crate into a
//!   pressure-test of vle-thermo's public API (PLAN §1).
//!
//! At M0 the adapter is a thin smoke path. The real `ThermoProvider`-style
//! struct (per-stage K, H, and derivative evaluation with a single
//! reference-state convention set once for the whole column) is built out from
//! M1/M5; see `PLAN.md` §7 for the need-vs-provides mapping and the FD-interim
//! → analytic-derivative upstream plan.

use vle_thermo::eos::{CubicEos, LiquidModel, VaporModel};
use vle_thermo::flash::SystemSpec;
use vle_thermo::flash::bubble::bubble_temperature;
use vle_thermo::mixing::MixingRule;
use vle_thermo::types::Component;

use crate::types::{Result, StagesError};

/// The bubble temperature of an equimolar methanol(1)/water(2) mixture at
/// 101.325 kPa, computed through vle-thermo with the Peng–Robinson EOS and
/// classical mixing.
///
/// This is the M0 end-to-end smoke path: it exists to prove that the
/// vle-thermo dependency links, that a `SystemSpec` round-trips through the
/// flash layer, and that the FFI boundary is wired both from Rust tests and
/// (via [`crate::py_bindings`]) from Python. It is **not** a validated
/// property — methanol/water is strongly non-ideal and a bare cubic EOS with
/// no activity model is only approximate here. The pinned literature bubble
/// point arrives at M1, once the component database and activity models are
/// wired (PLAN §9 validation ladder).
///
/// Returns the bubble temperature in **K**.
pub fn smoke_bubble_temperature() -> Result<f64> {
    // Standard critical constants (Tc in K, Pc in kPa, ω dimensionless).
    let components = [methanol(), water()];
    let spec = SystemSpec {
        components: &components,
        vapor: VaporModel::Cubic(CubicEos::PR1976),
        liquid: LiquidModel::Cubic(CubicEos::PR1976),
        mixing_rule: MixingRule::Classical,
        kij: &[],
        aij: &[],
        vl: &[],
        delta: &[],
        sat_models: &[],
        ge_model: None,
    };

    let res = bubble_temperature(&spec, 101.325, &[0.5, 0.5], 1e-9, 200)
        .map_err(|e| StagesError::Thermo(e.to_string()))?;
    Ok(res.value)
}

fn methanol() -> Component {
    Component {
        name: "methanol".into(),
        tc: 512.6,
        pc: 8090.0,
        omega: 0.559,
        ..Component::default()
    }
}

fn water() -> Component {
    Component {
        name: "water".into(),
        tc: 647.096,
        pc: 22064.0,
        omega: 0.3443,
        ..Component::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The smoke path evaluates end-to-end and returns a physically plausible
    /// bubble temperature. The window is deliberately wide — this asserts the
    /// dependency works, not a validated value (see the function docs).
    #[test]
    fn smoke_bubble_temperature_is_physical() {
        let t = smoke_bubble_temperature().expect("bubble_temperature should converge");
        assert!(
            (280.0..400.0).contains(&t),
            "methanol/water bubble T = {t} K is outside the plausible window [280, 400)"
        );
    }
}
