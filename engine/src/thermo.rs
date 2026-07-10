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
//! From M1 the adapter is a concrete [`ThermoSystem`]: an owned, reusable
//! description of a mixture + model selection that can answer bubble-point
//! questions. It stays a struct (not a trait) until a surrogate model or a
//! mock genuinely needs to slot in — no premature abstraction (PLAN §7).
//!
//! ## Reference-state discipline
//!
//! One convention for the whole column, set once here, never per-stage.
//! Enthalpy queries arrive at M2 (Ponchon–Savarit); when they do, the
//! reference state lives in this struct.
//!
//! ## Units
//!
//! Canonical engine units throughout: temperature **K**, pressure **kPa**
//! (absolute), mole fractions dimensionless.

use vle_thermo::activity::ActivityModel;
use vle_thermo::eos::{CubicEos, LiquidModel, PhaseId, VaporModel};
use vle_thermo::flash::bubble::{bubble_pressure, bubble_temperature};
use vle_thermo::flash::{SystemSpec, phase_enthalpy_entropy};
use vle_thermo::mixing::MixingRule;
use vle_thermo::types::Component;

use crate::types::{Result, StagesError};

/// Convergence tolerance passed to vle-thermo's bubble-point iterations.
const BUBBLE_TOL: f64 = 1e-9;
/// Iteration budget for vle-thermo's bubble-point iterations.
const BUBBLE_MAX_ITER: usize = 200;

/// Reference temperature for the enthalpy datum — **298.15 K** (25 °C).
///
/// One convention for the whole column (PLAN §7). Ponchon–Savarit works with
/// enthalpy *differences*, so the datum choice cancels; 298.15 K / 101.325 kPa
/// is the conventional thermochemical standard state and keeps the numbers
/// readable.
const T_REF: f64 = 298.15;
/// Reference pressure for the enthalpy datum — **101.325 kPa** (1 atm).
const P_REF: f64 = 101.325;

/// Which phase's molar enthalpy to evaluate — the adapter's own phase tag, so
/// no `vle_thermo` type crosses the module boundary (mirrors [`BubblePoint`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    /// Saturated (or subcooled) liquid.
    Liquid,
    /// Saturated (or superheated) vapor.
    Vapor,
}

impl Phase {
    /// Map to vle-thermo's internal phase id (kept private to this module).
    fn to_phase_id(self) -> PhaseId {
        match self {
            Phase::Liquid => PhaseId::Liquid,
            Phase::Vapor => PhaseId::Vapor,
        }
    }
}

/// A bubble-point evaluation result: the solved value plus the incipient
/// vapor composition and K-values, re-exported in adapter-owned form so no
/// vle-thermo type crosses the module boundary.
#[derive(Debug, Clone)]
pub struct BubblePoint {
    /// The solved quantity — bubble temperature in **K** (for
    /// [`ThermoSystem::bubble_temperature`]) or bubble pressure in **kPa**
    /// (for [`ThermoSystem::bubble_pressure`]).
    pub value: f64,
    /// Incipient vapor mole fractions `y*` in equilibrium with the liquid.
    pub y: Vec<f64>,
    /// K-values `K_i = y_i / x_i` at the bubble point.
    pub k: Vec<f64>,
}

/// Which thermodynamic route evaluates the liquid phase.
///
/// vle-thermo supports both classic formulations; stages-thermo exposes the
/// two used on the M1 pedagogical ladder:
///
/// - **φ-φ**: one cubic EOS for both phases (benzene–toluene with
///   Peng–Robinson — near-ideal hydrocarbon pairs).
/// - **γ-φ**: an activity-coefficient model for the liquid + ideal-gas vapor
///   (methanol–water with van Laar — the nonideal aqueous systems validated
///   in vle Chapter IV).
#[derive(Debug, Clone)]
enum ModelKind {
    /// Cubic EOS both phases; the inner matrix is the (symmetric) binary
    /// interaction parameter matrix `kij` (empty ⇒ all zero).
    PhiPhi { eos: CubicEos, kij: Vec<Vec<f64>> },
    /// Activity-coefficient liquid + ideal-gas vapor; `aij` is the activity
    /// model's binary parameter matrix (`aij[0][1] = A₁₂`, `aij[1][0] = A₂₁`,
    /// zero diagonal), and `alpha` is the NRTL non-randomness matrix — **empty
    /// for every model except NRTL**, which reads it (vle-thermo 0.11's option-B
    /// parallel `alpha` field on `SystemSpec`).
    GammaPhi {
        activity: ActivityModel,
        aij: Vec<Vec<f64>>,
        alpha: Vec<Vec<f64>>,
    },
}

/// An owned mixture + model selection that answers equilibrium questions.
///
/// `vle_thermo::flash::SystemSpec` borrows all its data, so this struct owns
/// the components and parameter matrices and assembles a fresh (cheap, `Copy`)
/// spec per call. Construct once per column problem, reuse for every
/// evaluation.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "python", pyo3::pyclass)]
pub struct ThermoSystem {
    components: Vec<Component>,
    model: ModelKind,
    /// Enthalpy-datum temperature in **K** — [`T_REF`], set once here, never
    /// per-stage (PLAN §7 reference-state discipline).
    t_ref: f64,
    /// Enthalpy-datum pressure in **kPa** — [`P_REF`].
    p_ref: f64,
}

impl ThermoSystem {
    /// A φ-φ system: Peng–Robinson (1976) for both phases, classical mixing,
    /// all `kij = 0`. Components are loaded from vle-thermo's built-in
    /// database by (case-insensitive) name.
    ///
    /// This is the right default for near-ideal hydrocarbon pairs like
    /// benzene–toluene.
    ///
    /// # Errors
    /// [`StagesError::Thermo`] if any name is not in the database.
    pub fn peng_robinson(names: &[&str]) -> Result<Self> {
        Ok(Self {
            components: load_components(names)?,
            model: ModelKind::PhiPhi {
                eos: CubicEos::PR1976,
                kij: Vec::new(),
            },
            t_ref: T_REF,
            p_ref: P_REF,
        })
    }

    /// A γ-φ system: van Laar activity coefficients for the liquid, ideal-gas
    /// vapor. Components are loaded from vle-thermo's built-in database by
    /// (case-insensitive) name.
    ///
    /// # Arguments
    /// * `a12`, `a21` — the dimensionless van Laar parameters, in the
    ///   component order given by `names` (`a12` multiplies the first
    ///   component's infinite-dilution log-γ: `ln γ₁^∞ = A₁₂`).
    ///
    /// For methanol(1)–water(2), vle Chapter IV (Orbey & Sandler Table 4.5)
    /// regressed `A₁₂ = 0.5853`, `A₂₁ = 0.3458`.
    ///
    /// # Errors
    /// [`StagesError::Thermo`] if any name is not in the database.
    pub fn van_laar(names: &[&str; 2], a12: f64, a21: f64) -> Result<Self> {
        Ok(Self {
            components: load_components(names)?,
            model: ModelKind::GammaPhi {
                activity: ActivityModel::VanLaar,
                aij: vec![vec![0.0, a12], vec![a21, 0.0]],
                // van Laar reads no non-randomness matrix — leave it empty so
                // `SystemSpec.alpha` is ignored.
                alpha: Vec::new(),
            },
            t_ref: T_REF,
            p_ref: P_REF,
        })
    }

    /// A γ-φ system: **NRTL** activity coefficients (Renon & Prausnitz, *AIChE
    /// J.* **1968**, 14, 135) for the liquid, ideal-gas vapor. Components are
    /// loaded from vle-thermo's built-in database by (case-insensitive) name.
    ///
    /// NRTL is the general aqueous-organic infrastructure the ladder needs from
    /// M2 on (Ponchon–Savarit's ammonia–water showcase, and the later
    /// extractive/azeotropic ternaries). Its binary form has three parameters:
    /// two interaction energies plus one shared non-randomness factor.
    ///
    /// # Arguments
    /// * `a12`, `a21` — the NRTL interaction energies `gᵢⱼ − gⱼⱼ` in
    ///   **kJ/kmol**, in the component order given by `names`. vle-thermo forms
    ///   `τᵢⱼ = aᵢⱼ / (R·T)` internally (R = 8.31451 kJ/(kmol·K)), so these are
    ///   the *energy* parameters, **not** the dimensionless τ.
    /// * `alpha12` — the (symmetric) non-randomness parameter α₁₂ = α₂₁,
    ///   dimensionless, typically 0.2–0.47.
    ///
    /// # Errors
    /// [`StagesError::Thermo`] if any name is not in the database.
    pub fn nrtl(names: &[&str; 2], a12: f64, a21: f64, alpha12: f64) -> Result<Self> {
        Ok(Self {
            components: load_components(names)?,
            model: ModelKind::GammaPhi {
                activity: ActivityModel::Nrtl,
                // aᵢⱼ = gᵢⱼ − gⱼⱼ, zero diagonal (a component with itself).
                aij: vec![vec![0.0, a12], vec![a21, 0.0]],
                // Symmetric non-randomness matrix; the diagonal is unused.
                alpha: vec![vec![0.0, alpha12], vec![alpha12, 0.0]],
            },
            t_ref: T_REF,
            p_ref: P_REF,
        })
    }

    /// The enthalpy-datum temperature in **K** ([`T_REF`]).
    pub fn t_ref(&self) -> f64 {
        self.t_ref
    }

    /// The enthalpy-datum pressure in **kPa** ([`P_REF`]).
    pub fn p_ref(&self) -> f64 {
        self.p_ref
    }

    /// Number of components in the mixture.
    pub fn n_components(&self) -> usize {
        self.components.len()
    }

    /// Component names, in the order compositions are indexed.
    pub fn component_names(&self) -> Vec<String> {
        self.components.iter().map(|c| c.name.clone()).collect()
    }

    /// Bubble temperature of a liquid of composition `x` at pressure `p`.
    ///
    /// # Arguments
    /// * `p` — pressure in **kPa** (absolute)
    /// * `x` — liquid mole fractions (must sum to 1)
    ///
    /// # Returns
    /// [`BubblePoint`] with `value` = bubble temperature in **K**.
    pub fn bubble_temperature(&self, p: f64, x: &[f64]) -> Result<BubblePoint> {
        self.check_composition(x)?;
        // `with_spec` hands the closure a borrowed SystemSpec assembled from
        // this struct's owned data — the borrow can't outlive the call, which
        // is exactly the lifetime discipline SystemSpec's design asks for.
        self.with_spec(|spec| {
            let r = bubble_temperature(spec, p, x, BUBBLE_TOL, BUBBLE_MAX_ITER)
                .map_err(|e| StagesError::Thermo(e.to_string()))?;
            Ok(BubblePoint {
                value: r.value,
                y: r.incipient,
                k: r.k,
            })
        })
    }

    /// Bubble pressure of a liquid of composition `x` at temperature `t`.
    ///
    /// # Arguments
    /// * `t` — temperature in **K**
    /// * `x` — liquid mole fractions (must sum to 1)
    ///
    /// # Returns
    /// [`BubblePoint`] with `value` = bubble pressure in **kPa** (absolute).
    pub fn bubble_pressure(&self, t: f64, x: &[f64]) -> Result<BubblePoint> {
        self.check_composition(x)?;
        self.with_spec(|spec| {
            let r = bubble_pressure(spec, t, x, BUBBLE_TOL, BUBBLE_MAX_ITER)
                .map_err(|e| StagesError::Thermo(e.to_string()))?;
            Ok(BubblePoint {
                value: r.value,
                y: r.incipient,
                k: r.k,
            })
        })
    }

    /// Molar enthalpy of `comp` in the given `phase` at `(t, p)`.
    ///
    /// This is the energy quantity Ponchon–Savarit needs: the H–x–y diagram is
    /// saturated-liquid and saturated-vapor enthalpy versus composition. The
    /// call routes to vle-thermo's `phase_enthalpy_entropy`, which for the γ-φ
    /// liquid assembles ideal-gas terms **minus** the per-component
    /// Clausius–Clapeyron condensation enthalpy **plus** the excess enthalpy
    /// Hᴱ — the γ-φ enthalpy route (see the M2 lesson in
    /// `docs/theory/ponchon-savarit.md`). Entropy is computed alongside but
    /// **deliberately discarded** — P–S never uses it; widen the return to a
    /// pair only when a later milestone needs S.
    ///
    /// # Arguments
    /// * `t` — temperature in **K**
    /// * `p` — pressure in **kPa** (absolute)
    /// * `comp` — phase mole fractions (must sum to 1, length = n components)
    /// * `phase` — [`Phase::Liquid`] or [`Phase::Vapor`]
    ///
    /// # Returns
    /// Molar enthalpy in **kJ/kmol**, relative to the system's `t_ref`/`p_ref`
    /// datum with per-component reference enthalpy zero (P–S uses differences).
    ///
    /// # Errors
    /// [`StagesError::Dimension`] on a length mismatch; [`StagesError::Thermo`]
    /// if the property evaluation fails.
    pub fn phase_enthalpy(&self, t: f64, p: f64, comp: &[f64], phase: Phase) -> Result<f64> {
        self.check_composition(comp)?;
        let phase_id = phase.to_phase_id();
        self.with_spec(|spec| {
            // `&[]` for h_ref/s_ref selects the zero datum. Upstream returns an
            // (H, S) tuple; `.0` takes H and drops S (see the doc comment).
            let (h, _s) = phase_enthalpy_entropy(
                spec,
                t,
                p,
                comp,
                phase_id,
                self.t_ref,
                self.p_ref,
                &[],
                &[],
            )
            .map_err(|e| StagesError::Thermo(e.to_string()))?;
            Ok(h)
        })
    }

    /// Assemble the borrowing `SystemSpec` and run `f` against it.
    fn with_spec<T>(&self, f: impl FnOnce(&SystemSpec<'_>) -> Result<T>) -> Result<T> {
        let spec = match &self.model {
            ModelKind::PhiPhi { eos, kij } => SystemSpec {
                components: &self.components,
                vapor: VaporModel::Cubic(*eos),
                liquid: LiquidModel::Cubic(*eos),
                mixing_rule: MixingRule::Classical,
                kij,
                aij: &[],
                alpha: &[],
                vl: &[],
                delta: &[],
                sat_models: &[],
                ge_model: None,
            },
            ModelKind::GammaPhi {
                activity,
                aij,
                alpha,
            } => SystemSpec {
                components: &self.components,
                vapor: VaporModel::IdealGas,
                liquid: LiquidModel::Activity(*activity),
                mixing_rule: MixingRule::Classical,
                kij: &[],
                aij,
                // Empty for every model except NRTL (van Laar leaves it `[]`).
                alpha,
                vl: &[],
                delta: &[],
                sat_models: &[],
                // `ge_model` couples an activity model into a GE-based cubic
                // mixing rule (Wong–Sandler etc.) — not what a plain γ-φ
                // liquid uses. The activity model rides in `liquid` above.
                ge_model: None,
            },
        };
        f(&spec)
    }

    fn check_composition(&self, x: &[f64]) -> Result<()> {
        if x.len() != self.components.len() {
            return Err(StagesError::Dimension(format!(
                "composition has {} entries but the system has {} components",
                x.len(),
                self.components.len()
            )));
        }
        Ok(())
    }
}

/// Load components from vle-thermo's built-in database, erroring with the
/// list of available names on a miss.
fn load_components(names: &[&str]) -> Result<Vec<Component>> {
    names
        .iter()
        .map(|name| {
            vle_thermo::db::component(name).ok_or_else(|| {
                StagesError::Thermo(format!(
                    "component {name:?} not in the vle-thermo database (available: {})",
                    vle_thermo::db::available().join(", ")
                ))
            })
        })
        .collect()
}

/// The bubble temperature of an equimolar methanol(1)/water(2) mixture at
/// 101.325 kPa, computed through vle-thermo with the Peng–Robinson EOS and
/// classical mixing.
///
/// This is the M0 end-to-end smoke path, kept as the cheapest cross-FFI
/// health check (see [`crate::py_bindings`]). Since M1 it runs through
/// [`ThermoSystem`] with database-loaded components. It is **not** a
/// validated property — methanol/water is strongly non-ideal and a bare
/// cubic EOS with no activity model is only approximate here; the validated
/// van Laar route is [`ThermoSystem::van_laar`].
///
/// Returns the bubble temperature in **K**.
pub fn smoke_bubble_temperature() -> Result<f64> {
    let system = ThermoSystem::peng_robinson(&["methanol", "water"])?;
    Ok(system.bubble_temperature(101.325, &[0.5, 0.5])?.value)
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

    /// Pure-component bubble T of a φ-φ system reproduces the normal boiling
    /// point from the database within a couple of kelvin (PR is not exact at
    /// low pressure, but it must be close).
    #[test]
    fn benzene_toluene_pure_endpoints_near_boiling_points() {
        let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
        // Benzene boils at 353.24 K, toluene at 383.78 K (1 atm).
        let t_b = sys.bubble_temperature(101.325, &[1.0, 0.0]).unwrap().value;
        let t_t = sys.bubble_temperature(101.325, &[0.0, 1.0]).unwrap().value;
        assert!((t_b - 353.24).abs() < 3.0, "benzene Tb = {t_b} K");
        assert!((t_t - 383.78).abs() < 3.0, "toluene Tb = {t_t} K");
    }

    /// The van Laar γ-φ route reproduces the vle Chapter IV methanol–water
    /// system: positive deviation from Raoult's law means the bubble pressure
    /// at 298.15 K exceeds the Raoult (linear) interpolation of the pure
    /// vapor pressures.
    #[test]
    fn methanol_water_van_laar_positive_deviation() {
        let sys = ThermoSystem::van_laar(&["methanol", "water"], 0.5853, 0.3458).unwrap();
        let p_pure_m = sys.bubble_pressure(298.15, &[1.0, 0.0]).unwrap().value;
        let p_pure_w = sys.bubble_pressure(298.15, &[0.0, 1.0]).unwrap().value;
        let p_mid = sys.bubble_pressure(298.15, &[0.5, 0.5]).unwrap().value;
        let raoult = 0.5 * p_pure_m + 0.5 * p_pure_w;
        assert!(
            p_mid > raoult,
            "van Laar should give positive deviation: P(0.5) = {p_mid} kPa ≤ Raoult {raoult} kPa"
        );
    }

    /// The NRTL γ-φ route runs end-to-end (exercising vle-thermo 0.11's new
    /// `alpha` matrix path) and gives the expected positive deviation for
    /// ethanol–water — a strongly non-ideal, azeotrope-forming pair.
    #[test]
    fn ethanol_water_nrtl_positive_deviation() {
        // Aspen-style NRTL for ethanol(1)/water(2): Δg₁₂ = −458.7,
        // Δg₂₁ = 5574 kJ/kmol, α = 0.303 (VLE-regressed positive-deviation
        // azeotrope). Exact reproduction isn't the point here — the sign is.
        let sys = ThermoSystem::nrtl(&["ethanol", "water"], -458.7, 5574.0, 0.303).unwrap();
        let p_pure_e = sys.bubble_pressure(298.15, &[1.0, 0.0]).unwrap().value;
        let p_pure_w = sys.bubble_pressure(298.15, &[0.0, 1.0]).unwrap().value;
        let p_mid = sys.bubble_pressure(298.15, &[0.5, 0.5]).unwrap().value;
        let raoult = 0.5 * p_pure_e + 0.5 * p_pure_w;
        assert!(
            p_mid > raoult,
            "NRTL should give positive deviation: P(0.5) = {p_mid} kPa ≤ Raoult {raoult} kPa"
        );
    }

    /// The γ-φ enthalpy path runs and respects the sign and rough magnitude of
    /// the latent heat: at a bubble point, saturated-vapor molar enthalpy
    /// exceeds saturated-liquid enthalpy by tens of MJ/kmol.
    #[test]
    fn methanol_water_vapor_enthalpy_exceeds_liquid() {
        let sys = ThermoSystem::van_laar(&["methanol", "water"], 0.5853, 0.3458).unwrap();
        let p = 101.325;
        let x = [0.4, 0.6];
        let bp = sys.bubble_temperature(p, &x).unwrap();
        let t = bp.value;
        let h_liq = sys.phase_enthalpy(t, p, &x, Phase::Liquid).unwrap();
        // Saturated vapor is the incipient vapor y*(x) at the same T, P.
        let h_vap = sys.phase_enthalpy(t, p, &bp.y, Phase::Vapor).unwrap();
        assert!(
            h_vap > h_liq,
            "latent heat must be positive: H_vap = {h_vap} ≤ H_liq = {h_liq} kJ/kmol"
        );
        let dh = h_vap - h_liq;
        assert!(
            (10_000.0..80_000.0).contains(&dh),
            "ΔH_vap = {dh} kJ/kmol is outside the plausible latent-heat window"
        );
    }

    /// The reference-state datum is the standard thermochemical state, set once.
    #[test]
    fn reference_state_is_standard() {
        let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
        assert_eq!(sys.t_ref(), 298.15);
        assert_eq!(sys.p_ref(), 101.325);
    }

    /// Unknown component names produce a Thermo error naming the culprit.
    #[test]
    fn unknown_component_is_a_thermo_error() {
        let err = ThermoSystem::peng_robinson(&["benzene", "unobtainium"]).unwrap_err();
        assert!(matches!(err, StagesError::Thermo(_)));
        assert!(err.to_string().contains("unobtainium"));
    }

    /// Composition length must match the component count.
    #[test]
    fn dimension_mismatch_is_caught() {
        let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
        let err = sys.bubble_temperature(101.325, &[1.0]).unwrap_err();
        assert!(matches!(err, StagesError::Dimension(_)));
    }
}
