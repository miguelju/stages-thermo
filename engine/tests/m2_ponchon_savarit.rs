//! M2 integration tests — the Ponchon–Savarit surface end-to-end, against
//! invariants and cross-checks with McCabe–Thiele.
//!
//! Unit-level coverage lives inside the modules; these exercise the public
//! crate API the way a downstream consumer (and the Python bindings) does, on
//! two systems: near-ideal benzene–toluene (φ-φ Peng–Robinson) and nonideal
//! methanol–water (γ-φ van Laar), plus the NRTL enthalpy path.

use stages_thermo::binary::equilibrium::EnthalpyCurve;
use stages_thermo::binary::mccabe_thiele::{McCabeThieleSpec, mccabe_thiele, rmin};
use stages_thermo::binary::ponchon_savarit::{PonchonSavaritSpec, ponchon_savarit};
use stages_thermo::column::CondenserKind;
use stages_thermo::thermo::{Phase, ThermoSystem};

fn ps_spec(reflux: f64) -> PonchonSavaritSpec {
    PonchonSavaritSpec {
        x_distillate: 0.95,
        x_bottoms: 0.05,
        z_feed: 0.50,
        q: 1.0,
        reflux,
        condenser: CondenserKind::Total,
    }
}

/// The saturated enthalpy curves are physically ordered and the phase enthalpy
/// path is self-consistent: `EnthalpyCurve` samples equal direct adapter calls.
#[test]
fn enthalpy_curve_matches_direct_adapter_calls() {
    let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
    let ec = EnthalpyCurve::from_thermo(&sys, 101.325, 101).unwrap();
    let eq = ec.equilibrium();
    // Pick an interior grid point and reproduce its enthalpies directly.
    let i = 40;
    let (x, y, t) = (eq.x_samples()[i], eq.y_samples()[i], eq.t_samples()[i]);
    let h_l = sys
        .phase_enthalpy(t, 101.325, &[x, 1.0 - x], Phase::Liquid)
        .unwrap();
    let h_v = sys
        .phase_enthalpy(t, 101.325, &[y, 1.0 - y], Phase::Vapor)
        .unwrap();
    assert!((ec.h_liq_samples()[i] - h_l).abs() < 1e-6);
    assert!((ec.h_vap_samples()[i] - h_v).abs() < 1e-6);
    assert!(h_v > h_l, "positive latent heat expected");
}

/// Near-ideal benzene–toluene: Ponchon–Savarit and McCabe–Thiele agree to
/// within about a stage (CMO is a good assumption here), and the M–T value
/// reproduces the M1-pinned N ≈ 12.22 at R = 1.5·R_min.
#[test]
fn benzene_toluene_consistency_with_mccabe_thiele() {
    let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
    let ec = EnthalpyCurve::from_thermo(&sys, 101.325, 401).unwrap();
    let r = 1.5 * rmin(ec.equilibrium(), 0.95, 0.05, 0.50, 1.0).unwrap().r_min;

    let ps = ponchon_savarit(&ec, ps_spec(r)).unwrap();
    let mt = mccabe_thiele(
        ec.equilibrium(),
        McCabeThieleSpec {
            x_distillate: 0.95,
            x_bottoms: 0.05,
            z_feed: 0.50,
            q: 1.0,
            reflux: r,
            murphree: 1.0,
            condenser: CondenserKind::Total,
        },
    )
    .unwrap();

    assert!(
        (mt.n_stages - 12.22).abs() < 0.1,
        "M–T N = {} should reproduce the M1-pinned 12.22",
        mt.n_stages
    );
    assert!(
        (ps.n_stages - mt.n_stages).abs() < 1.2,
        "P–S N = {} vs M–T N = {} should agree within ~1 stage",
        ps.n_stages,
        mt.n_stages
    );
    assert_eq!(ps.feed_stage, mt.feed_stage);
}

/// Nonideal methanol–water: the heat of mixing makes CMO imperfect, so
/// Ponchon–Savarit and McCabe–Thiele give *different* stage counts — the
/// CMO-error demonstration. Both remain finite and physical.
#[test]
fn methanol_water_cmo_error_gap() {
    let sys = ThermoSystem::van_laar(&["methanol", "water"], 0.5853, 0.3458).unwrap();
    let ec = EnthalpyCurve::from_thermo(&sys, 101.325, 401).unwrap();
    let r = 1.5 * rmin(ec.equilibrium(), 0.95, 0.05, 0.50, 1.0).unwrap().r_min;

    let ps = ponchon_savarit(&ec, ps_spec(r)).unwrap();
    let mt = mccabe_thiele(
        ec.equilibrium(),
        McCabeThieleSpec {
            x_distillate: 0.95,
            x_bottoms: 0.05,
            z_feed: 0.50,
            q: 1.0,
            reflux: r,
            murphree: 1.0,
            condenser: CondenserKind::Total,
        },
    )
    .unwrap();

    assert!(ps.n_stages.is_finite() && mt.n_stages.is_finite());
    // The two methods disagree here (unlike near-ideal benzene–toluene).
    assert!(
        (ps.n_stages - mt.n_stages).abs() > 0.1,
        "expected a CMO-error gap for methanol–water, got P–S {} vs M–T {}",
        ps.n_stages,
        mt.n_stages
    );
}

/// Mass + energy closure on a converged construction (an invariant on every
/// P–S result): the lever-rule product split and the difference-point energy
/// balance both close.
#[test]
fn mass_and_energy_closure() {
    let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
    let ec = EnthalpyCurve::from_thermo(&sys, 101.325, 201).unwrap();
    let res = ponchon_savarit(&ec, ps_spec(2.0)).unwrap();

    let (xd, xb, zf): (f64, f64, f64) = (0.95, 0.05, 0.50);
    let d = (zf - xb) / (xd - xb);
    let b = (xd - zf) / (xd - xb);
    // Component mass balance: D·x_D + B·x_B = F·z_F (per mole feed).
    assert!((d * xd + b * xb - zf).abs() < 1e-12);

    // Energy balance around the whole column: F·h_F + Q_R = D·h_D + B·h_B + Q_C.
    let h_f = res.feed_point.1;
    let h_d = ec.h_liquid_of_x(xd).unwrap();
    let h_b = ec.h_liquid_of_x(xb).unwrap();
    let lhs = h_f + res.q_reboiler;
    let rhs = d * h_d + b * h_b + res.q_condenser;
    assert!(
        (lhs - rhs).abs() < 1e-6 * (lhs.abs() + rhs.abs() + 1.0),
        "energy not closed: {lhs} vs {rhs}"
    );
    assert!(res.q_condenser > 0.0 && res.q_reboiler > 0.0);
}

/// The NRTL enthalpy path builds a valid H–x–y curve and a converging
/// construction (exercises vle-thermo 0.11's NRTL `alpha` route end-to-end).
#[test]
fn nrtl_enthalpy_curve_constructs() {
    // Ethanol–water NRTL (Aspen-style, kJ/kmol): Δg₁₂ = −458.7, Δg₂₁ = 5574,
    // α = 0.303.
    let sys = ThermoSystem::nrtl(&["ethanol", "water"], -458.7, 5574.0, 0.303).unwrap();
    let ec = EnthalpyCurve::from_thermo(&sys, 101.325, 201).unwrap();
    // Vapor above liquid everywhere in the interior.
    let (hl, hv) = (ec.h_liq_samples(), ec.h_vap_samples());
    for i in 1..hl.len() - 1 {
        assert!(hv[i] > hl[i]);
    }
    // A modest separation converges (ethanol–water azeotrope near x≈0.9, so
    // keep x_D well below it).
    let res = ponchon_savarit(
        &ec,
        PonchonSavaritSpec {
            x_distillate: 0.70,
            x_bottoms: 0.05,
            z_feed: 0.30,
            q: 1.0,
            reflux: 3.0,
            condenser: CondenserKind::Total,
        },
    )
    .unwrap();
    assert!(res.n_stages.is_finite() && res.n_stages > 1.0);
}
