//! M1 integration tests: the McCabe–Thiele pipeline against real vle-thermo
//! thermodynamics — the two notebook systems, benzene–toluene (φ-φ,
//! Peng–Robinson) and methanol–water (γ-φ, van Laar from vle Chapter IV).
//!
//! The literature-pinned values live in the notebook's assertion cells
//! (CLAUDE.md, "Validation Strategy" layer 2); these tests pin the physical
//! invariants and textbook-magnitude windows so `cargo test` alone catches a
//! broken pipeline.

use stages_thermo::binary::mccabe_thiele::McCabeThieleSpec;
use stages_thermo::binary::{EquilibriumCurve, mccabe_thiele, rmin, total_reflux};
use stages_thermo::column::{BinaryColumn, CondenserKind, Feed};
use stages_thermo::thermo::ThermoSystem;

fn benzene_toluene_curve() -> EquilibriumCurve {
    let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
    EquilibriumCurve::from_thermo(&sys, 101.325, 101).unwrap()
}

/// Benzene–toluene at 1 atm: the textbook near-ideal pair. The relative
/// volatility must sit in the classic window (α ≈ 2.3–2.6 through the
/// composition range — S&H Table 7.3 territory), the temperature profile
/// must fall monotonically from toluene's boiling point to benzene's, and
/// the curve must be everywhere above the diagonal (no azeotrope).
#[test]
fn benzene_toluene_curve_is_textbook_shaped() {
    let curve = benzene_toluene_curve();
    for &x in &[0.1, 0.3, 0.5, 0.7, 0.9] {
        let a = curve.relative_volatility(x).unwrap();
        assert!(
            (2.1..2.8).contains(&a),
            "α({x}) = {a} outside the textbook window"
        );
        let y = curve.y_of_x(x).unwrap();
        assert!(y > x, "no azeotrope expected: y({x}) = {y}");
    }
    let t = curve.t_samples();
    assert!(
        t.windows(2).all(|w| w[1] < w[0]),
        "T(x) must fall monotonically toward pure benzene"
    );
    assert!((t[0] - 383.78).abs() < 3.0, "T(x=0) should be ~toluene bp");
    assert!(
        (t[t.len() - 1] - 353.24).abs() < 3.0,
        "T(x=1) should be ~benzene bp"
    );
}

/// The classic design: equimolar saturated-liquid feed to 95/5 products.
/// R_min must agree with the constant-α Underwood estimate computed from the
/// curve's own α at the feed pinch, and the full construction at 1.5·R_min
/// must land in the textbook stage-count range.
#[test]
fn benzene_toluene_classic_design() {
    let curve = benzene_toluene_curve();
    let r = rmin(&curve, 0.95, 0.05, 0.50, 1.0).unwrap();
    // Underwood's binary closed form with α evaluated at the pinch — the
    // curve is mildly non-constant-α, so allow a generous 10%.
    let alpha = curve.relative_volatility(r.pinch.0).unwrap();
    let underwood = 1.0 / (alpha - 1.0) * (0.95 / 0.50 - alpha * 0.05 / 0.50);
    assert!(
        (r.r_min - underwood).abs() / underwood < 0.10,
        "R_min = {} vs Underwood-at-pinch {underwood}",
        r.r_min
    );
    assert!(
        !r.tangent,
        "benzene–toluene is concave — feed pinch expected"
    );

    let res = mccabe_thiele(
        &curve,
        McCabeThieleSpec {
            x_distillate: 0.95,
            x_bottoms: 0.05,
            z_feed: 0.50,
            q: 1.0,
            reflux: 1.5 * r.r_min,
            murphree: 1.0,
            condenser: CondenserKind::Total,
        },
    )
    .unwrap();
    let nmin = total_reflux(&curve, 0.95, 0.05, 1.0).unwrap().n_min;
    assert!(
        res.n_stages > nmin && res.n_stages < 20.0,
        "N = {} (N_min = {nmin}) outside the plausible design range",
        res.n_stages
    );
    assert!(res.feed_stage > 1 && res.feed_stage < res.stages.len());

    // Gilliland sanity: at R = 1.5 R_min the reduced stage count
    // (N − Nmin)/(N + 1) classically falls near 0.2–0.55.
    let reduced = (res.n_stages - nmin) / (res.n_stages + 1.0);
    assert!(
        (0.1..0.6).contains(&reduced),
        "Gilliland-reduced N = {reduced}"
    );
}

/// Methanol–water with the vle Chapter IV van Laar parameters: strongly
/// non-ideal but azeotrope-free — y > x across the whole range, and the
/// curve must sit well above the constant-α-from-endpoints idealization at
/// low methanol (activity coefficients lift the dilute end).
#[test]
fn methanol_water_van_laar_curve() {
    let sys = ThermoSystem::van_laar(&["methanol", "water"], 0.5853, 0.3458).unwrap();
    let curve = EquilibriumCurve::from_thermo(&sys, 101.325, 101).unwrap();
    for &x in &[0.05, 0.2, 0.5, 0.8, 0.95] {
        let y = curve.y_of_x(x).unwrap();
        assert!(y > x, "methanol–water has no azeotrope: y({x}) = {y}");
    }
    // Boiling endpoints: water 373.12 K, methanol 337.7 K.
    let t = curve.t_samples();
    assert!(
        (t[0] - 373.12).abs() < 2.0,
        "T(x=0) = {} not ~water bp",
        t[0]
    );
    assert!(
        (t[t.len() - 1] - 337.7).abs() < 2.0,
        "T(x=1) = {} not ~methanol bp",
        t[t.len() - 1]
    );
    // A column on this curve solves.
    let res = mccabe_thiele(
        &curve,
        McCabeThieleSpec {
            x_distillate: 0.95,
            x_bottoms: 0.04,
            z_feed: 0.40,
            q: 1.0,
            reflux: 1.5,
            murphree: 1.0,
            condenser: CondenserKind::Total,
        },
    )
    .unwrap();
    assert!(res.n_stages.is_finite() && res.n_stages < 30.0);
}

/// The binary column model's material balances tie the curve work to product
/// rates.
#[test]
fn column_material_balance() {
    let col = BinaryColumn {
        pressure: 101.325,
        condenser: CondenserKind::Total,
        feed: Feed {
            rate: 100.0,
            z: 0.5,
            q: 1.0,
        },
        x_distillate: 0.95,
        x_bottoms: 0.05,
    };
    let d = col.distillate_rate().unwrap();
    let b = col.bottoms_rate().unwrap();
    assert!((d + b - 100.0).abs() < 1e-9);
    assert!((d * 0.95 + b * 0.05 - 50.0).abs() < 1e-9);
}
