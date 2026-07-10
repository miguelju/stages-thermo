//! The Ponchon–Savarit graphical method for binary distillation — rung 2 of
//! the pedagogical ladder (PLAN §4), where the **energy balance enters** and
//! constant-molal-overflow (CMO) stops being assumed.
//!
//! ## Method and reference
//!
//! Ponchon, M. Étude graphique de la distillation fractionnée. *La Technique
//! Moderne* **1921**, 13, 20 and 55. Savarit, R. *Arts et Métiers* **1922**.
//! Implementation-level treatment: Seader, J. D.; Henley, E. J.; Roper, D. K.
//! *Separation Process Principles*, Ch. 7 (the enthalpy–composition
//! construction). Full derivation: `docs/theory/ponchon-savarit.md`.
//!
//! ## The idea, versus McCabe–Thiele
//!
//! McCabe–Thiele steps a staircase between the equilibrium curve and straight
//! *operating lines* on an x–y square, which is exact only under CMO (equal
//! molar latent heats, negligible heat of mixing). Ponchon–Savarit works on
//! the **enthalpy–composition (H–x–y) diagram** — saturated-liquid enthalpy
//! `h_L(x)` and saturated-vapor enthalpy `H_V(y)` versus composition — and
//! closes the energy balance *exactly* through two **difference points**
//! (poles):
//!
//! ```text
//! Δ_D = (x_D, Q'_D)   top pole, above the diagram; Q'_D = h_D + Q_C/D
//! Δ_B = (x_B, Q'_B)   bottom pole, below the diagram; Q'_B = h_B − Q_R/B
//! ```
//!
//! The reflux ratio sets the height of `Δ_D`:
//! `R = L/D = (Q'_D − H_V1)/(H_V1 − h_L0)` (ratio of vertical segments on the
//! H–x diagram; S&H Ch. 7), so `Q'_D = H_V1 + R·(H_V1 − h_L0)` with `H_V1` the
//! saturated-vapor enthalpy at `y = x_D` and `h_L0` the saturated-liquid
//! reflux enthalpy at `x_D`. The overall balance makes `Δ_D`, the feed point
//! `F = (z_F, h_F)`, and `Δ_B` **collinear**, which fixes `Δ_B`.
//!
//! ## Stepping (code ↔ diagram)
//!
//! The staircase mirrors [`super::mccabe_thiele`], with two substitutions:
//!
//! - the **horizontal step to the equilibrium curve** becomes a **tie line**:
//!   the liquid `x_n` leaving stage `n` is in equilibrium with the vapor `y_n`
//!   (`x_n = x*(y_n)` — the *same* inverse the McCabe–Thiele stepper uses;
//!   on the H–x–y diagram this is the tie line joining `(x_n, h_L(x_n))` to
//!   `(y_n, H_V(y_n))`);
//! - the **vertical step to the operating line** becomes a **pole line**: the
//!   passing streams `L_n` (liquid, `x_n`) and `V_{n+1}` (vapor from below)
//!   plus the section's pole `Δ` are collinear, so `y_{n+1}` is where the line
//!   through `Δ` and `(x_n, h_L(x_n))` cuts the saturated-vapor curve.
//!
//! Stages are stepped **top-down** (stage 1 at the top; repo convention). The
//! feed stage is the first stage stepped on the stripping pole `Δ_B`; the
//! switch happens when the stage liquid passes the feed composition `z_F`
//! (which, for a saturated-liquid feed, coincides with McCabe–Thiele's
//! optimal-feed rule). A partial reboiler is an equilibrium stage and is
//! counted.
//!
//! For a **near-ideal, equal-latent-heat** system (benzene–toluene) the pole
//! lines reproduce the CMO operating lines, so Ponchon–Savarit and
//! McCabe–Thiele return the same stage count — the M2 consistency check. Where
//! the heat of mixing is large (methanol–water, and above all ammonia–water)
//! the two diverge, and Ponchon–Savarit is the one that is right.
//!
//! ## Units
//!
//! Enthalpies **kJ/kmol**, compositions dimensionless, pressure **kPa**.
//! Duties are reported **per mole of feed** (kJ/kmol of feed).

use crate::binary::equilibrium::EnthalpyCurve;
use crate::binary::mccabe_thiele::StagePoint;
use crate::column::CondenserKind;
use crate::types::{Result, StagesError};

/// Hard cap on stepped stages: hitting it means the staircase pinched
/// (`R` at or below the Ponchon–Savarit minimum) rather than a real column.
const MAX_STAGES: usize = 500;
/// Minimum per-stage composition progress before the stepper declares a pinch.
const PINCH_PROGRESS: f64 = 1e-10;

/// Design specification for a Ponchon–Savarit construction.
///
/// The composition/feed/reflux fields mirror [`super::mccabe_thiele::McCabeThieleSpec`]
/// (minus Murphree efficiency — Ponchon–Savarit here counts theoretical
/// stages) so the two methods can be run against the *same* spec and their
/// stage counts compared directly.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct PonchonSavaritSpec {
    /// Distillate light-component mole fraction `x_D`.
    pub x_distillate: f64,
    /// Bottoms light-component mole fraction `x_B`.
    pub x_bottoms: f64,
    /// Feed light-component mole fraction `z_F`.
    pub z_feed: f64,
    /// Feed thermal condition `q` (1 = saturated liquid, 0 = saturated vapor).
    /// The feed enthalpy is `h_F = q·h_L(z_F) + (1−q)·H_V(z_F)`.
    pub q: f64,
    /// External reflux ratio `R = L/D`.
    pub reflux: f64,
    /// Condenser kind. **v1 supports [`CondenserKind::Total`]**; a partial
    /// condenser (an extra equilibrium stage above tray 1) is rejected — see
    /// the module docs.
    pub condenser: CondenserKind,
}

impl PonchonSavaritSpec {
    fn validate(&self) -> Result<()> {
        for (name, v) in [
            ("x_D", self.x_distillate),
            ("x_B", self.x_bottoms),
            ("z_F", self.z_feed),
        ] {
            if !(0.0 < v && v < 1.0) {
                return Err(StagesError::Dimension(format!(
                    "{name} must be strictly inside (0, 1), got {v}"
                )));
            }
        }
        if !(self.x_bottoms < self.z_feed && self.z_feed < self.x_distillate) {
            return Err(StagesError::Dimension(format!(
                "specs must satisfy x_B < z_F < x_D, got x_B = {}, z_F = {}, x_D = {}",
                self.x_bottoms, self.z_feed, self.x_distillate
            )));
        }
        // NaN-safe positivity guard (a bare `> 0.0` inverted with `!` trips
        // clippy's neg_cmp_op_on_partial_ord; the NaN check is explicit).
        if self.reflux.is_nan() || self.reflux <= 0.0 {
            return Err(StagesError::Dimension(format!(
                "reflux ratio must be positive, got {}",
                self.reflux
            )));
        }
        if matches!(self.condenser, CondenserKind::Partial) {
            return Err(StagesError::Infeasible(
                "Ponchon–Savarit v1 supports a total condenser only; a partial condenser adds an \
                 equilibrium stage above the top tray that is not yet modelled"
                    .into(),
            ));
        }
        Ok(())
    }
}

/// The full Ponchon–Savarit construction — every geometric object on the
/// H–x–y diagram plus the stage count and duties, so notebooks and the MCP
/// layer can query or redraw anything (PLAN §8: rich result objects, never
/// bare numbers).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct PonchonSavaritResult {
    /// The specification this construction was built from.
    pub spec: PonchonSavaritSpec,
    /// Number of theoretical stages, fractional (the last stage is counted by
    /// the fraction of its step needed to reach `x_B`). Includes the partial
    /// reboiler.
    pub n_stages: f64,
    /// The stage corners `(x_n, y_n)` — liquid and vapor leaving each stage,
    /// top-down. (Reuses [`StagePoint`]; the enthalpy coordinates live in
    /// [`Self::tie_lines`].)
    pub stages: Vec<StagePoint>,
    /// Optimal feed stage, 1-based from the top (first stage stepped on `Δ_B`).
    pub feed_stage: usize,
    /// Top difference point `Δ_D = (x_D, Q'_D)`; ordinate in **kJ/kmol**.
    pub delta_d: (f64, f64),
    /// Bottom difference point `Δ_B = (x_B, Q'_B)`; ordinate in **kJ/kmol**.
    pub delta_b: (f64, f64),
    /// Feed point `F = (z_F, h_F)`; ordinate in **kJ/kmol**.
    pub feed_point: (f64, f64),
    /// Condenser duty **per mole of feed**, `Q_C/F` in **kJ/kmol** (> 0 = heat
    /// removed).
    pub q_condenser: f64,
    /// Reboiler duty **per mole of feed**, `Q_R/F` in **kJ/kmol** (> 0 = heat
    /// added).
    pub q_reboiler: f64,
    /// The tie lines actually stepped: each is
    /// `((x_n, h_L(x_n)), (y_n, H_V(y_n)))`, ready to plot on the H–x–y frame.
    pub tie_lines: Vec<((f64, f64), (f64, f64))>,
}

/// Find the vapor composition `y` where the **pole line** — the straight line
/// through the pole `Δ` and the liquid point `(x_l, h_l)` — cuts the
/// saturated-vapor curve, searching the composition window `(x_l, y_upper]`.
///
/// Collinearity of `Δ = (x_p, h_p)`, `L = (x_l, h_l)` and a candidate vapor
/// point `V = (y, H_V(y))` is the 2-D cross product
/// `g(y) = (x_l − x_p)·(H_V(y) − h_p) − (h_l − h_p)·(y − x_p) = 0`.
/// `H_V` is piecewise-linear on the vapor grid, so `g` is piecewise-linear and
/// the root inside a bracketing segment is exact (mirrors the sign-change scan
/// in [`super::mccabe_thiele`]'s `q_line_curve_intersection`).
fn pole_line_vapor_intersection(
    curve: &EnthalpyCurve,
    pole: (f64, f64),
    liq: (f64, f64),
    y_upper: f64,
) -> Result<f64> {
    let (xp, hp) = pole;
    let (xl, hl) = liq;
    let g = |y: f64, hv: f64| (xl - xp) * (hv - hp) - (hl - hp) * (y - xp);

    let ys = curve.equilibrium().y_samples();
    let hvs = curve.h_vap_samples();
    // Anchor the scan just above the liquid composition; walk the vapor grid
    // upward to `y_upper` and take the first sign change of `g`.
    let mut prev_y = xl;
    let mut prev_g = g(xl, curve.h_vapor_of_y(xl)?);
    for j in 0..ys.len() {
        let yj = ys[j];
        if yj <= xl {
            continue;
        }
        if yj > y_upper + 1e-9 {
            break;
        }
        let gj = g(yj, hvs[j]);
        if (prev_g < 0.0 && gj >= 0.0) || (prev_g > 0.0 && gj <= 0.0) {
            // Linear root between prev_y and yj (both grid-exact).
            let denom = prev_g - gj;
            let t = if denom.abs() > 0.0 {
                prev_g / denom
            } else {
                0.0
            };
            return Ok(prev_y + t * (yj - prev_y));
        }
        prev_y = yj;
        prev_g = gj;
    }
    Err(StagesError::Convergence(format!(
        "pole line from Δ = ({xp:.4}, {hp:.1}) through liquid x = {xl:.4} does not cut the \
         saturated-vapor curve below y = {y_upper:.4} — R too close to the minimum (pinched)"
    )))
}

/// Run the full Ponchon–Savarit construction for a design specification.
///
/// # Arguments
/// * `curve` — the binary enthalpy–composition curve at column pressure
/// * `spec` — the design spec (compositions, `q`, `R`, condenser)
///
/// # Errors
/// - [`StagesError::Infeasible`] for an unsupported spec (partial condenser)
///   or an infeasible separation.
/// - [`StagesError::Convergence`] if `R` is at or below the Ponchon–Savarit
///   minimum reflux (the staircase pinches before reaching `x_B`).
pub fn ponchon_savarit(
    curve: &EnthalpyCurve,
    spec: PonchonSavaritSpec,
) -> Result<PonchonSavaritResult> {
    spec.validate()?;
    let (xd, xb, zf, r) = (spec.x_distillate, spec.x_bottoms, spec.z_feed, spec.reflux);

    // --- Difference points (poles). ---
    let h_l0 = curve.h_liquid_of_x(xd)?; // saturated-liquid reflux at x_D
    let h_v1 = curve.h_vapor_of_y(xd)?; // saturated vapor to the total condenser (y_1 = x_D)
    // Top pole ordinate from the reflux ratio (module docs).
    let q_prime_d = h_v1 + r * (h_v1 - h_l0);
    let delta_d = (xd, q_prime_d);

    // Feed point (z_F, h_F); h_F from the thermal condition q via the curve.
    let h_lf = curve.h_liquid_of_x(zf)?;
    let h_vf = curve.h_vapor_of_y(zf)?;
    let h_f = spec.q * h_lf + (1.0 - spec.q) * h_vf;
    let feed_point = (zf, h_f);

    // Δ_D, F, Δ_B collinear (overall balance): extrapolate the Δ_D–F line to
    // x = x_B to place the bottom pole.
    if (zf - xd).abs() < 1e-12 {
        return Err(StagesError::Infeasible(
            "z_F equals x_D — degenerate feed/product specification".into(),
        ));
    }
    let slope_df = (h_f - q_prime_d) / (zf - xd);
    let q_prime_b = q_prime_d + slope_df * (xb - xd);
    let delta_b = (xb, q_prime_b);

    // --- Duties per mole of feed (basis F = 1; lever-rule split). ---
    let d_over_f = (zf - xb) / (xd - xb);
    let b_over_f = (xd - zf) / (xd - xb);
    let h_d = h_l0; // total condenser: distillate is saturated liquid at x_D
    let h_b = curve.h_liquid_of_x(xb)?; // saturated-liquid bottoms
    let q_condenser = d_over_f * (q_prime_d - h_d);
    let q_reboiler = b_over_f * (h_b - q_prime_b);

    // --- Step the staircase top-down. ---
    let mut stages: Vec<StagePoint> = Vec::new();
    let mut tie_lines: Vec<((f64, f64), (f64, f64))> = Vec::new();
    let mut y_n = xd; // vapor to the total condenser: y_1 = x_D
    let mut x_prev = xd;
    let mut in_stripping = false;
    let mut feed_stage = 0usize;
    let mut n_stages = 0.0_f64;

    for n in 1..=MAX_STAGES {
        // Tie line (equilibrium): liquid leaving stage n is in equilibrium
        // with the vapor y_n leaving it — the same inverse McCabe–Thiele uses.
        let x_n = curve.equilibrium().x_of_y(y_n)?;
        let h_ln = curve.h_liquid_of_x(x_n)?;
        let h_vn = curve.h_vapor_of_y(y_n)?;
        tie_lines.push(((x_n, h_ln), (y_n, h_vn)));
        stages.push(StagePoint {
            index: n,
            x: x_n,
            y: y_n,
        });

        if x_n <= xb {
            // Fractional final stage: fraction of this step needed to hit x_B.
            let frac = if x_prev - x_n > 0.0 {
                ((x_prev - xb) / (x_prev - x_n)).clamp(0.0, 1.0)
            } else {
                1.0
            };
            n_stages = (n - 1) as f64 + frac;
            break;
        }
        if x_prev - x_n < PINCH_PROGRESS {
            return Err(StagesError::Convergence(format!(
                "staircase pinched at x = {x_n:.6} after {n} stages (R too close to the \
                 Ponchon–Savarit minimum reflux)"
            )));
        }

        // Optimal-feed switch: once the stage liquid passes the feed
        // composition, subsequent steps ride the stripping pole Δ_B.
        if !in_stripping && x_n <= zf {
            in_stripping = true;
            feed_stage = n;
        }
        let pole = if in_stripping { delta_b } else { delta_d };

        // Pole line (operating): vapor V_{n+1} entering stage n from below.
        let y_next = pole_line_vapor_intersection(curve, pole, (x_n, h_ln), y_n)?;
        x_prev = x_n;
        y_n = y_next;

        if n == MAX_STAGES {
            return Err(StagesError::Convergence(format!(
                "exceeded {MAX_STAGES} stages without reaching x_B = {xb} — spec pinched"
            )));
        }
    }
    if feed_stage == 0 {
        // Reached x_B without leaving the rectifying section — feed at the base.
        feed_stage = stages.len();
    }

    Ok(PonchonSavaritResult {
        spec,
        n_stages,
        stages,
        feed_stage,
        delta_d,
        delta_b,
        feed_point,
        q_condenser,
        q_reboiler,
        tie_lines,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary::mccabe_thiele::{McCabeThieleSpec, mccabe_thiele, rmin};
    use crate::thermo::ThermoSystem;

    fn benzene_toluene_curve() -> EnthalpyCurve {
        let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
        EnthalpyCurve::from_thermo(&sys, 101.325, 201).unwrap()
    }

    fn spec(reflux: f64) -> PonchonSavaritSpec {
        PonchonSavaritSpec {
            x_distillate: 0.95,
            x_bottoms: 0.05,
            z_feed: 0.50,
            q: 1.0,
            reflux,
            condenser: CondenserKind::Total,
        }
    }

    /// The construction runs end-to-end: a finite fractional stage count, a
    /// feed stage strictly inside the column, and the top/bottom poles bracket
    /// the diagram (Δ_D above the vapor curve, Δ_B below the liquid curve).
    #[test]
    fn construction_is_well_formed() {
        let curve = benzene_toluene_curve();
        let res = ponchon_savarit(&curve, spec(2.0)).unwrap();
        assert!(res.n_stages > 1.0 && res.n_stages < 40.0);
        assert!(res.feed_stage >= 1 && res.feed_stage <= res.stages.len());
        // Δ_D sits above the saturated-vapor enthalpy at x_D; Δ_B below the
        // saturated-liquid enthalpy at x_B.
        assert!(res.delta_d.1 > curve.h_vapor_of_y(0.95).unwrap());
        assert!(res.delta_b.1 < curve.h_liquid_of_x(0.05).unwrap());
    }

    /// Energy closure: the pole ordinates are consistent with the reported
    /// duties, and the overall balance `D·Q'_D − B·Q'_B` (per mole feed) equals
    /// `Q_C − Q_R + F·h_F` — the difference-point construction's balance.
    #[test]
    fn energy_balance_closes() {
        let curve = benzene_toluene_curve();
        let res = ponchon_savarit(&curve, spec(2.0)).unwrap();
        let d = (0.50 - 0.05) / (0.95 - 0.05);
        let b = (0.95 - 0.50) / (0.95 - 0.05);
        // Feed, distillate and bottoms enthalpy balance closed by the duties:
        // F·h_F + Q_R = D·h_D + B·h_B + Q_C  (per mole of feed, F = 1).
        let h_f = res.feed_point.1;
        let h_d = curve.h_liquid_of_x(0.95).unwrap();
        let h_b = curve.h_liquid_of_x(0.05).unwrap();
        let lhs = h_f + res.q_reboiler;
        let rhs = d * h_d + b * h_b + res.q_condenser;
        assert!(
            (lhs - rhs).abs() < 1e-6 * (lhs.abs() + rhs.abs() + 1.0),
            "energy balance not closed: LHS {lhs} vs RHS {rhs} kJ/kmol"
        );
        // Both duties physically signed.
        assert!(res.q_condenser > 0.0, "Q_C = {}", res.q_condenser);
        assert!(res.q_reboiler > 0.0, "Q_R = {}", res.q_reboiler);
    }

    /// **The M2 consistency check**: for near-ideal benzene–toluene (nearly
    /// equal molar latent heats, tiny heat of mixing), CMO is a good
    /// assumption, so Ponchon–Savarit and McCabe–Thiele agree to within about
    /// a stage on the same spec.
    #[test]
    fn agrees_with_mccabe_thiele_on_benzene_toluene() {
        let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
        let ec = EnthalpyCurve::from_thermo(&sys, 101.325, 401).unwrap();
        let r = 1.5 * rmin(ec.equilibrium(), 0.95, 0.05, 0.50, 1.0).unwrap().r_min;

        let ps = ponchon_savarit(&ec, spec(r)).unwrap();
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
            (ps.n_stages - mt.n_stages).abs() < 1.2,
            "P–S N = {} vs M–T N = {} should agree within ~1 stage for near-ideal \
             benzene–toluene",
            ps.n_stages,
            mt.n_stages
        );
    }

    /// More reflux ⇒ fewer stages (approaching N_min), same as McCabe–Thiele.
    #[test]
    fn stages_decrease_with_reflux() {
        let curve = benzene_toluene_curve();
        let n_low = ponchon_savarit(&curve, spec(1.5)).unwrap().n_stages;
        let n_high = ponchon_savarit(&curve, spec(4.0)).unwrap().n_stages;
        assert!(
            n_high < n_low,
            "N should fall with reflux: N(1.5) = {n_low}, N(4.0) = {n_high}"
        );
    }

    /// A partial condenser is rejected in v1 with a clear error.
    #[test]
    fn partial_condenser_rejected() {
        let curve = benzene_toluene_curve();
        let s = PonchonSavaritSpec {
            condenser: CondenserKind::Partial,
            ..spec(2.0)
        };
        assert!(matches!(
            ponchon_savarit(&curve, s),
            Err(StagesError::Infeasible(_))
        ));
    }
}
