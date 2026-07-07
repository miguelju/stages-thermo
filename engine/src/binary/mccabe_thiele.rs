//! The McCabe–Thiele graphical method for binary distillation — rung 1 of
//! the pedagogical ladder (PLAN §4).
//!
//! ## Method and reference
//!
//! McCabe, W. L.; Thiele, E. W. Graphical Design of Fractionating Columns.
//! *Ind. Eng. Chem.* **1925**, 17 (6), 605–611. Implementation-level
//! treatment: Seader, J. D.; Henley, E. J.; Roper, D. K. *Separation Process
//! Principles*, Ch. 7 (equation numbers below are S&H's).
//!
//! ## Symbol map (code ↔ textbook)
//!
//! ```text
//! spec.reflux        R = L/D          external reflux ratio
//! spec.x_distillate  x_D              distillate light-component fraction
//! spec.x_bottoms     x_B              bottoms light-component fraction
//! spec.z_feed        z_F              feed light-component fraction
//! spec.q             q                feed thermal condition (L̄ = L + qF)
//! spec.murphree      E_MV             Murphree vapor efficiency
//! ```
//!
//! The three construction lines (light-component balances around column
//! sections):
//!
//! ```text
//! rectifying (S&H eq. 7-9):   y = R/(R+1) · x + x_D/(R+1)
//! q-line     (S&H eq. 7-26):  y = q/(q−1) · x − z_F/(q−1)     (vertical at q = 1)
//! stripping  (S&H eq. 7-12):  through (x_B, x_B) and the rectifying ∩ q-line point
//! ```
//!
//! Stages are stepped **top-down** (stage 1 at the top; repo convention):
//! from a point on an operating line, move horizontally to the equilibrium
//! curve (one theoretical stage), then vertically back to the operating
//! line. The step where the staircase crosses the operating-line
//! intersection is the **optimal feed stage** (S&H §7.2.3). A partial
//! reboiler is an equilibrium stage and is included in the count; a total
//! condenser is not (a partial condenser adds one — see
//! [`CondenserKind`]).
//!
//! ## Minimum reflux and pinches
//!
//! R_min is found geometrically, tangent pinches included, by slope
//! extremization over the sampled curve rather than assuming the pinch sits
//! at the feed point:
//!
//! - **Rectifying side**: the operating line anchored at `(x_D, x_D)` must
//!   stay *below* the equilibrium curve on `[x_q*, x_D]` (where `(x_q*,
//!   y_q*)` is the q-line ∩ curve point). The limiting slope is
//!   `max (x_D − y_e)/(x_D − x_e)` over curve samples — the maximum is the
//!   feed pinch for concave curves and the tangent pinch when the curve
//!   bulges (e.g. ethanol–water).
//! - **Stripping side**: the line anchored at `(x_B, x_B)` must stay *below*
//!   the curve on `[x_B, x_q*]`; the limiting slope is
//!   `min (y_e − x_B)/(x_e − x_B)`, converted to an equivalent R through the
//!   feed-section balances `L̄ = L + qF`, `V̄ = V − (1−q)F`.
//!
//! R_min is the larger of the two; for a normal concave curve both give the
//! same feed pinch (asserted in the tests). R_min uses the *true*
//! equilibrium curve — the classical construction — even when stepping later
//! applies a Murphree efficiency.
//!
//! ## Murphree efficiency
//!
//! Stepping with `E_MV < 1` replaces the equilibrium curve by the
//! pseudo-curve `y_eff(x) = y_op(x) + E_MV · (y*(x) − y_op(x))` between the
//! operating line and the true curve (S&H §7.4). The pseudo-curve is applied
//! to every stage including the reboiler — a slightly conservative
//! simplification (a real partial reboiler is a true equilibrium stage).

use crate::binary::equilibrium::EquilibriumCurve;
use crate::column::CondenserKind;
use crate::types::{Result, StagesError};

/// Hard cap on stepped stages: hitting it means the staircase pinched
/// (R ≤ R_min or an unreachable spec) rather than a real 500-stage column.
const MAX_STAGES: usize = 500;
/// Minimum per-stage composition progress before the stepper declares a
/// pinch.
const PINCH_PROGRESS: f64 = 1e-10;

/// Design specification for a McCabe–Thiele construction.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct McCabeThieleSpec {
    /// Distillate light-component mole fraction `x_D`.
    pub x_distillate: f64,
    /// Bottoms light-component mole fraction `x_B`.
    pub x_bottoms: f64,
    /// Feed light-component mole fraction `z_F`.
    pub z_feed: f64,
    /// Feed thermal condition `q` (1 = saturated liquid, 0 = saturated
    /// vapor).
    pub q: f64,
    /// External reflux ratio `R = L/D`.
    pub reflux: f64,
    /// Murphree vapor efficiency `E_MV` in (0, 1]; 1 = theoretical stages.
    pub murphree: f64,
    /// Condenser kind — a partial condenser is itself an equilibrium stage.
    pub condenser: CondenserKind,
}

impl McCabeThieleSpec {
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
        // NaN-safe positivity guard (a plain `> 0.0` inverted with `!` trips
        // clippy's neg_cmp_op_on_partial_ord; the NaN check is explicit).
        if self.reflux.is_nan() || self.reflux <= 0.0 {
            return Err(StagesError::Dimension(format!(
                "reflux ratio must be positive, got {}",
                self.reflux
            )));
        }
        if !(0.0 < self.murphree && self.murphree <= 1.0) {
            return Err(StagesError::Dimension(format!(
                "Murphree efficiency must be in (0, 1], got {}",
                self.murphree
            )));
        }
        Ok(())
    }
}

/// A straight construction line `y = slope · x + intercept`.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct Line {
    /// Slope (dimensionless).
    pub slope: f64,
    /// Intercept at x = 0 (mole fraction).
    pub intercept: f64,
}

impl Line {
    fn y_at(&self, x: f64) -> f64 {
        self.slope * x + self.intercept
    }
}

/// One theoretical (or pseudo-, when `E_MV < 1`) stage of the staircase.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct StagePoint {
    /// Stage number, 1-based from the top.
    pub index: usize,
    /// Liquid light-component fraction leaving the stage.
    pub x: f64,
    /// Vapor light-component fraction leaving the stage.
    pub y: f64,
}

/// The result of a minimum-reflux analysis.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct RminResult {
    /// Minimum external reflux ratio `R_min`.
    pub r_min: f64,
    /// The controlling pinch point `(x, y)` on the equilibrium curve.
    pub pinch: (f64, f64),
    /// `true` if the controlling pinch is a tangent pinch (an interior curve
    /// point) rather than the feed (q-line) pinch.
    pub tangent: bool,
    /// The q-line ∩ equilibrium-curve point `(x_q*, y_q*)`.
    pub feed_point: (f64, f64),
}

/// The full McCabe–Thiele construction — every geometric object on the
/// diagram plus the stage count, so notebooks and the MCP layer can query or
/// redraw anything (PLAN §8: rich result objects, never bare numbers).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct McCabeThieleResult {
    /// The specification this construction was built from.
    pub spec: McCabeThieleSpec,
    /// Number of theoretical stages, fractional (the last stage is counted
    /// by the fraction of its step needed to reach `x_B`). Includes the
    /// partial reboiler, and the condenser too when it is partial.
    pub n_stages: f64,
    /// The stage corners `(x_n, y_n)` on the (pseudo-)equilibrium curve,
    /// top-down.
    pub stages: Vec<StagePoint>,
    /// Optimal feed stage, 1-based from the top (the stage whose step
    /// crosses the operating-line intersection).
    pub feed_stage: usize,
    /// Minimum-reflux analysis for this spec (computed on the true curve).
    pub rmin: RminResult,
    /// Rectifying operating line (S&H eq. 7-9).
    pub rectifying: Line,
    /// Stripping operating line (S&H eq. 7-12).
    pub stripping: Line,
    /// Intersection of the two operating lines (on the q-line).
    pub intersection: (f64, f64),
    /// The full staircase polyline, starting at `(x_D, x_D)`, alternating
    /// horizontal and vertical segments — ready to plot.
    pub staircase: Vec<(f64, f64)>,
}

/// The total-reflux construction: operating lines collapse onto the y = x
/// diagonal and the stage count is the minimum, `N_min` (S&H §7.2.4 —
/// Fenske's equation is the algebraic constant-α equivalent).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct TotalRefluxResult {
    /// Minimum number of theoretical stages `N_min` (fractional, reboiler
    /// included).
    pub n_min: f64,
    /// Stage corners, top-down.
    pub stages: Vec<StagePoint>,
    /// The staircase polyline against the diagonal.
    pub staircase: Vec<(f64, f64)>,
}

/// Intersection of the rectifying line and the q-line — the point the
/// stripping line must also pass through (S&H §7.2.2).
///
/// For `q = 1` the q-line is vertical at `z_F`, so the intersection is
/// `(z_F, rectifying(z_F))`; otherwise the two lines are solved directly.
fn operating_intersection(rect: Line, z_feed: f64, q: f64) -> Result<(f64, f64)> {
    if (q - 1.0).abs() < 1e-12 {
        return Ok((z_feed, rect.y_at(z_feed)));
    }
    let mq = q / (q - 1.0);
    // Solve slope·x + intercept = mq·(x − z_F) + z_F.
    let denom = rect.slope - mq;
    if denom.abs() < 1e-12 {
        return Err(StagesError::Infeasible(
            "q-line is parallel to the rectifying line (degenerate q/R combination)".into(),
        ));
    }
    let x = (z_feed - mq * z_feed - rect.intercept) / denom;
    Ok((x, rect.y_at(x)))
}

/// The q-line ∩ equilibrium-curve point `(x_q*, y_q*)` — the feed pinch
/// candidate.
///
/// Solved on the sampled curve: special-cased for the vertical (`q = 1`) and
/// horizontal (`q = 0`) lines, otherwise by scanning the grid outward from
/// `z_F` for the first sign change of `y*(x) − y_qline(x)` and solving the
/// linear–linear crossing inside that segment.
fn q_line_curve_intersection(curve: &EquilibriumCurve, z_feed: f64, q: f64) -> Result<(f64, f64)> {
    if (q - 1.0).abs() < 1e-12 {
        return Ok((z_feed, curve.y_of_x(z_feed)?));
    }
    if q.abs() < 1e-12 {
        return Ok((curve.x_of_y(z_feed)?, z_feed));
    }
    let mq = q / (q - 1.0);
    let line = |x: f64| mq * (x - z_feed) + z_feed;
    // Walk the grid away from z_F on the side where the q-line climbs above
    // the diagonal: right for a subcooled-liquid slope (mq > 1), left
    // otherwise (mq < 1 covers 0 < q < 1 and superheated vapor).
    let xs = curve.x_samples();
    let f = |x: f64| -> Result<f64> { Ok(curve.y_of_x(x)? - line(x)) };
    let mut prev_x = z_feed;
    let mut prev_f = f(prev_x)?; // > 0: the curve starts above the q-line at z_F
    let indices: Vec<usize> = if mq > 1.0 {
        (0..xs.len()).filter(|&i| xs[i] > z_feed).collect()
    } else {
        (0..xs.len()).rev().filter(|&i| xs[i] < z_feed).collect()
    };
    for i in indices {
        let xi = xs[i];
        let fi = f(xi)?;
        if prev_f >= 0.0 && fi < 0.0 {
            // Linear root between prev_x and xi.
            let t = prev_f / (prev_f - fi);
            let xr = prev_x + t * (xi - prev_x);
            return Ok((xr, line(xr)));
        }
        prev_x = xi;
        prev_f = fi;
    }
    Err(StagesError::Infeasible(format!(
        "the q-line (q = {q}) never meets the equilibrium curve — check z_F = {z_feed} and q"
    )))
}

/// Minimum reflux ratio by pinch detection on the true equilibrium curve,
/// tangent pinches included (module docs, "Minimum reflux and pinches").
///
/// # Arguments
/// * `curve` — the equilibrium curve
/// * `x_distillate`, `x_bottoms`, `z_feed` — product/feed light-component
///   mole fractions, `x_B < z_F < x_D`
/// * `q` — feed thermal condition
pub fn rmin(
    curve: &EquilibriumCurve,
    x_distillate: f64,
    x_bottoms: f64,
    z_feed: f64,
    q: f64,
) -> Result<RminResult> {
    let (xd, xb) = (x_distillate, x_bottoms);
    if !(xb < z_feed && z_feed < xd) {
        return Err(StagesError::Dimension(format!(
            "specs must satisfy x_B < z_F < x_D, got x_B = {xb}, z_F = {z_feed}, x_D = {xd}"
        )));
    }
    let feed_point = q_line_curve_intersection(curve, z_feed, q)?;
    let (xq, yq) = feed_point;

    // --- Rectifying side: max slope of the line (x_D, x_D) → (x_e, y_e). ---
    let mut m_rect = f64::NEG_INFINITY;
    let mut pinch_rect = feed_point;
    let candidates = curve
        .x_samples()
        .iter()
        .zip(curve.y_samples())
        .map(|(&x, &y)| (x, y))
        .filter(|&(x, _)| x >= xq && x < xd)
        .chain(std::iter::once(feed_point));
    for (xe, ye) in candidates {
        if xd - xe < 1e-12 {
            continue;
        }
        let m = (xd - ye) / (xd - xe);
        if m > m_rect {
            m_rect = m;
            pinch_rect = (xe, ye);
        }
    }
    if m_rect >= 1.0 {
        return Err(StagesError::Infeasible(format!(
            "x_D = {xd} is unreachable at any reflux (equilibrium curve at or below the diagonal — azeotrope?)"
        )));
    }
    let m_rect = m_rect.max(0.0); // a negative slope means even R = 0 clears the curve
    let r_rect = m_rect / (1.0 - m_rect);

    // --- Stripping side: min slope of the line (x_B, x_B) → (x_e, y_e),
    //     converted to an equivalent R through the feed-section balances. ---
    let mut s_strip = f64::INFINITY;
    let mut pinch_strip = feed_point;
    let candidates = curve
        .x_samples()
        .iter()
        .zip(curve.y_samples())
        .map(|(&x, &y)| (x, y))
        .filter(|&(x, _)| x > xb && x <= xq)
        .chain(std::iter::once(feed_point));
    for (xe, ye) in candidates {
        if xe - xb < 1e-12 {
            continue;
        }
        let s = (ye - xb) / (xe - xb);
        if s < s_strip {
            s_strip = s;
            pinch_strip = (xe, ye);
        }
    }
    // With d = D/F from the material balances, L̄/V̄ = s at the pinch solves
    // to R (module docs). s → 1 as x_B recovery → total (R → ∞).
    let d = (z_feed - xb) / (xd - xb);
    if s_strip <= 1.0 + 1e-12 {
        return Err(StagesError::Infeasible(format!(
            "x_B = {xb} is unreachable at any boilup (stripping line pinned to the diagonal)"
        )));
    }
    let r_strip = (q + s_strip * (1.0 - q) - s_strip * d) / (d * (s_strip - 1.0));

    let (r_min, pinch) = if r_rect >= r_strip {
        (r_rect, pinch_rect)
    } else {
        (r_strip, pinch_strip)
    };
    // The pinch is "tangent" when it sits strictly away from the feed point.
    let tangent = (pinch.0 - xq).abs() > 1e-9 || (pinch.1 - yq).abs() > 1e-9;
    Ok(RminResult {
        r_min,
        pinch,
        tangent,
        feed_point,
    })
}

/// Solve the horizontal step: the liquid composition `x` on the
/// (pseudo-)equilibrium curve at vapor composition `y_target`.
///
/// With `E_MV = 1` this is the direct inverse `x*(y)`. With `E_MV < 1` the
/// pseudo-curve `y_eff(x) = y_op(x) + E · (y*(x) − y_op(x))` is inverted by
/// bisection on `[0, x_hi]` — `y_eff` is monotone increasing, so the root is
/// unique.
fn step_to_curve(
    curve: &EquilibriumCurve,
    op_line: Line,
    murphree: f64,
    y_target: f64,
    x_hi: f64,
) -> Result<f64> {
    if (murphree - 1.0).abs() < 1e-12 {
        return curve.x_of_y(y_target);
    }
    let g = |x: f64| -> Result<f64> {
        let y_op = op_line.y_at(x);
        Ok(y_op + murphree * (curve.y_of_x(x)? - y_op) - y_target)
    };
    let (mut lo, mut hi) = (0.0_f64, x_hi);
    let g_lo = g(lo)?;
    let g_hi = g(hi)?;
    if g_lo > 0.0 || g_hi < 0.0 {
        // The pseudo-curve can't reach y_target within [0, x_hi] — the
        // staircase has pinched against the operating line.
        return Err(StagesError::Convergence(format!(
            "pseudo-equilibrium step cannot reach y = {y_target} (E_MV = {murphree}) — \
             the staircase has pinched; raise R or E_MV"
        )));
    }
    // 60 bisections shrink the bracket below 1e-18 — far past f64 precision.
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if g(mid)? < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// Run the full McCabe–Thiele construction for a design specification.
///
/// # Arguments
/// * `curve` — the binary equilibrium curve at column pressure
/// * `spec` — the design spec (compositions, `q`, `R`, `E_MV`, condenser)
///
/// # Errors
/// - [`StagesError::Infeasible`] if the spec is unreachable at any reflux.
/// - [`StagesError::Convergence`] if `R ≤ R_min` (the staircase pinches
///   before reaching `x_B`).
pub fn mccabe_thiele(
    curve: &EquilibriumCurve,
    spec: McCabeThieleSpec,
) -> Result<McCabeThieleResult> {
    spec.validate()?;
    let (xd, xb, r) = (spec.x_distillate, spec.x_bottoms, spec.reflux);

    // S&H eq. 7-9: the rectifying line through (x_D, x_D) with slope R/(R+1).
    let rectifying = Line {
        slope: r / (r + 1.0),
        intercept: xd / (r + 1.0),
    };
    let intersection = operating_intersection(rectifying, spec.z_feed, spec.q)?;
    let (xi, yi) = intersection;
    if yi < xi {
        return Err(StagesError::Convergence(format!(
            "operating-line intersection ({xi:.4}, {yi:.4}) fell below the diagonal — R = {r} \
             is below R_min for this q"
        )));
    }
    // S&H eq. 7-12 in two-point form: the stripping line joins (x_B, x_B) to
    // the intersection point.
    let stripping = Line {
        slope: (yi - xb) / (xi - xb),
        intercept: xb - (yi - xb) / (xi - xb) * xb,
    };

    let rmin_result = rmin(curve, xd, xb, spec.z_feed, spec.q)?;
    if r <= rmin_result.r_min {
        return Err(StagesError::Convergence(format!(
            "R = {r} ≤ R_min = {:.6}: the staircase pinches before reaching x_B",
            rmin_result.r_min
        )));
    }

    // --- Step the staircase top-down. ---
    let mut stages = Vec::new();
    let mut staircase = vec![(xd, xd)];
    let mut x_prev = xd; // liquid composition entering the step (x_{n-1})
    let mut y_n = xd; // vapor composition leaving stage n (y_1 = x_D)
    let mut in_stripping = false;
    let mut feed_stage = 0usize;
    let mut n_stages = 0.0_f64;

    for n in 1..=MAX_STAGES {
        let op_line = if in_stripping { stripping } else { rectifying };
        let x_n = step_to_curve(curve, op_line, spec.murphree, y_n, x_prev)?;
        staircase.push((x_n, y_n));
        stages.push(StagePoint {
            index: n,
            x: x_n,
            y: y_n,
        });

        if x_n <= xb {
            // Fractional final stage: the fraction of this step actually
            // needed to reach x_B (S&H §7.2.3).
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
                "staircase pinched at x = {x_n:.6} after {n} stages (R too close to R_min = {:.6})",
                rmin_result.r_min
            )));
        }

        // Optimal-feed switch: once the stage liquid passes the operating-
        // line intersection, subsequent steps ride the stripping line. The
        // first stage stepped on the stripping line is the feed stage.
        if !in_stripping && x_n <= xi {
            in_stripping = true;
            feed_stage = n;
        }

        let y_next = if in_stripping {
            stripping.y_at(x_n)
        } else {
            rectifying.y_at(x_n)
        };
        staircase.push((x_n, y_next));
        x_prev = x_n;
        y_n = y_next;

        if n == MAX_STAGES {
            return Err(StagesError::Convergence(format!(
                "exceeded {MAX_STAGES} stages without reaching x_B = {xb} — spec pinched"
            )));
        }
    }
    if feed_stage == 0 {
        // The staircase reached x_B without leaving the rectifying section —
        // the feed enters at the bottom.
        feed_stage = stages.len();
    }

    Ok(McCabeThieleResult {
        spec,
        n_stages,
        stages,
        feed_stage,
        rmin: rmin_result,
        rectifying,
        stripping,
        intersection,
        staircase,
    })
}

/// The total-reflux construction: stepping between the equilibrium curve and
/// the y = x diagonal gives the minimum stage count `N_min` (S&H §7.2.4).
///
/// # Arguments
/// * `curve` — the binary equilibrium curve
/// * `x_distillate`, `x_bottoms` — product specs, `x_B < x_D`, both in (0, 1)
/// * `murphree` — Murphree vapor efficiency in (0, 1]
pub fn total_reflux(
    curve: &EquilibriumCurve,
    x_distillate: f64,
    x_bottoms: f64,
    murphree: f64,
) -> Result<TotalRefluxResult> {
    let (xd, xb) = (x_distillate, x_bottoms);
    if !(0.0 < xb && xb < xd && xd < 1.0) {
        return Err(StagesError::Dimension(format!(
            "need 0 < x_B < x_D < 1, got x_B = {xb}, x_D = {xd}"
        )));
    }
    if !(0.0 < murphree && murphree <= 1.0) {
        return Err(StagesError::Dimension(format!(
            "Murphree efficiency must be in (0, 1], got {murphree}"
        )));
    }
    // At total reflux both operating lines collapse onto the diagonal.
    let diagonal = Line {
        slope: 1.0,
        intercept: 0.0,
    };
    let mut stages = Vec::new();
    let mut staircase = vec![(xd, xd)];
    let mut x_prev = xd;
    let mut y_n = xd;
    let mut n_min = 0.0_f64;

    for n in 1..=MAX_STAGES {
        let x_n = step_to_curve(curve, diagonal, murphree, y_n, x_prev)?;
        staircase.push((x_n, y_n));
        stages.push(StagePoint {
            index: n,
            x: x_n,
            y: y_n,
        });
        if x_n <= xb {
            let frac = if x_prev - x_n > 0.0 {
                ((x_prev - xb) / (x_prev - x_n)).clamp(0.0, 1.0)
            } else {
                1.0
            };
            n_min = (n - 1) as f64 + frac;
            break;
        }
        if x_prev - x_n < PINCH_PROGRESS {
            return Err(StagesError::Convergence(format!(
                "total-reflux staircase pinched at x = {x_n:.6} (azeotrope between x_B and x_D?)"
            )));
        }
        staircase.push((x_n, x_n)); // vertical drop to the diagonal
        x_prev = x_n;
        y_n = x_n;
        if n == MAX_STAGES {
            return Err(StagesError::Convergence(format!(
                "exceeded {MAX_STAGES} stages at total reflux — x_B = {xb} unreachable"
            )));
        }
    }
    Ok(TotalRefluxResult {
        n_min,
        stages,
        staircase,
    })
}

/// Stage count as a function of reflux ratio — the N(R) design curve.
///
/// Returns `(r, n)` pairs; `n` is NaN where the construction failed (R at or
/// below R_min), mirroring the batch NaN-on-fail contract (PLAN §8).
///
/// # Arguments
/// * `curve` — the binary equilibrium curve
/// * `spec` — a template spec; its `reflux` field is replaced by each entry
///   of `r_values`
/// * `r_values` — the reflux ratios to sweep
pub fn n_vs_r(
    curve: &EquilibriumCurve,
    spec: McCabeThieleSpec,
    r_values: &[f64],
) -> Vec<(f64, f64)> {
    r_values
        .iter()
        .map(|&r| {
            let s = McCabeThieleSpec { reflux: r, ..spec };
            match mccabe_thiele(curve, s) {
                Ok(res) => (r, res.n_stages),
                Err(_) => (r, f64::NAN),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::column::CondenserKind;

    fn spec(reflux: f64) -> McCabeThieleSpec {
        McCabeThieleSpec {
            x_distillate: 0.95,
            x_bottoms: 0.05,
            z_feed: 0.50,
            q: 1.0,
            reflux,
            murphree: 1.0,
            condenser: CondenserKind::Total,
        }
    }

    /// Underwood's binary closed form is exact on a constant-α curve with a
    /// saturated-liquid feed (S&H eq. 7-24):
    /// R_min = 1/(α−1) · [x_D/z_F − α(1−x_D)/(1−z_F)].
    #[test]
    fn rmin_matches_underwood_closed_form() {
        let alpha = 2.5;
        let curve = EquilibriumCurve::constant_alpha(alpha, 2001).unwrap();
        let r = rmin(&curve, 0.95, 0.05, 0.50, 1.0).unwrap();
        let exact = 1.0 / (alpha - 1.0) * (0.95 / 0.50 - alpha * 0.05 / 0.50);
        assert!(
            (r.r_min - exact).abs() < 2e-3,
            "R_min = {}, Underwood exact = {exact}",
            r.r_min
        );
        // Concave curve ⇒ feed pinch, not tangent.
        assert!(!r.tangent);
        assert!((r.pinch.0 - 0.50).abs() < 1e-6);
    }

    /// A saturated-vapor feed (q = 0) moves the pinch to y = z_F: the pinch
    /// composition satisfies y* = z_F on the curve.
    #[test]
    fn rmin_saturated_vapor_feed_pinch() {
        let curve = EquilibriumCurve::constant_alpha(2.5, 2001).unwrap();
        let r = rmin(&curve, 0.95, 0.05, 0.50, 0.0).unwrap();
        assert!(
            (r.feed_point.1 - 0.50).abs() < 1e-9,
            "y_q* should equal z_F"
        );
        // q = 0 needs more reflux than q = 1 at the same specs.
        let r_liq = rmin(&curve, 0.95, 0.05, 0.50, 1.0).unwrap();
        assert!(r.r_min > r_liq.r_min);
    }

    /// Fenske's equation is exact on a constant-α curve:
    /// N_min = ln[(x_D/(1−x_D)) · ((1−x_B)/x_B)] / ln α.
    #[test]
    fn total_reflux_matches_fenske() {
        let alpha = 2.5;
        let curve = EquilibriumCurve::constant_alpha(alpha, 2001).unwrap();
        let res = total_reflux(&curve, 0.95, 0.05, 1.0).unwrap();
        let fenske = ((0.95_f64 / 0.05) * (0.95 / 0.05)).ln() / alpha.ln();
        // The whole-stage counts must agree exactly. The fractional parts
        // differ by convention: the graphical construction measures the last
        // partial step linearly in x (S&H §7.2.3), while Fenske's closed form
        // is logarithmic in the composition ratio — a ~0.1-stage difference
        // here, never a whole stage.
        assert_eq!(res.n_min.ceil(), fenske.ceil(), "whole-stage counts differ");
        assert!(
            (res.n_min - fenske).abs() < 0.15,
            "N_min = {}, Fenske = {fenske}",
            res.n_min
        );
    }

    /// The full construction at R = 1.5·R_min: finite stage count above
    /// N_min, a feed stage strictly inside the column, and a staircase that
    /// starts at (x_D, x_D) and ends at or below x_B.
    #[test]
    fn construction_at_operating_reflux() {
        let curve = EquilibriumCurve::constant_alpha(2.5, 2001).unwrap();
        let r_min = rmin(&curve, 0.95, 0.05, 0.50, 1.0).unwrap().r_min;
        let res = mccabe_thiele(&curve, spec(1.5 * r_min)).unwrap();
        let nmin = total_reflux(&curve, 0.95, 0.05, 1.0).unwrap().n_min;
        assert!(res.n_stages > nmin, "N = {} ≤ N_min = {nmin}", res.n_stages);
        assert!(
            res.n_stages < 30.0,
            "N = {} implausibly large",
            res.n_stages
        );
        assert!(res.feed_stage > 1 && res.feed_stage < res.stages.len());
        assert_eq!(res.staircase[0], (0.95, 0.95));
        assert!(res.stages.last().unwrap().x <= 0.05 + 1e-12);
        // The intersection sits on both operating lines.
        let (xi, yi) = res.intersection;
        assert!((res.rectifying.y_at(xi) - yi).abs() < 1e-12);
        assert!((res.stripping.y_at(xi) - yi).abs() < 1e-9);
    }

    /// More reflux ⇒ fewer stages, monotonically, approaching N_min.
    #[test]
    fn n_decreases_with_reflux() {
        let curve = EquilibriumCurve::constant_alpha(2.5, 2001).unwrap();
        let pts = n_vs_r(&curve, spec(0.0), &[1.3, 1.8, 2.5, 4.0, 8.0]);
        let ns: Vec<f64> = pts.iter().map(|&(_, n)| n).collect();
        assert!(ns.iter().all(|n| n.is_finite()));
        assert!(
            ns.windows(2).all(|w| w[1] < w[0]),
            "N(R) not decreasing: {ns:?}"
        );
        let nmin = total_reflux(&curve, 0.95, 0.05, 1.0).unwrap().n_min;
        assert!(*ns.last().unwrap() > nmin);
    }

    /// R below R_min fails with a Convergence error, and n_vs_r maps that to
    /// NaN.
    #[test]
    fn below_rmin_is_caught() {
        let curve = EquilibriumCurve::constant_alpha(2.5, 2001).unwrap();
        let r_min = rmin(&curve, 0.95, 0.05, 0.50, 1.0).unwrap().r_min;
        assert!(matches!(
            mccabe_thiele(&curve, spec(0.9 * r_min)),
            Err(StagesError::Convergence(_))
        ));
        let pts = n_vs_r(&curve, spec(0.0), &[0.9 * r_min]);
        assert!(pts[0].1.is_nan());
    }

    /// Murphree efficiency below 1 increases the stage count.
    #[test]
    fn murphree_increases_stage_count() {
        let curve = EquilibriumCurve::constant_alpha(2.5, 2001).unwrap();
        let ideal = mccabe_thiele(&curve, spec(2.0)).unwrap();
        let real = mccabe_thiele(
            &curve,
            McCabeThieleSpec {
                murphree: 0.7,
                ..spec(2.0)
            },
        )
        .unwrap();
        assert!(
            real.n_stages > ideal.n_stages,
            "E_MV = 0.7 gave {} stages vs {} ideal",
            real.n_stages,
            ideal.n_stages
        );
    }

    /// A tangent pinch is detected on a curve with an inflection: a distorted
    /// curve that dips toward the diagonal above the feed point forces the
    /// controlling pinch off the feed.
    #[test]
    fn tangent_pinch_detected() {
        // Blend a high-α curve at low x with a near-diagonal stretch at high
        // x — qualitatively an ethanol–water shape (without the azeotrope).
        let n = 2001;
        let mut x = Vec::with_capacity(n);
        let mut y = Vec::with_capacity(n);
        for i in 0..n {
            let xi = i as f64 / (n - 1) as f64;
            // α decays from 6 at x=0 to 1.05 at x=1 — the curve flattens
            // toward the diagonal at high x, creating a tangent pinch.
            let alpha = 1.05 + 4.95 * (1.0 - xi).powi(2);
            y.push(alpha * xi / (1.0 + (alpha - 1.0) * xi));
            x.push(xi);
        }
        let curve = EquilibriumCurve::from_points(x, y, Vec::new(), None).unwrap();
        let r = rmin(&curve, 0.90, 0.05, 0.30, 1.0).unwrap();
        assert!(
            r.tangent,
            "expected a tangent pinch, got feed pinch at {:?}",
            r.pinch
        );
        assert!(
            r.pinch.0 > 0.30,
            "tangent pinch should sit above the feed composition"
        );
        // The construction still works above that R_min.
        let res = mccabe_thiele(
            &curve,
            McCabeThieleSpec {
                x_distillate: 0.90,
                x_bottoms: 0.05,
                z_feed: 0.30,
                q: 1.0,
                reflux: 1.3 * r.r_min,
                murphree: 1.0,
                condenser: CondenserKind::Total,
            },
        )
        .unwrap();
        assert!(res.n_stages.is_finite());
    }
}
