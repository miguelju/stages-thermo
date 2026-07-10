//! The binary vapor–liquid equilibrium curve `y*(x)` — the foundation every
//! graphical binary method stands on.
//!
//! ## Textbook anchoring
//!
//! For a binary mixture at fixed column pressure `P`, each liquid composition
//! `x` (light component) has a bubble point: the temperature `T(x)` at which
//! the first vapor bubble forms, and the composition `y*(x)` of that
//! incipient vapor. The locus `(x, y*(x))` is the **equilibrium curve** of
//! the McCabe–Thiele diagram (Seader, Henley & Roper, *Separation Process
//! Principles*, §7.1; the classic x–y diagram of McCabe & Thiele, *Ind. Eng.
//! Chem.* **1925**, 17, 605).
//!
//! Where the textbook assumes constant relative volatility
//! `y = αx / (1 + (α − 1)x)` (S&H eq. 7-13), this module computes the real
//! curve by sweeping bubble-point calculations through the
//! [`ThermoSystem`](crate::thermo::ThermoSystem) adapter — vle-thermo's
//! K-values replace the constant-α idealization. The constant-α form is kept
//! as [`EquilibriumCurve::constant_alpha`], both as a teaching device and as
//! the analytic oracle for the stepping kernels (Fenske's equation is exact
//! on it — see the tests in [`super::mccabe_thiele`]).
//!
//! ## Units
//!
//! Pressure **kPa** (absolute), temperature **K**, compositions are light-
//! component mole fractions (dimensionless).

use crate::thermo::{Phase, ThermoSystem};
use crate::types::{Result, StagesError};

/// A sampled binary equilibrium curve `y*(x)`, with the bubble-point
/// temperature profile `T(x)` when it came from a thermodynamic model.
///
/// The samples are stored on a strictly increasing `x` grid including both
/// endpoints; queries interpolate linearly between samples, so the default
/// grid (101 points) reproduces the underlying model to well below the line
/// width of any plotted diagram.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "python", pyo3::pyclass)]
pub struct EquilibriumCurve {
    /// Column pressure in **kPa** (absolute); `None` for synthetic curves
    /// (constant-α or raw points without a pressure).
    pressure: Option<f64>,
    /// Light-component liquid mole fractions — strictly increasing, spanning
    /// [0, 1].
    x: Vec<f64>,
    /// Equilibrium vapor mole fractions `y*(x_i)` — nondecreasing.
    y: Vec<f64>,
    /// Bubble temperatures `T(x_i)` in **K**; empty for synthetic curves.
    t: Vec<f64>,
}

impl EquilibriumCurve {
    /// Generate the curve from a thermodynamic system by sweeping
    /// bubble-temperature calculations over an even `x` grid.
    ///
    /// # Arguments
    /// * `system` — a **binary** [`ThermoSystem`], light component first
    /// * `pressure` — column pressure in **kPa** (absolute)
    /// * `n_points` — number of samples (≥ 5; 101 is a good default)
    ///
    /// # Errors
    /// - [`StagesError::Dimension`] if the system is not binary or
    ///   `n_points < 5`.
    /// - [`StagesError::Thermo`] if a bubble-point calculation fails, or if
    ///   the first component turns out to be the *heavier* one (`y* < x` at
    ///   the dilute end) — reorder the components.
    pub fn from_thermo(system: &ThermoSystem, pressure: f64, n_points: usize) -> Result<Self> {
        if system.n_components() != 2 {
            return Err(StagesError::Dimension(format!(
                "equilibrium curve needs a binary system, got {} components",
                system.n_components()
            )));
        }
        if n_points < 5 {
            return Err(StagesError::Dimension(format!(
                "n_points must be at least 5, got {n_points}"
            )));
        }

        let mut x = Vec::with_capacity(n_points);
        let mut y = Vec::with_capacity(n_points);
        let mut t = Vec::with_capacity(n_points);
        for i in 0..n_points {
            // Even grid over [0, 1]; integer arithmetic first so the
            // endpoints are exactly 0.0 and 1.0 (no floating-point drift).
            let xi = i as f64 / (n_points - 1) as f64;
            let bp = system.bubble_temperature(pressure, &[xi, 1.0 - xi])?;
            x.push(xi);
            // Pin the thermodynamic identities y*(0) = 0 and y*(1) = 1
            // exactly — pure components leave no room for a composition
            // difference, and downstream stepping relies on the endpoints.
            let yi = if i == 0 {
                0.0
            } else if i == n_points - 1 {
                1.0
            } else {
                bp.y[0]
            };
            y.push(yi);
            t.push(bp.value);
        }

        // The light component must be listed first: at dilute x the vapor
        // must be enriched (y* > x). Checking the first interior point keeps
        // azeotropic systems valid (they cross y = x at high x, not low x).
        if y[1] <= x[1] {
            return Err(StagesError::Thermo(format!(
                "component order looks inverted: y*({:.3}) = {:.4} ≤ x — list the more volatile component first",
                x[1], y[1]
            )));
        }

        Ok(Self {
            pressure: Some(pressure),
            x,
            y,
            t,
        })
    }

    /// The constant-relative-volatility idealization
    /// `y = αx / (1 + (α − 1)x)` (S&H eq. 7-13), sampled on an even grid.
    ///
    /// Used for teaching and as the analytic test oracle: Fenske's and
    /// Underwood's closed forms are exact on this curve.
    ///
    /// # Arguments
    /// * `alpha` — relative volatility α > 0 (α > 1 for a meaningful light
    ///   component)
    /// * `n_points` — number of samples (≥ 5)
    pub fn constant_alpha(alpha: f64, n_points: usize) -> Result<Self> {
        // `alpha > 0.0` written positively would let NaN slip through — the
        // explicit NaN check keeps the guard clippy-clean AND NaN-safe.
        if alpha.is_nan() || alpha <= 0.0 {
            return Err(StagesError::Dimension(format!(
                "alpha must be positive, got {alpha}"
            )));
        }
        if n_points < 5 {
            return Err(StagesError::Dimension(format!(
                "n_points must be at least 5, got {n_points}"
            )));
        }
        let mut x = Vec::with_capacity(n_points);
        let mut y = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let xi = i as f64 / (n_points - 1) as f64;
            x.push(xi);
            y.push(alpha * xi / (1.0 + (alpha - 1.0) * xi));
        }
        Ok(Self {
            pressure: None,
            x,
            y,
            t: Vec::new(),
        })
    }

    /// Build a curve from raw `(x, y)` samples — e.g. digitized literature
    /// data — with an optional temperature profile and pressure.
    ///
    /// # Arguments
    /// * `x` — strictly increasing light-component liquid fractions spanning
    ///   `[0, 1]` (first sample 0, last sample 1)
    /// * `y` — equilibrium vapor fractions, same length, nondecreasing
    /// * `t` — bubble temperatures in **K** (empty slice if unavailable)
    /// * `pressure` — the pressure in **kPa** the data was measured at, if
    ///   known
    pub fn from_points(
        x: Vec<f64>,
        y: Vec<f64>,
        t: Vec<f64>,
        pressure: Option<f64>,
    ) -> Result<Self> {
        if x.len() != y.len() || x.len() < 5 {
            return Err(StagesError::Dimension(format!(
                "x and y must have equal length ≥ 5, got {} and {}",
                x.len(),
                y.len()
            )));
        }
        if !t.is_empty() && t.len() != x.len() {
            return Err(StagesError::Dimension(format!(
                "t must be empty or match x's length {}, got {}",
                x.len(),
                t.len()
            )));
        }
        if x[0] != 0.0 || *x.last().unwrap() != 1.0 {
            return Err(StagesError::Dimension(
                "x must span [0, 1] with x[0] = 0 and x[last] = 1".into(),
            ));
        }
        if !x.windows(2).all(|w| w[1] > w[0]) {
            return Err(StagesError::Dimension(
                "x must be strictly increasing".into(),
            ));
        }
        if !y.windows(2).all(|w| w[1] >= w[0]) {
            return Err(StagesError::Dimension(
                "y must be nondecreasing (binary VLE curves are monotone in x)".into(),
            ));
        }
        if y.iter().any(|&v| !(0.0..=1.0).contains(&v)) {
            return Err(StagesError::Dimension("y values must lie in [0, 1]".into()));
        }
        Ok(Self { pressure, x, y, t })
    }

    /// Column pressure in **kPa** (absolute), if the curve has one.
    pub fn pressure(&self) -> Option<f64> {
        self.pressure
    }

    /// The sampled liquid compositions (strictly increasing, spans [0, 1]).
    pub fn x_samples(&self) -> &[f64] {
        &self.x
    }

    /// The sampled equilibrium vapor compositions.
    pub fn y_samples(&self) -> &[f64] {
        &self.y
    }

    /// The sampled bubble temperatures in **K** (empty for synthetic curves).
    pub fn t_samples(&self) -> &[f64] {
        &self.t
    }

    /// Interpolated equilibrium vapor fraction `y*(x)`.
    ///
    /// # Arguments
    /// * `x` — liquid light-component mole fraction in [0, 1]
    pub fn y_of_x(&self, x: f64) -> Result<f64> {
        interp(&self.x, &self.y, x)
    }

    /// Interpolated inverse `x*(y)` — the liquid composition in equilibrium
    /// with vapor fraction `y`. Well-defined because binary `y*(x)` is
    /// monotone (enforced at construction).
    ///
    /// # Arguments
    /// * `y` — vapor light-component mole fraction in [0, 1]
    pub fn x_of_y(&self, y: f64) -> Result<f64> {
        interp(&self.y, &self.x, y)
    }

    /// Interpolated bubble temperature `T(x)` in **K**.
    ///
    /// # Errors
    /// [`StagesError::Dimension`] if the curve has no temperature data
    /// (synthetic curves).
    pub fn temperature_of_x(&self, x: f64) -> Result<f64> {
        if self.t.is_empty() {
            return Err(StagesError::Dimension(
                "this curve has no temperature data (synthetic curve)".into(),
            ));
        }
        interp(&self.x, &self.t, x)
    }

    /// Point relative volatility `α(x) = [y/(1−y)] / [x/(1−x)]` from the
    /// interpolated curve, for `x` strictly inside (0, 1).
    pub fn relative_volatility(&self, x: f64) -> Result<f64> {
        if !(0.0 < x && x < 1.0) {
            return Err(StagesError::Dimension(format!(
                "relative volatility needs x strictly inside (0, 1), got {x}"
            )));
        }
        let y = self.y_of_x(x)?;
        Ok((y / (1.0 - y)) / (x / (1.0 - x)))
    }
}

/// A binary equilibrium curve **augmented with saturated-phase molar
/// enthalpies** — the enthalpy–composition (H–x–y) diagram that
/// Ponchon–Savarit (M2) constructs on.
///
/// It wraps an [`EquilibriumCurve`] and adds, on the same `x` grid:
///
/// - `h_liq[i]` — saturated-**liquid** molar enthalpy at `(x_i, T(x_i))`;
/// - `h_vap[i]` — saturated-**vapor** molar enthalpy of the equilibrium vapor
///   `y*(x_i)` at the *same* bubble temperature `T(x_i)`.
///
/// The **tie line** of the H–x–y diagram connects the liquid point
/// `(x_i, h_liq[i])` to its equilibrium vapor point `(y*(x_i), h_vap[i])` —
/// both at `T(x_i)`, because a liquid and its incipient vapor share a
/// temperature at the bubble point (Seader, Henley & Roper, Ch. 7).
///
/// Enthalpy lives in a **separate** type (rather than as extra fields on
/// `EquilibriumCurve`) so synthetic constant-α curves — which carry no
/// enthalpy model — aren't forced to hold it.
///
/// ## Units
/// Enthalpies in **kJ/kmol** relative to the [`ThermoSystem`] datum; pressure
/// **kPa** (absolute), temperature **K**, compositions dimensionless.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "python", pyo3::pyclass)]
pub struct EnthalpyCurve {
    curve: EquilibriumCurve,
    /// Saturated-liquid molar enthalpy at each `x_i`, **kJ/kmol**.
    h_liq: Vec<f64>,
    /// Saturated-vapor molar enthalpy at each `y*(x_i)`, **kJ/kmol**.
    h_vap: Vec<f64>,
}

impl EnthalpyCurve {
    /// **Route (a)** — the H–x–y chart *computed* from a thermodynamic model:
    /// sweep the bubble curve and evaluate saturated-liquid and
    /// saturated-vapor molar enthalpies at each grid point.
    ///
    /// # Arguments
    /// * `system` — a **binary** [`ThermoSystem`], light component first
    /// * `pressure` — column pressure in **kPa** (absolute)
    /// * `n_points` — number of samples (≥ 5; 101 is a good default)
    ///
    /// # Errors
    /// [`StagesError::Dimension`] / [`StagesError::Thermo`] as
    /// [`EquilibriumCurve::from_thermo`], plus any enthalpy-evaluation failure.
    pub fn from_thermo(system: &ThermoSystem, pressure: f64, n_points: usize) -> Result<Self> {
        let curve = EquilibriumCurve::from_thermo(system, pressure, n_points)?;
        let n = curve.x_samples().len();
        let mut h_liq = Vec::with_capacity(n);
        let mut h_vap = Vec::with_capacity(n);
        for i in 0..n {
            // Bubble temperature at this liquid composition; liquid `x_i` and
            // its incipient vapor `y*(x_i)` are both evaluated there.
            let (xi, yi, ti) = (
                curve.x_samples()[i],
                curve.y_samples()[i],
                curve.t_samples()[i],
            );
            h_liq.push(system.phase_enthalpy(ti, pressure, &[xi, 1.0 - xi], Phase::Liquid)?);
            h_vap.push(system.phase_enthalpy(ti, pressure, &[yi, 1.0 - yi], Phase::Vapor)?);
        }
        Ok(Self {
            curve,
            h_liq,
            h_vap,
        })
    }

    /// **Route (b)** — literature H–x–y data fed directly: digitized enthalpy +
    /// VLE points (e.g. Ibrahim–Klein / Tillner-Roth for ammonia–water). Runs
    /// the pure graphical construction on reference data, with no
    /// thermodynamic model — the "textbook, faithfully" path (M2 lesson).
    ///
    /// # Arguments
    /// * `x`, `y`, `t` — as [`EquilibriumCurve::from_points`] (strictly
    ///   increasing `x` spanning [0, 1]; nondecreasing `y`; `t` in **K** or
    ///   empty)
    /// * `h_liq` — saturated-liquid molar enthalpy at each `x_i`, **kJ/kmol**
    /// * `h_vap` — saturated-vapor molar enthalpy at each `y_i`, **kJ/kmol**
    /// * `pressure` — the pressure in **kPa** the data was measured at, if known
    pub fn from_points(
        x: Vec<f64>,
        y: Vec<f64>,
        t: Vec<f64>,
        h_liq: Vec<f64>,
        h_vap: Vec<f64>,
        pressure: Option<f64>,
    ) -> Result<Self> {
        let n = x.len();
        if h_liq.len() != n || h_vap.len() != n {
            return Err(StagesError::Dimension(format!(
                "h_liq and h_vap must match x's length {n}, got {} and {}",
                h_liq.len(),
                h_vap.len()
            )));
        }
        // Reuse EquilibriumCurve's shape validation (monotonicity, span, etc.).
        let curve = EquilibriumCurve::from_points(x, y, t, pressure)?;
        Ok(Self {
            curve,
            h_liq,
            h_vap,
        })
    }

    /// The underlying equilibrium curve (`x`, `y*`, `T`, pressure).
    pub fn equilibrium(&self) -> &EquilibriumCurve {
        &self.curve
    }

    /// Saturated-liquid molar-enthalpy samples, **kJ/kmol**.
    pub fn h_liq_samples(&self) -> &[f64] {
        &self.h_liq
    }

    /// Saturated-vapor molar-enthalpy samples, **kJ/kmol**.
    pub fn h_vap_samples(&self) -> &[f64] {
        &self.h_vap
    }

    /// Interpolated saturated-liquid molar enthalpy `h_L(x)` in **kJ/kmol**.
    ///
    /// # Arguments
    /// * `x` — liquid light-component mole fraction in [0, 1]
    pub fn h_liquid_of_x(&self, x: f64) -> Result<f64> {
        interp(self.curve.x_samples(), &self.h_liq, x)
    }

    /// Interpolated saturated-vapor molar enthalpy `H_V(y)` in **kJ/kmol**.
    ///
    /// Parameterized by vapor composition: the vapor grid is `y*(x_i)`, which
    /// is monotone (binary `y*(x)` is), so `interp` over `(y*, h_vap)` is
    /// well-defined.
    ///
    /// # Arguments
    /// * `y` — vapor light-component mole fraction in [0, 1]
    pub fn h_vapor_of_y(&self, y: f64) -> Result<f64> {
        interp(self.curve.y_samples(), &self.h_vap, y)
    }
}

/// Piecewise-linear interpolation of `ys` over the nondecreasing grid `xs`.
///
/// `xs` must be nondecreasing (the curve constructors enforce it); ties are
/// handled by taking the first bracketing segment with nonzero width.
fn interp(xs: &[f64], ys: &[f64], x: f64) -> Result<f64> {
    let (lo, hi) = (xs[0], *xs.last().unwrap());
    if x < lo - 1e-12 || x > hi + 1e-12 {
        return Err(StagesError::Dimension(format!(
            "query {x} outside the sampled range [{lo}, {hi}]"
        )));
    }
    let x = x.clamp(lo, hi);
    // `partition_point` is a binary search: it returns the index of the first
    // element for which the predicate is false, i.e. the first grid point
    // ≥ x. O(log n) per query.
    let idx = xs.partition_point(|&v| v < x);
    if idx == 0 {
        return Ok(ys[0]);
    }
    let (x0, x1) = (xs[idx - 1], xs[idx.min(xs.len() - 1)]);
    let (y0, y1) = (ys[idx - 1], ys[idx.min(ys.len() - 1)]);
    if x1 == x0 {
        return Ok(y0);
    }
    Ok(y0 + (y1 - y0) * (x - x0) / (x1 - x0))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The constant-α curve reproduces S&H eq. 7-13 exactly at the samples
    /// and interpolates tightly between them.
    #[test]
    fn constant_alpha_matches_closed_form() {
        let curve = EquilibriumCurve::constant_alpha(2.5, 201).unwrap();
        for &x in &[0.1, 0.25, 0.5, 0.75, 0.9] {
            let exact = 2.5 * x / (1.0 + 1.5 * x);
            let y = curve.y_of_x(x).unwrap();
            assert!((y - exact).abs() < 5e-5, "y({x}) = {y}, exact {exact}");
        }
    }

    /// x_of_y inverts y_of_x within interpolation error.
    #[test]
    fn inverse_roundtrip() {
        let curve = EquilibriumCurve::constant_alpha(2.5, 201).unwrap();
        for &x in &[0.05, 0.3, 0.6, 0.95] {
            let y = curve.y_of_x(x).unwrap();
            let x_back = curve.x_of_y(y).unwrap();
            assert!((x_back - x).abs() < 1e-9, "roundtrip {x} → {y} → {x_back}");
        }
    }

    /// Point relative volatility recovers α on a constant-α curve.
    #[test]
    fn relative_volatility_recovers_alpha() {
        let curve = EquilibriumCurve::constant_alpha(2.5, 401).unwrap();
        for &x in &[0.2, 0.5, 0.8] {
            let a = curve.relative_volatility(x).unwrap();
            assert!((a - 2.5).abs() < 5e-3, "α({x}) = {a}");
        }
    }

    /// from_points validates shape errors.
    #[test]
    fn from_points_validation() {
        // Too short.
        assert!(
            EquilibriumCurve::from_points(vec![0.0, 1.0], vec![0.0, 1.0], vec![], None).is_err()
        );
        // Doesn't span [0, 1].
        assert!(
            EquilibriumCurve::from_points(
                vec![0.1, 0.3, 0.5, 0.7, 1.0],
                vec![0.2, 0.5, 0.7, 0.85, 1.0],
                vec![],
                None
            )
            .is_err()
        );
        // Decreasing y.
        assert!(
            EquilibriumCurve::from_points(
                vec![0.0, 0.25, 0.5, 0.75, 1.0],
                vec![0.0, 0.6, 0.5, 0.9, 1.0],
                vec![],
                None
            )
            .is_err()
        );
    }

    /// The computed enthalpy curve orders the phases physically: saturated
    /// vapor lies above saturated liquid at every interior composition
    /// (positive latent heat), and interpolation round-trips a grid sample.
    #[test]
    fn enthalpy_curve_vapor_above_liquid() {
        use crate::thermo::ThermoSystem;
        let sys = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
        let ec = EnthalpyCurve::from_thermo(&sys, 101.325, 51).unwrap();
        let (hl, hv) = (ec.h_liq_samples(), ec.h_vap_samples());
        for i in 1..hl.len() - 1 {
            assert!(
                hv[i] > hl[i],
                "H_vap[{i}] = {} ≤ H_liq[{i}] = {} kJ/kmol",
                hv[i],
                hl[i]
            );
        }
        let x = ec.equilibrium().x_samples()[10];
        assert!((ec.h_liquid_of_x(x).unwrap() - hl[10]).abs() < 1e-6);
    }

    /// `EnthalpyCurve::from_points` rejects mismatched enthalpy-array lengths
    /// and accepts consistent literature data (the route-(b) path).
    #[test]
    fn enthalpy_from_points_validation() {
        let x = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let y = vec![0.0, 0.4, 0.6, 0.8, 1.0];
        // Wrong h_liq length.
        assert!(
            EnthalpyCurve::from_points(
                x.clone(),
                y.clone(),
                vec![],
                vec![0.0, 1.0, 2.0],
                vec![10.0, 11.0, 12.0, 13.0, 14.0],
                None
            )
            .is_err()
        );
        // Consistent data.
        assert!(
            EnthalpyCurve::from_points(
                x,
                y,
                vec![],
                vec![0.0, 1.0, 2.0, 3.0, 4.0],
                vec![10.0, 11.0, 12.0, 13.0, 14.0],
                None
            )
            .is_ok()
        );
    }
}
