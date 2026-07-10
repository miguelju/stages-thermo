//! Python bindings for the stages-thermo engine.
//!
//! Builds a Python extension module `stages._engine` when this crate is
//! compiled with the `python` feature. End-users import the higher-level
//! `stages` Python wrapper, which re-exports the pieces they need.
//!
//! ## M1 surface
//!
//! - [`ThermoSystem`] — the thermo adapter (φ-φ Peng–Robinson and γ-φ van
//!   Laar constructors).
//! - [`EquilibriumCurve`] — `from_thermo` / `constant_alpha` / `from_points`
//!   constructors plus the interpolation queries.
//! - `mccabe_thiele`, `rmin`, `total_reflux`, `n_vs_r` — the McCabe–Thiele
//!   construction, returning the rich result pyclasses.
//! - [`BinaryColumn`] / [`Feed`] / [`CondenserKind`] — the binary column
//!   model with its material balances.
//!
//! Units cross this boundary in the canonical engine set — K, kPa (absolute),
//! mole fractions — as plain floats. Unit-string ergonomics (pint) live in
//! the Python wrapper, not here.
//!
//! Per `CLAUDE.md`'s "PyO3 Bindings Rule", every milestone that adds public
//! Rust functionality also exposes it here in the same commit series, with at
//! least one round-trip test in `python/tests/`.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use crate::binary::equilibrium::{EnthalpyCurve, EquilibriumCurve};
use crate::binary::mccabe_thiele::{
    self as mt, Line, McCabeThieleResult, McCabeThieleSpec, RminResult, StagePoint,
    TotalRefluxResult,
};
use crate::binary::ponchon_savarit::{self as ps, PonchonSavaritResult, PonchonSavaritSpec};
use crate::column::{BinaryColumn, CondenserKind, Feed};
use crate::thermo::{Phase, ThermoSystem};
use crate::types::StagesError;

/// Map engine errors onto idiomatic Python exceptions: bad inputs and
/// impossible specs raise `ValueError`; thermodynamic and convergence
/// failures raise `RuntimeError`.
fn to_py_err(e: StagesError) -> PyErr {
    match e {
        StagesError::Dimension(_) | StagesError::Infeasible(_) => {
            PyValueError::new_err(e.to_string())
        }
        StagesError::Thermo(_) | StagesError::Convergence(_) => {
            PyRuntimeError::new_err(e.to_string())
        }
    }
}

/// Parse the Python-facing condenser string ("total" | "partial").
fn parse_condenser(kind: &str) -> PyResult<CondenserKind> {
    match kind.to_ascii_lowercase().as_str() {
        "total" => Ok(CondenserKind::Total),
        "partial" => Ok(CondenserKind::Partial),
        other => Err(PyValueError::new_err(format!(
            "condenser must be 'total' or 'partial', got {other:?}"
        ))),
    }
}

/// Parse the Python-facing phase string ("liquid" | "vapor").
fn parse_phase(phase: &str) -> PyResult<Phase> {
    match phase.to_ascii_lowercase().as_str() {
        "liquid" | "l" => Ok(Phase::Liquid),
        "vapor" | "vapour" | "v" => Ok(Phase::Vapor),
        other => Err(PyValueError::new_err(format!(
            "phase must be 'liquid' or 'vapor', got {other:?}"
        ))),
    }
}

#[pymethods]
impl ThermoSystem {
    /// φ-φ system: Peng–Robinson both phases, classical mixing, kij = 0.
    /// Components come from vle-thermo's built-in database by name.
    #[staticmethod]
    #[pyo3(name = "peng_robinson")]
    fn py_peng_robinson(names: Vec<String>) -> PyResult<Self> {
        let refs: Vec<&str> = names.iter().map(String::as_str).collect();
        Self::peng_robinson(&refs).map_err(to_py_err)
    }

    /// γ-φ system: van Laar liquid + ideal-gas vapor, for a binary pair.
    /// `a12`/`a21` are the dimensionless van Laar parameters in the given
    /// component order.
    #[staticmethod]
    #[pyo3(name = "van_laar")]
    fn py_van_laar(names: [String; 2], a12: f64, a21: f64) -> PyResult<Self> {
        Self::van_laar(&[names[0].as_str(), names[1].as_str()], a12, a21).map_err(to_py_err)
    }

    /// γ-φ system: NRTL liquid + ideal-gas vapor, for a binary pair.
    /// `a12`/`a21` are the NRTL interaction energies gᵢⱼ − gⱼⱼ in **kJ/kmol**
    /// (vle-thermo forms τᵢⱼ = aᵢⱼ/(R·T) internally); `alpha` is the shared
    /// non-randomness parameter (dimensionless, typically 0.2–0.47).
    #[staticmethod]
    #[pyo3(name = "nrtl", signature = (names, a12, a21, alpha = 0.3))]
    fn py_nrtl(names: [String; 2], a12: f64, a21: f64, alpha: f64) -> PyResult<Self> {
        Self::nrtl(&[names[0].as_str(), names[1].as_str()], a12, a21, alpha).map_err(to_py_err)
    }

    /// Component names, in composition-index order.
    #[getter]
    fn components(&self) -> Vec<String> {
        self.component_names()
    }

    /// Enthalpy-datum temperature [K] (set once for the whole column).
    #[getter(t_ref)]
    fn py_t_ref(&self) -> f64 {
        self.t_ref()
    }

    /// Enthalpy-datum pressure [kPa].
    #[getter(p_ref)]
    fn py_p_ref(&self) -> f64 {
        self.p_ref()
    }

    /// Molar enthalpy [kJ/kmol] of composition `comp` in `phase`
    /// ("liquid" | "vapor") at temperature `t` [K] and pressure `p` [kPa].
    #[pyo3(name = "phase_enthalpy")]
    fn py_phase_enthalpy(&self, t: f64, p: f64, comp: Vec<f64>, phase: &str) -> PyResult<f64> {
        let ph = parse_phase(phase)?;
        self.phase_enthalpy(t, p, &comp, ph).map_err(to_py_err)
    }

    /// Bubble temperature [K] at pressure `p` [kPa] for liquid composition
    /// `x`. Returns `(T, y, k)`.
    #[pyo3(name = "bubble_temperature")]
    fn py_bubble_temperature(&self, p: f64, x: Vec<f64>) -> PyResult<(f64, Vec<f64>, Vec<f64>)> {
        let bp = self.bubble_temperature(p, &x).map_err(to_py_err)?;
        Ok((bp.value, bp.y, bp.k))
    }

    /// Bubble pressure [kPa] at temperature `t` [K] for liquid composition
    /// `x`. Returns `(P, y, k)`.
    #[pyo3(name = "bubble_pressure")]
    fn py_bubble_pressure(&self, t: f64, x: Vec<f64>) -> PyResult<(f64, Vec<f64>, Vec<f64>)> {
        let bp = self.bubble_pressure(t, &x).map_err(to_py_err)?;
        Ok((bp.value, bp.y, bp.k))
    }

    fn __repr__(&self) -> String {
        format!("ThermoSystem({})", self.component_names().join(", "))
    }
}

#[pymethods]
impl EquilibriumCurve {
    /// Sweep bubble-temperature calculations over an even x grid.
    /// `pressure` in kPa (absolute); light component first in `system`.
    #[staticmethod]
    #[pyo3(name = "from_thermo", signature = (system, pressure, n_points = 101))]
    fn py_from_thermo(system: &ThermoSystem, pressure: f64, n_points: usize) -> PyResult<Self> {
        Self::from_thermo(system, pressure, n_points).map_err(to_py_err)
    }

    /// The constant-relative-volatility idealization y = αx/(1 + (α−1)x).
    #[staticmethod]
    #[pyo3(name = "constant_alpha", signature = (alpha, n_points = 101))]
    fn py_constant_alpha(alpha: f64, n_points: usize) -> PyResult<Self> {
        Self::constant_alpha(alpha, n_points).map_err(to_py_err)
    }

    /// Build from raw (x, y[, t]) samples, e.g. digitized literature data.
    #[staticmethod]
    #[pyo3(name = "from_points", signature = (x, y, t = None, pressure = None))]
    fn py_from_points(
        x: Vec<f64>,
        y: Vec<f64>,
        t: Option<Vec<f64>>,
        pressure: Option<f64>,
    ) -> PyResult<Self> {
        Self::from_points(x, y, t.unwrap_or_default(), pressure).map_err(to_py_err)
    }

    /// Sampled liquid compositions (strictly increasing, spans [0, 1]).
    #[getter]
    fn x(&self) -> Vec<f64> {
        self.x_samples().to_vec()
    }

    /// Sampled equilibrium vapor compositions y*(x).
    #[getter]
    fn y(&self) -> Vec<f64> {
        self.y_samples().to_vec()
    }

    /// Sampled bubble temperatures [K] (empty for synthetic curves).
    #[getter]
    fn t(&self) -> Vec<f64> {
        self.t_samples().to_vec()
    }

    /// Column pressure [kPa], if the curve has one.
    #[getter]
    #[pyo3(name = "pressure")]
    fn py_pressure(&self) -> Option<f64> {
        self.pressure()
    }

    /// Interpolated equilibrium vapor fraction y*(x).
    #[pyo3(name = "y_of_x")]
    fn py_y_of_x(&self, x: f64) -> PyResult<f64> {
        self.y_of_x(x).map_err(to_py_err)
    }

    /// Interpolated inverse x*(y).
    #[pyo3(name = "x_of_y")]
    fn py_x_of_y(&self, y: f64) -> PyResult<f64> {
        self.x_of_y(y).map_err(to_py_err)
    }

    /// Interpolated bubble temperature T(x) [K].
    #[pyo3(name = "temperature_of_x")]
    fn py_temperature_of_x(&self, x: f64) -> PyResult<f64> {
        self.temperature_of_x(x).map_err(to_py_err)
    }

    /// Point relative volatility α(x) from the curve.
    #[pyo3(name = "relative_volatility")]
    fn py_relative_volatility(&self, x: f64) -> PyResult<f64> {
        self.relative_volatility(x).map_err(to_py_err)
    }

    fn __repr__(&self) -> String {
        match self.pressure() {
            Some(p) => format!(
                "EquilibriumCurve({} points, P = {p} kPa)",
                self.x_samples().len()
            ),
            None => format!(
                "EquilibriumCurve({} points, synthetic)",
                self.x_samples().len()
            ),
        }
    }
}

#[pymethods]
impl EnthalpyCurve {
    /// Route (a): sweep the bubble curve and evaluate saturated-liquid and
    /// saturated-vapor molar enthalpies at each grid point. `pressure` in kPa
    /// (absolute); light component first in `system`.
    #[staticmethod]
    #[pyo3(name = "from_thermo", signature = (system, pressure, n_points = 101))]
    fn py_from_thermo(system: &ThermoSystem, pressure: f64, n_points: usize) -> PyResult<Self> {
        Self::from_thermo(system, pressure, n_points).map_err(to_py_err)
    }

    /// Route (b): literature H–x–y data fed directly. `x`, `y` (and optional
    /// `t` [K]) as `EquilibriumCurve.from_points`; `h_liq`/`h_vap` are the
    /// saturated-liquid/vapor molar enthalpies [kJ/kmol] at each sample.
    #[staticmethod]
    #[pyo3(name = "from_points", signature = (x, y, h_liq, h_vap, t = None, pressure = None))]
    fn py_from_points(
        x: Vec<f64>,
        y: Vec<f64>,
        h_liq: Vec<f64>,
        h_vap: Vec<f64>,
        t: Option<Vec<f64>>,
        pressure: Option<f64>,
    ) -> PyResult<Self> {
        Self::from_points(x, y, t.unwrap_or_default(), h_liq, h_vap, pressure).map_err(to_py_err)
    }

    /// Sampled liquid compositions.
    #[getter]
    fn x(&self) -> Vec<f64> {
        self.equilibrium().x_samples().to_vec()
    }

    /// Sampled equilibrium vapor compositions y*(x).
    #[getter]
    fn y(&self) -> Vec<f64> {
        self.equilibrium().y_samples().to_vec()
    }

    /// Sampled bubble temperatures [K] (empty for reference curves with no T).
    #[getter]
    fn t(&self) -> Vec<f64> {
        self.equilibrium().t_samples().to_vec()
    }

    /// Sampled saturated-liquid molar enthalpies [kJ/kmol].
    #[getter]
    fn h_liq(&self) -> Vec<f64> {
        self.h_liq_samples().to_vec()
    }

    /// Sampled saturated-vapor molar enthalpies [kJ/kmol].
    #[getter]
    fn h_vap(&self) -> Vec<f64> {
        self.h_vap_samples().to_vec()
    }

    /// Interpolated saturated-liquid molar enthalpy h_L(x) [kJ/kmol].
    #[pyo3(name = "h_liquid_of_x")]
    fn py_h_liquid_of_x(&self, x: f64) -> PyResult<f64> {
        self.h_liquid_of_x(x).map_err(to_py_err)
    }

    /// Interpolated saturated-vapor molar enthalpy H_V(y) [kJ/kmol].
    #[pyo3(name = "h_vapor_of_y")]
    fn py_h_vapor_of_y(&self, y: f64) -> PyResult<f64> {
        self.h_vapor_of_y(y).map_err(to_py_err)
    }

    /// Interpolated equilibrium vapor fraction y*(x).
    #[pyo3(name = "y_of_x")]
    fn py_y_of_x(&self, x: f64) -> PyResult<f64> {
        self.equilibrium().y_of_x(x).map_err(to_py_err)
    }

    /// Interpolated inverse x*(y).
    #[pyo3(name = "x_of_y")]
    fn py_x_of_y(&self, y: f64) -> PyResult<f64> {
        self.equilibrium().x_of_y(y).map_err(to_py_err)
    }

    fn __repr__(&self) -> String {
        format!(
            "EnthalpyCurve({} points)",
            self.equilibrium().x_samples().len()
        )
    }
}

#[pymethods]
impl Feed {
    /// Feed of `rate` kmol/h with light-component fraction `z` and thermal
    /// condition `q` (1 = saturated liquid, 0 = saturated vapor).
    #[new]
    #[pyo3(signature = (rate, z, q = 1.0))]
    fn py_new(rate: f64, z: f64, q: f64) -> Self {
        Self { rate, z, q }
    }

    fn __repr__(&self) -> String {
        format!(
            "Feed(rate={} kmol/h, z={}, q={})",
            self.rate, self.z, self.q
        )
    }
}

#[pymethods]
impl BinaryColumn {
    /// A binary two-product column at `pressure` kPa (absolute).
    #[new]
    #[pyo3(signature = (pressure, feed, x_distillate, x_bottoms, condenser = "total"))]
    fn py_new(
        pressure: f64,
        feed: Feed,
        x_distillate: f64,
        x_bottoms: f64,
        condenser: &str,
    ) -> PyResult<Self> {
        let col = Self {
            pressure,
            condenser: parse_condenser(condenser)?,
            feed,
            x_distillate,
            x_bottoms,
        };
        col.validate().map_err(to_py_err)?;
        Ok(col)
    }

    /// Distillate molar flow D [kmol/h] from the material balances.
    #[pyo3(name = "distillate_rate")]
    fn py_distillate_rate(&self) -> PyResult<f64> {
        self.distillate_rate().map_err(to_py_err)
    }

    /// Bottoms molar flow B = F − D [kmol/h].
    #[pyo3(name = "bottoms_rate")]
    fn py_bottoms_rate(&self) -> PyResult<f64> {
        self.bottoms_rate().map_err(to_py_err)
    }

    fn __repr__(&self) -> String {
        format!(
            "BinaryColumn(P={} kPa, zF={}, xD={}, xB={})",
            self.pressure, self.feed.z, self.x_distillate, self.x_bottoms
        )
    }
}

/// Minimum reflux by pinch detection (tangent pinches included).
///
/// Returns an `RminResult` with `r_min`, the controlling `pinch` point,
/// whether it is `tangent`, and the q-line/curve `feed_point`.
#[pyfunction]
#[pyo3(signature = (curve, x_distillate, x_bottoms, z_feed, q = 1.0))]
fn rmin(
    curve: &EquilibriumCurve,
    x_distillate: f64,
    x_bottoms: f64,
    z_feed: f64,
    q: f64,
) -> PyResult<RminResult> {
    mt::rmin(curve, x_distillate, x_bottoms, z_feed, q).map_err(to_py_err)
}

/// The full McCabe–Thiele construction. Compositions are light-component
/// mole fractions; `q` is the feed thermal condition; `murphree` the vapor
/// Murphree efficiency; `condenser` is "total" or "partial".
#[pyfunction]
#[pyo3(signature = (curve, x_distillate, x_bottoms, z_feed, reflux, q = 1.0, murphree = 1.0, condenser = "total"))]
#[allow(clippy::too_many_arguments)]
fn mccabe_thiele(
    curve: &EquilibriumCurve,
    x_distillate: f64,
    x_bottoms: f64,
    z_feed: f64,
    reflux: f64,
    q: f64,
    murphree: f64,
    condenser: &str,
) -> PyResult<McCabeThieleResult> {
    let spec = McCabeThieleSpec {
        x_distillate,
        x_bottoms,
        z_feed,
        q,
        reflux,
        murphree,
        condenser: parse_condenser(condenser)?,
    };
    mt::mccabe_thiele(curve, spec).map_err(to_py_err)
}

/// Total-reflux construction: the minimum stage count N_min.
#[pyfunction]
#[pyo3(signature = (curve, x_distillate, x_bottoms, murphree = 1.0))]
fn total_reflux(
    curve: &EquilibriumCurve,
    x_distillate: f64,
    x_bottoms: f64,
    murphree: f64,
) -> PyResult<TotalRefluxResult> {
    mt::total_reflux(curve, x_distillate, x_bottoms, murphree).map_err(to_py_err)
}

/// Stage count vs reflux ratio: returns (r, n) pairs, with n = NaN where the
/// construction failed (R at or below R_min) — the batch NaN-on-fail
/// contract.
#[pyfunction]
#[pyo3(signature = (curve, r_values, x_distillate, x_bottoms, z_feed, q = 1.0, murphree = 1.0, condenser = "total"))]
#[allow(clippy::too_many_arguments)]
fn n_vs_r(
    curve: &EquilibriumCurve,
    r_values: Vec<f64>,
    x_distillate: f64,
    x_bottoms: f64,
    z_feed: f64,
    q: f64,
    murphree: f64,
    condenser: &str,
) -> PyResult<Vec<(f64, f64)>> {
    let spec = McCabeThieleSpec {
        x_distillate,
        x_bottoms,
        z_feed,
        q,
        reflux: 1.0, // template value; replaced per sweep entry
        murphree,
        condenser: parse_condenser(condenser)?,
    };
    Ok(mt::n_vs_r(curve, spec, &r_values))
}

/// The full Ponchon–Savarit construction on an enthalpy–composition curve.
/// Compositions are light-component mole fractions; `q` is the feed thermal
/// condition; `condenser` is "total" (v1 supports total only). Returns a rich
/// `PonchonSavaritResult` (stages, poles, per-mole-feed duties, tie lines).
#[pyfunction]
#[pyo3(signature = (curve, x_distillate, x_bottoms, z_feed, reflux, q = 1.0, condenser = "total"))]
fn ponchon_savarit(
    curve: &EnthalpyCurve,
    x_distillate: f64,
    x_bottoms: f64,
    z_feed: f64,
    reflux: f64,
    q: f64,
    condenser: &str,
) -> PyResult<PonchonSavaritResult> {
    let spec = PonchonSavaritSpec {
        x_distillate,
        x_bottoms,
        z_feed,
        q,
        reflux,
        condenser: parse_condenser(condenser)?,
    };
    ps::ponchon_savarit(curve, spec).map_err(to_py_err)
}

/// Return the engine crate's version string (matches `Cargo.toml`).
#[pyfunction]
fn version() -> String {
    crate::version().to_string()
}

/// Bubble temperature (K) of equimolar methanol/water at 101.325 kPa — the M0
/// smoke path, kept as the cheapest cross-FFI health check. See
/// [`crate::thermo::smoke_bubble_temperature`].
#[pyfunction]
fn smoke_bubble_temperature() -> PyResult<f64> {
    crate::thermo::smoke_bubble_temperature().map_err(to_py_err)
}

/// The `stages._engine` native module. The function name must match the last
/// component of `module-name` in `python/pyproject.toml` (`stages._engine`).
#[pymodule]
fn _engine(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", crate::version())?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(smoke_bubble_temperature, m)?)?;
    m.add_function(wrap_pyfunction!(rmin, m)?)?;
    m.add_function(wrap_pyfunction!(mccabe_thiele, m)?)?;
    m.add_function(wrap_pyfunction!(total_reflux, m)?)?;
    m.add_function(wrap_pyfunction!(n_vs_r, m)?)?;
    m.add_function(wrap_pyfunction!(ponchon_savarit, m)?)?;
    m.add_class::<ThermoSystem>()?;
    m.add_class::<EquilibriumCurve>()?;
    m.add_class::<EnthalpyCurve>()?;
    m.add_class::<Feed>()?;
    m.add_class::<BinaryColumn>()?;
    m.add_class::<CondenserKind>()?;
    m.add_class::<McCabeThieleSpec>()?;
    m.add_class::<Line>()?;
    m.add_class::<StagePoint>()?;
    m.add_class::<RminResult>()?;
    m.add_class::<McCabeThieleResult>()?;
    m.add_class::<TotalRefluxResult>()?;
    m.add_class::<PonchonSavaritSpec>()?;
    m.add_class::<PonchonSavaritResult>()?;
    Ok(())
}
