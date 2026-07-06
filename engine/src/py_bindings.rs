//! Python bindings for the stages-thermo engine.
//!
//! Builds a Python extension module `stages._engine` when this crate is
//! compiled with the `python` feature. End-users import the higher-level
//! `stages` Python wrapper, which re-exports the pieces they need.
//!
//! ## What lives here at M0
//!
//! The bootstrap surface exists to prove the PyO3 boundary end-to-end (wheel
//! builds, installs, imports, and the smoke test passes on every CI platform):
//!
//! - [`version()`] — the crate's semver string.
//! - [`smoke_bubble_temperature()`] — the M0 vle-thermo smoke path, so the
//!   dependency is exercised across the FFI boundary from Python too.
//!
//! ## Adding bindings in M1+
//!
//! Per `CLAUDE.md`'s "PyO3 Bindings Rule", every milestone that adds public
//! Rust functionality also exposes it here in the same commit series, with at
//! least one round-trip test in `python/tests/`. CI runs the Python tests on
//! every wheel, so a missing binding is a hard failure, not a review oversight.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

/// Return the engine crate's version string (matches `Cargo.toml`).
#[pyfunction]
fn version() -> String {
    crate::version().to_string()
}

/// Bubble temperature (K) of equimolar methanol/water at 101.325 kPa — the M0
/// smoke path. See [`crate::thermo::smoke_bubble_temperature`].
#[pyfunction]
fn smoke_bubble_temperature() -> PyResult<f64> {
    crate::thermo::smoke_bubble_temperature().map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// The `stages._engine` native module. The function name must match the last
/// component of `module-name` in `python/pyproject.toml` (`stages._engine`).
#[pymodule]
fn _engine(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", crate::version())?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(smoke_bubble_temperature, m)?)?;
    Ok(())
}
