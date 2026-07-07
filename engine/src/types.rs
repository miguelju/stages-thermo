//! Core error and report types shared across the crate.
//!
//! This is intentionally thin at M0. As solvers land, `SolveReport`
//! (iteration count, residual-norm history, damping steps, homotopy path) and
//! the richer result objects described in `PLAN.md` §8 grow here alongside the
//! per-method result structs.

use thiserror::Error;

/// Errors surfaced by the stages-thermo engine.
///
/// Thermodynamic failures originating in `vle-thermo` are wrapped as
/// [`StagesError::Thermo`] with the upstream message preserved, so callers see
/// one error type regardless of which layer failed.
#[derive(Debug, Error)]
pub enum StagesError {
    /// A thermodynamic evaluation failed in the vle-thermo layer. The string
    /// carries the upstream error's `Display` text (vle-thermo's error enums
    /// are not part of this crate's public API surface).
    #[error("thermo evaluation failed: {0}")]
    Thermo(String),

    /// A dimension mismatch between, e.g., a composition vector and the
    /// component count.
    #[error("dimension mismatch: {0}")]
    Dimension(String),

    /// An iterative solver failed to converge within its iteration budget.
    #[error("solver did not converge: {0}")]
    Convergence(String),

    /// The specification cannot be met by any operating condition — e.g. a
    /// distillate purity beyond an azeotrope, reachable at no reflux ratio.
    /// Distinct from [`StagesError::Convergence`]: a convergence failure
    /// might yield to a better solver or seed; an infeasible spec never will.
    #[error("specification is infeasible: {0}")]
    Infeasible(String),
}

/// Convenience alias for fallible engine operations.
pub type Result<T> = std::result::Result<T, StagesError>;
