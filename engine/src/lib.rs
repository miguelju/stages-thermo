//! # stages-thermo — a staged-separation (distillation) column solver
//!
//! `stages-thermo` is a learning library **and** a fast steady-state column
//! engine for staged separations. It walks the full pedagogical ladder —
//! McCabe–Thiele → Ponchon–Savarit → Fenske–Underwood–Gilliland shortcut →
//! rigorous MESH (Wang–Henke, Burningham–Otto, Naphtali–Sandholm) — with every
//! method implemented from scratch and anchored to its textbook equations, and
//! it exposes a granular, batch-capable API ("numpy for distillation columns").
//!
//! The thermodynamics — K-values, enthalpies, and their derivatives — come
//! **entirely** from [`vle_thermo`](https://crates.io/crates/vle-thermo). This
//! crate adds no thermodynamics of its own; it is the first downstream consumer
//! of vle-thermo. See `PLAN.md` §1 for the full vision.
//!
//! ## Status
//!
//! Milestone 0 (repo bootstrap). The only computational surface today is the
//! [`thermo`] adapter smoke path, which proves the vle-thermo dependency links
//! and evaluates end-to-end. The method modules (`binary`, `shortcut`,
//! `rigorous`, `numerics`, `column`) land milestone by milestone per
//! `ROADMAP.md`; the target module tree is documented in `PLAN.md` §6.
//!
//! ## Design rule (PLAN §7)
//!
//! **Exactly one module talks to vle-thermo:** [`thermo`]. Everything else
//! consumes it through a provider interface, so a surrogate thermo model (the
//! inside-out inner loop, a stretch milestone) or a test mock can slot in
//! without touching the solvers.

pub mod thermo;
pub mod types;

// Python bindings are compiled only with the `python` feature (enabled by
// maturin when building the wheel). `cargo add stages-thermo` gets a pure-Rust
// crate with no PyO3 in its dependency closure.
#[cfg(feature = "python")]
mod py_bindings;

/// Return this crate's version string (matches `Cargo.toml`).
///
/// Mirrors vle-thermo's `version()` so downstream Python can introspect which
/// engine build is loaded.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_semver_shaped() {
        let v = version();
        let parts: Vec<&str> = v.split('.').collect();
        assert!(parts.len() >= 3, "version {v} is not semver-shaped");
        assert!(parts[0].parse::<u64>().is_ok(), "major of {v} not numeric");
    }
}
