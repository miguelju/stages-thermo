//! Binary column model — the M1 subset.
//!
//! Stages are numbered **top-down**: 1 = condenser, N = reboiler (the
//! repo-wide convention, CLAUDE.md "Domain Context"). Compositions refer to
//! the **light** (more volatile) component throughout the binary layer, the
//! standard McCabe–Thiele convention (Seader, Henley & Roper, *Separation
//! Process Principles*, Ch. 7).
//!
//! ## Textbook anchoring
//!
//! With `F`, `D`, `B` the feed/distillate/bottoms molar flows and `z_F`,
//! `x_D`, `x_B` their light-component mole fractions, the two steady-state
//! material balances are
//!
//! ```text
//! overall:  F = D + B                        (S&H eq. 7-2)
//! light:    F z_F = D x_D + B x_B           (S&H eq. 7-3)
//! ⇒ D = F (z_F − x_B) / (x_D − x_B)
//! ```

use crate::types::{Result, StagesError};

/// The condenser at the top of the column (stage 1).
///
/// A **total** condenser condenses the whole overhead vapor, so the reflux
/// and distillate share the vapor's composition and the condenser is *not* an
/// equilibrium stage — stepping starts on the diagonal at `x_D`. A
/// **partial** condenser condenses only the reflux and takes the distillate
/// as vapor; the vapor–liquid split makes it one extra equilibrium stage
/// (S&H §7.2).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[cfg_attr(feature = "python", pyo3::pyclass(eq, eq_int))]
pub enum CondenserKind {
    /// Condense everything; distillate is liquid at composition `x_D`.
    #[default]
    Total,
    /// Condense only the reflux; distillate is vapor at composition `x_D`,
    /// and the condenser counts as an equilibrium stage.
    Partial,
}

/// A single feed to a binary column.
///
/// `q` is the **thermal condition** — the fraction of the feed that joins
/// the liquid flowing down the column (S&H eq. 7-18):
///
/// - `q = 1`: saturated liquid (bubble point)
/// - `q = 0`: saturated vapor (dew point)
/// - `0 < q < 1`: two-phase feed
/// - `q > 1`: subcooled liquid; `q < 0`: superheated vapor
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct Feed {
    /// Feed molar flow `F` in **kmol/h**.
    pub rate: f64,
    /// Light-component mole fraction `z_F` of the feed (0–1).
    pub z: f64,
    /// Thermal condition `q` (dimensionless, see type docs).
    pub q: f64,
}

/// A binary two-product column: one feed, a distillate and a bottoms product,
/// a single column pressure — everything McCabe–Thiele needs.
///
/// This is the M1 "binary-sufficient subset" of the full column model
/// (PLAN §6); the multicomponent `Column` with multi-feed, side draws,
/// duties, and pressure profiles lands at M5.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "python", pyo3::pyclass(get_all))]
pub struct BinaryColumn {
    /// Column pressure in **kPa** (absolute), taken uniform across stages at
    /// this rung of the ladder.
    pub pressure: f64,
    /// Condenser kind (total condensers are the McCabe–Thiele default).
    pub condenser: CondenserKind,
    /// The single feed.
    pub feed: Feed,
    /// Light-component mole fraction of the distillate, `x_D`.
    pub x_distillate: f64,
    /// Light-component mole fraction of the bottoms, `x_B`.
    pub x_bottoms: f64,
}

impl BinaryColumn {
    /// Validate the specification: fractions in (0, 1), the separation
    /// ordered `x_B < z_F < x_D`, and a positive feed rate.
    pub fn validate(&self) -> Result<()> {
        let Feed { rate, z, q: _ } = self.feed;
        // NaN-safe positivity guard (a plain `> 0.0` inverted with `!` trips
        // clippy's neg_cmp_op_on_partial_ord; the NaN check is explicit).
        if rate.is_nan() || rate <= 0.0 {
            return Err(StagesError::Dimension(format!(
                "feed rate must be positive, got {rate} kmol/h"
            )));
        }
        for (name, v) in [
            ("z_F", z),
            ("x_D", self.x_distillate),
            ("x_B", self.x_bottoms),
        ] {
            if !(0.0..=1.0).contains(&v) {
                return Err(StagesError::Dimension(format!(
                    "{name} must be a mole fraction in [0, 1], got {v}"
                )));
            }
        }
        if !(self.x_bottoms < z && z < self.x_distillate) {
            return Err(StagesError::Dimension(format!(
                "specs must satisfy x_B < z_F < x_D, got x_B = {}, z_F = {z}, x_D = {}",
                self.x_bottoms, self.x_distillate
            )));
        }
        Ok(())
    }

    /// Distillate molar flow `D` in **kmol/h** from the material balances
    /// (S&H eqs. 7-2/7-3): `D = F (z_F − x_B) / (x_D − x_B)`.
    pub fn distillate_rate(&self) -> Result<f64> {
        self.validate()?;
        Ok(self.feed.rate * (self.feed.z - self.x_bottoms) / (self.x_distillate - self.x_bottoms))
    }

    /// Bottoms molar flow `B = F − D` in **kmol/h**.
    pub fn bottoms_rate(&self) -> Result<f64> {
        Ok(self.feed.rate - self.distillate_rate()?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn column() -> BinaryColumn {
        BinaryColumn {
            pressure: 101.325,
            condenser: CondenserKind::Total,
            feed: Feed {
                rate: 100.0,
                z: 0.5,
                q: 1.0,
            },
            x_distillate: 0.95,
            x_bottoms: 0.05,
        }
    }

    /// The symmetric 50/50 → 95/5 split puts half the feed in each product,
    /// and D + B recovers F exactly.
    #[test]
    fn material_balance_closes() {
        let c = column();
        let d = c.distillate_rate().unwrap();
        let b = c.bottoms_rate().unwrap();
        assert!((d - 50.0).abs() < 1e-12, "D = {d}");
        assert!((d + b - c.feed.rate).abs() < 1e-12);
        // Light-component balance: F z = D xD + B xB.
        let light = d * c.x_distillate + b * c.x_bottoms;
        assert!((light - c.feed.rate * c.feed.z).abs() < 1e-12);
    }

    /// Mis-ordered specs are rejected with a Dimension error.
    #[test]
    fn misordered_specs_rejected() {
        let mut c = column();
        c.x_bottoms = 0.6; // above z_F
        assert!(matches!(
            c.distillate_rate(),
            Err(StagesError::Dimension(_))
        ));
    }
}
