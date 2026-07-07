//! The column domain model — the objects every method shares.
//!
//! At M1 this holds the **binary-sufficient subset** (PLAN §6): a feed, a
//! condenser kind, and a two-product binary column with its overall and
//! light-component material balances. The full multicomponent model —
//! multi-feed, side draws, per-stage duties, pressure profiles, the spec
//! system — arrives at M5 (MESH infrastructure) and grows in this module.

pub mod model;

pub use model::{BinaryColumn, CondenserKind, Feed};
