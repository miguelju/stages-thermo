//! Integration smoke test — exercises the crate through its public API only,
//! the same way a downstream `cargo add stages-thermo` consumer would.
//!
//! This is the M0 correctness anchor: it proves the vle-thermo dependency
//! resolves from crates.io and evaluates end-to-end. Real validation (the
//! Seader/Henley literature tables, cross-simulator checks) lands with each
//! method milestone; see `PLAN.md` §9.

use stages_thermo::thermo::smoke_bubble_temperature;
use stages_thermo::version;

#[test]
fn version_matches_cargo_pkg() {
    assert_eq!(version(), env!("CARGO_PKG_VERSION"));
}

#[test]
fn methanol_water_bubble_temperature_evaluates() {
    // vle-thermo is reached through the adapter; a converged, physically
    // plausible bubble temperature proves the whole dependency path works.
    let t = smoke_bubble_temperature().expect("bubble temperature should converge");
    assert!(
        (280.0..400.0).contains(&t),
        "methanol/water bubble T = {t} K outside the plausible window [280, 400)"
    );
}
