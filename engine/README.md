# stages-thermo (engine)

**A staged-separation (distillation) learning library and fast steady-state
column solver, built on [`vle-thermo`](https://crates.io/crates/vle-thermo).**

This is the Rust core crate. It walks the full pedagogical ladder of column
methods — McCabe–Thiele → Ponchon–Savarit → Fenske–Underwood–Gilliland
shortcut → rigorous MESH (Wang–Henke, Burningham–Otto, Naphtali–Sandholm) —
each implemented from scratch and anchored to its textbook equations, and
exposes a granular, batch-capable API ("numpy for distillation columns").

All thermodynamics (K-values, enthalpies, derivatives) come from `vle-thermo`;
this crate adds none of its own.

```toml
[dependencies]
stages-thermo = "0.1"
```

> **Status:** `0.2.x` ships the first two rungs of the ladder — the binary
> **McCabe–Thiele** (Milestone 1) and **Ponchon–Savarit** (Milestone 2) layers:
> equilibrium and enthalpy–composition (H–x–y) curves from real thermodynamics,
> minimum reflux by geometric pinch detection (tangent pinches included), stage
> stepping with Murphree efficiency, total reflux, N(R), the energy-exact
> difference-point construction (with condenser/reboiler duties), the NRTL γ-φ
> model, per-phase enthalpies, and the binary column material balances.
> Multicomponent and rigorous MESH solvers land milestone by milestone — see the
> repo's `ROADMAP.md`. The API may still move before `1.0`.

```rust
use stages_thermo::binary::mccabe_thiele::McCabeThieleSpec;
use stages_thermo::binary::{EquilibriumCurve, mccabe_thiele, rmin};
use stages_thermo::column::CondenserKind;
use stages_thermo::thermo::ThermoSystem;

// Peng–Robinson benzene–toluene at 1 atm, components from vle-thermo's
// built-in database (light component first). Units: K, kPa absolute.
let system = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
let curve = EquilibriumCurve::from_thermo(&system, 101.325, 101).unwrap();

// Minimum reflux by pinch detection, then the full construction at 1.5·R_min.
let r = rmin(&curve, 0.95, 0.05, 0.50, 1.0).unwrap();
let design = mccabe_thiele(
    &curve,
    McCabeThieleSpec {
        x_distillate: 0.95,
        x_bottoms: 0.05,
        z_feed: 0.50,
        q: 1.0,
        reflux: 1.5 * r.r_min,
        murphree: 1.0,
        condenser: CondenserKind::Total,
    },
)
.unwrap();

// Rich result objects, never bare numbers: stages, staircase polyline,
// operating lines, pinch analysis — everything on the diagram is queryable.
println!(
    "N = {:.2} stages, feed stage {}, R_min = {:.3} (tangent pinch: {})",
    design.n_stages, design.feed_stage, r.r_min, r.tangent
);
```

Rung 2 — **Ponchon–Savarit** — steps the same-shaped staircase on the
enthalpy–composition (H–x–y) diagram, closing the energy balance through two
difference points and returning the condenser/reboiler duties:

```rust
use stages_thermo::binary::equilibrium::EnthalpyCurve;
use stages_thermo::binary::ponchon_savarit::{PonchonSavaritSpec, ponchon_savarit};
use stages_thermo::column::CondenserKind;
use stages_thermo::thermo::ThermoSystem;

let system = ThermoSystem::peng_robinson(&["benzene", "toluene"]).unwrap();
// The H–x–y curve: saturated-liquid and -vapor enthalpies alongside y*(x).
let ec = EnthalpyCurve::from_thermo(&system, 101.325, 201).unwrap();
let ps = ponchon_savarit(
    &ec,
    PonchonSavaritSpec {
        x_distillate: 0.95,
        x_bottoms: 0.05,
        z_feed: 0.50,
        q: 1.0,
        reflux: 1.75,
        condenser: CondenserKind::Total,
    },
)
.unwrap();
println!(
    "N = {:.2} stages; Q_C/F = {:.0}, Q_R/F = {:.0} kJ/kmol feed",
    ps.n_stages, ps.q_condenser, ps.q_reboiler
);
```

The same constructions run on a γ-φ activity-model curve (`ThermoSystem::van_laar`
or the NRTL `ThermoSystem::nrtl`), the constant-α idealization
(`EquilibriumCurve::constant_alpha` — on which Fenske's and Underwood's closed
forms are exact, which is how the stepping kernels are unit-tested), or raw
literature data (`EquilibriumCurve::from_points` / `EnthalpyCurve::from_points`).

The Python bindings (PyO3) are gated behind the `python` feature, which maturin
enables when building the wheel published to PyPI as `stages-thermo`
(`import stages`). `cargo add stages-thermo` gets a pure-Rust crate with no
PyO3 in its dependency closure.

## License

MIT © Miguel Roberto Jackson Ugueto
