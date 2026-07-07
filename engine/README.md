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

> **Status:** `0.1.x` ships the first rung of the ladder — the binary
> McCabe–Thiele layer (Milestone 1): equilibrium curves from real
> thermodynamics, minimum reflux by geometric pinch detection (tangent pinches
> included), stage stepping with Murphree efficiency, total reflux, N(R), and
> the binary column material balances. Multicomponent and rigorous MESH
> solvers land milestone by milestone — see the repo's `ROADMAP.md`. The API
> may still move before `1.0`.

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

The same construction runs on a γ-φ activity-model curve
(`ThermoSystem::van_laar`), the constant-α idealization
(`EquilibriumCurve::constant_alpha` — on which Fenske's and Underwood's closed
forms are exact, which is how the stepping kernels are unit-tested), or raw
literature data (`EquilibriumCurve::from_points`).

The Python bindings (PyO3) are gated behind the `python` feature, which maturin
enables when building the wheel published to PyPI as `stages-thermo`
(`import stages`). `cargo add stages-thermo` gets a pure-Rust crate with no
PyO3 in its dependency closure.

## License

MIT © Miguel Roberto Jackson Ugueto
