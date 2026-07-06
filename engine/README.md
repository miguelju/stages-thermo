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
stages-thermo = "0.0"
```

> **Status:** `0.0.x` is a name-holding stub (Milestone 0, repo bootstrap). The
> only computational surface today is the vle-thermo adapter smoke path. Method
> modules land milestone by milestone — see the repo's `ROADMAP.md`. Do not
> depend on this for real work until `1.0`.

```rust
use stages_thermo::thermo::smoke_bubble_temperature;

// M0 smoke path: equimolar methanol/water bubble T at 101.325 kPa, via
// vle-thermo (Peng–Robinson, classical mixing). Approximate — proves the
// dependency links, not a validated value.
let t_bubble = smoke_bubble_temperature().unwrap();
```

The Python bindings (PyO3) are gated behind the `python` feature, which maturin
enables when building the wheel published to PyPI as `stages-thermo`
(`import stages`). `cargo add stages-thermo` gets a pure-Rust crate with no
PyO3 in its dependency closure.

## License

MIT © Miguel Roberto Jackson Ugueto
