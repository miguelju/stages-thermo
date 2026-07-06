# stages-thermo

**A staged-separation (distillation) learning library and fast steady-state
column solver, built on [`vle-thermo`](https://pypi.org/project/vle-thermo/).**

`import stages` — walk the full pedagogical ladder of column methods
(McCabe–Thiele → Ponchon–Savarit → Fenske–Underwood–Gilliland shortcut →
rigorous MESH), each implemented from scratch and anchored to its textbook
equations, with a granular, batch-capable API ("numpy for distillation
columns"). All thermodynamics come from `vle-thermo`; this package adds none of
its own.

```sh
pip install stages-thermo      # import name is `stages`
```

> **Status:** `0.0.x` is a name-holding stub (Milestone 0, repo bootstrap). The
> only surface today is the vle-thermo smoke path. Column methods land
> milestone by milestone — see the repo's `ROADMAP.md`. Not usable for real
> work until `1.0`.

```python
import stages

stages.__version__
# '0.0.1'

# M0 smoke path: equimolar methanol/water bubble T [K] at 101.325 kPa, computed
# through vle-thermo. Approximate — proves the dependency links end-to-end,
# not a validated value.
stages.smoke_bubble_temperature()
```

The native core is a Rust crate (`stages-thermo` on crates.io) with PyO3
bindings; wheels are abi3 (`cp310-abi3-*`), so one wheel per (OS, arch) covers
CPython 3.10+.

## License

MIT © Miguel Roberto Jackson Ugueto
