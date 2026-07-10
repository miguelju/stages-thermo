"""M2 binding tests: the Ponchon–Savarit surface through the built wheel.

Every new M2 PyO3 binding is exercised at least once (CLAUDE.md "PyO3 Bindings
Rule"): ThermoSystem.nrtl / .phase_enthalpy / .t_ref / .p_ref, EnthalpyCurve
(both constructors + accessors + queries), ponchon_savarit, and the
PonchonSavaritResult / PonchonSavaritSpec fields. Deeper numerical validation
lives in the Rust suite and the notebook's pinned assertion cells.
"""

import math

import pytest

import stages


# --- ThermoSystem: NRTL + enthalpy ----------------------------------------


def test_nrtl_constructor_positive_deviation() -> None:
    sys = stages.ThermoSystem.nrtl(["ethanol", "water"], -458.7, 5574.0, 0.303)
    assert sys.components == ["ethanol", "water"]
    p_mid, _, _ = sys.bubble_pressure(298.15, [0.5, 0.5])
    p_e, _, _ = sys.bubble_pressure(298.15, [1.0, 0.0])
    p_w, _, _ = sys.bubble_pressure(298.15, [0.0, 1.0])
    assert p_mid > 0.5 * p_e + 0.5 * p_w  # positive deviation from Raoult


def test_reference_state_getters() -> None:
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    assert sys.t_ref == pytest.approx(298.15)
    assert sys.p_ref == pytest.approx(101.325)


def test_phase_enthalpy_latent_heat_positive() -> None:
    sys = stages.ThermoSystem.van_laar(["methanol", "water"], 0.5853, 0.3458)
    x = [0.4, 0.6]
    t, y, _ = sys.bubble_temperature(101.325, x)
    h_liq = sys.phase_enthalpy(t, 101.325, x, "liquid")
    h_vap = sys.phase_enthalpy(t, 101.325, y, "vapor")
    assert h_vap > h_liq
    assert 10_000.0 < h_vap - h_liq < 80_000.0  # kJ/kmol latent-heat window


def test_phase_enthalpy_bad_phase_raises() -> None:
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    with pytest.raises(ValueError, match="phase"):
        sys.phase_enthalpy(360.0, 101.325, [0.5, 0.5], "plasma")


# --- EnthalpyCurve ---------------------------------------------------------


def test_enthalpy_curve_from_thermo() -> None:
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    ec = stages.EnthalpyCurve.from_thermo(sys, 101.325, n_points=51)
    assert len(ec.x) == len(ec.y) == len(ec.t) == 51
    assert len(ec.h_liq) == len(ec.h_vap) == 51
    # Saturated vapor lies above saturated liquid in the interior.
    for i in range(1, 50):
        assert ec.h_vap[i] > ec.h_liq[i]
    # Interpolation queries round-trip a grid sample.
    xi = ec.x[10]
    assert ec.h_liquid_of_x(xi) == pytest.approx(ec.h_liq[10], rel=1e-9)
    assert ec.y_of_x(0.5) > 0.5  # benzene is the light component


def test_enthalpy_curve_from_points() -> None:
    x = [0.0, 0.25, 0.5, 0.75, 1.0]
    y = [0.0, 0.4, 0.6, 0.8, 1.0]
    h_liq = [0.0, 1000.0, 2000.0, 3000.0, 4000.0]
    h_vap = [30000.0, 31000.0, 32000.0, 33000.0, 34000.0]
    ec = stages.EnthalpyCurve.from_points(x, y, h_liq, h_vap)
    assert ec.h_liq == h_liq
    assert ec.h_vap == h_vap
    assert ec.h_liquid_of_x(0.5) == pytest.approx(2000.0)
    # Mismatched enthalpy length raises.
    with pytest.raises(ValueError):
        stages.EnthalpyCurve.from_points(x, y, [0.0, 1.0], h_vap)


# --- ponchon_savarit -------------------------------------------------------


def test_ponchon_savarit_construction() -> None:
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    ec = stages.EnthalpyCurve.from_thermo(sys, 101.325, n_points=201)
    res = stages.ponchon_savarit(
        ec, x_distillate=0.95, x_bottoms=0.05, z_feed=0.5, reflux=2.0
    )
    # Rich result object.
    assert 1.0 < res.n_stages < 40.0
    assert 1 <= res.feed_stage <= len(res.stages)
    assert res.spec.reflux == 2.0
    assert res.spec.x_distillate == 0.95
    # Poles bracket the diagram.
    assert res.delta_d[0] == pytest.approx(0.95)
    assert res.delta_b[0] == pytest.approx(0.05)
    assert res.delta_d[1] > ec.h_vapor_of_y(0.95)
    assert res.delta_b[1] < ec.h_liquid_of_x(0.05)
    # Duties physically signed, tie lines one per stage.
    assert res.q_condenser > 0.0
    assert res.q_reboiler > 0.0
    assert len(res.tie_lines) == len(res.stages)
    # Feed point at z_F.
    assert res.feed_point[0] == pytest.approx(0.5)


def test_ponchon_savarit_energy_closure() -> None:
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    ec = stages.EnthalpyCurve.from_thermo(sys, 101.325, n_points=201)
    res = stages.ponchon_savarit(
        ec, x_distillate=0.95, x_bottoms=0.05, z_feed=0.5, reflux=2.0
    )
    xd, xb, zf = 0.95, 0.05, 0.5
    d = (zf - xb) / (xd - xb)
    b = (xd - zf) / (xd - xb)
    h_f = res.feed_point[1]
    lhs = h_f + res.q_reboiler
    rhs = d * ec.h_liquid_of_x(xd) + b * ec.h_liquid_of_x(xb) + res.q_condenser
    assert math.isclose(lhs, rhs, rel_tol=1e-6, abs_tol=1e-3)


def test_ponchon_savarit_partial_condenser_rejected() -> None:
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    ec = stages.EnthalpyCurve.from_thermo(sys, 101.325, n_points=51)
    with pytest.raises(ValueError, match="total condenser"):
        stages.ponchon_savarit(
            ec,
            x_distillate=0.95,
            x_bottoms=0.05,
            z_feed=0.5,
            reflux=2.0,
            condenser="partial",
        )


def _plain_curve(ec: "stages.EnthalpyCurve") -> "stages.EquilibriumCurve":
    """Rebuild the plain equilibrium curve from an enthalpy curve's samples so
    the McCabe–Thiele functions (which take an EquilibriumCurve) share the exact
    same x/y/T grid."""
    return stages.EquilibriumCurve.from_points(
        list(ec.x), list(ec.y), list(ec.t), None
    )


def test_ponchon_savarit_vs_mccabe_thiele_benzene_toluene() -> None:
    """Near-ideal: the two methods agree within ~1 stage on the same spec."""
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    ec = stages.EnthalpyCurve.from_thermo(sys, 101.325, n_points=401)
    eq = _plain_curve(ec)
    r = 1.5 * stages.rmin(eq, 0.95, 0.05, 0.5).r_min
    ps = stages.ponchon_savarit(
        ec, x_distillate=0.95, x_bottoms=0.05, z_feed=0.5, reflux=r
    )
    mt = stages.mccabe_thiele(
        eq, x_distillate=0.95, x_bottoms=0.05, z_feed=0.5, reflux=r
    )
    assert abs(ps.n_stages - mt.n_stages) < 1.2
