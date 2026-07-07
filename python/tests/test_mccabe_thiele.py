"""M1 binding tests: the McCabe–Thiele surface through the built wheel.

Every M1 PyO3 binding is exercised at least once (CLAUDE.md "PyO3 Bindings
Rule"): ThermoSystem, EquilibriumCurve (all three constructors + queries),
mccabe_thiele, rmin, total_reflux, n_vs_r, Feed, BinaryColumn, and the rich
result objects' fields. Deeper numerical validation lives in the Rust test
suite and the notebook's pinned assertion cells; these tests pin the physics
lightly and the API contract firmly.
"""

import math

import pytest

import stages


# --- ThermoSystem ---------------------------------------------------------


def test_thermo_system_peng_robinson() -> None:
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    assert sys.components == ["benzene", "toluene"]
    t, y, k = sys.bubble_temperature(101.325, [0.5, 0.5])
    # Between the two boiling points, benzene enriched in the vapor.
    assert 353.0 < t < 384.0
    assert y[0] > 0.5
    assert k[0] > 1.0 > k[1]


def test_thermo_system_van_laar_positive_deviation() -> None:
    sys = stages.ThermoSystem.van_laar(["methanol", "water"], 0.5853, 0.3458)
    p_mid, _, _ = sys.bubble_pressure(298.15, [0.5, 0.5])
    p_m, _, _ = sys.bubble_pressure(298.15, [1.0, 0.0])
    p_w, _, _ = sys.bubble_pressure(298.15, [0.0, 1.0])
    assert p_mid > 0.5 * p_m + 0.5 * p_w  # positive deviation from Raoult


def test_unknown_component_raises() -> None:
    with pytest.raises(RuntimeError, match="unobtainium"):
        stages.ThermoSystem.peng_robinson(["benzene", "unobtainium"])


# --- EquilibriumCurve ------------------------------------------------------


def test_curve_from_thermo() -> None:
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    curve = stages.EquilibriumCurve.from_thermo(sys, 101.325, n_points=51)
    assert curve.pressure == 101.325
    assert len(curve.x) == len(curve.y) == len(curve.t) == 51
    assert curve.x[0] == 0.0 and curve.x[-1] == 1.0
    assert curve.y[0] == 0.0 and curve.y[-1] == 1.0
    # Benzene is the light component: curve above the diagonal.
    assert curve.y_of_x(0.5) > 0.5
    # Temperature falls from toluene's bp toward benzene's.
    assert curve.t[0] > curve.t[-1]
    assert 2.0 < curve.relative_volatility(0.5) < 3.0


def test_curve_component_order_check() -> None:
    sys = stages.ThermoSystem.peng_robinson(["toluene", "benzene"])
    with pytest.raises(RuntimeError, match="more volatile"):
        stages.EquilibriumCurve.from_thermo(sys, 101.325)


def test_curve_constant_alpha_and_inverse() -> None:
    curve = stages.EquilibriumCurve.constant_alpha(2.5, n_points=201)
    assert curve.pressure is None
    x = 0.4
    y = curve.y_of_x(x)
    assert math.isclose(y, 2.5 * x / (1 + 1.5 * x), abs_tol=5e-5)
    assert math.isclose(curve.x_of_y(y), x, abs_tol=1e-9)
    with pytest.raises(ValueError):
        curve.temperature_of_x(0.5)  # synthetic curve has no T data


def test_curve_from_points() -> None:
    xs = [0.0, 0.25, 0.5, 0.75, 1.0]
    ys = [0.0, 0.45, 0.7, 0.88, 1.0]
    curve = stages.EquilibriumCurve.from_points(xs, ys, pressure=101.325)
    assert curve.y_of_x(0.25) == 0.45
    with pytest.raises(ValueError):
        stages.EquilibriumCurve.from_points(xs, list(reversed(ys)))


# --- McCabe–Thiele ---------------------------------------------------------


@pytest.fixture(scope="module")
def bt_curve() -> "stages.EquilibriumCurve":
    sys = stages.ThermoSystem.peng_robinson(["benzene", "toluene"])
    return stages.EquilibriumCurve.from_thermo(sys, 101.325, n_points=101)


def test_rmin(bt_curve) -> None:
    r = stages.rmin(bt_curve, x_distillate=0.95, x_bottoms=0.05, z_feed=0.5, q=1.0)
    assert 0.8 < r.r_min < 1.5  # textbook ballpark for benzene–toluene 95/5
    assert not r.tangent  # concave curve → feed pinch
    px, py = r.pinch
    assert math.isclose(px, 0.5, abs_tol=1e-6)  # saturated liquid: pinch at z_F
    assert py > px


def test_mccabe_thiele_construction(bt_curve) -> None:
    r_min = stages.rmin(bt_curve, 0.95, 0.05, 0.5).r_min
    res = stages.mccabe_thiele(
        bt_curve,
        x_distillate=0.95,
        x_bottoms=0.05,
        z_feed=0.5,
        reflux=1.5 * r_min,
        q=1.0,
    )
    n_min = stages.total_reflux(bt_curve, 0.95, 0.05).n_min
    assert n_min < res.n_stages < 25.0
    assert 1 < res.feed_stage < len(res.stages)
    # Rich result object: geometry present and consistent.
    assert res.staircase[0] == (0.95, 0.95)
    assert res.stages[-1].x <= 0.05 + 1e-12
    xi, yi = res.intersection
    assert math.isclose(res.rectifying.slope * xi + res.rectifying.intercept, yi, abs_tol=1e-12)
    assert math.isclose(res.stripping.slope * xi + res.stripping.intercept, yi, abs_tol=1e-9)
    assert res.spec.reflux == pytest.approx(1.5 * r_min)
    # Murphree < 1 needs more stages.
    res_e = stages.mccabe_thiele(
        bt_curve, 0.95, 0.05, 0.5, reflux=1.5 * r_min, murphree=0.7
    )
    assert res_e.n_stages > res.n_stages


def test_below_rmin_raises(bt_curve) -> None:
    r_min = stages.rmin(bt_curve, 0.95, 0.05, 0.5).r_min
    with pytest.raises(RuntimeError, match="R_min"):
        stages.mccabe_thiele(bt_curve, 0.95, 0.05, 0.5, reflux=0.8 * r_min)


def test_n_vs_r_nan_contract(bt_curve) -> None:
    r_min = stages.rmin(bt_curve, 0.95, 0.05, 0.5).r_min
    pairs = stages.n_vs_r(
        bt_curve, [0.5 * r_min, 1.5 * r_min, 3.0 * r_min], 0.95, 0.05, 0.5
    )
    assert len(pairs) == 3
    assert math.isnan(pairs[0][1])  # below R_min → NaN, not an exception
    assert pairs[1][1] > pairs[2][1] > 0  # N falls with R


def test_condenser_string_validation(bt_curve) -> None:
    with pytest.raises(ValueError, match="condenser"):
        stages.mccabe_thiele(bt_curve, 0.95, 0.05, 0.5, reflux=2.0, condenser="magic")


# --- Column model ----------------------------------------------------------


def test_binary_column_balances() -> None:
    col = stages.BinaryColumn(
        pressure=101.325,
        feed=stages.Feed(rate=100.0, z=0.5, q=1.0),
        x_distillate=0.95,
        x_bottoms=0.05,
    )
    d, b = col.distillate_rate(), col.bottoms_rate()
    assert math.isclose(d + b, 100.0, abs_tol=1e-9)
    assert math.isclose(0.95 * d + 0.05 * b, 50.0, abs_tol=1e-9)


def test_binary_column_validates() -> None:
    with pytest.raises(ValueError):
        stages.BinaryColumn(
            pressure=101.325,
            feed=stages.Feed(rate=100.0, z=0.5),
            x_distillate=0.4,  # below z_F
            x_bottoms=0.05,
        )


# --- Plotting (headless smoke) ---------------------------------------------


def test_plotting_smoke(bt_curve) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from stages import plotting

    res = stages.mccabe_thiele(bt_curve, 0.95, 0.05, 0.5, reflux=2.0)
    ax = plotting.plot_mccabe_thiele(res, bt_curve, show_rmin=True)
    assert ax.get_title().startswith("McCabe–Thiele")
    tr = stages.total_reflux(bt_curve, 0.95, 0.05)
    ax2 = plotting.plot_total_reflux(tr, bt_curve)
    assert "Total reflux" in ax2.get_title()
    import matplotlib.pyplot as plt

    plt.close("all")
