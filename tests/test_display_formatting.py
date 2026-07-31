from webapp.display_formatting import lcoe_y_axis_settings


def test_lcoe_axis_remains_linear_for_ordinary_values():
    scale, ymin, ymax = lcoe_y_axis_settings([100, 250, 400])

    assert scale == "linear"
    assert ymin == 0
    assert abs(ymax - 460) < 1e-9


def test_lcoe_axis_switches_to_log_and_includes_extreme_curve():
    scale, ymin, ymax = lcoe_y_axis_settings([300, 900, 200_000])

    assert scale == "log"
    assert ymin == 10
    assert abs(ymax - 230_000) < 1e-9


def test_lcoe_axis_ignores_nonfinite_and_nonpositive_values():
    scale, ymin, ymax = lcoe_y_axis_settings(
        [float("nan"), float("inf"), -50, 0],
    )

    assert (scale, ymin, ymax) == ("linear", 0, 410)
