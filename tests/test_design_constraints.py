from webapp.design_constraints import safe_height_interval


def test_safe_height_interval_enforces_exclusive_thirty_year_limit():
    interval = safe_height_interval(
        lambda height: height * 365,
        1,
        40,
    )

    assert interval == (1, 29)


def test_safe_height_interval_skips_subcritical_and_short_lifetimes():
    interval = safe_height_interval(
        lambda height: height * 30,
        1,
        10,
    )

    assert interval == (3, 10)


def test_safe_height_interval_returns_none_without_a_usable_height():
    assert safe_height_interval(lambda _height: 0, 10, 20) is None


def test_safe_height_interval_uses_a_contiguous_slider_safe_run():
    lifetimes = {
        1: 100,
        2: 200,
        3: 20_000,
        4: 300,
        5: 400,
        6: 500,
    }

    assert safe_height_interval(lifetimes.__getitem__, 1, 6) == (4, 6)
