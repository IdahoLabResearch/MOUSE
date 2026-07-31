"""Display-only formatting helpers for MOUSE economic results."""

from __future__ import annotations

import math
from typing import Optional


def round_cost_for_display(value: float) -> float:
    """Round a cost so about half of its integer digits are trailing zeros.

    The underlying economic calculations retain full precision.  This helper
    only reduces the precision shown to users.  Examples:

    * 1,234,345 -> 1,234,000
    * 5,467 -> 5,500
    """
    number = float(value)
    if not math.isfinite(number) or number == 0.0:
        return number
    integer_digits = max(1, int(math.floor(math.log10(abs(number)))) + 1)
    trailing_zero_digits = integer_digits // 2
    return float(round(number, -trailing_zero_digits))


def format_cost_for_display(
    value: Optional[float],
    *,
    currency_symbol: str = "$",
) -> str:
    """Format a cost using the reduced-precision display convention."""
    if value is None:
        return "N/A"
    number = float(value)
    if not math.isfinite(number):
        return "N/A"
    rounded = round_cost_for_display(number)
    return f"{currency_symbol}{rounded:,.0f}"


def lcoe_y_axis_settings(
    upper_band_values,
    *,
    linear_floor: float = 410.0,
    linear_cap: float = 800.0,
) -> tuple[str, float, float]:
    """Choose a readable LCOE axis, switching to log for extreme values."""
    finite_values = [
        float(value)
        for value in upper_band_values
        if math.isfinite(float(value)) and float(value) > 0
    ]
    if not finite_values:
        return "linear", 0.0, linear_floor

    maximum = max(finite_values)
    sorted_values = sorted(finite_values)
    percentile_index = round(0.9 * (len(sorted_values) - 1))
    percentile_90 = sorted_values[percentile_index]

    if maximum <= linear_cap:
        ymax = min(linear_cap, max(linear_floor, percentile_90 * 1.15))
        return "linear", 0.0, ymax

    padded_maximum = maximum * 1.15
    if not math.isfinite(padded_maximum):
        padded_maximum = maximum
    # The lowest market benchmark is $29/MWh. A $10 lower bound keeps
    # every benchmark visible while providing a clean logarithmic decade.
    return "log", 10.0, max(linear_floor, padded_maximum)
