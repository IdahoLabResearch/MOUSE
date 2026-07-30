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
