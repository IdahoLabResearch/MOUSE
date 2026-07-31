"""Input constraints derived from the reactor fuel-lifetime models."""

from __future__ import annotations

from collections.abc import Callable


MIN_USEFUL_LIFETIME_DAYS = 90
MAX_FUEL_LIFETIME_DAYS = 30 * 365


def safe_height_interval(
    estimate_lifetime: Callable[[int], int],
    minimum_height: int,
    maximum_height: int,
    *,
    minimum_lifetime_days: int = MIN_USEFUL_LIFETIME_DAYS,
    maximum_lifetime_days: int = MAX_FUEL_LIFETIME_DAYS,
) -> tuple[int, int] | None:
    """Return the longest contiguous integer-height interval with a safe lifetime.

    A contiguous interval is required because the height slider exposes every
    integer between its endpoints. The upper lifetime limit is exclusive, so a
    30-year estimate itself is not selectable.
    """
    runs: list[tuple[int, int]] = []
    run_start = None

    for height in range(int(minimum_height), int(maximum_height) + 1):
        lifetime = estimate_lifetime(height)
        is_safe = minimum_lifetime_days <= lifetime < maximum_lifetime_days
        if is_safe and run_start is None:
            run_start = height
        elif not is_safe and run_start is not None:
            runs.append((run_start, height - 1))
            run_start = None

    if run_start is not None:
        runs.append((run_start, int(maximum_height)))
    if not runs:
        return None

    return max(runs, key=lambda interval: interval[1] - interval[0])
