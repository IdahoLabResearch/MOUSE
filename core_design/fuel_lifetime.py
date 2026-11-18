# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
"""Fuel-lifetime interpretation shared by depletion and cost workflows.

This module deliberately has no OpenMC dependency. A predicted crossing is
reported for review, but only a simulated, bracketed crossing is cost eligible.
"""

from dataclasses import dataclass
import math
from typing import Mapping, MutableMapping, Sequence


INTERPOLATED = 'Interpolated k=1 crossing'
EXTRAPOLATED = 'Extrapolated beyond depletion schedule'
LOWER_BOUND = 'Simulated endpoint lower bound'
SUBCRITICAL = 'BOL subcritical'


@dataclass(frozen=True)
class FuelLifetimeResult:
    """A lifetime estimate and the evidence supporting its interpretation."""

    lifetime_days: float
    status: str
    extrapolated: bool
    simulated_endpoint_days: float


class InvalidFuelLifetimeError(ValueError):
    """The reported lifetime must not be used to calculate costs."""


def estimate_fuel_lifetime(
    time_days: Sequence[float], keff_values: Sequence[float]
) -> FuelLifetimeResult:
    """Interpret a depletion history without confusing a prediction with EOL.

    Return the first downward k=1 crossing when one is simulated. Otherwise,
    extrapolate the final two points only when they have a declining slope and
    yield a finite future crossing. A flat/rising or single-point history gives
    a lower bound at the simulated endpoint. Initial k <= 1 gives zero operating
    lifetime, even if later points rise above one.
    """
    try:
        times = [float(value) for value in time_days]
        keff = [float(value) for value in keff_values]
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError('Depletion times and keff values must be numeric.') from error
    if not times or len(times) != len(keff):
        raise ValueError('Depletion times and keff values must have equal, nonzero lengths.')
    if not all(math.isfinite(value) and value >= 0 for value in times):
        raise ValueError('Depletion times must be finite and nonnegative.')
    if any(later <= earlier for earlier, later in zip(times, times[1:])):
        raise ValueError('Depletion times must be strictly increasing.')
    if not all(math.isfinite(value) and value > 0 for value in keff):
        raise ValueError('Depletion keff values must be finite and positive.')

    endpoint = times[-1]
    if keff[0] <= 1.0:
        return FuelLifetimeResult(0.0, SUBCRITICAL, False, endpoint)

    for index in range(1, len(keff)):
        k1, k2 = keff[index - 1], keff[index]
        if k1 >= 1.0 and k2 <= 1.0:
            fraction = (k1 - 1.0) / (k1 - k2)
            crossing = times[index - 1] + fraction * (times[index] - times[index - 1])
            return FuelLifetimeResult(crossing, INTERPOLATED, False, endpoint)

    if len(keff) >= 2:
        slope = (keff[-1] - keff[-2]) / (times[-1] - times[-2])
        if math.isfinite(slope) and slope < 0:
            predicted = endpoint + (keff[-1] - 1.0) / -slope
            if math.isfinite(predicted) and predicted > endpoint:
                return FuelLifetimeResult(predicted, EXTRAPOLATED, True, endpoint)

    return FuelLifetimeResult(endpoint, LOWER_BOUND, False, endpoint)


def record_fuel_lifetime(
    params: MutableMapping[str, object], result: FuelLifetimeResult, *, round_days: bool = True
) -> None:
    """Store/export the interpretation and print a clear user-facing message.

    The corrected-2D workflow keeps its existing nearest-day convention.
    Explicit-3D callers use ``round_days=False`` to retain their existing raw
    precision. Neither path rounds to fuel-cycle years or study bins.
    """
    params['Fuel Lifetime'] = (
        round(result.lifetime_days, 0) if round_days else result.lifetime_days
    )
    params['Fuel Lifetime Status'] = result.status
    params['Fuel Lifetime Extrapolated'] = result.extrapolated
    params['Fuel Lifetime Simulated Endpoint (days)'] = result.simulated_endpoint_days
    if result.status == SUBCRITICAL:
        print(
            'ERROR: Reactor is subcritical or exactly critical at the starting '
            'point (k_eff <= 1). Operating fuel lifetime: 0 days. '
            'Cost estimation will not run.'
        )
    elif result.status == INTERPOLATED:
        print(
            f"Estimated fuel cycle length: {params['Fuel Lifetime']:g} days "
            f'[{result.status}]. Simulated endpoint: '
            f'{result.simulated_endpoint_days:g} days.'
        )
    elif result.status == EXTRAPOLATED:
        print(
            f'WARNING: Provisional extrapolated fuel lifetime: '
            f'{result.lifetime_days:g} days; simulated endpoint: '
            f'{result.simulated_endpoint_days:g} days. Extend the depletion '
            'schedule until k_eff reaches 1 before lifecycle or cost estimation.'
        )
    else:
        print(
            f'WARNING: Fuel lifetime is only a lower bound at the simulated '
            f'endpoint of {result.simulated_endpoint_days:g} days. No finite '
            'future crossing can be inferred from the final trend. Extend the '
            'depletion schedule until k_eff reaches 1 before lifecycle or cost estimation.'
        )


def require_cost_eligible_lifetime(params: Mapping[str, object]) -> None:
    """Reject incomplete lifetime interpretations before costs or snapshots.

    Positive legacy/app lifetimes without interpretation metadata remain valid.
    Callers that historically omit lifetime entirely retain their existing
    parameter-validation behavior. Explicit metadata must describe a simulated
    crossing, and any supplied lifetime must be finite and positive.
    """
    status = params.get('Fuel Lifetime Status')
    if status == SUBCRITICAL:
        raise InvalidFuelLifetimeError(
            'Reactor is subcritical or exactly critical at the starting point '
            '(k_eff <= 1). Operating fuel lifetime is 0 days. Cost estimation '
            'will not run.'
        )
    if status in (EXTRAPOLATED, LOWER_BOUND) or params.get('Fuel Lifetime Extrapolated', False):
        raise InvalidFuelLifetimeError(
            f'Fuel lifetime is provisional ({status or EXTRAPOLATED}); '
            f"reported lifetime: {params.get('Fuel Lifetime', 'unavailable')} days; "
            'simulated endpoint: '
            f"{params.get('Fuel Lifetime Simulated Endpoint (days)', 'unavailable')} days. "
            'Extend the depletion schedule until a simulated k_eff=1 crossing '
            'is obtained. Lifecycle and cost estimation will not run.'
        )
    if 'Fuel Lifetime Status' in params and status != INTERPOLATED:
        raise InvalidFuelLifetimeError(
            f'Unrecognized fuel-lifetime status {status!r}. Cost estimation will not run.'
        )
    if 'Fuel Lifetime' not in params:
        if any(key in params for key in (
            'Fuel Lifetime Status', 'Fuel Lifetime Extrapolated',
            'Fuel Lifetime Simulated Endpoint (days)',
        )):
            raise InvalidFuelLifetimeError('Fuel lifetime is missing. Cost estimation will not run.')
        return
    try:
        lifetime = float(params['Fuel Lifetime'])
    except (TypeError, ValueError, OverflowError) as error:
        raise InvalidFuelLifetimeError(
            'Fuel lifetime must be finite and greater than zero. Cost estimation will not run.'
        ) from error
    if not math.isfinite(lifetime) or lifetime <= 0:
        raise InvalidFuelLifetimeError(
            f'Fuel lifetime must be finite and greater than zero; received {lifetime:g} days. '
            'Cost estimation will not run.'
        )
