"""Shared parameter helpers for Fleet Mode campus sizing."""

import math


def servicing_facility_allocation(
    fleet_size,
    servicing_rate,
    max_reactors_per_facility=1000,
):
    """Return an integer, balanced allocation of reactors among facilities.

    The design capacity is the largest assigned facility load. Cost estimates
    use that capacity for every replicated facility so differently loaded
    facilities still share one standard design.
    """
    if fleet_size <= 0:
        raise ValueError("'Fleet' must be greater than zero.")
    if int(fleet_size) != fleet_size:
        raise ValueError("'Fleet' must be a whole number of reactors.")
    if servicing_rate < 0:
        raise ValueError("'Servicing Rate' cannot be negative.")
    if max_reactors_per_facility <= 0:
        raise ValueError("'Max Reactors Per Servicing Facility' must be greater than zero.")
    if int(max_reactors_per_facility) != max_reactors_per_facility:
        raise ValueError("'Max Reactors Per Servicing Facility' must be a whole number.")

    fleet_size = int(fleet_size)
    max_reactors_per_facility = int(max_reactors_per_facility)
    facility_count = math.ceil(fleet_size / max_reactors_per_facility)
    minimum_load, facilities_with_extra_reactor = divmod(
        fleet_size, facility_count
    )
    reactor_counts = (
        (minimum_load + 1,) * facilities_with_extra_reactor
        + (minimum_load,) * (facility_count - facilities_with_extra_reactor)
    )
    design_capacity = max(reactor_counts)
    return (
        facility_count,
        reactor_counts,
        design_capacity,
        servicing_rate / facility_count,
    )


def servicing_facility_occ_learning_multipliers(
    facility_count,
    learning_rate=0.30,
    learning_cap=5,
):
    """Return OCC multipliers for sequentially built servicing facilities."""
    if facility_count < 1 or int(facility_count) != facility_count:
        raise ValueError("'Servicing Facility Count' must be a positive whole number.")
    if not 0 <= learning_rate < 1:
        raise ValueError("'Servicing Facility OCC Learning Rate' must be in [0, 1).")
    if learning_cap < 1 or int(learning_cap) != learning_cap:
        raise ValueError("'Servicing Facility Learning Cap' must be a positive whole number.")

    exponent = math.log2(1 - learning_rate)
    return tuple(
        min(facility_number, int(learning_cap)) ** exponent
        for facility_number in range(1, int(facility_count) + 1)
    )


def cumulative_unit_learning_multiplier(
    unit_count,
    learning_rate=0.15,
    learning_cap=100,
):
    """Return the sum of learned unit-cost multipliers for ``unit_count``.

    The conventional unit learning curve is applied to each sequential unit
    and the results are added. Fractional annual demand receives the
    multiplier of the next unit. Units beyond ``learning_cap`` retain the
    capped unit's multiplier.
    """
    if unit_count < 0:
        raise ValueError("'unit_count' cannot be negative.")
    if not 0 <= learning_rate < 1:
        raise ValueError("'Cask Learning Rate' must be in [0, 1).")
    if learning_cap < 1 or int(learning_cap) != learning_cap:
        raise ValueError("'Cask Learning Cap' must be a positive whole number.")
    if unit_count == 0:
        return 0.0

    learning_cap = int(learning_cap)
    exponent = math.log2(1 - learning_rate)
    whole_units = int(math.floor(unit_count))
    fractional_unit = unit_count - whole_units
    aggregate_multiplier = sum(
        min(unit_number, learning_cap) ** exponent
        for unit_number in range(1, whole_units + 1)
    )
    if fractional_unit:
        aggregate_multiplier += (
            fractional_unit
            * min(whole_units + 1, learning_cap) ** exponent
        )
    return aggregate_multiplier
