# Component replacement periods

Inner vessel structure, moderator, reflector, reactor control devices, and
moderator booster replacement costs accept either calendar years or fuel-cycle
counts. Each component can use a different unit. The outer vessel retains its
existing calendar-year input.

For example, to replace the moderator on a ten-year cost basis:

```python
params['A75: Moderator Replacement Period (years)'] = 10
```

To use two complete fuel cycles instead:

```python
params['A75: Moderator Replacement Period (cycles)'] = 2
```

The same naming pattern applies to `Inner Vessel Structure`, `Reflector`,
`Reactor Control Devices`, and `Moderator Booster`.

An explicit `(years)` input takes precedence when both forms are supplied.
Otherwise, the `(cycles)` input is multiplied by the complete cycle duration:

```text
cycle years = (Fuel Lifetime + Refueling Period
               + Startup Duration after Refueling) / 365
replacement years = replacement cycles * cycle years
```

These input durations are in days. The existing 365-day conversion is retained.
For a four-year complete cycle, two cycles mean eight replacement years. An
explicit ten-year input remains ten years, regardless of the cycle duration.

Existing example and web-app defaults are unchanged. Cycle-based settings remain
supported, including non-one-cycle moderator settings. If neither input is
supplied for the moderator or moderator booster, the existing one-cycle default
is retained. The inner vessel, reflector, and control devices still require an
explicit period in one of the two units.

Calendar-year inputs must be finite and positive. Invalid year inputs are
reported rather than silently falling back to cycles. Cycle counts must be
finite and nonnegative, with a positive full-cycle duration. The existing
zero-cycle convention is preserved: zero means no scheduled replacement cost,
not immediate replacement. This includes zeros produced by existing
`floor(10 / cycle_years)` defaults for cycles longer than ten years. Supply a
positive calendar interval or cycle count to include replacements in such a
case. An explicit positive year input also overrides a legacy zero-cycle value.

The resolved interval feeds the existing capital recovery factor in the FOAK
and NOAK annualized replacement accounts: 752 (inner vessel), 753 (moderator),
754 (reflector), 755 (control devices), and 757 (moderator booster). This feature
does not schedule replacement outages, align calendar intervals to refueling,
or add downtime and restart costs. It does not change the component's original
capital cost, the maintenance calculation, or the fuel-replacement interval.
