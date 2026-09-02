# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import copy


# Define placeholders for the variables
fuel_pin = 'FUEL'
moderator_pin = 'MODERATOR'
shutdown_pin = 'SHUTDOWN'

SUPPORTED_LTMR_SHUTDOWN_ROD_COUNTS = (6, 12)

# Your list structure with placeholders
_LTMR_PINS_ARRANGEMENT_6 = [

[moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin,\
 fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] * 6,

[fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin,\
 moderator_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] *6,

[moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin,\
 fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] *6,
 
[fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin,\
 moderator_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin] *6, 

[moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin,\
 fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] *6, 

[fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin,\
 moderator_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin] *6,
 
[moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin,\
 fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] *6, 
[fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin,\
 moderator_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin] *6,
[moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin,\
 fuel_pin, fuel_pin, fuel_pin, fuel_pin] *6, 

[fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin,\
 moderator_pin, fuel_pin, moderator_pin] *6,

[moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin,\
 fuel_pin, fuel_pin] *6,

[fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin,\
 moderator_pin] * 6,
    



[moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin] * 6,
[fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin] * 6,
[moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] * 6,

[fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin] * 6,
[shutdown_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] * 6,
[fuel_pin, fuel_pin, moderator_pin, fuel_pin, moderator_pin, fuel_pin] * 6,

[moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] * 6,
[fuel_pin, fuel_pin, moderator_pin, fuel_pin] * 6,
[moderator_pin, fuel_pin, fuel_pin] * 6,

[fuel_pin, moderator_pin] * 6,
[fuel_pin] * 6,
[moderator_pin]
]


# The existing six shutdown channels occupy one sixfold-symmetric orbit on
# the 42-position ring.  The predefined twelve-rod option adds a second
# sixfold-symmetric orbit on that same ring by replacing six fuel positions.
_LTMR_SHUTDOWN_RING_INDEX = -8
_LTMR_PINS_ARRANGEMENT_12 = copy.deepcopy(_LTMR_PINS_ARRANGEMENT_6)
_LTMR_PINS_ARRANGEMENT_12[_LTMR_SHUTDOWN_RING_INDEX] = [
    shutdown_pin,
    fuel_pin,
    fuel_pin,
    shutdown_pin,
    fuel_pin,
    fuel_pin,
    fuel_pin,
] * 6

_LTMR_PIN_ARRANGEMENTS = {
    6: _LTMR_PINS_ARRANGEMENT_6,
    12: _LTMR_PINS_ARRANGEMENT_12,
}


def get_ltmr_pins_arrangement(number_of_shutdown_rods=6):
    """Return a copy of the predefined LTMR map for 6 or 12 rods."""
    if number_of_shutdown_rods not in SUPPORTED_LTMR_SHUTDOWN_ROD_COUNTS:
        raise ValueError(
            "Number of Shutdown Rods must be one of "
            f"{list(SUPPORTED_LTMR_SHUTDOWN_ROD_COUNTS)}, got "
            f"{number_of_shutdown_rods!r}."
        )

    return copy.deepcopy(_LTMR_PIN_ARRANGEMENTS[number_of_shutdown_rods])


# Backward-compatible default for existing LTMR examples that import the
# original module-level arrangement directly.
LTMR_pins_arrangement = get_ltmr_pins_arrangement(6)
