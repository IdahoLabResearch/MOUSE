# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import numpy as np 


REMOTE_MONITORED_OPERATION_MODES = frozenset({
    "Remotely Monitored",
    "On-Site Staffed and Remotely Monitored",
})

ONSITE_STAFFED_OPERATION_MODES = frozenset({
    "On-Site Staffed",
    "On-Site Staffed and Remotely Monitored",
})

VALID_OPERATION_MODES = (
    REMOTE_MONITORED_OPERATION_MODES | ONSITE_STAFFED_OPERATION_MODES
)


def operation_mode_includes_remote_monitoring(operation_mode):
    """Return whether *operation_mode* requires 24/7 remote monitoring."""
    return operation_mode in REMOTE_MONITORED_OPERATION_MODES

def reactor_operation(params):
    # Backward compatibility for external inputs created before the parameter
    # was renamed. Maintained MOUSE inputs use the clearer canonical name.
    if ('Number of On-Site Operators per Shift' not in params
            and 'Number of Operators' in params):
        params['Number of On-Site Operators per Shift'] = params['Number of Operators']
    
    # Refueling
    # how many times you add the fuel over the entire reactor lifetime
    add_fuel_num = int(np.floor(365*params['Levelization Period']/ 
                                (params['Refueling Period'] + params['Fuel Lifetime'])))

    num_of_refuel_days_per_year = params['Refueling Period'] *\
        add_fuel_num/params['Levelization Period']
    
    #how many FTES per operator per year (for refueling)
    FTEs_per_operator_per_year_for_refueling =  num_of_refuel_days_per_year * params['Work Hours Per Shift']/ params['Hours Per FTE']
    params['FTEs Per Operator Per Year Per Refueling'] = FTEs_per_operator_per_year_for_refueling
    
    #how many days to startup after refueling (per year)
    num_startup_days_after_refuel_per_year =  add_fuel_num *\
        params['Startup Duration after Refueling']/params['Levelization Period']

    #how many FTES per operator per year (for startup after refueling)
    FTEs_per_operator_per_year_for_startup_after_refueling =  num_startup_days_after_refuel_per_year * params['Work Hours Per Shift']/ params['Hours Per FTE']

    #how many days to startup after emergency shutdown (per year)
    num_startup_days_after_shutdown_per_year = params['Startup Duration after Emergency Shutdown'] *\
        params['Emergency Shutdowns Per Year']

    #how many FTES per operator per year (for startup after emergency shutdown)
    FTEs_per_operator_per_year_for_startup_after_emergency_shutdown =    num_startup_days_after_shutdown_per_year  * params['Work Hours Per Shift']/ params['Hours Per FTE']
         
    Capacity_factor  = 1 - ((num_of_refuel_days_per_year +num_startup_days_after_refuel_per_year + num_startup_days_after_shutdown_per_year )/365)
    params['Capacity Factor'] = Capacity_factor 
    params['Annual Electricity Production'] = Capacity_factor * params['Power MWe'] * 365 * 24 # MWe.hour
    operation_mode = params['Operation Mode']
    if operation_mode == "Remotely Monitored":
        params['FTEs Per Onsite Operator Per Year'] =   FTEs_per_operator_per_year_for_startup_after_refueling + FTEs_per_operator_per_year_for_startup_after_emergency_shutdown
    elif operation_mode in ONSITE_STAFFED_OPERATION_MODES:
        params['FTEs Per Onsite Operator Per Year'] =  params['FTEs Per Onsite Operator (24/7)']
    else:
        valid_modes = ", ".join(sorted(VALID_OPERATION_MODES))
        raise ValueError(
            f"Unsupported Operation Mode {operation_mode!r}. "
            f"Choose one of: {valid_modes}."
        )
