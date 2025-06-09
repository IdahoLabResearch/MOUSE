
import numpy as np 

def reactor_operation(params):
    
    # Refueling
    # how many times you add the fuel over the entire reactor lifetime
    add_fuel_num = int(np.floor(365*params['Levelization Period']/ 
                                (params['Refueling Period'] + params['fuel_lifetime_days'])))

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
    params['Annual Electricity Production'] = Capacity_factor * params['Power MWe'] * 354 * 24 # MWe.hour
    if params['Operation Mode'] == "Autonomous":
        params['FTEs Per Onsite Operator Per Year'] =   FTEs_per_operator_per_year_for_startup_after_refueling + FTEs_per_operator_per_year_for_startup_after_emergency_shutdown
    elif params['Operation Mode'] == "Non-Autonomous":
        params['FTEs Per Onsite Operator Per Year'] =  params['FTEs Per Onsite Operator (24/7)']

