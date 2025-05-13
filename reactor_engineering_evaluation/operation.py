
import numpy as np 

def reactor_operation(params):
    
    # Refueling
    # how many times you add the fuel over the entire reactor lifetime
    add_fuel_num = int(np.floor(365*params['levelization_period_years']/ 
                                (params['refueling_period_days'] + params['fuel_lifetime_days'])))

    num_of_refuel_days_per_year = params['refueling_period_days'] *\
        add_fuel_num/params['levelization_period_years']
    
    #how many people multiplied by how many days of refueling per year
    people_by_days_refueling = params['num_people_required_per_refueling'] * num_of_refuel_days_per_year
    
    
    #how many days to startup after refueling (per year)
    num_startup_days_after_refuel_per_year =  add_fuel_num *\
        params['duration_to_startup_after_refueling_days']/params['levelization_period_years']

    people_by_days_startup_after_refueling = params['num_people_required_per_startup'] * num_startup_days_after_refuel_per_year

     #how many days to startup after emergency shutdown (per year)
    num_startup_days_after_shutdown_per_year = params['duration_to_startup_after_shutdown_days'] *\
        params['number_of_unanticipated_shutdowns_per_year']
    people_by_days_startup_after_shutdown =   params['num_people_required_per_startup'] * num_startup_days_after_shutdown_per_year  
     
    total_people_by_days_for_startup  = people_by_days_startup_after_refueling + people_by_days_startup_after_shutdown 
    
    Capacity_factor  = 1 - ((num_of_refuel_days_per_year +num_startup_days_after_refuel_per_year + num_startup_days_after_shutdown_per_year )/365)

    return people_by_days_refueling, total_people_by_days_for_startup ,Capacity_factor    
