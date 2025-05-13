from cost.baseline_costs import cost_dictionary_data
from cost.cost_utils import update_high_level_costs, energy_cost_levelized

#Cost Scaling methods
def cost_scale_general(db, account_number, old_var, new_var, scaling_exponent, old_cost, constant_cost):
    if old_var == None and new_var == None and old_cost == None:
        new_cost = constant_cost
    else:    
        new_cost == old_cost * ((new_var/old_var)**scaling_exponent) + constant_cost
       #add new cost to the table
    db.loc[db.Account == account_number, 'Estimated Cost (2023 $)'] =  new_cost      
    
    
def scale_unit_costs(db, account_number, unit_cost, new_var):
    new_cost =   unit_cost * new_var
    db.loc[db.Account == account_number, 'Estimated Cost (2023 $)'] =  new_cost  


def add_customized_cost_estimate(db, account_number, customized_cost_estimate):
    db.loc[db.Account == account_number, 'Estimated Cost (2023 $)'] =  customized_cost_estimate 
    # print ( f"The {db.loc[db.Account == account_number, 'Account Title'].values[0] } account has been updated assuming {scaling_method} scaling \n")
    
   

def cost_estimate(code_of_account , params) :
    power_kWe = 1000 * params['power_MW_th'] * params['thermal_efficiency'] 

      #Create new database
    cost_new_reactor = code_of_account[['Account', 'Account Title', 'Estimated Cost (2023 $)']].copy()  
    
    # scaling based on unit costs: accounts 11, 
    scale_unit_costs(cost_new_reactor, 11, cost_dictionary_data['Land Unit Cost'], params['land_area_acres'])

    # scaling based on the general scaling method: accounts 13, 15
    cost_scale_general(cost_new_reactor, 13, None, None, None, None,  cost_dictionary_data['Licensing Cost'])
    cost_scale_general(cost_new_reactor, 15, None, None, None, None,  cost_dictionary_data['Plant Studies'])
    cost_scale_general(cost_new_reactor, '212A', None, None, None, None,  cost_dictionary_data['Pit Preparation'])
    scale_unit_costs(cost_new_reactor, '212B', cost_dictionary_data['Building Construction Cost'], params['building_area_feet_squared'])
    scale_unit_costs(cost_new_reactor, '212C', cost_dictionary_data['Building Construction Cost'], params['building_area_feet_squared'])

    db_updated = update_high_level_costs(cost_new_reactor, power_kWe) 
    
    (db_updated.loc[db_updated.Account == 'Annualized Cost ($/MWh)', 'Estimated Cost (2023 $)'])   =\
        (db_updated.loc[db_updated.Account == 'Annualized Cost ($/ year)', 'Estimated Cost (2023 $)'])\
            .values[0] / (power_kWe/ 1000 * params['capacity_factor'] * 365 * 24)
        
    cap_cost = (db_updated.loc[db_updated.Account == 'Total Capital Investment (TCI)', 'Estimated Cost (2023 $)']).values[0]
    ann_cost = (db_updated.loc[db_updated.Account == 'Annualized Cost ($/ year)', 'Estimated Cost (2023 $)']).values[0]
    ann_cost_levelized = (db_updated.loc[db_updated.Account == 'Annualized Cost ($/MWh)', 'Estimated Cost (2023 $)']).values[0]

    lcoe1 = energy_cost_levelized(params['levelization_period_years'], cap_cost, ann_cost,\
        params['interest_rate'], power_kWe/1000, params['capacity_factor']  )

    (db_updated.loc[db_updated.Account == 'LCOE', 'Estimated Cost (2023 $)'])   = lcoe1 
    # pretty_table = prettify(db_updated, f"{params['reactor type']} {power_kWe/ 1000} MWe" )
    return db_updated, cap_cost, ann_cost, ann_cost_levelized, lcoe1