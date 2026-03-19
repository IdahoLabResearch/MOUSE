# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import numpy as np
import pandas as pd
from cost.sampling import sampler

def non_standard_cost_scale(account, unit_cost, scaling_variable_value, exponent, params):
    # pumps
    if account == 222.11 or account == 222.12:
        cost_multiplier = (0.2 / (1 - params['Pump Isentropic Efficiency'])) + 1
        cost = cost_multiplier * unit_cost * pow(scaling_variable_value,exponent)
    
    # compressors
    elif account == 222.13:
        if 'Primary Loop Count' in params.keys():
            # Account for multiple primary loops and their individual rated load
            ## PR1: Updated cost correlation based on ANL/NSE-20/28 in-place of default
            ## due to inherent uncertainty of compressor pressure ratio between GCMR designs
            cost_multiplier = (((params['Primary Loop Outlet Temperature'] - 273.15)/650)**1.29 *
                                (params['Primary Loop Compressor Power']/1e6/2.6)**0.74)
            cost = cost_multiplier * unit_cost
        else:
            # Old Correlation kept as backup
            cost_multiplier = (1 / (0.95 - params['Compressor Isentropic Efficiency'])) * params['Compressor Pressure Ratio'] * np.log(params['Compressor Pressure Ratio'])
            cost = cost_multiplier * unit_cost * pow(scaling_variable_value,exponent)
    
    elif account == 253:
        params['Enrichment Category'] = 'HALEU' if params['Enrichment'] >= 0.1 else 'LEU'
        if params['Enrichment'] >= 0.2:
            print("\033[91m ERROR: Enrichment is too high \033[0m")
            raise ValueError("Enrichment is too high")
        
        if params['Enrichment Category']== 'LEU':
            # LEU: standard commercial SWU pricing
            cost = unit_cost * pow(scaling_variable_value,exponent)   # unit_cost = $184.2/SWU from DB
        else:
            # HALEU: premium until sufficient market volume
            noak_units = params['NOAK Unit Number']
            haleu_threshold = params.get('HALEU NOAK Threshold', 100)
            effective_unit_cost = 184.2 if noak_units >= haleu_threshold else unit_cost  # unit_cost = $1000/SWU from DB
            cost = effective_unit_cost * pow(scaling_variable_value,exponent) 

    elif account == 711:
        cost_multiplier = params['FTEs Per Onsite Operator Per Year'] 
        cost = cost_multiplier * unit_cost * pow(scaling_variable_value,exponent)
    elif account == 712:
        cost_multiplier = params['FTEs Per Offsite Operator (24/7)']
        cost = cost_multiplier * unit_cost * pow(1 / scaling_variable_value, exponent) 
    elif account == 713:
        cost_multiplier = params['FTEs Per Security Staff (24/7)']
        cost = cost_multiplier * unit_cost * pow(scaling_variable_value,exponent)       
    elif account == 721:
        cost_multiplier = params['Annual Coolant Supply Frequency']
        cost = cost_multiplier * unit_cost * scaling_variable_value
 
    elif account == 81:
        cost_multiplier =  params['FTEs Per Operator Per Year Per Refueling'] 
        cost = cost_multiplier * unit_cost * pow(scaling_variable_value, exponent)
    return cost



def scale_redundant_BOP_and_primary_loop(df, params):
    # Scales special cases to handle redundant / multiple coolant, BoP loops
    escalation_year = params['Escalation Year']
    cost_col = f'FOAK Estimated Cost (${escalation_year })'

    if 'Primary Loop Count' in params.keys():
        df.loc[df['Account'].astype(str).str.startswith('222'), cost_col] *= params['Primary Loop Count']
    if 'BoP Count' in params.keys():
        # Balance of Plant
        df.loc[df['Account'].astype(str).str.startswith('232'), cost_col] *= params['BoP Count']
        # Balance of Plant Building - Assumed to be High 40' CONEX Container with 20 cm wall thickness (including conex wall)
        df.loc[df['Account'].astype(str).str.startswith('213.1'), cost_col] *= params['BoP Count']
    if 'Primary Loop Purification' in params.keys():
        df.loc[df['Account'] == 226, cost_col] *= int(params['Primary Loop Purification'])

    return df



def scale_cost(initial_database, params):
    scaled_cost = initial_database[['Account', 'Level', 'Account Title', 'FOAK to NOAK Multiplier Type',\
                                    "Fixed Cost Low End", "Fixed Cost High End", "Fixed Cost Distribution",\
                                    "Unit Cost Low End", "Unit Cost High End", "Unit Cost Distribution",\
                                    "Exponent std",  "Exponent Max", "Exponent Min", "Exponent Distribution"]]
    
    escalation_year = params['Escalation Year']
    

    # Iterate through each row in the DataFrame
    for index, row in initial_database.iterrows():
        
        # Check if cost data are available (fixed or unit cost)
        if row['Fixed Cost ($)'] > 0 or	row['Unit Cost'] > 0:
            
            scaling_variable_value = params[row['Scaling Variable']] if pd.notna(row['Scaling Variable']) else 0
            
            # Calculate the 'Estimated Cost
            fixed_cost_0 = row['Adjusted Fixed Cost ($)'] 
            fixed_cost_lo = row['Adjusted Fixed Cost Low End ($)'] 
            fixed_cost_hi = row['Adjusted Fixed Cost High End ($)'] 
            fixed_cost_dist = row['Fixed Cost Distribution']

            if pd.notna(row['Fixed Cost ($)']):
                if params['Number of Samples'] > 1:
                    if fixed_cost_dist == 'Lognormal':
                        fixed_cost = sampler("Lognormal", low_cost=fixed_cost_lo, high_cost=fixed_cost_hi, class3_cost=fixed_cost_0)
                    elif fixed_cost_dist == 'Uniform': 
                        fixed_cost = sampler('Uniform', low=fixed_cost_lo, high=fixed_cost_hi)
                    else:
                        fixed_cost = fixed_cost_0
                else:
                    fixed_cost = fixed_cost_0
            else:
                fixed_cost = 0    
            
            unit_cost_0 = row['Adjusted Unit Cost ($)'] 
            unit_cost_lo = row['Adjusted Unit Cost Low End ($)'] 
            unit_cost_hi = row['Adjusted Unit Cost High End ($)'] 
            unit_cost_dist = row['Unit Cost Distribution']

            if pd.notna(row['Unit Cost']):
                if params['Number of Samples'] > 1:
                    if unit_cost_dist == 'Lognormal':
                        unit_cost = sampler("Lognormal", low_cost=unit_cost_lo, high_cost=unit_cost_hi, class3_cost=unit_cost_0)
                    elif unit_cost_dist == 'Uniform': 
                        unit_cost = sampler('Uniform', low=unit_cost_lo, high=unit_cost_hi)
                    else:
                        unit_cost = unit_cost_0
                else:
                    unit_cost =unit_cost_0

            else:
                unit_cost = 0  

            scaling_variable_ref_value  = row['Scaling Variable Ref Value']
            exponent_0 = row['Exponent']
            exponent_min = row['Exponent Min']
            exponent_max = row['Exponent Max']
            exponent_std = row['Exponent std']
            exponent_dist = row['Exponent Distribution']

            if pd.notna(row['Exponent']):
                if params['Number of Samples'] > 1:
                    if exponent_dist == 'Truncated Normal':
                        exponent = sampler("Truncated Normal", mean=exponent_0, std=exponent_std, lower_bound=exponent_min, upper_bound=exponent_max)
                    else:
                        exponent = exponent_0
                else:
                    exponent = exponent_0
            
            if row['Standard Cost Equation?'] == 'standard' :
                
                if pd.notna(row['Scaling Variable']) and scaling_variable_value == 0:
                    estimated_cost = 0
                
                else:     
                    # Check if there is a ref value for the scaling variable or just use the unit cost
                    if row['Scaling Variable Ref Value'] > 0:
                        estimated_cost = fixed_cost +\
                        unit_cost * pow(scaling_variable_value,exponent) /(pow(scaling_variable_ref_value,exponent-1))

                    else:
                        # Calculate the 'Estimated Cost
                        estimated_cost = fixed_cost + unit_cost * scaling_variable_value
            
            elif row['Standard Cost Equation?'] == 'nonstandard':
                if pd.notna(row['Scaling Variable']) and scaling_variable_value == 0:
                    estimated_cost = 0
                else:    
                    estimated_cost = non_standard_cost_scale(row['Account'],\
                    unit_cost, scaling_variable_value, exponent, params)


            # Assign the calculated value to the corresponding row in the new DataFrame
            scaled_cost.at[index, f'FOAK Estimated Cost (${escalation_year })'] = estimated_cost
    return scaled_cost


def scale_central_facility_cost(initial_database, params):
    """
    Scale costs for central facility accounts.
    Similar to scale_cost() but includes Count Scaling Variable support.
    """
    scaled_cost = initial_database[['Account', 'Level', 'Account Title', 'FOAK to NOAK Multiplier Type',
                                    "Fixed Cost Low End", "Fixed Cost High End", "Fixed Cost Distribution",
                                    "Unit Cost Low End", "Unit Cost High End", "Unit Cost Distribution",
                                    "Exponent std", "Exponent Max", "Exponent Min", "Exponent Distribution"]]

    escalation_year = params['Escalation Year']
    params['Constant'] = 1

    for index, row in initial_database.iterrows():

        if row['Fixed Cost ($)'] > 0 or row['Unit Cost'] > 0:

            scaling_variable_value = params[row['Scaling Variable']] if pd.notna(row['Scaling Variable']) else 0
            count_variable_value = (params[row['Count Scaling Variable']] * row['Count per Variable']
                                    if pd.notna(row['Count Scaling Variable']) else 0)

            fixed_cost_0 = row['Adjusted Fixed Cost ($)']
            fixed_cost_lo = row['Adjusted Fixed Cost Low End ($)']
            fixed_cost_hi = row['Adjusted Fixed Cost High End ($)']
            fixed_cost_dist = row['Fixed Cost Distribution']

            if pd.notna(row['Fixed Cost ($)']):
                if params['Number of Samples'] > 1:
                    if fixed_cost_dist == 'Lognormal':
                        fixed_cost = sampler("Lognormal", low_cost=fixed_cost_lo, high_cost=fixed_cost_hi, class3_cost=fixed_cost_0)
                    elif fixed_cost_dist == 'Uniform':
                        fixed_cost = sampler('Uniform', low=fixed_cost_lo, high=fixed_cost_hi)
                    else:
                        fixed_cost = fixed_cost_0
                else:
                    fixed_cost = fixed_cost_0
            else:
                fixed_cost = 0

            unit_cost_0 = row['Adjusted Unit Cost ($)']
            unit_cost_lo = row['Adjusted Unit Cost Low End ($)']
            unit_cost_hi = row['Adjusted Unit Cost High End ($)']
            unit_cost_dist = row['Unit Cost Distribution']

            if pd.notna(row['Unit Cost']):
                if params['Number of Samples'] > 1:
                    if unit_cost_dist == 'Lognormal':
                        unit_cost = sampler("Lognormal", low_cost=unit_cost_lo, high_cost=unit_cost_hi, class3_cost=unit_cost_0)
                    elif unit_cost_dist == 'Uniform':
                        unit_cost = sampler('Uniform', low=unit_cost_lo, high=unit_cost_hi)
                    else:
                        unit_cost = unit_cost_0
                else:
                    unit_cost = unit_cost_0
            else:
                unit_cost = 0

            scaling_variable_ref_value = row['Scaling Variable Ref Value']
            exponent_0 = row['Exponent']
            exponent_min = row['Exponent Min']
            exponent_max = row['Exponent Max']
            exponent_std = row['Exponent std']
            exponent_dist = row['Exponent Distribution']

            if pd.notna(row['Exponent']):
                if params['Number of Samples'] > 1:
                    if exponent_dist == 'Truncated Normal':
                        exponent = sampler("Truncated Normal", mean=exponent_0, std=exponent_std, lower_bound=exponent_min, upper_bound=exponent_max)
                    else:
                        exponent = exponent_0
                else:
                    exponent = exponent_0

            if row['Standard Cost Equation?'] == 'standard':

                if pd.notna(row['Scaling Variable']) and scaling_variable_value == 0:
                    estimated_cost = 0
                else:
                    if row['Scaling Variable Ref Value'] > 0:
                        estimated_cost = (fixed_cost
                                          + unit_cost * pow(scaling_variable_value, exponent)
                                          / (pow(scaling_variable_ref_value, exponent - 1)))
                    else:
                        estimated_cost = fixed_cost + unit_cost * scaling_variable_value

                # Apply count scaling if specified
                if pd.notna(row['Count Scaling Variable']):
                    if count_variable_value == 0:
                        estimated_cost = 0
                    else:
                        estimated_cost = estimated_cost * count_variable_value

            elif row['Standard Cost Equation?'] == 'nonstandard':
                if pd.notna(row['Scaling Variable']) and scaling_variable_value == 0:
                    estimated_cost = 0
                else:
                    estimated_cost = non_standard_cost_scale(row['Account'],
                                                             unit_cost, scaling_variable_value, exponent, params)

            scaled_cost.at[index, f'FOAK Estimated Cost (${escalation_year })'] = estimated_cost
        else:
                # Explicitly assign 0 so these rows don't remain NaN
                # and cause parent account aggregation to fail
                scaled_cost.at[index, f'FOAK Estimated Cost (${escalation_year })'] = 0
    return scaled_cost
