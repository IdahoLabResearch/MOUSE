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
        if params['Enrichment'] < 0.1:
            cost_premium = 1
        elif  0.1 <= params['Enrichment'] < 0.2:
            cost_premium = 1.15
        elif  0.2 <params['Enrichment']:
            print("\033[91m ERROR: Enrichment is too high \033[0m")
        cost = cost_premium * unit_cost *pow(scaling_variable_value,exponent) 
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



def scale_GCMR_accounts(df, params):
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
