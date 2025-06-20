import numpy as np
from cost.code_of_account_processing import get_estimated_cost_column

def calculate_accounts_31_32_75_82_cost( df, params):
    
    # Find the column name that starts with 'Estimated Cosoption == "other costs"t'
    estimated_cost_col_F = get_estimated_cost_column(df, 'F')
    estimated_cost_col_N = get_estimated_cost_column(df, 'N')

    for estimated_cost_col in [estimated_cost_col_F, estimated_cost_col_N ]:
        # Filter the DataFrame for accounts 21, 22, and 23
        filtered_df = df[df['Account'].isin([21, 22, 23])]

        # Sum the values in the 'Estimated Cost' column for the filtered accounts
        tot_field_direct_cost = filtered_df[estimated_cost_col].sum()

        acct_31_cost = params['indirect to direct field-related cost'] * tot_field_direct_cost # This ratio is based on MARVEL
        df.loc[df['Account'] == 31, estimated_cost_col] = acct_31_cost

        # To calculate the cost of factory and construction supervision (Account 32), 
        # the ratio of the factory and field indirect costs (Account 31) to the reactor systems cost (account 22) 
        # is calculated and multiplied by the cost of structures and improvements (Account 21)
        df.loc[df['Account'] == 32, estimated_cost_col] = df.loc[df['Account'] == 21, estimated_cost_col].values[0] * (df.loc[df['Account'] == 31, estimated_cost_col].values[0] / df.loc[df['Account'] == 22, estimated_cost_col].values[0])
        
        # decommisioning cost
        df.loc[df['Account'] == 75, estimated_cost_col] = df.loc[df['Account'] == 20, estimated_cost_col].values[0] * params['Maintenance to Direct Cost Ratio']

        cycle_length = params['Fuel Lifetime'] / 365.25
        df.loc[df['Account'] == 82, estimated_cost_col] = df.loc[df['Account'] == 25, estimated_cost_col].values[0] /  cycle_length

    return df


def calculate_interest_cost(params, OCC):
    interest_rate = params['Interest Rate']
    construction_duration = params['Construction Duration']
    debt_to_equity_ratio = params['Debt To Equity Ratio'] 
    # Interest rate from this equation (from Levi)
    B =(1+ np.exp((np.log(1+ interest_rate)) * construction_duration/12))
    C  =((np.log(1+ interest_rate)*(construction_duration/12)/3.14)**2+1)
    Interest_expenses = debt_to_equity_ratio*OCC*((0.5*B/C)-1)
    return Interest_expenses


def calculate_high_level_capital_costs(df, params):
    power_kWe = 1000 * params['Power MWe']
     # List of accounts to sum
    accounts_to_sum = [10, 20, 30, 40, 50]

    # Create the OCC account "OCC" with the new total cost
    df = df._append({'Account': 'OCC','Account Title' : 'Overnight Capital Cost'}, ignore_index=True)
    df = df._append({'Account': 'OCC per kW','Account Title' : 'Overnight Capital Cost per kW' }, ignore_index=True)
    df = df._append({'Account': 'OCC excl. fuel','Account Title' : 'Overnight Capital Cost Excluding Fuel'}, ignore_index=True)
    df = df._append({'Account': 'OCC excl. fuel per kW','Account Title' : 'Overnight Capital Cost Excluding Fuel per kW'}, ignore_index=True)

    # Find the column that starts with "Estimated Cost"
    cost_column_F = get_estimated_cost_column(df, 'F')
    cost_column_N = get_estimated_cost_column(df, 'N')

    for cost_column in [cost_column_F, cost_column_N]:
        # Calculate the sum of costs for the specified accounts
        occ_cost = df[df['Account'].isin(accounts_to_sum)][cost_column].sum()
        df.loc[df['Account'] == 'OCC', cost_column] = occ_cost
        df.loc[df['Account'] == 'OCC per kW', cost_column] = occ_cost/ power_kWe
        
        #OCC excluding the fuel
        occ_excl_fuel = occ_cost - (df.loc[df['Account'] == 25, cost_column].values[0])
        df.loc[df['Account'] == 'OCC excl. fuel', cost_column] = occ_excl_fuel
        df.loc[df['Account'] == 'OCC excl. fuel per kW', cost_column] = occ_excl_fuel/ power_kWe

        df.loc[df['Account'] == 62, cost_column] =  calculate_interest_cost(params, occ_cost)
    return df


def calculate_TCI(df, params):
    power_kWe = 1000 * params['Power MWe']

    df = df._append({'Account': 'TCI','Account Title' : 'Total Capital Investment'}, ignore_index=True)
    df = df._append({'Account': 'TCI per kW','Account Title' : 'Total Capital Investment per kW'}, ignore_index=True)

    # List of accounts to sum
    accounts_to_sum = ['OCC', 60]
    # Find the column that starts with "Estimated Cost"
    cost_column_F = get_estimated_cost_column(df, 'F')
    cost_column_N = get_estimated_cost_column(df, 'N')
    for cost_column in [cost_column_F , cost_column_N]:
        # Calculate the sum of costs for the specified accounts
        tci_cost = df[df['Account'].isin(accounts_to_sum)][cost_column].sum()
        # Update the existing account "OCC" with the new total cost
        df.loc[df['Account'] == 'TCI', cost_column] = tci_cost
        df.loc[df['Account'] == 'TCI per kW', cost_column] = tci_cost/power_kWe

    return df


def energy_cost_levelized(params, df):

    df = df._append({'Account': 'AC','Account Title' : 'Annualized Cost'}, ignore_index=True)
    df = df._append({'Account': 'AC per MWh','Account Title' : 'Annualized Cost per MWh'}, ignore_index=True)
    df = df._append({'Account': 'LCOE','Account Title' : 'Levelized Cost Of Energy ($/MWh)'}, ignore_index=True)
  
    
    plant_lifetime_years = params['Levelization Period']
    discount_rate = params['Interest Rate']
    power_MWe = params['Power MWe']
    capacity_factor = params['Capacity Factor']
    estimated_cost_col_F = get_estimated_cost_column(df, 'F')
    estimated_cost_col_N = get_estimated_cost_column(df, 'N')

    for estimated_cost_col in [estimated_cost_col_F, estimated_cost_col_N]:

        cap_cost = df.loc[df['Account'] == 'TCI', estimated_cost_col].values[0]
        ann_cost = df.loc[df['Account'] == 70, estimated_cost_col].values[0]  + df.loc[df['Account'] == 80, estimated_cost_col].values[0] 
        levelized_ann_cost = ann_cost / params['Annual Electricity Production'] 
        df.loc[df['Account'] == 'AC', estimated_cost_col] = ann_cost
        df.loc[df['Account'] == 'AC per MWh', estimated_cost_col] = levelized_ann_cost
        sum_cost = 0 # initialization 
        sum_elec = 0
        
        for i in range(  1 + plant_lifetime_years) :
            
            if i == 0:
            # assuming that the cap cost is split between the cons years
                cap_cost_per_year = cap_cost
                annual_cost = 0
                elec_gen = 0
                
            elif i >0:
                cap_cost_per_year  = 0
                annual_cost = ann_cost
                elec_gen = power_MWe *capacity_factor * 365 * 24       # MW hour. 
            sum_cost +=  (cap_cost_per_year + annual_cost)/ ((1+ discount_rate)**i) 
            sum_elec += elec_gen/ ((1 + discount_rate)**i) 
        lcoe =  sum_cost/ sum_elec
        
        df.loc[df['Account'] == 'LCOE', estimated_cost_col] = lcoe    
    return df