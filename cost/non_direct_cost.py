import numpy as np
import pandas as pd
from cost.code_of_account_processing import get_estimated_cost_column

def _crf(rate, period):
    # Returns the Capital Recovery Factor  based on the discount rate and period
    numer = rate * (1 + rate)**period
    denum = (1 + rate)**period - 1
    factor = numer/denum 
    
    ## If components are not set for replacement (i.e. period == 0) return 0
    if np.array(factor).size > 1:
        factor[factor == np.inf] = 0
    return factor



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
        
        # A75: Annualized Capital Expenditures
        # Check if replacement period is specified in params
        refueling_period = params['Fuel Lifetime'] + params['Refueling Period'] + params['Startup Duration after Refueling']
        refueling_period_yr = refueling_period / 365
        params_df = pd.DataFrame(params.items(), columns=['keys', 'values'])
        if params_df.loc[params_df['keys'].str.contains('replacement', case=False), 'keys'].size > 0:
            # Input Case includes period replacement of internals (e.g. GCMS)
            # Replacements are assumed to match with refueling so #cycles are used instead of #years
            A20_replacement_period = refueling_period_yr * np.array([params['A75: RPV Replacement Period (cycles)'],
                                                                     params['A75: Core Barrel Replacement Period (cycles)'],
                                                                     1, # Moderator Block Replacement Period (cycles)
                                                                     params['A75: Reflector Replacement Period (cycles)'],
                                                                     params['A75: Drum Replacement Period (cycles)'],
                                                                     params.get('A75: Integrated HX Replacement Period (cycles)', 0),])
            ## Keep the same ordering as `A20_replacement_period`
            A20_capital_cost = np.array([df.loc[df['Account'] == 221.12, estimated_cost_col].values.sum(), 
                                         df.loc[df['Account'] == 221.13,  estimated_cost_col].values.sum(), 
                                         df.loc[df['Account'] == 221.33,  estimated_cost_col].values.sum(),
                                         df.loc[df['Account'] == 221.31,  estimated_cost_col].values.sum(),
                                         df.loc[df['Account'] == 221.2,   estimated_cost_col].values.sum(),
                                         df.loc[df['Account'].isin([222.1, 222.2, 222.3, 222.61]), estimated_cost_col].values.sum()])
            annualized_replacement_cost = (A20_capital_cost*_crf(params['Interest Rate'], A20_replacement_period))
            A20_other_cost = df.loc[df['Account'] == 20, estimated_cost_col].values[0] - A20_capital_cost.sum()
            annualized_other_cost = A20_other_cost * params['Maintenance to Direct Cost Ratio']
            # For non-specified CAPEX components, use the old method of saving 
            # `params['Mainenance to Direct Cost Ratio']` * CAPEX annually
            df.loc[df['Account'] == 751, estimated_cost_col] = annualized_replacement_cost[0]
            df.loc[df['Account'] == 752, estimated_cost_col] = annualized_replacement_cost[1]
            df.loc[df['Account'] == 753, estimated_cost_col] = annualized_replacement_cost[2]
            df.loc[df['Account'] == 754, estimated_cost_col] = annualized_replacement_cost[3]
            df.loc[df['Account'] == 755, estimated_cost_col] = annualized_replacement_cost[4]
            df.loc[df['Account'] == 756, estimated_cost_col] = annualized_replacement_cost[5]
            df.loc[df['Account'] == 759, estimated_cost_col] = annualized_other_cost
        else:
            # If no A75's specified in `params`, rely on
            # `params['Mainenance to Direct Cost Ratio']` * CAPEX annually
            df.loc[df['Account'] == 75, estimated_cost_col] = df.loc[df['Account'] == 20, estimated_cost_col].values[0] * params['Maintenance to Direct Cost Ratio']

        # A82: Annualized Fuel Cost
        lump_fuel_cost = df.loc[df['Account'] == 25, estimated_cost_col].values[0]
        annualized_fuel_cost = lump_fuel_cost*_crf(params['Interest Rate'], refueling_period_yr)
        df.loc[df['Account'] == 82, estimated_cost_col] = annualized_fuel_cost

    return df



def calculate_decommissioning_cost(df, params):
    # A78: Annualized Decommissioning Cost
    # Find the column name that starts with 'Estimated Cost option == "other costs"'
    estimated_cost_col_F = get_estimated_cost_column(df, 'F')
    estimated_cost_col_N = get_estimated_cost_column(df, 'N')

    for estimated_cost_col in [estimated_cost_col_F, estimated_cost_col_N ]:
        capex = df.loc[df['Account'].isin([10, 20]), estimated_cost_col].sum()
        AR = params['Annual Return']
        LP = params ['Levelization Period']
        
        if 'A78: CAPEX to Decommissioning Cost Ratio' not in params.keys():
            # Key is not specified. Use the default recommeneded value
            # PR#1: Chosen over previous unit_cost based estiamte
            # Estimating A78 based on a fraction of CAPEX suggested by 
            # Venneri, (2023) (15%) and INL/EXT-21-63067 (9%)
            params['A78: CAPEX to Decommissioning Cost Ratio'] = 0.15

        decommissioning_fv_cost = capex * params['A78: CAPEX to Decommissioning Cost Ratio']
        fv_to_pv_of_annuity = -AR/(1- pow(1+AR, LP))
        annualized_decommisioning_cost = decommissioning_fv_cost * fv_to_pv_of_annuity     
        df.loc[df['Account'] == 78, estimated_cost_col] = annualized_decommisioning_cost

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