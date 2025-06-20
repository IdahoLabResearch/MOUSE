
import pandas as pd
import numpy as np
from reactor_engineering_evaluation.operation import *

# **************************************************************************************************************************
#                                                Sec. 0 :Inflation
# **************************************************************************************************************************


def calculate_inflation_multiplier(file_path, base_dollar_year, cost_type, escalation_year):
 
    # Read the Excel file
    df = pd.read_excel(file_path, sheet_name = "Inflation Adjustment")

     # Check if the target year exists in the dataframe
    if base_dollar_year not in df['Year'].values:
        print(f"\033[91mBase Year : {base_dollar_year} not found in the Excel file.\033[0m")
        

    if escalation_year not in df['Year'].values:
        print(f"\033[91mEscalation Year:  {escalation_year} not found in the Excel file.\033[0m")



    if cost_type =='NA':
        multiplier = 1
    else:    
        # Get the multiplier for the target year
        multiplier = df.loc[df['Year'] == base_dollar_year, cost_type].values[0] / df.loc[df['Year'] == escalation_year, cost_type].values[0]

    return multiplier    

            

# # **************************************************************************************************************************
# #                                                Sec. 1 : Baseline Costs (dollars)
# # **************************************************************************************************************************



def escalate_cost_database(file_name, escalation_year, params):
    """
    Reads an Excel file with a specified sheet name into a Pandas DataFrame.
    
    Parameters:
    file_name (str): The name of the Excel file.
    sheet_name (str): The name of the sheet to read from the Excel file.
    
    Returns:
    pd.DataFrame: A DataFrame containing the data from the specified sheet.
    """
    # Read the Excel file into a Pandas DataFrame
    df = pd.read_excel(file_name, sheet_name= "Cost Database")

        
    # Initialize an empty list to store the results (inflation multipliers)
    inflation_multipliers = []

    # Iterate through each row in the dataframe
    for idx, row in df.iterrows():
        # Check if the value in 'Fixed Cost ($)' is non-NaN
        if not pd.isna(row['Fixed Cost ($)']) or not pd.isna(row['Unit Cost']):
                # Calculate the inflation multiplier for the current row
                multiplier = calculate_inflation_multiplier(file_name, row['Dollar Year'], row['Type'], escalation_year)
                
        else:
                multiplier = 0    
        
        # Append the result to the list
        inflation_multipliers.append(multiplier)    


    # Assign the list of results to the new column in the dataframe
    df['inflation_multiplier'] = inflation_multipliers


    # Add the inflation-adjusted columns
    df['Adjusted Fixed Cost ($)'] = df['Fixed Cost ($)'] * df['inflation_multiplier']
    df['Adjusted Unit Cost ($)'] = df['Unit Cost'] * df['inflation_multiplier']
    

    # also read extra economic parameters that do not need to go through the inflation adjustment
    df_extra_params = pd.read_excel(file_name, sheet_name="Economics Parameters")
    # Create the dictionary
    extra_economic_parameters = dict(zip(df_extra_params["Parameter"], df_extra_params["Value"]))
    # Add the extra economic parameters to the params dictionary using a for loop
    for parameter, value in extra_economic_parameters.items():
        params[parameter] = value
    return df



def remove_irrelevant_account(df, params):
    # Filter the dataframe
    indices_to_drop = []
    
    # Iterate over the rows and print the "Optional Variable" if it's not NaN
    for index, row in df.iterrows():
        if not pd.isna(row['Optional Variable']):
            if row['Optional Variable'] in params and params[row['Optional Variable']] == row['Optional Value']:
                    print("\n")
                    print(f"For the cost of the Account {row['Account']}: {row['Account Name']}, the {row['Optional Variable']} is selected to be {row['Optional Value']}")
            else:
                indices_to_drop.append(index)   
    
    # Drop the rows with collected indices (irrelevant accounts)
    df.drop(indices_to_drop, inplace=True)

    return df 

def non_standard_cost_scale(account, unit_cost, scaling_variable_value, exponent, params):
    

    # pumps
    if account == 222.11 or account == 222.12:
        cost_multiplier = (0.2 / (1 - params['Pump Isentropic Efficiency'])) + 1
        cost = cost_multiplier * unit_cost * pow(scaling_variable_value,exponent)
    
    # compressors
    elif account == 222.13:
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
    
    elif account == 78:
        AR = params['Annual Return']
        LP = params ['Levelization Period']
        cost_multiplier = - AR/(1- pow(1+AR, LP))
        cost = cost_multiplier * unit_cost * pow(scaling_variable_value,exponent)       
    elif account == 81:
        cost_multiplier =  params['FTEs Per Operator Per Year Per Refueling'] 
        cost = cost_multiplier * unit_cost * pow(scaling_variable_value, exponent)
    return cost
    


def scale_cost(initial_database, params):
    scaled_cost = initial_database[['Account', 'Level','Account Title', 'FOAK to NOAK Multiplier Type']]
    
    escalation_year = params['Escalation Year']
    

    # Iterate through each row in the DataFrame
    for index, row in initial_database.iterrows():
        
        # Check if cost data are available (fixed or unit cost)
        if row['Fixed Cost ($)'] > 0 or	row['Unit Cost'] > 0:
            
            scaling_variable_value = params[row['Scaling Variable']] if pd.notna(row['Scaling Variable']) else 0
            
            # Calculate the 'Estimated Cost
            fixed_cost = row['Adjusted Fixed Cost ($)'] if pd.notna(row['Fixed Cost ($)']) else 0
            unit_cost = row['Adjusted Unit Cost ($)'] if pd.notna(row['Unit Cost']) else 0
            scaling_variable_ref_value  = row['Scaling Variable Ref Value']
            exponent = row['Exponent']

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
            scaled_cost.at[index, f'FOAK Estimated Cost (${escalation_year } $)'] = estimated_cost
    return scaled_cost  

def find_children_accounts(df):
    # Find the column name that starts with "Estimated Cost"
    estimated_cost_column = [col for col in df.columns if col.startswith("FOAK Estimated Cost")][0]

    # Initialize a list for children accounts
    children_accounts = [None] * len(df)
    
    for target_level in range(4, -1, -1):
        source_level = target_level + 1
        # Iterate over the dataframe
        for i in range(len(df)):
            if df.iloc[i]['Level'] == target_level and pd.isna(df.iloc[i][estimated_cost_column]):
                children = []
                for j in range(i + 1, len(df)):
                    if df.iloc[j]['Level'] == source_level:
                        children.append(str(df.iloc[j]['Account']))  # Convert to string
                    elif df.iloc[j]['Level'] < source_level:
                        break
                # Convert the list to a comma-separated string
                children_str = ','.join(children) if children else None
                children_accounts[i] = children_str

    # Assign the list to the DataFrame
    df['Children Accounts'] = children_accounts
    return df



def get_estimated_cost_column(df, option):
    if option == 'F':
        for col in df.columns:
            if col.startswith("FOAK Estimated Cost"):
                return col
    elif option == 'N'   :
        for col in df.columns:
            if col.startswith("NOAK Estimated Cost"):
                return col         
    return None




def calculate_high_level_accounts_cost(df, target_level, option, FOAK_or_NOAK):
    cost_column = get_estimated_cost_column(df, FOAK_or_NOAK)
    print(f"Updating costs of the level {target_level} accounts for the {cost_column}")

    # Determine the prefix condition based on the option parameter
    if option == "base":
        valid_prefixes = ('1', '2')
    elif option == "other":
        valid_prefixes = ('3', '4', '5')
    elif option == "finance": 
        valid_prefixes = ('6')  
    elif option == "annual": 
        valid_prefixes = ('7', '8')      
    else:
        raise ValueError("Invalid option. Choose 'base' or 'other' or 'finance' or 'annual'.")

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Check if the account starts with the valid prefixes
        if str(row["Account"]).startswith(valid_prefixes):
            if row["Level"] == target_level and pd.isna(row[cost_column]):
                print(f"Updating Account {row['Account']}")
                children_accounts = row["Children Accounts"]
                
                if not pd.isna(children_accounts):
                    children_accounts_list = children_accounts.split(",")

                    # Initialize the sum
                    total_sum = 0
                    # Iterate through each account in the children_accounts_list
                    for account in children_accounts_list:
                        # Convert the account to a float
                        account_value = float(account)
                        # Add the corresponding value from the DataFrame to the total_sum
                        total_sum += df[df["Account"] == account_value][cost_column].values[0]
                    df.at[index, cost_column] = total_sum

    return df



def update_high_level_costs(scaled_cost, option):
    # input is the scaled cost
    df_with_children_accounts =  find_children_accounts(scaled_cost)
    for level in range(4, -1, -1):
        df_updated = calculate_high_level_accounts_cost(df_with_children_accounts, level, option, 'F') # FOAK
        df_updated_2 = calculate_high_level_accounts_cost(df_updated, level, option, 'N')     # NOAK
    return df_updated_2


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

def save_params_to_excel_file(excel_file, params):
    # Convert the Parameters object to a dictionary
    params_dict = dict(params)
    # Create a DataFrame from the dictionary
    df = pd.DataFrame(list(params_dict.items()), columns=['Parameter', 'Value'])
    # Write the DataFrame to an Excel file with a specified sheet name
    df.to_excel(excel_file, sheet_name='Parameters', index=False)

def transform_dataframe(df):
    """
    Divides all values in the specified column by one million, except the last two rows,
    rounds to one non-zero digit after the decimal point, and appends 'M'. Converts the last two rows to integers.

    Parameters:
    df (pd.DataFrame): The dataframe containing the data.

    Returns:
    pd.DataFrame: The modified dataframe.
    """
    df = df.drop(columns=['Children Accounts', 'Level', 'FOAK to NOAK Multiplier Type', 'Multiplier'])

    column_name_F = get_estimated_cost_column(df, 'F')
    column_name_N = get_estimated_cost_column(df, 'N')
    for column_name in [column_name_F , column_name_N]:
  
    
        # Dividing all the elements in 'column_name' by a million, except the last two rows
        df.iloc[:-2, df.columns.get_loc(column_name)] = df.iloc[:-2, df.columns.get_loc(column_name)] / 1000000
        
        # Rounding to one non-zero digit after the decimal point
        df.iloc[:-2, df.columns.get_loc(column_name)] = df.iloc[:-2, df.columns.get_loc(column_name)].apply(lambda x: round(x, -int(np.floor(np.log10(abs(x)))) + 1) if x != 0 else 0)
        
        # Appending 'M' to the modified values
        df.iloc[:-2, df.columns.get_loc(column_name)] = df.iloc[:-2, df.columns.get_loc(column_name)].astype(str) + ' M'
        
        # Converting the last two rows to integers
        df.iloc[-2:, df.columns.get_loc(column_name)] = df.iloc[-2:, df.columns.get_loc(column_name)].astype(int)
    
    return df

def learning_rate_muliplier(learning_rate, number_of_units):
    return pow(1 -learning_rate ,np.log2(number_of_units))

def FOAK_to_NOAK(df, params):
    # Function to map the multiplier type to the respective value
    def get_multiplier(multiplier_type):
        if multiplier_type in ['No Learning', 'Licensing Learning', 'Factory Primary Structure', 'Factory Drums', 'Factory Other']:
            return params[multiplier_type]
        elif multiplier_type == 'Onsite Learning' :
            multiplier =  learning_rate_muliplier(params['Assumed Onsite Learning Rate'], params['Assumed Number Of Units For Onsite Learning'])
            return multiplier
        else:
            return np.nan
    
    # Apply the function to create the new column
    df['Multiplier'] = df['FOAK to NOAK Multiplier Type'].apply(get_multiplier)
    foak_col = get_estimated_cost_column(df, 'F')
    noak_column = foak_col.replace("FOAK", "NOAK")
    df[noak_column ] = df['Multiplier'] * df[foak_col ]
    return df



def bottom_up_cost_estimate(cost_database_filename, params, output_filename):
    escalated_cost = escalate_cost_database(cost_database_filename, params['Escalation Year'], params)
    escalated_cost_cleaned = remove_irrelevant_account(escalated_cost, params)
    reactor_operation(params)
    scaled_cost = scale_cost(escalated_cost_cleaned, params)
    NOAK_COA = FOAK_to_NOAK(scaled_cost, params)

    updated_cost = update_high_level_costs(scaled_cost, 'base' )
    updated_cost_with_indirect_cost = calculate_accounts_31_32_75_82_cost(updated_cost, params)
    updated_accounts_10_40 = update_high_level_costs(updated_cost_with_indirect_cost, 'other' )
    high_Level_capital_cost = calculate_high_level_capital_costs(updated_accounts_10_40, params)
    
    updated_accounts_10_60 = update_high_level_costs(high_Level_capital_cost , 'finance' )
    TCI = calculate_TCI(updated_accounts_10_60, params )
    updated_accounts_70_80 = update_high_level_costs(TCI , 'annual' )
    Final_COA = energy_cost_levelized(params, updated_accounts_70_80)
    presented_COA = transform_dataframe(Final_COA  )

    # Create an ExcelWriter object
    with pd.ExcelWriter(output_filename) as writer:
        # Save the presented_COA DataFrame to the first sheet
        presented_COA.to_excel(writer, sheet_name="cost estimate", index=False)
        # Save the parameters to the second sheet
        save_params_to_excel_file(writer, params)
        
    print(f"\n\nThe cost estimate and all the paramters are saved at {output_filename}\n\n")

