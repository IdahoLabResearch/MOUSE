
import pandas as pd
from cost.cost_utils import *
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



def escalate_cost_database(file_name, escalation_year):
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
        if params['enrichment'] < 0.1:
            cost_premium = 1
        elif  0.1 <= params['enrichment'] < 0.2:
            cost_premium = 1.15
        elif  0.2 <params['enrichment']:
            print("\033[91m ERROR: Enrichment is too high \033[0m")
        cost = cost_premium * unit_cost *pow(scaling_variable_value,exponent)    
    return cost
    


def scale_cost(initial_database, params):
    scaled_cost = initial_database[['Account', 'Level','Account Title']]
    
    escalation_year = params['escalation_year']
    

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
                # Check if there is a ref value for the scaling variable or just use the unit cost
                if row['Scaling Variable Ref Value'] > 0:
                    estimated_cost = fixed_cost +\
                    unit_cost * pow(scaling_variable_value,exponent) /(pow(scaling_variable_ref_value,exponent-1))

                else:
                    # Calculate the 'Estimated Cost
                    estimated_cost = fixed_cost + unit_cost * scaling_variable_value
            
            elif row['Standard Cost Equation?'] == 'nonstandard':
                estimated_cost = non_standard_cost_scale(row['Account'],\
                 unit_cost, scaling_variable_value, exponent, params)


            # Assign the calculated value to the corresponding row in the new DataFrame
            scaled_cost.at[index, f'Estimated Cost (${escalation_year } $))'] = estimated_cost
    return scaled_cost  

def find_children_accounts(df):
    # Find the column name that starts with "Estimated Cost"
    estimated_cost_column = [col for col in df.columns if col.startswith("Estimated Cost")][0]

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
                print("children now", children_str, "parent now", df.iloc[i]['Account'])

    # Assign the list to the DataFrame
    df['Children Accounts'] = children_accounts
    return df



def get_estimated_cost_column(df):
    for col in df.columns:
        if col.startswith("Estimated Cost"):
            return col
    return None

def calculate_high_level_accounts_cost(df, target_level, option):
    cost_column = get_estimated_cost_column(df)
    print(f"Updating costs of the level {target_level} accounts")

    # Determine the prefix condition based on the option parameter
    if option == "base":
        valid_prefixes = ('1', '2')
    elif option == "other":
        valid_prefixes = ('3', '4', '5')
    else:
        raise ValueError("Invalid option. Choose 'base' or 'other'.")

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Check if the account starts with the valid prefixes
        if str(row["Account"]).startswith(valid_prefixes):
            if row["Level"] == target_level and pd.isna(row[cost_column]):
                print(f"Updating Account {row['Account']}")
                children_accounts = row["Children Accounts"]
                print(children_accounts)
                
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
        df_updated = calculate_high_level_accounts_cost(df_with_children_accounts, level, option)
    return df_updated




def calculate_accounts_31_32_cost( df):
    
    # Find the column name that starts with 'Estimated Cosoption == "other costs"t'
    estimated_cost_col = get_estimated_cost_column(df)

    # Filter the DataFrame for accounts 21, 22, and 23
    filtered_df = df[df['Account'].isin([21, 22, 23])]

    # Sum the values in the 'Estimated Cost' column for the filtered accounts
    tot_field_direct_cost = filtered_df[estimated_cost_col].sum()

    acct_31_cost = 0.07 * tot_field_direct_cost # This ratio is based on MARVEL
    df.loc[df['Account'] == 31, estimated_cost_col] = acct_31_cost

    # To calculate the cost of factory and construction supervision (Account 32), 
    # the ratio of the factory and field indirect costs (Account 31) to the reactor systems cost (account 22) 
    # is calculated and multiplied by the cost of structures and improvements (Account 21)
    df.loc[df['Account'] == 32, estimated_cost_col] = df.loc[df['Account'] == 21, estimated_cost_col].values[0] * (df.loc[df['Account'] == 31, estimated_cost_col].values[0] / df.loc[df['Account'] == 22, estimated_cost_col].values[0])
    return df

# def update_accounts_30_40(df):
#     # Ensure 'Account' is treated as a string for proper filtering
#     df['Account'] = df['Account'].astype(str)
    
#     # Identify the column that starts with 'Estimated Cost'
#     cost_column = get_estimated_cost_column(df)

#     # Calculate the sum for account 30
#     accounts_30s = [str(i) for i in range(31, 40)]
#     df_30s = df[df['Account'].isin(accounts_30s)]
#     sum_30s = df_30s[cost_column].sum()

#     # Calculate the sum for account 40
#     accounts_40s = [str(i) for i in range(41, 50)]
#     df_40s = df[df['Account'].isin(accounts_40s)]
#     sum_40s = df_40s[cost_column].sum()

#     # Update or add the rows for accounts 30 and 40
#     if '30' in df['Account'].values:
#         df.loc[df['Account'] == '30', cost_column] = sum_30s
#     else:
#         df = df.append({'Account': '30', cost_column: sum_30s}, ignore_index=True)
        
#     if '40' in df['Account'].values:
#         df.loc[df['Account'] == '40', cost_column] = sum_40s
#     else:
#         df = df.append({'Account': '40', cost_column: sum_40s}, ignore_index=True)

#     return df
   
def calculate_high_level_capital_costs(df):
     # List of accounts to sum
    accounts_to_sum = [10, 20, 30, 40]

    # Find the column that starts with "Estimated Cost"
    cost_column = get_estimated_cost_column(df)
    
    # Calculate the sum of costs for the specified accounts
    occ_cost = df[df['Account'].isin(accounts_to_sum)][cost_column].sum()
    print(occ_cost)
    print(df[df['Account'].isin(accounts_to_sum)][cost_column], 'OCC')
    # Update the existing account "OCC" with the new total cost
    df = df._append({'Account': 'OCC','Account Title' : 'Overnight Capital Cost' ,cost_column: occ_cost}, ignore_index=True)
    return df


def bottom_up_cost_estimate(cost_database_filename, params):
    escalated_cost = escalate_cost_database(cost_database_filename, params['escalation_year'])
    escalated_cost_cleaned = remove_irrelevant_account(escalated_cost, params)
    scaled_cost = scale_cost(escalated_cost_cleaned, params)
    updated_cost = update_high_level_costs(scaled_cost, 'base' )
    updated_cost_with_indirect_cost = calculate_accounts_31_32_cost(updated_cost)
    updated_accounts_10_40 = update_high_level_costs(updated_cost_with_indirect_cost, 'other' )
    high_Level_catpial_cost = calculate_high_level_capital_costs(updated_accounts_10_40)
    return high_Level_catpial_cost

