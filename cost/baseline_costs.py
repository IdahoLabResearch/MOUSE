
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


    # df['inflation_multiplier'] = df.apply(lambda row: calculate_inflation_multiplier(file_name, row['Dollar Year'], row['Type'], escalation_year), axis=1)

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



def update_high_level_costs(df1):
    # Find the column that starts with 'Estimated Cost'
    estimated_cost_column = [col for col in df1.columns if col.startswith('Estimated Cost')][0]

    # Function to update costs for a given level based on the next level
    def update_level_costs(current_level, next_level):
        for index, row in df1.iterrows():
            account_starts_with = str(row['Account'])[0]
            condition = account_starts_with in ('1', '2')  
            
            if row['Level'] == current_level and pd.isna(row[estimated_cost_column]) and condition:
                sum_cost = 0
                for i in range(index + 1, len(df1)):
                    if df1.iloc[i]['Level'] == next_level:
                        sum_cost += df1.iloc[i][estimated_cost_column]
                    elif df1.iloc[i]['Level'] < next_level:
                        break
                df1.at[index, estimated_cost_column] = sum_cost

    # Update levels starting from the highest to the lowest
    for current_level in range(4, -1, -1):
        update_level_costs(current_level, current_level + 1)
    
    return df1

# def update_high_level_costs(df, option):
#     # Find the column name that starts with 'Estimated Cost'
#     estimated_cost_col = [col for col in df.columns if col.startswith('Estimated Cost')][0]
    
#     # Define the levels to process
#     levels = [4, 3, 2, 1, 0]

#     # Define the account digit conditions based on the option
#     if option == "base cost":
#         valid_first_digits = ['1', '2']
#     elif option == "other costs":
#         valid_first_digits = ['3', '4', '5']
#     else:
#         raise ValueError("Invalid option. Please choose 'base cost' or 'other costs'.")

#     for level in levels:

#         # Iterate over the rows
#         for i in range(len(df)):
#             if pd.isna(df.iloc[i, df.columns.get_loc(estimated_cost_col)]):


#                     if df.iloc[i, df.columns.get_loc('Level')] == level:   
#                         total_cost = 0
#                         for j in range(i + 1, len(df)):

#                                 if df.iloc[j, df.columns.get_loc('Level')] <= level:
#                                     break
#                                 if df.iloc[j, df.columns.get_loc('Level')] == level + 1:
#                                     account_num = str(df.iloc[j, df.columns.get_loc('Account')])
                                  
#                                     if account_num[0] in valid_first_digits:
#                                         total_cost += df.iloc[j, df.columns.get_loc(estimated_cost_col)] 
#                         df.iloc[i, df.columns.get_loc(estimated_cost_col)] = total_cost
    
#     return df

def calculate_accounts_31_32_cost( df):
    
    # Find the column name that starts with 'Estimated Cosoption == "other costs"t'
    estimated_cost_col = [col for col in df.columns if col.startswith('Estimated Cost')][0]

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

def update_accounts_30_40(df):
    # Ensure 'Account' is treated as a string for proper filtering
    df['Account'] = df['Account'].astype(str)
    
    # Identify the column that starts with 'Estimated Cost'
    cost_column = [col for col in df.columns if col.startswith('Estimated Cost')][0]

    # Calculate the sum for account 30
    accounts_30s = [str(i) for i in range(31, 40)]
    df_30s = df[df['Account'].isin(accounts_30s)]
    sum_30s = df_30s[cost_column].sum()

    # Calculate the sum for account 40
    accounts_40s = [str(i) for i in range(41, 50)]
    df_40s = df[df['Account'].isin(accounts_40s)]
    sum_40s = df_40s[cost_column].sum()

    # Update or add the rows for accounts 30 and 40
    if '30' in df['Account'].values:
        df.loc[df['Account'] == '30', cost_column] = sum_30s
    else:
        df = df.append({'Account': '30', cost_column: sum_30s}, ignore_index=True)
        
    if '40' in df['Account'].values:
        df.loc[df['Account'] == '40', cost_column] = sum_40s
    else:
        df = df.append({'Account': '40', cost_column: sum_40s}, ignore_index=True)

    return df
   


def bottom_up_cost_estimate(cost_database_filename, params):
    escalated_cost = escalate_cost_database(cost_database_filename, params['escalation_year'])
    escalated_cost_cleaned = remove_irrelevant_account(escalated_cost, params)
    scaled_cost = scale_cost(escalated_cost_cleaned, params)
    updated_cost = update_high_level_costs(scaled_cost) #, option = "base cost") # update base cost: pre-construction and direct costs
    updated_cost_with_indirect_cost = calculate_accounts_31_32_cost(updated_cost)
    updated_accounts_10_40 = update_accounts_30_40(updated_cost_with_indirect_cost)
    return updated_accounts_10_40

