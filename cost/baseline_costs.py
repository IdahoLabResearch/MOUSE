
import pandas as pd
from cost.cost_utils import *
# **************************************************************************************************************************
#                                                Sec. 0 :Inflation
# **************************************************************************************************************************


def calculate_inflation_multiplier(file_path, base_dollar_year, cost_type, escalation_year):
    # Read the Excel file
    df = pd.read_excel(file_path, sheet_name = "Inflation Adjustment")

    # Check if the target year exists in the dataframe
    if base_dollar_year in df['Year'].values:
        if  escalation_year in df['Year'].values:
            # Get the multiplier for the target year
            multiplier = df.loc[df['Year'] == base_dollar_year, cost_type].values[0] / df.loc[df['Year'] == escalation_year, cost_type].values[0]
            
            return multiplier
        else:  
            return f"Year {escalation_year} not found in the Excel file."  
    else:
        return f"Year {base_dollar_year} not found in the Excel file."

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
    df['inflation_multiplier'] = df.apply(lambda row: calculate_inflation_multiplier(file_name, row['Dollar Year'], row['Type'], escalation_year), axis=1)

    # Add the inflation-adjusted columns
    df['Adjusted Fixed Cost ($)'] = df['Fixed Cost ($)'] * df['inflation_multiplier']
    df['Adjusted Unit Cost ($)'] = df['Unit Cost ($)'] * df['inflation_multiplier']
    
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
        if row['Fixed Cost ($)'] > 0 or	row['Unit Cost ($)'] > 0:
            
            scaling_variable_value = params[row['Scaling Variable']] if pd.notna(row['Scaling Variable']) else 0
            
            # Calculate the 'Estimated Cost
            fixed_cost = row['Adjusted Fixed Cost ($)'] if pd.notna(row['Fixed Cost ($)']) else 0
            unit_cost = row['Adjusted Unit Cost ($)'] if pd.notna(row['Unit Cost ($)']) else 0
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
            scaled_cost.at[index, f'Estimated Cost (${escalation_year } M$))'] = np.round(estimated_cost/ 1000000,2)
    return scaled_cost  




def update_high_level_costs(df, option):
    # Find the column name that starts with 'Estimated Cost'
    estimated_cost_col = [col for col in df.columns if col.startswith('Estimated Cost')][0]
    
    # Define the levels to process
    levels = [4, 3, 2, 1, 0]

    # Define the account digit conditions based on the option
    if option == "base cost":
        valid_first_digits = ['1', '2']
    elif option == "other costs":
        valid_first_digits = ['3', '4', '5']
    else:
        raise ValueError("Invalid option. Please choose 'base cost' or 'other costs'.")

    for level in levels:
        # Iterate over the rows
        for i in range(len(df)):
            if pd.isna(df.iloc[i, df.columns.get_loc(estimated_cost_col)]):
                if df.iloc[i, df.columns.get_loc('Level')] == level:   
                    total_cost = 0
                    for j in range(i + 1, len(df)):
                        if df.iloc[j, df.columns.get_loc('Level')] <= level:
                            break
                        if df.iloc[j, df.columns.get_loc('Level')] == level + 1:
                            account_num = str(df.iloc[j, df.columns.get_loc('Account')])
                            if account_num[0] in valid_first_digits:
                                total_cost += df.iloc[j, df.columns.get_loc(estimated_cost_col)] 
                    df.iloc[i, df.columns.get_loc(estimated_cost_col)] = total_cost
    df = df.drop(columns=['Level'])
    return df

def add_indirect_cost(df):
         

def bottom_up_cost_estimate(cost_database_filename, params):
    escalated_cost = escalate_cost_database(cost_database_filename, params['escalation_year'])
    escalated_cost_cleaned = remove_irrelevant_account(escalated_cost, params)
    scaled_cost = scale_cost(escalated_cost_cleaned, params)
    updated_cost = update_high_level_costs(scaled_cost, option = "base cost") # update base cost: pre-construction and direct costs
    return updated_cost