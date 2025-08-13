# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import pandas as pd

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
    df['Adjusted Fixed Cost Low End ($)'] = df['Fixed Cost Low End'] * df['inflation_multiplier']
    df['Adjusted Fixed Cost High End ($)'] = df['Fixed Cost High End'] * df['inflation_multiplier']
    
    df['Adjusted Unit Cost ($)'] = df['Unit Cost'] * df['inflation_multiplier']
    df['Adjusted Unit Cost Low End ($)'] = df['Unit Cost Low End'] * df['inflation_multiplier']
    df['Adjusted Unit Cost High End ($)'] = df['Unit Cost High End'] * df['inflation_multiplier']


    

    # also read extra economic parameters that do not need to go through the inflation adjustment
    df_extra_params = pd.read_excel(file_name, sheet_name="Economics Parameters")
    # Create the dictionary
    extra_economic_parameters = dict(zip(df_extra_params["Parameter"], df_extra_params["Value"]))
    # Add the extra economic parameters to the params dictionary using a for loop
    for parameter, value in extra_economic_parameters.items():
        params[parameter] = value
    return df



