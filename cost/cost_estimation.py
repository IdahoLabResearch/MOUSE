import os
import pandas as pd
import numpy as np
import csv
from cost.cost_escalation import escalate_cost_database
from cost.code_of_account_processing import remove_irrelevant_account, get_estimated_cost_column, find_children_accounts, create_cost_dictionary
from cost.cost_scaling import scale_cost, scale_GCMR_accounts
from cost.non_direct_cost import calculate_accounts_31_32_75_82_cost, calculate_decommissioning_cost, calculate_high_level_capital_costs, calculate_TCI, energy_cost_levelized
from reactor_engineering_evaluation.operation import reactor_operation


def calculate_high_level_accounts_cost(df, target_level, option, FOAK_or_NOAK):
    cost_column = get_estimated_cost_column(df, FOAK_or_NOAK)
    # print(f"Updating costs of the level {target_level} accounts for the {cost_column}")

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
                # print(f"Updating Account {row['Account']}")
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

    # Select all numerical columns
    numerical_columns = df.select_dtypes(include=[np.number]).columns
    # Convert all values in the numerical columns to integers
    df[numerical_columns] = df[numerical_columns].astype(int)

    # Remove rows where all values in numerical columns are zero
    df = df.loc[~(df[numerical_columns] == 0).all(axis=1)]
    return df



def learning_rate_multiplier(learning_rate, number_of_units):
    # Assumes that learning plateau's after 100th unit! - Jack
    # This is an arbitrary stopping point.
    return pow(1-learning_rate, np.log2(min(100, number_of_units)))



def FOAK_to_NOAK(df, params):
    # Additional Cost Scaling Based on Assumed Learning Rate
    # Learning Rate and Cost multiplier is based on 
    # DOI: 10.1080/00295450.2023.2206779
    # Cost Multiplier is capped after the 100th Unit for any component
    if 'NOAK Unit Number' not in params.keys():
        # Assume the default value if no `NOAK Unit Number` is specified.
        params['NOAK Unit Number'] = 10
        # Custom Check to see if input specifies which Nth-of-a-Kind
        # Default is ~10th with 20(2*10)-units assumed for Onsite Learning
    params['Assumed Number Of Units For Onsite Learning'] = params['NOAK Unit Number'] * 2
    
    for multiplier_type in ['No Learning', 
                            'Licensing Learning', 
                            'Factory Primary Structure', 
                            'Factory Drums',
                            'Factory Other', 
                            'Factory Be',
                            'Factory BeO',
                            'Non-nuclear off-the-shelf']:
        params[f"{multiplier_type} Cost Multiplier"] = learning_rate_multiplier(params[f'{multiplier_type}'], 
                                                                                params['NOAK Unit Number'])
    params['Onsite Learning Cost Multiplier'] = learning_rate_multiplier(params['Onsite Learning'], 
                                                                         params['Assumed Number Of Units For Onsite Learning'])

    # Function to map the multiplier type to the respective value
    def get_multiplier(multiplier_type):
        if multiplier_type in ['No Learning', 
                               'Licensing Learning', 
                               'Factory Primary Structure', 
                               'Factory Drums',
                               'Factory Other', 
                               'Factory Be',
                               'Factory BeO',
                               'Onsite Learning',
                               'Non-nuclear off-the-shelf']:
            return params[f"{multiplier_type} Cost Multiplier"]
        else:
            return np.nan
    
    # Apply the function to create the new column
    df['Multiplier'] = df['FOAK to NOAK Multiplier Type'].apply(get_multiplier)
    foak_col = get_estimated_cost_column(df, 'F')
    noak_column = foak_col.replace("FOAK", "NOAK")
    df[noak_column] = df['Multiplier'] * df[foak_col]
    return df



def reorder_dataframe(df):
    # List the desired order of the first two columns
    first_columns = ['Account', 'Account Title']

    # List all other columns in the dataframe
    other_columns = [col for col in df.columns if col not in first_columns]

    # Create the new column order
    new_column_order = first_columns + other_columns

    # Reorder the dataframe columns
    df = df[new_column_order]
    return df



def bottom_up_cost_estimate(cost_database_filename, params):
    escalated_cost = escalate_cost_database(cost_database_filename, params['Escalation Year'], params)
    escalated_cost_cleaned = remove_irrelevant_account(escalated_cost, params)
    reactor_operation(params)

    COA_list = []
    # sampling
    for i in range(params['Number of Samples']):
        if (i + 1) % 100 == 0:
            print(f"\n\nSample # {i+1}")

        scaled_cost = scale_cost(escalated_cost_cleaned, params)
        scaled_cost = scale_GCMR_accounts(scaled_cost, params)
        NOAK_COA = FOAK_to_NOAK(scaled_cost, params)

        updated_cost = update_high_level_costs(scaled_cost, 'base' )
        updated_cost_with_indirect_cost = calculate_accounts_31_32_75_82_cost(updated_cost, params)
        cost_with_decommissioning = calculate_decommissioning_cost(updated_cost_with_indirect_cost, params)
        updated_accounts_10_40 = update_high_level_costs(cost_with_decommissioning, 'other' )
        high_Level_capital_cost = calculate_high_level_capital_costs(updated_accounts_10_40, params)
        
        updated_accounts_10_60 = update_high_level_costs(high_Level_capital_cost, 'finance' )
        TCI = calculate_TCI(updated_accounts_10_60, params )
        updated_accounts_70_80 = update_high_level_costs(TCI , 'annual' )
        Final_COA = energy_cost_levelized(params, updated_accounts_70_80)
        FOAK_column = get_estimated_cost_column(Final_COA, 'F')
        NOAK_column = get_estimated_cost_column(Final_COA, 'N')
        Final_COA = Final_COA[['Account', 'Account Title',  FOAK_column  ,NOAK_column]]
        
        COA_list.append(Final_COA)
    
    # Concatenate all dataframes along the rows
    concatenated_df = pd.concat(COA_list)
    # Select only the numeric columns
    numeric_columns = concatenated_df.select_dtypes(include='number').columns

    # Calculate the mean for each numeric column grouped by the index
    mean_df = concatenated_df[numeric_columns].groupby(concatenated_df.index).mean()

    # Calculate the standard deviation for each numeric column grouped by the index
    if params["Number of Samples"] > 1:
        std_df = concatenated_df[numeric_columns].groupby(concatenated_df.index).std()
    else:
        std_df = concatenated_df[numeric_columns].groupby(concatenated_df.index).std(ddof=0)

    # Add the standard deviation columns to the mean dataframe
    mean_df[FOAK_column.replace('Cost', 'Cost std')] = std_df[FOAK_column]
    mean_df[NOAK_column.replace('Cost', 'Cost std')] = std_df[NOAK_column]

    # If you need to retain the non-numeric columns, you can join them back to the mean_df
    non_numeric_columns = concatenated_df.select_dtypes(exclude='number').groupby(concatenated_df.index).first()
    result_df = mean_df.join(non_numeric_columns)
    reordered_df = reorder_dataframe(result_df)
    return reordered_df



def parametric_studies(cost_database_filename, params, tracked_params_list, output_csv_filename):
    detatiled_cost_table = bottom_up_cost_estimate(cost_database_filename, params)
    tracked_costs = create_cost_dictionary(detatiled_cost_table, params, tracked_params_list)
    
    file_exists = os.path.isfile(output_csv_filename)
    
    with open(output_csv_filename, 'a', newline='') as csvfile:
        fieldnames = tracked_costs.keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        if not file_exists or os.stat(output_csv_filename).st_size == 0:
            writer.writeheader()
        
        writer.writerow(tracked_costs)
        print(f"Results are being saved on {output_csv_filename}")



def detailed_bottom_up_cost_estimate(cost_database_filename, params, output_filename):
    detatiled_cost_table = bottom_up_cost_estimate(cost_database_filename, params)
    pretty_df = transform_dataframe(detatiled_cost_table)
    # Create an ExcelWriter object
    with pd.ExcelWriter(output_filename) as writer:
        # Save the presented_COA DataFrame to the first sheet
        pretty_df.to_excel(writer, sheet_name="cost estimate", index=False)
        # Save the parameters to the second sheet
        save_params_to_excel_file(writer, params)
        
    print(f"\n\nThe cost estimate and all the paramters are saved at {output_filename}\n\n")

    # Return final DataFrame for post processing
    return detatiled_cost_table