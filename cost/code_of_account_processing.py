import pandas as pd

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
            if col.startswith("FOAK Estimated Cost ("):
                return col
    elif option == 'N'   :
        for col in df.columns:
            if col.startswith("NOAK Estimated Cost ("):
                return col       
    elif option == 'F std'   :
        for col in df.columns:
            if col.startswith("FOAK Estimated Cost std ("):
                return col  
    elif option == 'N std'   :
        for col in df.columns:
            if col.startswith("NOAK Estimated Cost std ("):
                return col                              
    return None



def create_cost_dictionary(df, params, tracked_params_list):
    # create a dictionary of costs we are interested in tracking
    
    # start with params we are tracking
    filtered_params = {key: params[key] for key in tracked_params_list if key in params}

    # Initialize the dictionary with all required accounts and default values
    accounts = ['OCC', 'OCC per kW', 'OCC excl. fuel', 'OCC excl. fuel per kW', 'TCI', 'TCI per kW', 'AC', 'AC per MWh', 'LCOE']
    cost_dict = {}
    
    for account in accounts:
        cost_dict[f"{account}_FOAK Estimated Cost"] = None
        cost_dict[f"{account}_NOAK Estimated Cost"] = None
        cost_dict[f"{account}_FOAK Estimated Cost std"] = None
        cost_dict[f"{account}_NOAK Estimated Cost std"] = None
    
    # Populate the dictionary with values from the dataframe
    for _, row in df.iterrows():
        account = row['Account']
        if account in accounts:
            cost_dict[f"{account}_FOAK Estimated Cost"] =     row[get_estimated_cost_column(df, 'F')]
            cost_dict[f"{account}_NOAK Estimated Cost"] =     row[get_estimated_cost_column(df, 'N')]
            cost_dict[f"{account}_FOAK Estimated Cost std"] = row[get_estimated_cost_column(df, 'F std')]
            cost_dict[f"{account}_NOAK Estimated Cost std"] = row[get_estimated_cost_column(df, 'N std')]  
    
    
    filtered_params.update(cost_dict)

    return filtered_params

