# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import pandas as pd

def remove_irrelevant_account(df, params):
    indices_to_drop = []
    
    has_sec_optional = 'Sec Optional Variable' in df.columns  # ← add this check

    for index, row in df.iterrows():
        def _optional_matches(param_val, expected_val):
            """Return True if param_val equals expected_val, or if param_val is a list that contains expected_val."""
            def _normalized(value):
                return value.strip() if isinstance(value, str) else value

            if isinstance(param_val, list):
                return _normalized(expected_val) in [_normalized(value) for value in param_val]
            return _normalized(param_val) == _normalized(expected_val)

        # Check for 'Optional Variable'
        if not pd.isna(row['Optional Variable']):
            if row['Optional Variable'] in params and _optional_matches(params[row['Optional Variable']], row['Optional Value']):
                print("\n")
                print(f"For the cost of the Account {row['Account']}: {row['Account Name']}, the {row['Optional Variable']} is selected to be {row['Optional Value']}")
                # Append the selected optional value to Account Title for clarity in the output
                df.at[index, 'Account Title'] = str(row['Account Title']) + ' - ' + str(row['Optional Value'])
            else:
                indices_to_drop.append(index)
                continue

        # Check for 'Sec Optional Variable' only if column exists
        if has_sec_optional and not pd.isna(row['Sec Optional Variable']):
            if row['Sec Optional Variable'] in params and _optional_matches(params[row['Sec Optional Variable']], row['Sec Optional Value']):
                print("\n")
                print(f"For the cost of the Account {row['Account']}: {row['Account Name']}, the {row['Sec Optional Variable']} is selected to be {row['Sec Optional Value']}")
                # Also append the sec optional value
                df.at[index, 'Account Title'] = str(df.at[index, 'Account Title']) + ' - ' + str(row['Sec Optional Value'])
            else:
                indices_to_drop.append(index)
                continue

    df.drop(indices_to_drop, inplace=True)
    return df
    return df



def find_children_accounts(df):
    # Find the column name that starts with "Estimated Cost"
    estimated_cost_column = [col for col in df.columns if col.startswith("FOAK Estimated Cost")][0]

    # Pull columns once into numpy arrays. df.iloc[i][col] in a nested loop
    # is dominated by per-call overhead — numpy indexing on a flat array
    # is 10–50× faster for the same logic.
    levels = df['Level'].to_numpy()
    is_nan_cost = pd.isna(df[estimated_cost_column].to_numpy())
    index_strs = [str(idx) for idx in df.index]
    n = len(df)

    children_accounts = [None] * n

    for target_level in range(4, -1, -1):
        source_level = target_level + 1
        for i in range(n):
            if levels[i] == target_level and is_nan_cost[i]:
                children = []
                for j in range(i + 1, n):
                    lvl_j = levels[j]
                    if lvl_j == source_level:
                        children.append(index_strs[j])
                    elif lvl_j < source_level:
                        break
                children_accounts[i] = ','.join(children) if children else None

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

    # Base accounts that are always tracked regardless of tax credit selection
    base_accounts = [
        'OCC', 'OCC per kW',
        'OCC excl. fuel', 'OCC excl. fuel per kW',
        'TCI', 'TCI per kW',
        'AC', 'AC per MWh',
        'LCOE'
    ]

    # Physics safety metrics — tracked from params directly (not from the cost dataframe)
    # These are always included if present in params; set to nan if not calculated
    # (e.g. when Shutdown Margin Calc or Isothermal Temperature Coefficients are False)
    physics_metrics = ['Temp Coeff 3D (2D corrected)', 'SDM 3D (2D corrected)']
    for metric in physics_metrics:
        if metric in params.keys():
            filtered_params[metric] = params[metric]

    # ITC-related accounts — only present if user provided 'ITC credit level' in params
    itc_accounts = [
        'OCC (ITC-adjusted)', 'OCC (ITC-adjusted) per kW',
        'TCI (ITC-adjusted)', 'TCI (ITC-adjusted) per kW',
        'LCOE (ITC-adjusted)'
    ] if 'ITC credit level' in params.keys() else []

    # PTC-related accounts — only present if user provided 'PTC credit value' in params
    ptc_accounts = [
        'LCOE with PTC'
    ] if 'PTC credit value' in params.keys() else []

    # Combine all accounts to track
    accounts = base_accounts + itc_accounts + ptc_accounts

    cost_dict = {}
    
    for account in accounts:
        cost_dict[f"{account}_FOAK Estimated Cost"] = None
        cost_dict[f"{account}_NOAK Estimated Cost"] = None
        cost_dict[f"{account}_FOAK Estimated Cost std"] = None
        cost_dict[f"{account}_NOAK Estimated Cost std"] = None
    
    # Populate the dictionary with values from the dataframe
    # If an account doesn't exist in the dataframe (e.g. ITC/PTC not used), it stays None
    for _, row in df.iterrows():
        account = row['Account']
        if account in accounts:
            cost_dict[f"{account}_FOAK Estimated Cost"] =     row[get_estimated_cost_column(df, 'F')]
            cost_dict[f"{account}_NOAK Estimated Cost"] =     row[get_estimated_cost_column(df, 'N')]
            cost_dict[f"{account}_FOAK Estimated Cost std"] = row[get_estimated_cost_column(df, 'F std')]
            cost_dict[f"{account}_NOAK Estimated Cost std"] = row[get_estimated_cost_column(df, 'N std')]  
    
    filtered_params.update(cost_dict)

    return filtered_params
