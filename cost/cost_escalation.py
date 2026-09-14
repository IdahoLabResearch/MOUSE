# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import pandas as pd
import numpy as np

from cost.cost_database_csv import read_cost_database_sheet

# **************************************************************************************************************************
#                                                Sec. 0 :Inflation
# **************************************************************************************************************************


def calculate_inflation_multiplier(file_path, base_dollar_year, cost_type, escalation_year):
    
    base_dollar_year = int(base_dollar_year)
    escalation_year  = int(escalation_year)
    
    df = read_cost_database_sheet(file_path, "Inflation Adjustment")
    # print("Shape:", df.shape)
    # print("First 5 rows raw:")
    # print(df.head(5))
    # print("\nYear column raw values:")
    # print(df['Year'].tolist())

    df = df.dropna(subset=['Year'])
    df['Year'] = df['Year'].astype(int)

    # --- DEBUG: remove after issue is resolved ---
    # print("Year column dtype:", df['Year'].dtype)
    # print("Year values:", df['Year'].values)
    # print("Looking for:", base_dollar_year, type(base_dollar_year))

    if base_dollar_year not in df['Year'].values:
        print(f"\033[91mBase Year : {base_dollar_year} not found in the Excel file.\033[0m")

    if escalation_year not in df['Year'].values:
        print(f"\033[91mEscalation Year:  {escalation_year} not found in the Excel file.\033[0m")

    if cost_type == 'NA':
        multiplier = 1
    else:    
        multiplier = df.loc[df['Year'] == base_dollar_year, cost_type].values[0] / \
                     df.loc[df['Year'] == escalation_year, cost_type].values[0]

    return multiplier


# # **************************************************************************************************************************
# #                                                Sec. 1 : Baseline Costs (dollars)
# # **************************************************************************************************************************


def resolve_value(val, params):
    """
    If val is numeric, return it.
    If val is a string, look it up in params.
    Otherwise return NaN.
    """
    if pd.isna(val):
        return np.nan
    if isinstance(val, (int, float, np.number)):
        return float(val)
    if isinstance(val, str):
        if val in params:
            return float(params[val])
        else:
            raise KeyError(f"Parameter '{val}' not found in params.")
    return np.nan


def escalate_cost_database(file_name, escalation_year, params, sheet_name="Cost Database"):
    """
    Reads an Excel file with a specified sheet name into a Pandas DataFrame.
    Escalates fixed and unit costs, allowing cost fields to reference params.
    """

    # Regenerate the Git-reviewable CSV mirror, then read the selected sheet
    # from that mirror so the runtime input exactly matches the tracked data.
    df = read_cost_database_sheet(file_name, sheet_name)

    # Helper function to resolve numeric or parameter-referenced values
    def resolve_value(val, params):
        if pd.isna(val):
            return np.nan
        if isinstance(val, (int, float, np.number)):
            return float(val)
        if isinstance(val, str):
            if val in params:
                return float(params[val])
            else:
                raise KeyError(f"Parameter '{val}' not found in params.")
        return np.nan

    # Resolve all cost columns to numeric values
    cost_columns = [
        'Fixed Cost ($)',
        'Fixed Cost Low End',
        'Fixed Cost High End',
        'Unit Cost',
        'Unit Cost Low End',
        'Unit Cost High End'
    ]

    for col in cost_columns:
        df[col] = df[col].apply(lambda x: resolve_value(x, params))

    # Initialize an empty list to store inflation multipliers
    inflation_multipliers = []

    # Iterate through each row in the DataFrame
    for _, row in df.iterrows():
        if not pd.isna(row['Fixed Cost ($)']) or not pd.isna(row['Unit Cost']):
            multiplier = calculate_inflation_multiplier(
                file_name,
                row['Dollar Year'],
                row['Type'],
                escalation_year
            )
        else:
            multiplier = 0

        inflation_multipliers.append(multiplier)

    # Assign the inflation multiplier column
    df['inflation_multiplier'] = inflation_multipliers

    # Inflation-adjusted columns
    df['Adjusted Fixed Cost ($)'] = df['Fixed Cost ($)'] * df['inflation_multiplier']
    df['Adjusted Fixed Cost Low End ($)'] = df['Fixed Cost Low End'] * df['inflation_multiplier']
    df['Adjusted Fixed Cost High End ($)'] = df['Fixed Cost High End'] * df['inflation_multiplier']

    df['Adjusted Unit Cost ($)'] = df['Unit Cost'] * df['inflation_multiplier']
    df['Adjusted Unit Cost Low End ($)'] = df['Unit Cost Low End'] * df['inflation_multiplier']
    df['Adjusted Unit Cost High End ($)'] = df['Unit Cost High End'] * df['inflation_multiplier']

    # Read extra economic parameters (no escalation)
    df_extra_params = read_cost_database_sheet(file_name, "Economics Parameters")
    extra_economic_parameters = dict(zip(df_extra_params["Parameter"], df_extra_params["Value"]))

    for parameter, value in extra_economic_parameters.items():
        params[parameter] = value

    return df
