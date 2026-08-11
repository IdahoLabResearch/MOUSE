# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import pandas as pd
import numpy as np

# **************************************************************************************************************************
#                                                Sec. 0 :Inflation
# **************************************************************************************************************************


def calculate_inflation_multiplier(file_path, base_dollar_year, cost_type, escalation_year):
    if pd.isna(base_dollar_year):
        raise ValueError("Dollar Year is missing.")
    if pd.isna(cost_type) or (isinstance(cost_type, str) and not cost_type.strip()):
        raise ValueError("Inflation Type is missing.")

    base_dollar_year = int(base_dollar_year)
    escalation_year  = int(escalation_year)
    
    df = pd.read_excel(file_path, sheet_name="Inflation Adjustment")
    # print("Shape:", df.shape)
    # print("First 5 rows raw:")
    # print(df.head(5))
    # print("\nYear column raw values:")
    # print(df['Year'].tolist())

    df = df.dropna(subset=['Year'])
    df['Year'] = df['Year'].astype(int)

    valid_cost_types = [column for column in df.columns if column != 'Year']
    if not isinstance(cost_type, str) or cost_type not in valid_cost_types:
        raise ValueError(
            f"Inflation Type '{cost_type}' is invalid. "
            f"Expected exactly one of: {', '.join(valid_cost_types)}."
        )

    # --- DEBUG: remove after issue is resolved ---
    # print("Year column dtype:", df['Year'].dtype)
    # print("Year values:", df['Year'].values)
    # print("Looking for:", base_dollar_year, type(base_dollar_year))

    if base_dollar_year not in df['Year'].values:
        raise ValueError(
            f"Dollar Year {base_dollar_year} is not available in the Inflation Adjustment sheet."
        )

    if escalation_year not in df['Year'].values:
        raise ValueError(
            f"Escalation Year {escalation_year} is not available in the Inflation Adjustment sheet."
        )

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

    # Read the Excel file into a Pandas DataFrame
    df = pd.read_excel(file_name, sheet_name=sheet_name)

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
            try:
                multiplier = calculate_inflation_multiplier(
                    file_name,
                    row['Dollar Year'],
                    row['Type'],
                    escalation_year
                )
            except ValueError as error:
                raise ValueError(
                    f"{sheet_name} account {row['Account']}: {error}"
                ) from error
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
    df_extra_params = pd.read_excel(file_name, sheet_name="Economics Parameters")
    extra_economic_parameters = dict(zip(df_extra_params["Parameter"], df_extra_params["Value"]))

    for parameter, value in extra_economic_parameters.items():
        params[parameter] = value

    return df
