# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import os
import copy
import math
from contextlib import contextmanager
import pandas as pd
import numpy as np
import csv
from cost.cost_escalation import escalate_cost_database
from cost.code_of_account_processing import remove_irrelevant_account, get_estimated_cost_column, find_children_accounts, create_cost_dictionary
from cost.cost_scaling import (scale_cost, scale_redundant_BOP_and_primary_loop,
                               scale_campus_cost, scale_central_facility_cost)
from cost.non_direct_cost import (validate_tax_credit_params, calculate_accounts_31_32_75_82_cost,
                                   calculate_decommissioning_cost, calculate_high_level_capital_costs,
                                   calculate_TCI, energy_cost_levelized,
                                   calculate_accounts_31_32_75_central_facility_cost,
                                   calculate_high_level_capital_costs_central_facility, calculate_TCI_central,
                                   calculate_servicing_campus_derived_costs,
                                   calculate_servicing_campus_capital_costs,
                                   calculate_servicing_campus_TCI,
                                   energy_cost_levelized_servicing_campus,
                                   calculate_manufacturing_campus_derived_costs,
                                   calculate_manufacturing_campus_ratio_costs,
                                   calculate_manufacturing_campus_capital_costs,
                                   calculate_manufacturing_campus_TCI,
                                   energy_cost_levelized_manufacturing_campus)
from cost.params_registry import PARAMS_REGISTRY, GROUP_ORDER
from reactor_engineering_evaluation.operation import reactor_operation
from cost.cost_drivers import cost_drivers_estimate


FLEET_MODE_REACTOR_EXCLUDED_ACCOUNTS = (
    213.2,   # Control Building
    214.11,  # Refueling Building
    214.12,  # Spent Fuel Building
    221.211, # Reactivity Control System Fabrication
    221.212, # Reactivity Control System Installation
    712,     # Remote Monitoring Technicians
    721,     # Coolant
    83,      # Spent Fuel Management
)


@contextmanager
def _numpy_sample_seed(seed):
    """Temporarily seed NumPy so one Monte Carlo realization is reproducible."""
    if seed is None:
        yield
        return
    previous_state = np.random.get_state()
    np.random.seed(seed)
    try:
        yield
    finally:
        np.random.set_state(previous_state)


def remove_fleet_mode_reactor_accounts(df):
    """Remove reactor accounts supplied by the shared Fleet Mode campuses.

    Parent accounts are removed together with all of their descendants.  The
    remaining reactor hierarchy is subsequently re-aggregated, so reactor OCC,
    financing, TCI, annual cost, and LCOE all use the Fleet Mode scope.
    """
    indices_to_drop = set()
    account_numbers = pd.to_numeric(df['Account'], errors='coerce')
    levels = df['Level'].to_numpy()

    for excluded_account in FLEET_MODE_REACTOR_EXCLUDED_ACCOUNTS:
        matching_positions = np.flatnonzero(
            np.isclose(account_numbers.to_numpy(dtype=float), excluded_account, equal_nan=False)
        )
        for position in matching_positions:
            parent_level = levels[position]
            indices_to_drop.add(df.index[position])
            for child_position in range(position + 1, len(df)):
                if levels[child_position] <= parent_level:
                    break
                indices_to_drop.add(df.index[child_position])

    return df.drop(index=list(indices_to_drop)).copy()



def calculate_high_level_accounts_cost(df, target_level, option, FOAK_or_NOAK):
    cost_column = get_estimated_cost_column(df, FOAK_or_NOAK)

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

    # Replace per-row iterrows with a vectorized eligibility mask. Only
    # the rows that actually need updating (typically a handful) enter
    # the Python loop — and even there the children sum uses pandas
    # vectorized .sum() instead of a per-child Python add.
    accounts_str = df['Account'].astype(str)
    mask = (
        accounts_str.str.startswith(valid_prefixes)
        & (df['Level'] == target_level)
        & df[cost_column].isna()
        & df['Children Accounts'].notna()
    )

    if not mask.any():
        return df

    children_col = df['Children Accounts']
    for index in df.index[mask]:
        children_idxs = [int(x) for x in children_col.at[index].split(',')]
        df.at[index, cost_column] = df.loc[children_idxs, cost_column].sum()

    return df


def update_high_level_costs(scaled_cost, option, sample):
    df_with_children_accounts = find_children_accounts(scaled_cost)
    no_subaccounts_list = []

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

    # Hoist column-name lookups out of the inner row loop — they were
    # being called 4× per row × ~100 rows × 5 levels = ~2000 redundant
    # column-name searches per update_high_level_costs call.
    foak_col = get_estimated_cost_column(df_with_children_accounts, 'F')
    noak_col = get_estimated_cost_column(df_with_children_accounts, 'N')

    for level in range(4, -1, -1):
        df_updated = calculate_high_level_accounts_cost(df_with_children_accounts, level, option, 'F')
        df_updated_2 = calculate_high_level_accounts_cost(df_updated, level, option, 'N')

        # Vectorized replacement for the per-row iterrows loop. Build masks
        # for FOAK and NOAK separately, then use bulk df.loc[mask, col] = 0
        # writes. Both masks read the same pre-mutation state so order is
        # irrelevant.
        accounts_str = df_updated_2['Account'].astype(str)
        common_mask = (
            accounts_str.str.startswith(valid_prefixes)
            & (df_updated_2['Level'] == level)
            & df_updated_2['Children Accounts'].isna()
        )
        foak_mask = common_mask & df_updated_2[foak_col].isna()
        noak_mask = common_mask & df_updated_2[noak_col].isna()

        if foak_mask.any():
            df_updated_2.loc[foak_mask, foak_col] = 0
            no_subaccounts_list.extend(df_updated_2.loc[foak_mask, 'Account'].tolist())
        if noak_mask.any():
            df_updated_2.loc[noak_mask, noak_col] = 0
            no_subaccounts_list.extend(df_updated_2.loc[noak_mask, 'Account'].tolist())
    
    if sample == 0:
        if no_subaccounts_list:
            print(f"Warning: The following accounts do not have any subaccounts: {', '.join(map(str, set(no_subaccounts_list))) }")
    return df_updated_2


def save_params_to_excel_file(excel_file, params):
    """
    Saves the params dictionary to the 'Parameters' sheet of the output Excel file.
    Parameters are organized into labeled groups, sorted alphabetically within each group,
    with units, descriptions, and source (User Input vs Calculated) for each parameter.
    Array parameters are summarized (BOL, EOL, min, max) rather than shown as raw lists.
    Parameters not found in the registry are placed in an 'Uncategorized' group with a warning.
    """

    def format_value(val):
        """
        Format a single scalar value for display.
        Converts numpy scalar types to native Python types to prevent
        Excel file corruption when openpyxl serializes the values.
        """
        # Handle complex types that openpyxl can't serialize
        if isinstance(val, dict):
            return str(val)
        if isinstance(val, np.ndarray):
            return str(val.tolist())
        # Handle numpy scalars first (before float check, since np.float64 is a subclass of float)
        if isinstance(val, np.floating):
            if np.isnan(val):
                return 'N/A'
            return float(val)
        if isinstance(val, np.integer):
            return int(val)
        if isinstance(val, np.bool_):
            return str(bool(val)).upper()
        # Handle native Python types
        if isinstance(val, float) and np.isnan(val):
            return 'N/A'
        if isinstance(val, bool):
            return str(val).upper()
        return val

    def handle_array(name, val, mode, units, description, source):
        """
        Expand an array parameter into multiple display rows based on mode:
          'summary' → BOL, EOL, min, max
          'steps'   → first step, last step, number of steps
          'as_is'   → single row with the list as a string
        Returns a list of (display_name, value, units, description, source) tuples.
        """
        rows = []
        if not isinstance(val, (list, tuple)) or len(val) == 0:
            rows.append((name, format_value(val), units, description, source))
            return rows

        if mode == 'summary':
            rows.append((f'{name} (BOL)',   float(round(val[0], 6)),   units, f'{description} — beginning of life', source))
            rows.append((f'{name} (EOL)',   float(round(val[-1], 6)),  units, f'{description} — end of life',       source))
            rows.append((f'{name} (min)',   float(round(min(val), 6)), units, f'{description} — minimum value',     source))
            rows.append((f'{name} (max)',   float(round(max(val), 6)), units, f'{description} — maximum value',     source))
        elif mode == 'steps':
            rows.append((f'{name} (first)', format_value(val[0]),  units, f'{description} — first step',     source))
            rows.append((f'{name} (last)',  format_value(val[-1]), units, f'{description} — last step',      source))
            rows.append((f'{name} (count)', len(val),              '',    f'{description} — number of steps', source))
        elif mode == 'as_is':
            rows.append((name, str(val), units, description, source))
        else:
            rows.append((name, str(val), units, description, source))
        return rows

    # ---------------------------------------------------------------
    # Build grouped rows from params using the registry
    # ---------------------------------------------------------------
    groups = {g: [] for g in GROUP_ORDER}
    params_dict = dict(params)

    for param_name, value in sorted(params_dict.items()):  # alphabetical within each group
        entry = PARAMS_REGISTRY.get(param_name)

        if entry is None:
            # Not in registry — place in Uncategorized with a warning marker
            if isinstance(value, (list, tuple)) and len(value) > 10:
                display_value = f'[list of {len(value)} items — see input file]'
            else:
                display_value = format_value(value)
            groups['Uncategorized'].append((
                param_name,
                display_value,
                '',
                '--- Not in params registry. Please add to cost/params_registry.py ---',
                'Unknown'
            ))
            continue

        # Skip hidden parameters
        if entry.get('hidden', False):
            continue

        # Tax Rate is only relevant when PTC is used (needed for gross-up calculation)
        # Skip it if PTC is not defined in params
        if param_name == 'Tax Rate' and 'PTC credit value' not in params_dict:
            continue

        units       = entry.get('units', '')
        description = entry.get('description', '')
        source      = entry.get('source', '')
        array_mode  = entry.get('array_mode', None)
        group       = entry.get('group', 'Uncategorized')

        if group not in groups:
            group = 'Uncategorized'

        if array_mode is not None and isinstance(value, (list, tuple)):
            rows = handle_array(param_name, value, array_mode, units, description, source)
            groups[group].extend(rows)
        else:
            groups[group].append((param_name, format_value(value), units, description, source))

    # ---------------------------------------------------------------
    # Build the final list of rows with group headers and separators
    # ---------------------------------------------------------------
    all_rows = []
    columns = ['Group', 'Parameter', 'Value', 'Units', 'Description', 'Source']

    for group_name in GROUP_ORDER:
        rows = groups.get(group_name, [])
        if not rows:
            continue  # skip empty groups

        # Group header row
        all_rows.append([f'--- {group_name.upper()} ---', '', '', '', '', ''])

        for (pname, pval, punits, pdesc, psource) in rows:
            all_rows.append([group_name, pname, pval, punits, pdesc, psource])

        # Blank separator row between groups
        all_rows.append(['', '', '', '', '', ''])

    # ---------------------------------------------------------------
    # Write to the Parameters sheet using the existing ExcelWriter
    # ---------------------------------------------------------------
    df = pd.DataFrame(all_rows, columns=columns)
    df.to_excel(excel_file, sheet_name='Parameters', index=False)

    total_params = sum(len(rows) for rows in groups.values())
    active_groups = sum(1 for g in GROUP_ORDER if groups.get(g))
    print(f"\n\nParameters saved — {total_params} entries across {active_groups} groups.\n\n")



def transform_dataframe(df):
    numerical_columns = df.select_dtypes(include=[np.number]).columns
    df = df.loc[~(df[numerical_columns] == 0).all(axis=1)]
    for col in numerical_columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        df[col] = df[col].fillna('-')
        df[col] = df[col].apply(lambda x: int(x) if x != '-' else x)
    return df


def format_campus_output(df, params, annual_cost_period):
    """Prepare a shared-campus table for its Excel output sheet.

    The estimator retains floating-point precision internally.  This formatter
    adds per-account contributions to fleet LCOE, gives the cost columns
    escalation-year-specific names, and converts displayed numerical results to
    integers without changing the values used by later calculations.
    """
    output = df.copy()
    escalation_year = params['Escalation Year']
    cost_column = f'Estimated Cost ($ {escalation_year})'
    cost_std_column = f'Estimated Cost std ($ {escalation_year})'
    lcoe_column = 'LCOE Contribution ($/MWh)'

    lifetime = int(params['Levelization Period'])
    annual_cost_period = int(annual_cost_period)
    discount_rate = params['Discount Rate']
    discount_sum = sum(
        (1 + discount_rate) ** -year
        for year in range(1, lifetime + 1)
    )
    annual_cost_discount_sum = sum(
        (1 + discount_rate) ** -year
        for year in range(1, annual_cost_period + 1)
    )
    fleet_annual_generation = (
        params['Fleet'] * params['Annual Electricity Production']
    )
    discounted_generation = fleet_annual_generation * discount_sum
    if discounted_generation <= 0:
        raise ValueError("Fleet discounted electricity generation must be greater than zero.")

    account_prefix = output['Account'].astype(str).str[0]
    capital_mask = account_prefix.isin(['1', '2', '3', '4', '5', '6'])
    annual_mask = account_prefix.isin(['7', '8'])

    output[lcoe_column] = np.nan
    output.loc[capital_mask, lcoe_column] = (
        output.loc[capital_mask, 'Value'] / discounted_generation
    )
    output.loc[annual_mask, lcoe_column] = (
        output.loc[annual_mask, 'Value'] * annual_cost_discount_sum / discounted_generation
    )

    output.rename(
        columns={
            'Value': cost_column,
            'Value std': cost_std_column,
        },
        inplace=True,
    )

    # Use pandas' nullable integer type so blank LCOE cells remain blank in
    # Excel while every displayed numerical result has no fractional part.
    for column in (cost_column, cost_std_column, lcoe_column):
        output[column] = output[column].apply(
            lambda value: int(value) if pd.notna(value) else pd.NA
        ).astype('Int64')

    return output


def format_servicing_campus_output(df, params):
    return format_campus_output(df, params, params['Levelization Period'])


def format_manufacturing_campus_output(df, params):
    return format_campus_output(df, params, params['Deployment Period'])


def learning_rate_multiplier(learning_rate, number_of_units):
    return pow(1-learning_rate, np.log2(min(100, number_of_units)))


def FOAK_to_NOAK(df, params):
    # Additional cost scaling based on an assumed learning rate.
    # Learning rate and cost multiplier are based on
    # DOI: 10.1080/00295450.2023.2206779.
    # Cost multiplier is capped at the 100th unit for any component.
    if 'NOAK Unit Number' not in params.keys():
        # Use the default value if 'NOAK Unit Number' is not specified.
        params['NOAK Unit Number'] = 10
        # Defaults to approximately the 10th unit, with 20 (2×10) units assumed for onsite learning.
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
    
    df['Multiplier'] = df['FOAK to NOAK Multiplier Type'].apply(get_multiplier)
    foak_col = get_estimated_cost_column(df, 'F')
    noak_column = foak_col.replace("FOAK", "NOAK")
    df[noak_column] = df['Multiplier'] * df[foak_col]
    return df


def reorder_dataframe(df):
    first_columns = ['Account', 'Account Title']
    other_columns = [col for col in df.columns if col not in first_columns]
    new_column_order = first_columns + other_columns
    df = df[new_column_order]
    return df


def bottom_up_cost_estimate(
    cost_database_filename,
    params,
    return_samples=False,
    sample_seeds=None,
):
    # Validate tax credit params early — before any simulation or cost calculation runs.
    # This catches cases where a user accidentally defines both ITC and PTC,
    # which are mutually exclusive under the IRA.
    validate_tax_credit_params(params)

    escalated_cost = escalate_cost_database(cost_database_filename, params['Escalation Year'], params)
    escalated_cost_cleaned = remove_irrelevant_account(escalated_cost, params)
    if params.get('Fleet Mode', False):
        escalated_cost_cleaned = remove_fleet_mode_reactor_accounts(
            escalated_cost_cleaned
        )
    reactor_operation(params)

    COA_list = []
    for i in range(params['Number of Samples']):
        if (i + 1) % 100 == 0:
            print(f"\n\nSample # {i+1}")

        seed = sample_seeds[i] if sample_seeds is not None else None
        with _numpy_sample_seed(seed):
            scaled_cost = scale_cost(escalated_cost_cleaned, params)
            scaled_cost = scale_redundant_BOP_and_primary_loop(scaled_cost, params)
            NOAK_COA = FOAK_to_NOAK(scaled_cost, params)

            updated_cost = update_high_level_costs(scaled_cost, 'base', i)
            updated_cost_with_indirect_cost = calculate_accounts_31_32_75_82_cost(updated_cost, params)
            cost_with_decommissioning = calculate_decommissioning_cost(updated_cost_with_indirect_cost, params)
            updated_accounts_10_40 = update_high_level_costs(cost_with_decommissioning, 'other', i)
            high_Level_capital_cost = calculate_high_level_capital_costs(updated_accounts_10_40, params)

            updated_accounts_10_60 = update_high_level_costs(high_Level_capital_cost, 'finance', i)
            TCI = calculate_TCI(updated_accounts_10_60, params)
            updated_accounts_70_80 = update_high_level_costs(TCI, 'annual', i)
            Final_COA = energy_cost_levelized(params, updated_accounts_70_80)
            FOAK_column = get_estimated_cost_column(Final_COA, 'F')
            NOAK_column = get_estimated_cost_column(Final_COA, 'N')
            Final_COA = Final_COA[['Account', 'Account Title', FOAK_column, NOAK_column]]

        COA_list.append(Final_COA)

    concatenated_df = pd.concat(COA_list)
    numeric_columns = concatenated_df.select_dtypes(include='number').columns
    mean_df = concatenated_df[numeric_columns].groupby(concatenated_df.index).mean()

    if params["Number of Samples"] > 1:
        std_df = concatenated_df[numeric_columns].groupby(concatenated_df.index).std()
    else:
        std_df = concatenated_df[numeric_columns].groupby(concatenated_df.index).std(ddof=0)

    mean_df[FOAK_column.replace('Cost', 'Cost std')] = std_df[FOAK_column]
    mean_df[NOAK_column.replace('Cost', 'Cost std')] = std_df[NOAK_column]

    non_numeric_columns = concatenated_df.select_dtypes(exclude='number').groupby(concatenated_df.index).first()
    result_df = mean_df.join(non_numeric_columns)
    reordered_df = reorder_dataframe(result_df)
    if return_samples:
        return reordered_df, COA_list
    return reordered_df


def bottom_up_cost_estimate_central(cost_database_filename, params):
    """
    Bottom-up cost estimate for central facility.
    Only runs if params['Estimate Central Facility'] is True.
    """
    get_central_facility_cost = params.get('Estimate Central Facility', False)

    if not get_central_facility_cost:
        return None

    escalated_central = escalate_cost_database(cost_database_filename,
                                                params['Escalation Year'],
                                                params,
                                                sheet_name='Central Facility Database')
    escalated_central_cleaned = remove_irrelevant_account(escalated_central, params)

    COA_list = []
    for i in range(params['Number of Samples']):
        if (i + 1) % 100 == 0:
            print(f"\n\nSample # {i+1}")

        scaled_cost = scale_central_facility_cost(escalated_central_cleaned, params)
        NOAK_COA = FOAK_to_NOAK(scaled_cost, params)

        updated_cost = update_high_level_costs(scaled_cost, 'base', i)
        updated_cost_with_indirect_cost = calculate_accounts_31_32_75_central_facility_cost(updated_cost, params)
        cost_with_decommissioning = calculate_decommissioning_cost(updated_cost_with_indirect_cost, params)
        updated_accounts_10_40 = update_high_level_costs(cost_with_decommissioning, 'other', i)
        high_Level_capital_cost = calculate_high_level_capital_costs_central_facility(updated_accounts_10_40, params)

        updated_accounts_10_60 = update_high_level_costs(high_Level_capital_cost, 'finance', i)
        TCI = calculate_TCI_central(updated_accounts_10_60, params)
        updated_accounts_70_80 = update_high_level_costs(TCI, 'annual', i)

        FOAK_column = get_estimated_cost_column(updated_accounts_70_80, 'F')
        NOAK_column = get_estimated_cost_column(updated_accounts_70_80, 'N')
        Final_COA = updated_accounts_70_80[['Account', 'Account Title', FOAK_column, NOAK_column]]

        COA_list.append(Final_COA)

    concatenated_df = pd.concat(COA_list)
    numeric_columns = concatenated_df.select_dtypes(include='number').columns
    mean_df = concatenated_df[numeric_columns].groupby(concatenated_df.index).mean()

    if params["Number of Samples"] > 1:
        std_df = concatenated_df[numeric_columns].groupby(concatenated_df.index).std()
    else:
        std_df = concatenated_df[numeric_columns].groupby(concatenated_df.index).std(ddof=0)

    mean_df[FOAK_column.replace('Cost', 'Cost std')] = std_df[FOAK_column]
    mean_df[NOAK_column.replace('Cost', 'Cost std')] = std_df[NOAK_column]

    non_numeric_columns = concatenated_df.select_dtypes(exclude='number').groupby(concatenated_df.index).first()
    result_df = mean_df.join(non_numeric_columns)
    reordered_df = reorder_dataframe(result_df)
    return reordered_df


def bottom_up_cost_estimate_servicing_campus(
    cost_database_filename,
    params,
    return_samples=False,
    sample_seeds=None,
):
    """Estimate the servicing campus when Fleet Mode is enabled.

    The servicing database assumes no learning, so the public result contains a
    single Value column rather than FOAK and NOAK columns. Two identical internal
    columns are retained only to reuse the established account aggregation code.
    """
    if not params.get('Fleet Mode', False):
        return None

    if 'Annual Electricity Production' not in params:
        raise KeyError(
            "'Annual Electricity Production' is missing. The reactor operation "
            "calculation must run before the servicing-campus estimate."
        )

    escalated_cost = escalate_cost_database(
        cost_database_filename,
        params['Escalation Year'],
        params,
        sheet_name='Servicing Campus Database',
    )
    cleaned_cost = remove_irrelevant_account(escalated_cost, params)

    samples = []
    for sample in range(params['Number of Samples']):
        if (sample + 1) % 100 == 0:
            print(f"\n\nServicing campus sample #{sample + 1}")

        seed = sample_seeds[sample] if sample_seeds is not None else None
        with _numpy_sample_seed(seed):
            scaled_cost = scale_campus_cost(cleaned_cost, params)
            value_column = get_estimated_cost_column(scaled_cost, 'F')
            internal_column = value_column.replace('FOAK', 'NOAK')
            scaled_cost[internal_column] = scaled_cost[value_column]

            # The servicing database intentionally contains empty placeholder
            # accounts (for example 75 and 78), so suppress the generic missing-
            # children warning while still aggregating all populated hierarchies.
            aggregation_sample = sample + 1
            updated_cost = update_high_level_costs(scaled_cost, 'base', aggregation_sample)
            updated_cost = calculate_servicing_campus_derived_costs(updated_cost, params)
            updated_cost = update_high_level_costs(updated_cost, 'other', aggregation_sample)
            updated_cost = calculate_servicing_campus_capital_costs(updated_cost, params)
            updated_cost = update_high_level_costs(updated_cost, 'finance', aggregation_sample)
            updated_cost = calculate_servicing_campus_TCI(updated_cost, params)
            updated_cost = update_high_level_costs(updated_cost, 'annual', aggregation_sample)
            final_cost = energy_cost_levelized_servicing_campus(params, updated_cost)

            sample_result = final_cost[['Account', 'Account Title', value_column]].copy()
            sample_result.rename(columns={value_column: 'Value'}, inplace=True)
        samples.append(sample_result)

    concatenated = pd.concat(samples)
    mean_values = concatenated[['Value']].groupby(concatenated.index).mean()
    if params['Number of Samples'] > 1:
        std_values = concatenated[['Value']].groupby(concatenated.index).std()
    else:
        std_values = concatenated[['Value']].groupby(concatenated.index).std(ddof=0)

    mean_values['Value std'] = std_values['Value']
    labels = concatenated[['Account', 'Account Title']].groupby(concatenated.index).first()
    result = reorder_dataframe(mean_values.join(labels))
    nonzero_mask = (result[['Value', 'Value std']] != 0).any(axis=1)
    result = result.loc[nonzero_mask].reset_index(drop=True)
    if return_samples:
        return result, samples
    return result


def bottom_up_cost_estimate_manufacturing_campus(
    cost_database_filename,
    params,
    reactor_cost_table,
    return_samples=False,
    sample_seeds=None,
    reactor_cost_samples=None,
):
    """Estimate the manufacturing campus when Fleet Mode is enabled."""
    if not params.get('Fleet Mode', False):
        return None
    if 'Annual Electricity Production' not in params:
        raise KeyError(
            "'Annual Electricity Production' is missing. The reactor operation "
            "calculation must run before the manufacturing-campus estimate."
        )
    if reactor_cost_table is None:
        raise ValueError(
            "The reactor cost table is required to calculate manufacturing account 223."
        )

    reactor_value_column = get_estimated_cost_column(reactor_cost_table, 'F')
    mean_reactor_occ_rows = reactor_cost_table.loc[reactor_cost_table['Account'] == 'OCC']
    if mean_reactor_occ_rows.empty:
        raise KeyError("The reactor cost table is missing the OCC summary account.")
    mean_reactor_occ = mean_reactor_occ_rows.iloc[0][reactor_value_column]
    if reactor_cost_samples is not None and len(reactor_cost_samples) != params['Number of Samples']:
        raise ValueError(
            "reactor_cost_samples must contain one reactor result per Monte Carlo sample."
        )

    escalated_cost = escalate_cost_database(
        cost_database_filename,
        params['Escalation Year'],
        params,
        sheet_name='Manufacturing Campus Database',
    )
    cleaned_cost = remove_irrelevant_account(escalated_cost, params)

    samples = []
    for sample in range(params['Number of Samples']):
        if (sample + 1) % 100 == 0:
            print(f"\n\nManufacturing campus sample #{sample + 1}")

        if reactor_cost_samples is None:
            reactor_occ = mean_reactor_occ
        else:
            reactor_sample = reactor_cost_samples[sample]
            sample_reactor_column = get_estimated_cost_column(reactor_sample, 'F')
            reactor_occ = _result_value(
                reactor_sample, 'OCC', sample_reactor_column
            )

        seed = sample_seeds[sample] if sample_seeds is not None else None
        with _numpy_sample_seed(seed):
            scaled_cost = scale_campus_cost(cleaned_cost, params, campus_type='manufacturing')
            value_column = get_estimated_cost_column(scaled_cost, 'F')
            internal_column = value_column.replace('FOAK', 'NOAK')
            scaled_cost[internal_column] = scaled_cost[value_column]

            scaled_cost = calculate_manufacturing_campus_derived_costs(
                scaled_cost, params, reactor_occ
            )
            aggregation_sample = sample + 1
            updated_cost = update_high_level_costs(scaled_cost, 'base', aggregation_sample)
            updated_cost = calculate_manufacturing_campus_ratio_costs(updated_cost, params)
            updated_cost = update_high_level_costs(updated_cost, 'other', aggregation_sample)
            updated_cost = calculate_manufacturing_campus_capital_costs(updated_cost, params)
            updated_cost = update_high_level_costs(updated_cost, 'finance', aggregation_sample)
            updated_cost = calculate_manufacturing_campus_TCI(updated_cost, params)
            updated_cost = update_high_level_costs(updated_cost, 'annual', aggregation_sample)
            final_cost = energy_cost_levelized_manufacturing_campus(params, updated_cost)

            sample_result = final_cost[['Account', 'Account Title', value_column]].copy()
            sample_result.rename(columns={value_column: 'Value'}, inplace=True)
        samples.append(sample_result)

    concatenated = pd.concat(samples)
    mean_values = concatenated[['Value']].groupby(concatenated.index).mean()
    if params['Number of Samples'] > 1:
        std_values = concatenated[['Value']].groupby(concatenated.index).std()
    else:
        std_values = concatenated[['Value']].groupby(concatenated.index).std(ddof=0)

    mean_values['Value std'] = std_values['Value']
    labels = concatenated[['Account', 'Account Title']].groupby(concatenated.index).first()
    result = reorder_dataframe(mean_values.join(labels))
    nonzero_mask = (result[['Value', 'Value std']] != 0).any(axis=1)
    result = result.loc[nonzero_mask].reset_index(drop=True)
    if return_samples:
        return result, samples
    return result


def create_campus_cost_dictionary(df, campus_name):
    """Extract shared-campus summary and high-level accounts for a CSV row."""
    accounts = {
        10: 'Capitalized Pre-Construction Costs',
        20: 'Capitalized Direct Costs',
        30: 'Capitalized Indirect Services Cost',
        40: 'Capitalized Training Costs',
        60: 'Capitalized Financial Costs',
        70: 'Annualized O&M Cost',
        'OCC': 'OCC',
        'OCC per reactor': 'OCC per reactor',
        'TCI': 'TCI',
        'TCI per reactor': 'TCI per reactor',
        'Annual Cost': 'Annual Cost',
        'Annual Cost per reactor': 'Annual Cost per reactor',
        'LCOE': 'LCOE',
    }
    tracked_costs = {}

    for account, label in accounts.items():
        matching_rows = df.loc[df['Account'] == account]
        if matching_rows.empty:
            raise KeyError(
                f"{campus_name} result is missing required tracked account '{account}'."
            )
        row = matching_rows.iloc[0]
        tracked_costs[f'{campus_name} {label}'] = row['Value']
        tracked_costs[f'{campus_name} {label} std'] = row['Value std']

    return tracked_costs


def create_servicing_campus_cost_dictionary(df):
    return create_campus_cost_dictionary(df, 'Servicing Campus')


def create_manufacturing_campus_cost_dictionary(df):
    return create_campus_cost_dictionary(df, 'Manufacturing Campus')


def _result_value(df, account, column='Value'):
    rows = df.loc[df['Account'] == account]
    if rows.empty:
        raise KeyError(f"Cost result is missing required account '{account}'.")
    return float(rows.iloc[0][column])


def _servicing_operating_state_params(params, operating_reactors, service_events):
    """Return campus parameters for one fleet operating year.

    Capital sizing remains at the final-fleet design stored in ``params``. Only
    variables that drive annual staffing, consumables, transportation, and
    servicing activity are updated for the current operating fleet and the
    cohort-based number of service events.
    """
    state = copy.deepcopy(params)
    state['Number of Samples'] = 1
    state['Generating Sites Count'] = operating_reactors
    state['Servicing Rate'] = service_events

    service_scale = np.rint(service_events / 3)
    state['SER Switchyard Average Power'] = 6 * service_scale ** 0.698970
    state['Servicing Hot Cell Count'] = np.ceil(
        service_events / state['Servicing Hot Cell Annual Rate']
    )
    state['He Gas Replenishment'] = (
        state['Servicing Hot Cell Count'] * state['He Gas Replenishment Per Hot Cell']
        + state['Radioactive Waste Processing Hot Cell Count']
        * state['He Gas Replenishment Per Hot Cell']
        + state['CoolantInventoryRPV_Mass'] * service_events
    )
    state['SER Number of Operators Per Shift'] = np.ceil(
        5.625 * service_scale ** 0.426
    )
    state['SER Engineering Headcount'] = np.ceil(
        20.0 * service_scale ** 0.301
    )

    reactor_trip_time = (
        state['Roundtrip Time Reactor Transport']
        + state['Dwell Time Reactor Transport GenSite']
        + state['Dwell Time Reactor Transport Serv']
    )
    supply_trip_time = (
        state['Roundtrip Time']
        + state['Dwell Time GenSite']
        + state['Dwell Time Serv']
    )
    state['Reactor Transport Vehicle Count'] = np.ceil(
        service_events * reactor_trip_time + 1
    )
    state['Helium Transport Truck Count'] = np.ceil(
        operating_reactors
        * state['Annual Coolant Supply Frequency']
        * supply_trip_time
        * 1.05
    )
    state['Water Tanker Truck Count'] = np.ceil(
        operating_reactors
        * state['Water Supply Frequency']
        * supply_trip_time
        * 1.05
    )
    state['Maintenance Truck Count'] = np.ceil(
        operating_reactors
        * state['Maintenance Visit Frequency']
        * supply_trip_time
        * 1.05
    )

    state['Reactor Transport Cask Count'] = round(
        service_events * (4 / 12) * 1.5, -1
    )
    state['Annual Used Fuel Cask Consumption'] = service_events
    state['Annual Reactor Cask Replacement'] = np.ceil(
        0.05 * state['Reactor Transport Cask Count']
    )
    state['Annual Radwaste Cask Consumption'] = 0.5 * service_events
    return state


def _fleet_cohort_schedule(params):
    """Build annual reactor manufacture, operation, retirement, and service counts."""
    fleet_size = int(params['Fleet'])
    production_rate = int(params['Production Rate'])
    deployment_years = int(params['Deployment Period'])
    operating_lifetime = int(params['Levelization Period'])
    cycle_length = float(params['Cycle Length'])
    commissioning_lag = int(params.get('Fleet Reactor Commissioning Lag', 1))

    if min(fleet_size, production_rate, deployment_years, operating_lifetime) <= 0:
        raise ValueError(
            "Fleet, Production Rate, Deployment Period, and Levelization Period "
            "must all be greater than zero."
        )
    if cycle_length <= 0:
        raise ValueError("'Cycle Length' must be greater than zero.")
    if commissioning_lag < 0:
        raise ValueError("'Fleet Reactor Commissioning Lag' cannot be negative.")

    cohorts = []
    remaining = fleet_size
    for manufacture_year in range(1, deployment_years + 1):
        manufactured = min(production_rate, remaining)
        if manufactured <= 0:
            break
        cohorts.append({
            'manufacture_year': manufacture_year,
            'operation_start': manufacture_year + commissioning_lag,
            'reactors': manufactured,
        })
        remaining -= manufactured
    if remaining:
        raise ValueError(
            "Production Rate multiplied by Deployment Period is insufficient "
            f"to manufacture the requested fleet; {remaining} reactors remain."
        )

    horizon = max(
        cohort['operation_start'] + operating_lifetime - 1
        for cohort in cohorts
    )
    service_events_by_year = {year: 0 for year in range(horizon + 1)}
    for cohort in cohorts:
        service_number = 1
        retirement_time = cohort['operation_start'] + operating_lifetime
        while True:
            event_time = cohort['operation_start'] + service_number * cycle_length
            if event_time >= retirement_time:
                break
            event_year = int(math.ceil(event_time - 1e-12))
            service_events_by_year[event_year] += cohort['reactors']
            service_number += 1

    rows = []
    for year in range(horizon + 1):
        manufactured = sum(
            cohort['reactors']
            for cohort in cohorts
            if cohort['manufacture_year'] == year
        )
        operating = sum(
            cohort['reactors']
            for cohort in cohorts
            if cohort['operation_start'] <= year
            < cohort['operation_start'] + operating_lifetime
        )
        retiring = sum(
            cohort['reactors']
            for cohort in cohorts
            if cohort['operation_start'] + operating_lifetime - 1 == year
        )
        rows.append({
            'Year': year,
            'Reactors Manufactured': manufactured,
            'Operating Reactors': operating,
            'Reactors Retiring at Year End': retiring,
            'Servicing Events': service_events_by_year[year],
        })
    return pd.DataFrame(rows)


def build_fleet_annual_cash_flow(
    cost_database_filename,
    params,
    reactor_cost_table,
    manufacturing_cost_table,
    servicing_cost_table,
    servicing_annual_costs=None,
):
    """Create the realistic cohort-based annual fleet cash-flow table."""
    schedule = _fleet_cohort_schedule(params)
    reactor_noak_column = get_estimated_cost_column(reactor_cost_table, 'N')
    reactor_tci = _result_value(reactor_cost_table, 'TCI', reactor_noak_column)
    reactor_annual = (
        _result_value(reactor_cost_table, 70, reactor_noak_column)
        + _result_value(reactor_cost_table, 80, reactor_noak_column)
    )
    manufacturing_tci = _result_value(manufacturing_cost_table, 'TCI')
    manufacturing_annual = _result_value(manufacturing_cost_table, 'Annual Cost')
    servicing_tci = _result_value(servicing_cost_table, 'TCI')
    reported_full_fleet_servicing_annual = _result_value(
        servicing_cost_table, 'Annual Cost'
    )

    # The annual servicing profile is evaluated deterministically for each
    # unique operating state, then normalized to the Monte Carlo mean reported
    # by the standalone servicing-campus estimate at the final-fleet state.
    annual_servicing_cache = {}

    def deterministic_servicing_annual(operating_reactors, service_events):
        if operating_reactors <= 0:
            return 0.0
        cache_key = (int(operating_reactors), float(service_events))
        if cache_key not in annual_servicing_cache:
            state_params = _servicing_operating_state_params(
                params, operating_reactors, service_events
            )
            state_table = bottom_up_cost_estimate_servicing_campus(
                cost_database_filename, state_params
            )
            annual_servicing_cache[cache_key] = _result_value(
                state_table, 'Annual Cost'
            )
        return annual_servicing_cache[cache_key]

    if servicing_annual_costs is None:
        deterministic_full_fleet_annual = deterministic_servicing_annual(
            int(params['Fleet']), float(params['Servicing Rate'])
        )
        servicing_normalization = (
            reported_full_fleet_servicing_annual / deterministic_full_fleet_annual
            if deterministic_full_fleet_annual > 0
            else 1.0
        )

    schedule['Fleet Electricity Generation (MWh)'] = (
        schedule['Operating Reactors'] * params['Annual Electricity Production']
    )
    schedule['Reactor TCI Deployed'] = (
        schedule['Reactors Manufactured'] * reactor_tci
    )
    schedule['Reactor Annual Cost'] = (
        schedule['Operating Reactors'] * reactor_annual
    )
    schedule['Manufacturing Campus TCI'] = 0.0
    schedule.loc[schedule['Year'] == 0, 'Manufacturing Campus TCI'] = (
        manufacturing_tci
    )
    schedule['Manufacturing Campus Annual Cost'] = (
        manufacturing_annual
        * schedule['Reactors Manufactured']
        / params['Production Rate']
    )
    schedule['Servicing Campus TCI'] = 0.0
    schedule.loc[schedule['Year'] == 0, 'Servicing Campus TCI'] = servicing_tci
    if servicing_annual_costs is None:
        schedule['Servicing Campus Annual Cost'] = schedule.apply(
            lambda row: deterministic_servicing_annual(
                row['Operating Reactors'], row['Servicing Events']
            ) * servicing_normalization,
            axis=1,
        )
    else:
        if len(servicing_annual_costs) != len(schedule):
            raise ValueError(
                "servicing_annual_costs must contain one value per cash-flow year."
            )
        schedule['Servicing Campus Annual Cost'] = np.asarray(
            servicing_annual_costs, dtype=float
        )

    component_cost_columns = [
        'Reactor TCI Deployed',
        'Reactor Annual Cost',
        'Manufacturing Campus TCI',
        'Manufacturing Campus Annual Cost',
        'Servicing Campus TCI',
        'Servicing Campus Annual Cost',
    ]
    schedule['Total Fleet Cost'] = schedule[component_cost_columns].sum(axis=1)
    schedule['Discount Factor'] = (
        1 + params['Discount Rate']
    ) ** -schedule['Year']
    schedule['Discounted Fleet Cost'] = (
        schedule['Total Fleet Cost'] * schedule['Discount Factor']
    )
    schedule['Discounted Fleet Generation (MWh)'] = (
        schedule['Fleet Electricity Generation (MWh)']
        * schedule['Discount Factor']
    )
    return schedule


def build_servicing_annual_cost_samples(
    cost_database_filename,
    params,
    sample_seeds,
):
    """Calculate each sample's servicing cost in every fleet operating year.

    The same seed is reused for all operating states within one sample. Thus a
    sampled wage, equipment price, or material price remains the same over that
    sample's complete project life while fleet activity changes by year.
    """
    schedule = _fleet_cohort_schedule(params)
    sample_count = len(sample_seeds)
    annual_cost_samples = np.zeros((sample_count, len(schedule)))
    state_cache = {}

    for row_index, row in schedule.iterrows():
        operating_reactors = int(row['Operating Reactors'])
        service_events = float(row['Servicing Events'])
        if operating_reactors <= 0:
            continue
        state_key = (operating_reactors, service_events)
        if state_key not in state_cache:
            state_params = _servicing_operating_state_params(
                params, operating_reactors, service_events
            )
            state_params['Number of Samples'] = sample_count
            _, state_samples = bottom_up_cost_estimate_servicing_campus(
                cost_database_filename,
                state_params,
                return_samples=True,
                sample_seeds=sample_seeds,
            )
            state_cache[state_key] = np.asarray([
                _result_value(sample_table, 'Annual Cost')
                for sample_table in state_samples
            ])
        annual_cost_samples[:, row_index] = state_cache[state_key]

    return annual_cost_samples


def build_fleet_total_cost_estimate(
    params,
    reactor_cost_table,
    manufacturing_cost_table,
    servicing_cost_table,
    annual_cash_flow,
):
    """Build the auditable high-level total cost table for Fleet Mode."""
    fleet_size = int(params['Fleet'])
    reactor_column = get_estimated_cost_column(reactor_cost_table, 'N')

    def reactor_value(account):
        return _result_value(reactor_cost_table, account, reactor_column)

    def campus_value(table, account, default=0.0):
        rows = table.loc[table['Account'] == account]
        return default if rows.empty else float(rows.iloc[0]['Value'])

    titles = {
        10: 'Capitalized Pre-Construction Costs',
        20: 'Capitalized Direct Costs',
        30: 'Capitalized Indirect Services Cost',
        40: 'Capitalized Training Costs',
        60: 'Capitalized Financial Costs',
        70: 'Annualized O&M Cost (Full-Fleet Deployment-Year Basis)',
        80: 'Annualized Fuel Cost',
    }
    rows = []

    def add_row(account, title, units, reactor, manufacturing, servicing):
        rows.append({
            'Account': account,
            'Account Title': title,
            'Units': units,
            'Reactor Fleet': reactor,
            'Manufacturing Campus': manufacturing,
            'Servicing Campus': servicing,
            'Total Fleet': reactor + manufacturing + servicing,
        })

    for account, title in titles.items():
        add_row(
            account,
            title,
            f"$ {params['Escalation Year']}",
            fleet_size * reactor_value(account),
            campus_value(manufacturing_cost_table, account),
            campus_value(servicing_cost_table, account),
        )

    reactor_occ = fleet_size * reactor_value('OCC')
    manufacturing_occ = campus_value(manufacturing_cost_table, 'OCC')
    servicing_occ = campus_value(servicing_cost_table, 'OCC')
    add_row('OCC', 'Overnight Capital Cost', f"$ {params['Escalation Year']}",
            reactor_occ, manufacturing_occ, servicing_occ)
    add_row('OCC per reactor', 'Overnight Capital Cost per Reactor',
            f"$ {params['Escalation Year']}/reactor",
            reactor_occ / fleet_size, manufacturing_occ / fleet_size,
            servicing_occ / fleet_size)

    reactor_tci = fleet_size * reactor_value('TCI')
    manufacturing_tci = campus_value(manufacturing_cost_table, 'TCI')
    servicing_tci = campus_value(servicing_cost_table, 'TCI')
    add_row('TCI', 'Total Capital Investment', f"$ {params['Escalation Year']}",
            reactor_tci, manufacturing_tci, servicing_tci)
    add_row('TCI per reactor', 'Total Capital Investment per Reactor',
            f"$ {params['Escalation Year']}/reactor",
            reactor_tci / fleet_size, manufacturing_tci / fleet_size,
            servicing_tci / fleet_size)

    full_reactor_annual = fleet_size * (reactor_value(70) + reactor_value(80))
    full_manufacturing_annual = campus_value(
        manufacturing_cost_table, 'Annual Cost'
    )
    full_servicing_annual = campus_value(servicing_cost_table, 'Annual Cost')
    add_row(
        'Annual Cost During Deployment',
        'Full-Fleet Annual Cost During Deployment',
        f"$ {params['Escalation Year']}/year",
        full_reactor_annual,
        full_manufacturing_annual,
        full_servicing_annual,
    )
    add_row(
        'Annual Cost per reactor During Deployment',
        'Full-Fleet Annual Cost per Reactor During Deployment',
        f"$ {params['Escalation Year']}/reactor-year",
        full_reactor_annual / fleet_size,
        full_manufacturing_annual / fleet_size,
        full_servicing_annual / fleet_size,
    )
    add_row(
        'Steady-State Annual Cost',
        'Full-Fleet Steady-State Annual Cost After Deployment',
        f"$ {params['Escalation Year']}/year",
        full_reactor_annual,
        0.0,
        full_servicing_annual,
    )
    add_row(
        'Steady-State Annual Cost per reactor',
        'Full-Fleet Steady-State Annual Cost per Reactor After Deployment',
        f"$ {params['Escalation Year']}/reactor-year",
        full_reactor_annual / fleet_size,
        0.0,
        full_servicing_annual / fleet_size,
    )

    discount = annual_cash_flow['Discount Factor']
    discounted_generation = annual_cash_flow[
        'Discounted Fleet Generation (MWh)'
    ].sum()
    if discounted_generation <= 0:
        raise ValueError("Discounted fleet electricity generation must be positive.")
    reactor_lcoe = (
        (annual_cash_flow['Reactor TCI Deployed']
         + annual_cash_flow['Reactor Annual Cost']) * discount
    ).sum() / discounted_generation
    manufacturing_lcoe = (
        (annual_cash_flow['Manufacturing Campus TCI']
         + annual_cash_flow['Manufacturing Campus Annual Cost']) * discount
    ).sum() / discounted_generation
    servicing_lcoe = (
        (annual_cash_flow['Servicing Campus TCI']
         + annual_cash_flow['Servicing Campus Annual Cost']) * discount
    ).sum() / discounted_generation
    add_row(
        'LCOE',
        'Total Fleet Levelized Cost of Electricity',
        '$/MWh',
        reactor_lcoe,
        manufacturing_lcoe,
        servicing_lcoe,
    )
    return pd.DataFrame(rows)


def _aggregate_fleet_total_samples(sample_tables):
    """Calculate fleet result means/stds only after complete sample aggregation."""
    if not sample_tables:
        raise ValueError("At least one fleet total sample is required.")
    result = sample_tables[0][['Account', 'Account Title', 'Units']].copy()
    value_columns = [
        'Reactor Fleet',
        'Manufacturing Campus',
        'Servicing Campus',
        'Total Fleet',
    ]
    ddof = 1 if len(sample_tables) > 1 else 0
    for column in value_columns:
        values = np.vstack([
            table[column].to_numpy(dtype=float) for table in sample_tables
        ])
        result[column] = values.mean(axis=0)
        result[f'{column} std'] = values.std(axis=0, ddof=ddof)

    ordered_columns = ['Account', 'Account Title', 'Units']
    for column in value_columns:
        ordered_columns.extend([column, f'{column} std'])
    return result[ordered_columns]


def _aggregate_fleet_cash_flow_samples(sample_tables):
    """Aggregate the annual fleet cash flows and retain annual uncertainty."""
    if not sample_tables:
        raise ValueError("At least one fleet cash-flow sample is required.")
    deterministic_columns = [
        'Year',
        'Reactors Manufactured',
        'Operating Reactors',
        'Reactors Retiring at Year End',
        'Servicing Events',
        'Fleet Electricity Generation (MWh)',
        'Discount Factor',
        'Discounted Fleet Generation (MWh)',
    ]
    cost_columns = [
        'Reactor TCI Deployed',
        'Reactor Annual Cost',
        'Manufacturing Campus TCI',
        'Manufacturing Campus Annual Cost',
        'Servicing Campus TCI',
        'Servicing Campus Annual Cost',
        'Total Fleet Cost',
        'Discounted Fleet Cost',
    ]
    result = sample_tables[0][deterministic_columns].copy()
    ddof = 1 if len(sample_tables) > 1 else 0
    for column in cost_columns:
        values = np.vstack([
            table[column].to_numpy(dtype=float) for table in sample_tables
        ])
        result[column] = values.mean(axis=0)
        result[f'{column} std'] = values.std(axis=0, ddof=ddof)

    ordered_columns = deterministic_columns[:6]
    for column in cost_columns[:-1]:
        ordered_columns.extend([column, f'{column} std'])
    ordered_columns.append('Discount Factor')
    ordered_columns.extend([
        'Discounted Fleet Cost',
        'Discounted Fleet Cost std',
        'Discounted Fleet Generation (MWh)',
    ])
    return result[ordered_columns]


def estimate_fleet_mode_costs(cost_database_filename, params):
    """Estimate paired reactor/campus samples and combine them at fleet level."""
    sample_count = int(params['Number of Samples'])
    reactor_seeds = [100_000 + sample for sample in range(sample_count)]
    servicing_seeds = [200_000 + sample for sample in range(sample_count)]
    manufacturing_seeds = [300_000 + sample for sample in range(sample_count)]

    reactor_result, reactor_samples = bottom_up_cost_estimate(
        cost_database_filename,
        params,
        return_samples=True,
        sample_seeds=reactor_seeds,
    )
    servicing_result, servicing_samples = bottom_up_cost_estimate_servicing_campus(
        cost_database_filename,
        params,
        return_samples=True,
        sample_seeds=servicing_seeds,
    )
    manufacturing_result, manufacturing_samples = (
        bottom_up_cost_estimate_manufacturing_campus(
            cost_database_filename,
            params,
            reactor_result,
            return_samples=True,
            sample_seeds=manufacturing_seeds,
            reactor_cost_samples=reactor_samples,
        )
    )

    servicing_annual_samples = build_servicing_annual_cost_samples(
        cost_database_filename,
        params,
        servicing_seeds,
    )
    cash_flow_samples = []
    total_cost_samples = []
    for sample in range(sample_count):
        cash_flow = build_fleet_annual_cash_flow(
            cost_database_filename,
            params,
            reactor_samples[sample],
            manufacturing_samples[sample],
            servicing_samples[sample],
            servicing_annual_costs=servicing_annual_samples[sample],
        )
        total_cost = build_fleet_total_cost_estimate(
            params,
            reactor_samples[sample],
            manufacturing_samples[sample],
            servicing_samples[sample],
            cash_flow,
        )
        cash_flow_samples.append(cash_flow)
        total_cost_samples.append(total_cost)

    return {
        'reactor': reactor_result,
        'servicing': servicing_result,
        'manufacturing': manufacturing_result,
        'cash_flow': _aggregate_fleet_cash_flow_samples(cash_flow_samples),
        'total': _aggregate_fleet_total_samples(total_cost_samples),
    }


def create_fleet_total_cost_dictionary(df):
    """Extract total-fleet value/std columns for one parametric-study row."""
    labels = {
        10: 'Capitalized Pre-Construction Costs',
        20: 'Capitalized Direct Costs',
        30: 'Capitalized Indirect Services Cost',
        40: 'Capitalized Training Costs',
        60: 'Capitalized Financial Costs',
        70: 'Annualized O&M Cost',
        80: 'Annualized Fuel Cost',
        'OCC': 'OCC',
        'OCC per reactor': 'OCC per reactor',
        'TCI': 'TCI',
        'TCI per reactor': 'TCI per reactor',
        'Annual Cost During Deployment': 'Annual Cost During Deployment',
        'Annual Cost per reactor During Deployment': (
            'Annual Cost per reactor During Deployment'
        ),
        'Steady-State Annual Cost': 'Steady-State Annual Cost',
        'Steady-State Annual Cost per reactor': (
            'Steady-State Annual Cost per reactor'
        ),
        'LCOE': 'LCOE',
    }
    tracked_costs = {}
    for account, label in labels.items():
        rows = df.loc[df['Account'] == account]
        if rows.empty:
            raise KeyError(
                f"Total fleet result is missing tracked account '{account}'."
            )
        row = rows.iloc[0]
        tracked_costs[f'Total Fleet {label}'] = row['Total Fleet']
        tracked_costs[f'Total Fleet {label} std'] = row['Total Fleet std']
    return tracked_costs


def parametric_studies(cost_database_filename, tracked_params_list):
    import inspect

    # Grab params and the calling script's path from the caller's frame automatically
    caller_frame = inspect.stack()[1][0]
    params = caller_frame.f_locals.get('params')
    if params is None:
        raise RuntimeError(
            "parametric_studies could not find 'params' in the calling scope. "
            "Make sure a variable named 'params' exists in the script that calls this function."
        )
    caller_file = caller_frame.f_globals.get('__file__', 'output')
    output_csv_filename = os.path.splitext(os.path.abspath(caller_file))[0] + '_output.csv'

    fleet_results = None
    if params.get('Fleet Mode', False):
        fleet_results = estimate_fleet_mode_costs(cost_database_filename, params)
        detailed_cost_table = fleet_results['reactor']
    else:
        detailed_cost_table = bottom_up_cost_estimate(cost_database_filename, params)
    tracked_costs = create_cost_dictionary(detailed_cost_table, params, tracked_params_list)

    # Fleet Mode adds its shared-campus results automatically, just as the
    # existing workflow automatically adds the reactor summary accounts.  The
    # user's tracked_params_list remains reserved for input/design parameters.
    if params.get('Fleet Mode', False):
        servicing_cost_table = fleet_results['servicing']
        tracked_costs.update(
            create_servicing_campus_cost_dictionary(servicing_cost_table)
        )
        manufacturing_cost_table = fleet_results['manufacturing']
        tracked_costs.update(
            create_manufacturing_campus_cost_dictionary(manufacturing_cost_table)
        )
        tracked_costs.update(
            create_fleet_total_cost_dictionary(fleet_results['total'])
        )

    file_exists = os.path.isfile(output_csv_filename)

    with open(output_csv_filename, 'a', newline='') as csvfile:
        fieldnames = tracked_costs.keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        if not file_exists or os.stat(output_csv_filename).st_size == 0:
            writer.writeheader()

        writer.writerow(tracked_costs)
        print(f"Results are being saved on {output_csv_filename}")


def detailed_bottom_up_cost_estimate(cost_database_filename):
    import inspect
    import os

    # Grab params and the calling script's path from the caller's frame automatically
    caller_frame = inspect.stack()[1][0]
    params = caller_frame.f_locals.get('params')
    if params is None:
        raise RuntimeError(
            "detailed_bottom_up_cost_estimate could not find 'params' in the calling scope. "
            "Make sure a variable named 'params' exists in the script that calls this function."
        )
    caller_file = caller_frame.f_globals.get('__file__', 'output')
    output_filename = os.path.splitext(os.path.abspath(caller_file))[0] + '_output.xlsx'

    if params.get('Fleet Mode', False):
        fleet_results = estimate_fleet_mode_costs(cost_database_filename, params)
        detailed_cost_table = fleet_results['reactor']
        detailed_servicing_cost_table = fleet_results['servicing']
        detailed_manufacturing_cost_table = fleet_results['manufacturing']
        fleet_annual_cash_flow = fleet_results['cash_flow']
        fleet_total_cost_table = fleet_results['total']
    else:
        detailed_cost_table = bottom_up_cost_estimate(cost_database_filename, params)
        detailed_servicing_cost_table = None
        detailed_manufacturing_cost_table = None
        fleet_annual_cash_flow = None
        fleet_total_cost_table = None
    detailed_central_cost_table = bottom_up_cost_estimate_central(cost_database_filename, params)
    pretty_df = transform_dataframe(detailed_cost_table)

    with pd.ExcelWriter(output_filename) as writer:
        pretty_df.to_excel(writer, sheet_name="cost estimate", index=False)

        if detailed_central_cost_table is not None:
            numerical_columns = detailed_central_cost_table.select_dtypes(include=[np.number]).columns
            nan_mask = detailed_central_cost_table[numerical_columns].isna().any(axis=1)
            if nan_mask.any():
                print("WARNING: NaN values in central facility accounts:")
                print(detailed_central_cost_table[nan_mask][['Account', 'Account Title'] + list(numerical_columns)])
            pretty_central_df = transform_dataframe(detailed_central_cost_table)
            pretty_central_df.to_excel(writer, sheet_name="central facility cost estimate", index=False)

        if detailed_servicing_cost_table is not None:
            numerical_columns = detailed_servicing_cost_table.select_dtypes(include=[np.number]).columns
            nan_mask = detailed_servicing_cost_table[numerical_columns].isna().any(axis=1)
            if nan_mask.any():
                print("WARNING: NaN values in servicing campus accounts:")
                print(detailed_servicing_cost_table.loc[nan_mask, ['Account', 'Account Title'] + list(numerical_columns)])
            servicing_output = format_servicing_campus_output(
                detailed_servicing_cost_table,
                params,
            )
            servicing_output.to_excel(
                writer,
                sheet_name="servicing campus cost estimate",
                index=False,
            )

        if detailed_manufacturing_cost_table is not None:
            numerical_columns = detailed_manufacturing_cost_table.select_dtypes(include=[np.number]).columns
            nan_mask = detailed_manufacturing_cost_table[numerical_columns].isna().any(axis=1)
            if nan_mask.any():
                print("WARNING: NaN values in manufacturing campus accounts:")
                print(detailed_manufacturing_cost_table.loc[nan_mask, ['Account', 'Account Title'] + list(numerical_columns)])
            manufacturing_output = format_manufacturing_campus_output(
                detailed_manufacturing_cost_table,
                params,
            )
            manufacturing_output.to_excel(
                writer,
                sheet_name="manufacturing campus cost estimate",
                index=False,
            )

        if fleet_total_cost_table is not None:
            fleet_total_output = fleet_total_cost_table.copy()
            cost_columns = [
                column for column in fleet_total_output.columns
                if column not in ['Account', 'Account Title', 'Units']
            ]
            lcoe_mask = fleet_total_output['Account'] == 'LCOE'
            fleet_total_output.loc[~lcoe_mask, cost_columns] = (
                fleet_total_output.loc[~lcoe_mask, cost_columns].round(0)
            )
            fleet_total_output.loc[lcoe_mask, cost_columns] = (
                fleet_total_output.loc[lcoe_mask, cost_columns].round(2)
            )
            fleet_total_output.to_excel(
                writer,
                sheet_name="fleet total cost estimate",
                index=False,
            )

        if fleet_annual_cash_flow is not None:
            fleet_cash_flow_output = fleet_annual_cash_flow.copy()
            integer_columns = [
                column for column in fleet_cash_flow_output.columns
                if column != 'Discount Factor'
            ]
            fleet_cash_flow_output[integer_columns] = (
                fleet_cash_flow_output[integer_columns].round(0)
            )
            fleet_cash_flow_output['Discount Factor'] = (
                fleet_cash_flow_output['Discount Factor'].round(8)
            )
            fleet_cash_flow_output.to_excel(
                writer,
                sheet_name="fleet annual cash flow",
                index=False,
            )

        save_params_to_excel_file(writer, params)

    # Always compute per-account LCOE contributions so they appear in Excel.
    # The PNG plot is only generated if params['plotting'] == "Y" —
    # that gate lives inside cost_drivers_estimate.
    lcoe_enriched_table, _ = cost_drivers_estimate(detailed_cost_table, params)

    if lcoe_enriched_table is not None:
        pretty_lcoe_df = transform_dataframe(lcoe_enriched_table)
        with pd.ExcelWriter(output_filename, mode='a', if_sheet_exists='replace') as writer:
            pretty_lcoe_df.to_excel(writer, sheet_name="cost estimate", index=False)

    print(f"\n\nThe cost estimate and all the parameters are saved at {output_filename}\n\n")
    return detailed_cost_table
