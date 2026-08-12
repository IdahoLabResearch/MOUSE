# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import os
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
                                   energy_cost_levelized_servicing_campus)
from cost.params_registry import PARAMS_REGISTRY, GROUP_ORDER
from reactor_engineering_evaluation.operation import reactor_operation
from cost.cost_drivers import cost_drivers_estimate



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


def format_servicing_campus_output(df, params):
    """Prepare the servicing-campus table for its Excel output sheet.

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
    discount_rate = params['Discount Rate']
    discount_sum = sum(
        (1 + discount_rate) ** -year
        for year in range(1, lifetime + 1)
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
        output.loc[annual_mask, 'Value'] * discount_sum / discounted_generation
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


def bottom_up_cost_estimate(cost_database_filename, params):
    # Validate tax credit params early — before any simulation or cost calculation runs.
    # This catches cases where a user accidentally defines both ITC and PTC,
    # which are mutually exclusive under the IRA.
    validate_tax_credit_params(params)

    escalated_cost = escalate_cost_database(cost_database_filename, params['Escalation Year'], params)
    escalated_cost_cleaned = remove_irrelevant_account(escalated_cost, params)
    reactor_operation(params)

    COA_list = []
    for i in range(params['Number of Samples']):
        if (i + 1) % 100 == 0:
            print(f"\n\nSample # {i+1}")

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


def bottom_up_cost_estimate_servicing_campus(cost_database_filename, params):
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
    return result.loc[nonzero_mask].reset_index(drop=True)


def create_servicing_campus_cost_dictionary(df):
    """Extract Fleet Mode summary metrics for a parametric-study CSV row."""
    accounts = {
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
                f"Servicing-campus result is missing required summary account '{account}'."
            )
        row = matching_rows.iloc[0]
        tracked_costs[f'Servicing Campus {label}'] = row['Value']
        tracked_costs[f'Servicing Campus {label} std'] = row['Value std']

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

    detailed_cost_table = bottom_up_cost_estimate(cost_database_filename, params)
    tracked_costs = create_cost_dictionary(detailed_cost_table, params, tracked_params_list)

    # Fleet Mode adds its shared-campus results automatically, just as the
    # existing workflow automatically adds the reactor summary accounts.  The
    # user's tracked_params_list remains reserved for input/design parameters.
    if params.get('Fleet Mode', False):
        servicing_cost_table = bottom_up_cost_estimate_servicing_campus(
            cost_database_filename,
            params,
        )
        tracked_costs.update(
            create_servicing_campus_cost_dictionary(servicing_cost_table)
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

    detailed_cost_table = bottom_up_cost_estimate(cost_database_filename, params)
    detailed_central_cost_table = bottom_up_cost_estimate_central(cost_database_filename, params)
    detailed_servicing_cost_table = bottom_up_cost_estimate_servicing_campus(cost_database_filename, params)
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
