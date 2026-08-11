# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import numpy as np
import pandas as pd
from cost.code_of_account_processing import get_estimated_cost_column

def validate_tax_credit_params(params):
    """
    Validates that the user has not selected both ITC and PTC simultaneously.
    These are mutually exclusive under the IRA — a project must choose one or neither.
    This should be called at the very start of the cost estimation workflow,
    before any simulation or cost calculation runs, to catch input errors early.

    @ In, params, dict, the user-defined parameters dictionary
    @ Out, None — raises ValueError if both ITC and PTC are defined
    """
    if 'ITC credit level' in params.keys() and 'PTC credit value' in params.keys():
        raise ValueError(
            "\n\n--- INPUT ERROR ---\n"
            "Both 'ITC credit level' and 'PTC credit value' are defined in params.\n"
            "Under the IRA, ITC and PTC are mutually exclusive — you must choose one or neither.\n"
            "Please remove one of them from your params and rerun.\n"
        )


def _crf(rate, period):
    # Returns the Capital Recovery Factor (CRF) based on the discount rate and period.
    # CRF converts a present value into a series of equal annual payments.
    # Formula: CRF = rate * (1 + rate)^period / ((1 + rate)^period - 1)
    # Special case: when rate == 0, CRF = 1/period (straight-line amortization, no interest)
    period = np.asarray(period, dtype=float)
    rate   = np.asarray(rate,   dtype=float)
    if rate.ndim == 0 and rate == 0:
        # scalar zero rate: straight-line amortization
        factor = np.where(period == 0, 0.0, 1.0 / period)
        return float(factor) if factor.ndim == 0 else factor
    numer = rate * (1 + rate)**period
    denum = (1 + rate)**period - 1
    factor = np.where(denum == 0, 0.0, numer / denum)
    return factor


def calculate_accounts_31_32_75_82_cost(df, params):
    estimated_cost_col_F = get_estimated_cost_column(df, 'F')
    estimated_cost_col_N = get_estimated_cost_column(df, 'N')

    for estimated_cost_col in [estimated_cost_col_F, estimated_cost_col_N]:
        filtered_df = df[df['Account'].isin([21, 22, 23])]
        tot_field_direct_cost = filtered_df[estimated_cost_col].sum()

        acct_31_cost = params['indirect to direct field-related cost'] * tot_field_direct_cost
        df.loc[df['Account'] == 31, estimated_cost_col] = acct_31_cost

        df.loc[df['Account'] == 32, estimated_cost_col] = df.loc[df['Account'] == 21, estimated_cost_col].values[0] * (df.loc[df['Account'] == 31, estimated_cost_col].values[0] / df.loc[df['Account'] == 22, estimated_cost_col].values[0])
        
        refueling_period = params['Fuel Lifetime'] + params['Refueling Period'] + params['Startup Duration after Refueling']
        refueling_period_yr = refueling_period / 365
        params_df = pd.DataFrame(params.items(), columns=['keys', 'values'])
        if params_df.loc[params_df['keys'].str.contains('replacement', case=False), 'keys'].size > 0:
            A20_replacement_period = refueling_period_yr * np.array([params['A75: Vessel Replacement Period (cycles)'],
                                                                    params['A75: Core Barrel Replacement Period (cycles)'],
                                                                     1,
                                                                     params['A75: Reflector Replacement Period (cycles)'],
                                                                     params['A75: Drum Replacement Period (cycles)'],
                                                                     params.get('A75: Integrated HX Replacement Period (cycles)', 0),])
            A20_capital_cost = np.array([df.loc[df['Account'] == 221.12, estimated_cost_col].values.sum(), 
                                         df.loc[df['Account'] == 221.13,  estimated_cost_col].values.sum(), 
                                         df.loc[df['Account'] == 221.33,  estimated_cost_col].values.sum(),
                                         df.loc[df['Account'] == 221.31,  estimated_cost_col].values.sum(),
                                         df.loc[df['Account'] == 221.2,   estimated_cost_col].values.sum(),
                                         df.loc[df['Account'].isin([222.1, 222.2, 222.3, 222.61]), estimated_cost_col].values.sum()])
            annualized_replacement_cost = (A20_capital_cost*_crf(params['Discount Rate'], A20_replacement_period))
            A20_other_cost = df.loc[df['Account'] == 20, estimated_cost_col].values[0] - A20_capital_cost.sum()
            annualized_other_cost = A20_other_cost * params['Maintenance to Direct Cost Ratio']
            df.loc[df['Account'] == 751, estimated_cost_col] = annualized_replacement_cost[0]
            df.loc[df['Account'] == 752, estimated_cost_col] = annualized_replacement_cost[1]
            df.loc[df['Account'] == 753, estimated_cost_col] = annualized_replacement_cost[2]
            df.loc[df['Account'] == 754, estimated_cost_col] = annualized_replacement_cost[3]
            df.loc[df['Account'] == 755, estimated_cost_col] = annualized_replacement_cost[4]
            df.loc[df['Account'] == 756, estimated_cost_col] = annualized_replacement_cost[5]
            df.loc[df['Account'] == 759, estimated_cost_col] = annualized_other_cost
        else:
            df.loc[df['Account'] == 75, estimated_cost_col] = df.loc[df['Account'] == 20, estimated_cost_col].values[0] * params['Maintenance to Direct Cost Ratio']

        lump_fuel_cost = df.loc[df['Account'] == 25, estimated_cost_col].values[0]
        annualized_fuel_cost = lump_fuel_cost*_crf(params['Discount Rate'], refueling_period_yr)
        df.loc[df['Account'] == 82, estimated_cost_col] = annualized_fuel_cost

    return df


def calculate_accounts_31_32_75_central_facility_cost(df, params):
    """
    Calculate indirect costs for central facility accounts (31, 32, 75).
    Similar to calculate_accounts_31_32_75_82_cost but for central facility.
    """
    estimated_cost_col_F = get_estimated_cost_column(df, 'F')
    estimated_cost_col_N = get_estimated_cost_column(df, 'N')

    for estimated_cost_col in [estimated_cost_col_F, estimated_cost_col_N]:
        # Filter for accounts 21, 22, 23, 24, 25, 27
        filtered_df = df[df['Account'].isin([21, 22, 23, 24, 25, 27])]
        tot_field_direct_cost = filtered_df[estimated_cost_col].sum()

        acct_31_cost = params['indirect to direct field-related cost'] * tot_field_direct_cost
        df.loc[df['Account'] == 31, estimated_cost_col] = acct_31_cost

        # Calculate Account 32 using ratio of Account 31 to reactor systems cost
        df.loc[df['Account'] == 32, estimated_cost_col] = (
            df.loc[df['Account'] == 21, estimated_cost_col].values[0]
            * (df.loc[df['Account'] == 31, estimated_cost_col].values[0]
               / df.loc[df['Account'].isin([22, 23, 24, 25]), estimated_cost_col].values[0])
        )

        # A75: Maintenance costs as percentage of direct costs
        df.loc[df['Account'] == 75, estimated_cost_col] = (
            df.loc[df['Account'] == 20, estimated_cost_col].values[0] * params['Maintenance to Direct Cost Ratio']
        )

    return df


def calculate_servicing_campus_derived_costs(df, params):
    """Populate servicing-campus accounts whose costs are defined by ratios."""
    cost_columns = [get_estimated_cost_column(df, option) for option in ('F', 'N')]

    for cost_column in cost_columns:
        direct_cost = df.loc[df['Account'] == 20, cost_column].sum()
        staffing_cost = df.loc[df['Account'].isin([711, 712, 713, 714, 715]), cost_column].sum()

        df.loc[df['Account'] == 30, cost_column] = (
            params['SER Account 30 to Account 20 Ratio'] * direct_cost
        )
        df.loc[df['Account'] == 716, cost_column] = (
            params['SER Account 716 to Accounts 711-715 Ratio'] * staffing_cost
        )
        df.loc[df['Account'] == 717, cost_column] = (
            params['SER Account 717 to Accounts 711-715 Ratio'] * staffing_cost
        )
        df.loc[df['Account'] == 741, cost_column] = (
            params['SER Account 741 to Account 20 Ratio'] * direct_cost
        )
        df.loc[df['Account'] == 742, cost_column] = (
            params['SER Account 742 to Account 20 Ratio'] * direct_cost
        )
        df.loc[df['Account'] == 743, cost_column] = (
            params['SER Account 743 to Account 20 Ratio'] * direct_cost
        )
        df.loc[df['Account'] == 747.1, cost_column] = (
            params['SER Account 747.1 to Accounts 741-743 Ratio']
            * df.loc[df['Account'].isin([741, 742, 743]), cost_column].sum()
        )
        df.loc[df['Account'] == 747.4, cost_column] = (
            params['SER Account 747.4 to Account 20 Ratio'] * direct_cost
        )

    return df


def calculate_interest_cost_servicing_campus(params, OCC):
    """Calculate servicing-campus interest during construction."""
    interest_rate = params['Interest Rate']
    construction_duration = params['SER Construction Duration']
    debt_to_equity_ratio = params['Debt To Equity Ratio']
    debt_fraction = debt_to_equity_ratio / (1 + debt_to_equity_ratio)
    B = 1 + np.exp(np.log(1 + interest_rate) * construction_duration / 12)
    C = (np.log(1 + interest_rate) * (construction_duration / 12) / 3.14) ** 2 + 1
    return debt_fraction * OCC * ((0.5 * B / C) - 1)


def calculate_servicing_campus_capital_costs(df, params):
    """Add OCC rows and Account 62 financing for the servicing campus."""
    generating_sites = params['Generating Sites Count']
    if generating_sites <= 0:
        raise ValueError("'Generating Sites Count' must be greater than zero.")

    df = pd.concat([
        df,
        pd.DataFrame([
            {'Account': 'OCC', 'Account Title': 'Overnight Capital Cost'},
            {'Account': 'OCC per reactor', 'Account Title': 'Overnight Capital Cost per Operating Reactor'},
        ]),
    ], ignore_index=True)

    accounts_to_sum = [10, 20, 30, 40, 50]
    for cost_column in [get_estimated_cost_column(df, option) for option in ('F', 'N')]:
        occ_cost = df.loc[df['Account'].isin(accounts_to_sum), cost_column].sum()
        df.loc[df['Account'] == 'OCC', cost_column] = occ_cost
        df.loc[df['Account'] == 'OCC per reactor', cost_column] = occ_cost / generating_sites
        df.loc[df['Account'] == 62, cost_column] = calculate_interest_cost_servicing_campus(params, occ_cost)

    return df


def calculate_servicing_campus_TCI(df, params):
    """Add total capital investment rows for the servicing campus."""
    generating_sites = params['Generating Sites Count']
    df = pd.concat([
        df,
        pd.DataFrame([
            {'Account': 'TCI', 'Account Title': 'Total Capital Investment'},
            {'Account': 'TCI per reactor', 'Account Title': 'Total Capital Investment per Operating Reactor'},
        ]),
    ], ignore_index=True)

    for cost_column in [get_estimated_cost_column(df, option) for option in ('F', 'N')]:
        tci_cost = df.loc[df['Account'].isin(['OCC', 60]), cost_column].sum()
        df.loc[df['Account'] == 'TCI', cost_column] = tci_cost
        df.loc[df['Account'] == 'TCI per reactor', cost_column] = tci_cost / generating_sites

    return df


def energy_cost_levelized_servicing_campus(params, df):
    """Add annual-cost summaries and the servicing contribution to fleet LCOE."""
    generating_sites = params['Generating Sites Count']
    fleet_size = params['Fleet']
    annual_generation_per_reactor = params['Annual Electricity Production']
    fleet_annual_generation = fleet_size * annual_generation_per_reactor
    if fleet_annual_generation <= 0:
        raise ValueError("Fleet annual electricity generation must be greater than zero.")

    df = pd.concat([
        df,
        pd.DataFrame([
            {'Account': 'Annual Cost', 'Account Title': 'Annual Servicing Campus Cost'},
            {'Account': 'Annual Cost per reactor', 'Account Title': 'Annual Servicing Campus Cost per Operating Reactor'},
            {'Account': 'LCOE', 'Account Title': 'Servicing Campus Contribution to Fleet LCOE ($/MWh)'},
        ]),
    ], ignore_index=True)

    lifetime = int(params['Levelization Period'])
    discount_rate = params['Discount Rate']
    discount_factors = np.array([(1 + discount_rate) ** -year for year in range(1, lifetime + 1)])
    discounted_generation = fleet_annual_generation * discount_factors.sum()

    for cost_column in [get_estimated_cost_column(df, option) for option in ('F', 'N')]:
        annual_cost = df.loc[df['Account'] == 70, cost_column].sum()
        capital_cost = df.loc[df['Account'] == 'TCI', cost_column].sum()
        discounted_cost = capital_cost + annual_cost * discount_factors.sum()

        df.loc[df['Account'] == 'Annual Cost', cost_column] = annual_cost
        df.loc[df['Account'] == 'Annual Cost per reactor', cost_column] = annual_cost / generating_sites
        df.loc[df['Account'] == 'LCOE', cost_column] = discounted_cost / discounted_generation

    return df


def calculate_decommissioning_cost(df, params):
    estimated_cost_col_F = get_estimated_cost_column(df, 'F')
    estimated_cost_col_N = get_estimated_cost_column(df, 'N')

    for estimated_cost_col in [estimated_cost_col_F, estimated_cost_col_N]:
        capex = df.loc[df['Account'].isin([10, 20]), estimated_cost_col].sum()
        AR = params['Annual Return']
        LP = params['Levelization Period']
        
        if 'A78: CAPEX to Decommissioning Cost Ratio' not in params.keys():
            params['A78: CAPEX to Decommissioning Cost Ratio'] = 0.15

        decommissioning_fv_cost = capex * params['A78: CAPEX to Decommissioning Cost Ratio']
        fv_to_pv_of_annuity = -AR/(1- pow(1+AR, LP))
        annualized_decommisioning_cost = decommissioning_fv_cost * fv_to_pv_of_annuity     
        df.loc[df['Account'] == 78, estimated_cost_col] = annualized_decommisioning_cost

    return df


def calculate_interest_cost(params, OCC):
    interest_rate = params['Interest Rate']
    construction_duration = params['Construction Duration']
    debt_to_equity_ratio = params['Debt To Equity Ratio']
    # Convert D:E ratio to debt fraction for the calculation.
    # e.g. D:E = 1.0 (1:1) → debt_fraction = 1/(1+1) = 0.5 (50% financed by debt)
    # e.g. D:E = 2.33      → debt_fraction = 2.33/3.33 ≈ 0.7 (70% financed by debt)
    debt_fraction = debt_to_equity_ratio / (1 + debt_to_equity_ratio)
    B = (1 + np.exp((np.log(1 + interest_rate)) * construction_duration / 12))
    C = ((np.log(1 + interest_rate) * (construction_duration / 12) / 3.14)**2 + 1)
    Interest_expenses = debt_fraction * OCC * ((0.5 * B / C) - 1)
    return Interest_expenses


def calculate_interest_cost_central(params, OCC):
    """Calculate interest cost for central facility using its construction duration."""
    interest_rate = params['Interest Rate']
    construction_duration = params['Central Facility Construction Duration']
    debt_to_equity_ratio = params['Debt To Equity Ratio']
    # Convert D:E ratio to debt fraction for the calculation (same logic as above).
    debt_fraction = debt_to_equity_ratio / (1 + debt_to_equity_ratio)
    B = (1 + np.exp((np.log(1 + interest_rate)) * construction_duration / 12))
    C = ((np.log(1 + interest_rate) * (construction_duration / 12) / 3.14)**2 + 1)
    Interest_expenses = debt_fraction * OCC * ((0.5 * B / C) - 1)
    return Interest_expenses


def calculate_high_level_capital_costs(df, params):
    power_kWe = 1000 * params['Power MWe']
    accounts_to_sum = [10, 20, 30, 40, 50]

    df = pd.concat([df, pd.DataFrame([{'Account': 'OCC','Account Title' : 'Overnight Capital Cost'}])], ignore_index=True)
    df = pd.concat([df, pd.DataFrame([{'Account': 'OCC per kW','Account Title' : 'Overnight Capital Cost per kW' }])], ignore_index=True)
    df = pd.concat([df, pd.DataFrame([{'Account': 'OCC excl. fuel','Account Title' : 'Overnight Capital Cost Excluding Fuel'}])], ignore_index=True)
    df = pd.concat([df, pd.DataFrame([{'Account': 'OCC excl. fuel per kW','Account Title' : 'Overnight Capital Cost Excluding Fuel per kW'}])], ignore_index=True)

    cost_column_F = get_estimated_cost_column(df, 'F')
    cost_column_N = get_estimated_cost_column(df, 'N')

    for cost_column in [cost_column_F, cost_column_N]:
        occ_cost = df[df['Account'].isin(accounts_to_sum)][cost_column].sum()
        df.loc[df['Account'] == 'OCC', cost_column] = occ_cost
        df.loc[df['Account'] == 'OCC per kW', cost_column] = occ_cost/ power_kWe
        
        occ_excl_fuel = occ_cost - (df.loc[df['Account'] == 25, cost_column].values[0])
        df.loc[df['Account'] == 'OCC excl. fuel', cost_column] = occ_excl_fuel
        df.loc[df['Account'] == 'OCC excl. fuel per kW', cost_column] = occ_excl_fuel/ power_kWe

        df.loc[df['Account'] == 62, cost_column] = calculate_interest_cost(params, occ_cost)
    return df


def calculate_high_level_capital_costs_central_facility(df, params):
    """Calculate OCC and interest costs for central facility."""
    power_kWe = 1000 * params['Power MWe'] * params['Maximum Number of Operating Reactors']
    accounts_to_sum = [10, 20, 30, 40, 50]

    df = pd.concat([df, pd.DataFrame([{'Account': 'OCC', 'Account Title': 'Overnight Capital Cost'}])], ignore_index=True)

    cost_column_F = get_estimated_cost_column(df, 'F')
    cost_column_N = get_estimated_cost_column(df, 'N')

    for cost_column in [cost_column_F, cost_column_N]:
        occ_cost = df[df['Account'].isin(accounts_to_sum)][cost_column].sum()
        df.loc[df['Account'] == 'OCC', cost_column] = occ_cost
        df.loc[df['Account'] == 'OCC per kW', cost_column] = occ_cost / power_kWe
        df.loc[df['Account'] == 62, cost_column] = calculate_interest_cost_central(params, occ_cost)
    return df


def calculate_TCI_central(df, params):
    """Calculate Total Capital Investment for central facility."""
    power_kWe = 1000 * params['Power MWe'] * params['Maximum Number of Operating Reactors']

    df = pd.concat([df, pd.DataFrame([{'Account': 'TCI', 'Account Title': 'Total Capital Investment'}])], ignore_index=True)
    df = pd.concat([df, pd.DataFrame([{'Account': 'TCI per kW', 'Account Title': 'Total Capital Investment per kW'}])], ignore_index=True)

    accounts_to_sum = ['OCC', 60]
    cost_column_F = get_estimated_cost_column(df, 'F')
    cost_column_N = get_estimated_cost_column(df, 'N')

    for cost_column in [cost_column_F, cost_column_N]:
        tci_cost = df[df['Account'].isin(accounts_to_sum)][cost_column].sum()
        df.loc[df['Account'] == 'TCI', cost_column] = tci_cost
        df.loc[df['Account'] == 'TCI per kW', cost_column] = tci_cost / power_kWe

    return df


# -----------------------------------------------------------------------------------------
# ITC (Investment Tax Credit) helper function
# -----------------------------------------------------------------------------------------
# The ITC is a one-time credit applied to the capital cost (OCC) of the plant.
# Under the IRA, the ITC level can be 6%, 30%, 40%, or 50% depending on whether
# the project meets prevailing wage, domestic content, and energy community requirements.
#
# This function returns a COST REDUCTION FACTOR (not the credit itself).
# The factor represents what fraction of the original OCC remains after the ITC is applied.
# Example: a 30% ITC reduces OCC to 73% of its original value → factor = 0.73
#
# itc_level is expressed as a fraction (e.g. 0.30 for 30%)
# Interpolation is used for ITC levels between the defined breakpoints.
# -----------------------------------------------------------------------------------------
def ITC_reduction_factor(itc_level):
    itc_values    = [0,    0.06,  0.3,   0.4,   0.5 ]  # ITC credit levels (fractions)
    reduction_factors = [1, 0.95,  0.73,  0.63,  0.53]  # corresponding OCC reduction factors
    # renamed from ITC_reduction_factor to avoid shadowing the function name
    return np.interp(itc_level, itc_values, reduction_factors)


def calculate_TCI(df, params):
    # -----------------------------------------------------------------------------------------
    # Total Capital Investment (TCI) = OCC + Account 60 (financing/interest costs)
    #
    # If an ITC credit level is provided in params, a second version of TCI is calculated
    # that reflects the reduced OCC after the ITC subsidy is applied:
    #   - OCC with ITC     = OCC × ITC_reduction_factor(itc_level)
    #   - TCI with ITC     = OCC with ITC + Account 60
    #
    # Note: Account 60 (financing costs) is NOT reduced by the ITC — only the OCC is.
    # This is consistent with how the ITC works in practice: it offsets capital investment,
    # not the financing charges on top of it.
    #
    # Output rows added to the cost dataframe:
    #   - 'OCC with ITC'         : reduced overnight capital cost
    #   - 'OCC with ITC per kW'  : reduced OCC normalized by plant capacity
    #   - 'TCI with ITC'         : reduced total capital investment
    #   - 'TCI with ITC per kW'  : reduced TCI normalized by plant capacity
    # -----------------------------------------------------------------------------------------
    power_kWe = 1000 * params['Power MWe']

    df = pd.concat([df, pd.DataFrame([{'Account': 'TCI','Account Title' : 'Total Capital Investment'}])], ignore_index=True)
    df = pd.concat([df, pd.DataFrame([{'Account': 'TCI per kW','Account Title' : 'Total Capital Investment per kW'}])], ignore_index=True)

    if 'ITC credit level' in params.keys():
        # Add ITC-adjusted output rows to the dataframe
        df = pd.concat([df, pd.DataFrame([{'Account': 'OCC (ITC-adjusted)',        'Account Title': 'Overnight Capital Cost Adjusted for the Investment Tax Credit'}])], ignore_index=True)
        df = pd.concat([df, pd.DataFrame([{'Account': 'OCC (ITC-adjusted) per kW', 'Account Title': 'Overnight Capital Cost Adjusted for the Investment Tax Credit per kW'}])], ignore_index=True)
        df = pd.concat([df, pd.DataFrame([{'Account': 'TCI (ITC-adjusted)',        'Account Title': 'Total Capital Investment Adjusted for the Investment Tax Credit'}])], ignore_index=True)
        df = pd.concat([df, pd.DataFrame([{'Account': 'TCI (ITC-adjusted) per kW', 'Account Title': 'Total Capital Investment Adjusted for the Investment Tax Credit per kW'}])], ignore_index=True)
        # note: ITC_cost_reduction_factor is computed inside the loop below for each cost column (FOAK and NOAK)

    accounts_to_sum = ['OCC', 60]
    cost_column_F = get_estimated_cost_column(df, 'F')
    cost_column_N = get_estimated_cost_column(df, 'N')

    # Per-column eligibility for ITC. FOAK column = unit 1; NOAK column = unit
    # 'NOAK Unit Number'. A unit qualifies only if its position is <= the cutoff
    # ('Number of Units Claiming ITC/PTC'). Default cutoff is effectively
    # infinite so omitting the param preserves the pre-existing behavior.
    n_credit = params.get('Number of Units Claiming ITC/PTC', 10**9)
    noak_unit = params.get('NOAK Unit Number', 10)
    eligible_by_column = {cost_column_F: 1 <= n_credit,
                          cost_column_N: noak_unit <= n_credit}

    for cost_column in [cost_column_F, cost_column_N]:
        # --- Baseline TCI (no ITC) ---
        tci_cost = df[df['Account'].isin(accounts_to_sum)][cost_column].sum()
        df.loc[df['Account'] == 'TCI', cost_column] = tci_cost
        df.loc[df['Account'] == 'TCI per kW', cost_column] = tci_cost / power_kWe

        if 'ITC credit level' in params.keys():
            if eligible_by_column[cost_column]:
                # --- ITC-adjusted TCI ---
                # Step 1: Get the reduction factor for the given ITC level
                ITC_cost_reduction_factor = ITC_reduction_factor(params['ITC credit level'])
                # Step 2: Apply the reduction factor to OCC to get the ITC-adjusted OCC
                # OCC_after_ITC is the reduced OCC value (not the savings amount)
                OCC_after_ITC = df.loc[df['Account'] == 'OCC', cost_column].values[0] * ITC_cost_reduction_factor
                # Step 3: Add financing costs (Account 60) to get TCI adjusted for ITC
                # Note: Account 60 is not reduced by the ITC
                tci_cost_with_itc = df.loc[df['Account'] == 60, cost_column].values[0] + OCC_after_ITC
            else:
                # This unit is past the IRA sunset cutoff — fall back to the
                # un-subsidized OCC/TCI so the ITC-adjusted columns show the
                # un-subsidized cost rather than blank/NaN.
                OCC_after_ITC = df.loc[df['Account'] == 'OCC', cost_column].values[0]
                tci_cost_with_itc = tci_cost
            df.loc[df['Account'] == 'OCC (ITC-adjusted)', cost_column] = OCC_after_ITC
            df.loc[df['Account'] == 'OCC (ITC-adjusted) per kW', cost_column] = OCC_after_ITC / power_kWe
            df.loc[df['Account'] == 'TCI (ITC-adjusted)', cost_column] = tci_cost_with_itc
            df.loc[df['Account'] == 'TCI (ITC-adjusted) per kW', cost_column] = tci_cost_with_itc / power_kWe

    return df


def energy_cost_levelized(params, df):
    # -----------------------------------------------------------------------------------------
    # LCOE (Levelized Cost of Energy) Calculation
    # ... (existing docstring unchanged)
    # -----------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------
    # Heat application cost reduction factors (hardcoded, based on process heat study)
    # For heat applications, the OCC is lower because no power conversion system is needed
    # (e.g. no turbine, generator, condenser). The annual O&M cost is also slightly reduced.
    # Source: [add your reference here]
    # -----------------------------------------------------------------------------------------
    HEAT_OCC_FACTOR         = 0.795  # OCC for heat = OCC_electric × 0.795
    HEAT_ANNUAL_COST_FACTOR = 0.966  # Annual O&M+fuel cost for heat = baseline × 0.966

    df = pd.concat([df, pd.DataFrame([{'Account': 'AC',         'Account Title': 'Annualized Cost'}])], ignore_index=True)
    df = pd.concat([df, pd.DataFrame([{'Account': 'AC per MWh', 'Account Title': 'Annualized Cost per MWh'}])], ignore_index=True)
    df = pd.concat([df, pd.DataFrame([{'Account': 'LCOE',       'Account Title': 'Levelized Cost Of Energy ($/MWh)'}])], ignore_index=True)

    if 'PTC credit value' in params.keys():
        df = pd.concat([df, pd.DataFrame([{'Account': 'LCOE with PTC', 'Account Title': 'Levelized Cost Of Energy with PTC ($/MWh)'}])], ignore_index=True)

    if 'ITC credit level' in params.keys():
        assert 'PTC credit value' not in params.keys(), '--error: Only PTC or ITC or None must be selected not both.'
        df = pd.concat([df, pd.DataFrame([{'Account': 'LCOE (ITC-adjusted)', 'Account Title': 'Levelized Cost Of Energy Adjusted for the Investment Tax Credit ($/MWh)'}])], ignore_index=True)

    df = pd.concat([df, pd.DataFrame([{'Account': 'LCOH',       'Account Title': 'Levelized Cost Of Heat ($/MWth)'}])], ignore_index=True)

    params.setdefault('Tax Rate', 0.21)

    plant_lifetime_years = params['Levelization Period']
    discount_rate        = params['Discount Rate']
    power_MWe            = params['Power MWe']
    capacity_factor      = params['Capacity Factor']
    thermal_efficiency   = params['Thermal Efficiency']
    estimated_cost_col_F = get_estimated_cost_column(df, 'F')
    estimated_cost_col_N = get_estimated_cost_column(df, 'N')

    # Per-column eligibility for ITC/PTC under the IRA sunset cutoff.
    # FOAK column = unit 1; NOAK column = unit 'NOAK Unit Number'. A unit
    # qualifies only if its position <= 'Number of Units Claiming ITC/PTC'.
    # Default cutoff is effectively infinite so omitting the param preserves
    # the pre-existing behavior. When ineligible, the credit-adjusted accounts
    # fall back to the un-subsidized LCOE.
    n_credit = params.get('Number of Units Claiming ITC/PTC', 10**9)
    noak_unit = params.get('NOAK Unit Number', 10)
    eligible_by_column = {estimated_cost_col_F: 1 <= n_credit,
                          estimated_cost_col_N: noak_unit <= n_credit}

    for estimated_cost_col in [estimated_cost_col_F, estimated_cost_col_N]:

        # -----------------------------------------------------------------------------------------
        # Baseline LCOE calculation (no tax credits) — unchanged
        # -----------------------------------------------------------------------------------------
        cap_cost          = df.loc[df['Account'] == 'TCI', estimated_cost_col].values[0]
        ann_cost          = df.loc[df['Account'] == 70, estimated_cost_col].values[0] + df.loc[df['Account'] == 80, estimated_cost_col].values[0]
        levelized_ann_cost = ann_cost / params['Annual Electricity Production']
        df.loc[df['Account'] == 'AC',        estimated_cost_col] = ann_cost
        df.loc[df['Account'] == 'AC per MWh', estimated_cost_col] = levelized_ann_cost

        sum_cost = 0
        sum_elec = 0
        for i in range(1 + plant_lifetime_years):
            if i == 0:
                cap_cost_per_year = cap_cost
                annual_cost       = 0
                elec_gen          = 0
            else:
                cap_cost_per_year = 0
                annual_cost       = ann_cost
                elec_gen          = power_MWe * capacity_factor * 365 * 24
            sum_cost += (cap_cost_per_year + annual_cost) / ((1 + discount_rate)**i)
            sum_elec += elec_gen / ((1 + discount_rate)**i)

        lcoe = sum_cost / sum_elec
        df.loc[df['Account'] == 'LCOE', estimated_cost_col] = lcoe

        # -----------------------------------------------------------------------------------------
        # LCOH (Levelized Cost of Heat) calculation
        #
        # For heat applications, the plant does not need a power conversion system,
        # so both capital and O&M costs are reduced by the factors defined above.
        #
        # The full heat cost chain (all intermediate values are behind the scenes):
        #   1. OCC_heat      = OCC × HEAT_OCC_FACTOR
        #   2. Interest_heat = calculate_interest_cost(params, OCC_heat)
        #   3. TCI_heat      = OCC_heat + Interest_heat
        #   4. ann_cost_heat = ann_cost × HEAT_ANNUAL_COST_FACTOR
        #   5. LCOE_heat     = PV(costs with TCI_heat, ann_cost_heat) / PV(electricity)
        #   6. LCOH          = LCOE_heat × Thermal Efficiency
        #
        # The ×Thermal Efficiency step converts from $/MWhe to $/MWth:
        # e.g. η=0.33 → 1 MWhe = 3 MWth → LCOH = LCOE_heat × 0.33 → cheaper per MWth
        # -----------------------------------------------------------------------------------------
        OCC           = df.loc[df['Account'] == 'OCC', estimated_cost_col].values[0]
        OCC_heat      = OCC * HEAT_OCC_FACTOR
        Interest_heat = calculate_interest_cost(params, OCC_heat)
        TCI_heat      = OCC_heat + Interest_heat
        ann_cost_heat = ann_cost * HEAT_ANNUAL_COST_FACTOR

        sum_cost_heat = 0
        sum_elec_heat = 0
        for i in range(1 + plant_lifetime_years):
            if i == 0:
                cap_cost_per_year = TCI_heat
                annual_cost_heat  = 0
                elec_gen          = 0
            else:
                cap_cost_per_year = 0
                annual_cost_heat  = ann_cost_heat
                elec_gen          = power_MWe * capacity_factor * 365 * 24
            sum_cost_heat += (cap_cost_per_year + annual_cost_heat) / ((1 + discount_rate)**i)
            sum_elec_heat += elec_gen / ((1 + discount_rate)**i)

            lcoe = sum_cost / sum_elec
            df.loc[df['Account'] == 'LCOE (ITC-adjusted)', estimated_cost_col] = lcoe

        # -----------------------------------------------------------------------------------------
        # LCOH — always computed last so it appears as the final row in the output table
        # -----------------------------------------------------------------------------------------
        OCC           = df.loc[df['Account'] == 'OCC', estimated_cost_col].values[0]
        OCC_heat      = OCC * HEAT_OCC_FACTOR
        Interest_heat = calculate_interest_cost(params, OCC_heat)
        TCI_heat      = OCC_heat + Interest_heat
        ann_cost_heat = ann_cost * HEAT_ANNUAL_COST_FACTOR

        sum_cost_heat = 0
        sum_elec_heat = 0
        for i in range(1 + plant_lifetime_years):
            if i == 0:
                cap_cost_per_year = TCI_heat
                annual_cost_heat  = 0
                elec_gen          = 0
            else:
                cap_cost_per_year = 0
                annual_cost_heat  = ann_cost_heat
                elec_gen          = power_MWe * capacity_factor * 365 * 24
            sum_cost_heat += (cap_cost_per_year + annual_cost_heat) / ((1 + discount_rate)**i)
            sum_elec_heat += elec_gen / ((1 + discount_rate)**i)

        lcoe_heat = sum_cost_heat / sum_elec_heat
        lcoh      = lcoe_heat * thermal_efficiency
        df.loc[df['Account'] == 'LCOH', estimated_cost_col] = lcoh
        # -----------------------------------------------------------------------------------------
        if 'PTC credit value' in params.keys():
            if eligible_by_column[estimated_cost_col]:
                sum_elec = 0
                sum_ptc  = 0
                assert 'PTC credit period' in params.keys(), 'error: If a PTC credit value is provided, a corresponding PTC credit period must be given as well.'
                try:
                    bonus_multiplier = 1.0 + params['domestic_content_bonus'] + params['energy_community_bonus']
                except:
                    print('--- warning: Assume no extra percentage on the credit')
                    bonus_multiplier = 1.0

                for i in range(1 + plant_lifetime_years):
                    if i == 0:
                        elec_gen = 0
                        ptc_gen  = 0
                    else:
                        elec_gen = power_MWe * capacity_factor * 365 * 24
                        if i > params['PTC credit period']:
                            ptc_gen = 0
                        else:
                            ptc_gen = elec_gen * (params['PTC credit value'] * bonus_multiplier) / (1 - params['Tax Rate'])
                    sum_ptc  += ptc_gen  / ((1 + discount_rate)**i)
                    sum_elec += elec_gen / ((1 + discount_rate)**i)

                estimated_ptc = sum_ptc / sum_elec
                df.loc[df['Account'] == 'LCOE with PTC', estimated_cost_col] = lcoe - estimated_ptc
            else:
                # This unit is past the IRA sunset cutoff — fall back to the
                # un-subsidized LCOE so the 'LCOE with PTC' cell shows the
                # un-subsidized cost rather than blank/NaN.
                df.loc[df['Account'] == 'LCOE with PTC', estimated_cost_col] = lcoe

        # -----------------------------------------------------------------------------------------
        # ITC adjustment — unchanged
        # -----------------------------------------------------------------------------------------
        if 'ITC credit level' in params.keys():
            cap_cost      = df.loc[df['Account'] == 'TCI (ITC-adjusted)', estimated_cost_col].values[0]
            ann_cost      = df.loc[df['Account'] == 70, estimated_cost_col].values[0] + df.loc[df['Account'] == 80, estimated_cost_col].values[0]
            levelized_ann_cost = ann_cost / params['Annual Electricity Production']
            df.loc[df['Account'] == 'AC',         estimated_cost_col] = ann_cost
            df.loc[df['Account'] == 'AC per MWh', estimated_cost_col] = levelized_ann_cost
            sum_cost = 0
            sum_elec = 0

            for i in range(1 + plant_lifetime_years):
                if i == 0:
                    cap_cost_per_year = cap_cost
                    annual_cost       = 0
                    elec_gen          = 0
                else:
                    cap_cost_per_year = 0
                    annual_cost       = ann_cost
                    elec_gen          = power_MWe * capacity_factor * 365 * 24
                sum_cost += (cap_cost_per_year + annual_cost) / ((1 + discount_rate)**i)
                sum_elec += elec_gen / ((1 + discount_rate)**i)

            lcoe = sum_cost / sum_elec
            df.loc[df['Account'] == 'LCOE (ITC-adjusted)', estimated_cost_col] = lcoe

    return df

    for estimated_cost_col in [estimated_cost_col_F, estimated_cost_col_N]:

        # -----------------------------------------------------------------------------------------
        # Baseline LCOE calculation (no tax credits)
        # Capital cost is paid upfront at year 0; O&M and fuel costs are paid annually from year 1
        # -----------------------------------------------------------------------------------------
        cap_cost = df.loc[df['Account'] == 'TCI', estimated_cost_col].values[0]
        ann_cost = df.loc[df['Account'] == 70, estimated_cost_col].values[0] + df.loc[df['Account'] == 80, estimated_cost_col].values[0]
        levelized_ann_cost = ann_cost / params['Annual Electricity Production']
        df.loc[df['Account'] == 'AC', estimated_cost_col] = ann_cost
        df.loc[df['Account'] == 'AC per MWh', estimated_cost_col] = levelized_ann_cost
        sum_cost = 0
        sum_elec = 0

        for i in range(1 + plant_lifetime_years):
            if i == 0:
                # Year 0: capital cost is paid, no electricity is produced
                cap_cost_per_year = cap_cost
                annual_cost = 0
                elec_gen = 0
            elif i > 0:
                # Years 1 to plant_lifetime: O&M and fuel costs are paid, electricity is produced
                cap_cost_per_year = 0
                annual_cost = ann_cost
                elec_gen = power_MWe * capacity_factor * 365 * 24  # MWh/year
            sum_cost += (cap_cost_per_year + annual_cost) / ((1 + discount_rate)**i)
            sum_elec += elec_gen / ((1 + discount_rate)**i)

        # Divide PV of costs by PV of electricity to get levelized $/MWh
        lcoe = sum_cost / sum_elec
        df.loc[df['Account'] == 'LCOE', estimated_cost_col] = lcoe

        # -----------------------------------------------------------------------------------------
        # PTC (Production Tax Credit) adjustment
        # -----------------------------------------------------------------------------------------
        # The PTC is a per-MWh tax credit earned for every MWh produced during the credit period.
        # To subtract it from the before-tax LCOE consistently, the PTC must be grossed up
        # by (1 - tax_rate) to convert it to its before-tax revenue equivalent.
        #
        # Without gross-up: $15/MWh credit would be subtracted directly → underestimates benefit
        # With gross-up:    $15/MWh / (1 - 0.21) = $18.99/MWh → correct before-tax equivalent
        #
        # The PTC is only earned during the credit period (e.g. 10 years).
        # Both PTC revenue and electricity production are discounted to present value,
        # then divided to give a levelized $/MWh equivalent that accounts for the fact
        # that credits only come in the first N years but electricity spans the full plant life.
        #
        # Bonus multiplier accounts for IRA stackable bonuses:
        #   - domestic_content_bonus : extra % for using US-made materials
        #   - energy_community_bonus : extra % for siting in an energy community
        # -----------------------------------------------------------------------------------------
        if 'PTC credit value' in params.keys():
            sum_elec = 0
            sum_ptc = 0
            assert 'PTC credit period' in params.keys(), 'error: If a PTC credit value is provided, a corresponding PTC credit period must be given as well.'
            try:
                # Apply stackable IRA bonus multipliers if provided
                bonus_multiplier = 1.0 + params['domestic_content_bonus'] + params['energy_community_bonus']
            except:
                print('--- warning: Assume no extra percentage on the credit')
                bonus_multiplier = 1.0

            for i in range(1 + plant_lifetime_years):
                if i == 0:
                    # Year 0: construction year, no electricity or PTC
                    elec_gen = 0
                    ptc_gen = 0
                elif i > 0:
                    elec_gen = power_MWe * capacity_factor * 365 * 24  # MWh/year
                    if i > params['PTC credit period']:
                        # Credit period has expired — no more PTC
                        ptc_gen = 0
                    else:
                        # PTC earned this year, grossed up to before-tax equivalent
                        # Gross-up formula: PTC_before_tax = PTC_credit / (1 - tax_rate)
                        ptc_gen = elec_gen * (params['PTC credit value'] * bonus_multiplier) / (1 - params['Tax Rate'])

                sum_ptc += ptc_gen / ((1 + discount_rate)**i)
                sum_elec += elec_gen / ((1 + discount_rate)**i)

            # Levelized PTC = PV(PTC revenue) / PV(electricity produced)
            # This gives the effective $/MWh reduction in LCOE due to the PTC
            estimated_ptc = sum_ptc / sum_elec
            df.loc[df['Account'] == 'LCOE with PTC', estimated_cost_col] = lcoe - estimated_ptc

        # -----------------------------------------------------------------------------------------
        # ITC (Investment Tax Credit) LCOE adjustment
        # -----------------------------------------------------------------------------------------
        # The ITC reduces the capital cost (OCC) upfront — already computed in calculate_TCI.
        # Here, the LCOE is recalculated using the reduced capital cost (TCI with ITC)
        # instead of the baseline TCI. O&M and fuel costs are unchanged.
        #
        # This produces two additional output metrics:
        #   - LCOE with ITC      : full LCOE using reduced capital cost
        #   - LCOE_cap_withitc   : capital-only component of LCOE with ITC
        # -----------------------------------------------------------------------------------------
        if 'ITC credit level' in params.keys():
            # Use the ITC-adjusted capital cost computed in calculate_TCI
            cap_cost = df.loc[df['Account'] == 'TCI (ITC-adjusted)', estimated_cost_col].values[0]
            ann_cost = df.loc[df['Account'] == 70, estimated_cost_col].values[0] + df.loc[df['Account'] == 80, estimated_cost_col].values[0]
            levelized_ann_cost = ann_cost / params['Annual Electricity Production']
            df.loc[df['Account'] == 'AC', estimated_cost_col] = ann_cost
            df.loc[df['Account'] == 'AC per MWh', estimated_cost_col] = levelized_ann_cost
            sum_cost = 0
            sum_elec = 0

            for i in range(1 + plant_lifetime_years):
                if i == 0:
                    # Year 0: ITC-adjusted capital cost is paid upfront
                    cap_cost_per_year = cap_cost
                    annual_cost = 0
                    elec_gen = 0
                elif i > 0:
                    cap_cost_per_year = 0
                    annual_cost = ann_cost
                    elec_gen = power_MWe * capacity_factor * 365 * 24  # MWh/year
                sum_cost += (cap_cost_per_year + annual_cost) / ((1 + discount_rate)**i)
                sum_elec += elec_gen / ((1 + discount_rate)**i)

            lcoe = sum_cost / sum_elec
            df.loc[df['Account'] == 'LCOE (ITC-adjusted)', estimated_cost_col] = lcoe

    return df
