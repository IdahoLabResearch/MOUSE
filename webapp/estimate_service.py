# Copyright 2025 Battelle Energy Alliance, LLC
# Released under the MIT License.
"""App-specific estimate execution for the MOUSE Streamlit UI.

This module intentionally stays free of Streamlit calls. The Streamlit app owns
UI state, cache decorators, and memory monitoring; this service owns the cost
engine run path used by the app.
"""

from __future__ import annotations

from contextlib import redirect_stdout
from dataclasses import dataclass
import io
from typing import Any, Optional, Tuple

import pandas as pd

from reactor_config import ESCALATION_YEAR, build_params
from cost.cost_drivers import cost_drivers_estimate
from cost.cost_estimation import bottom_up_cost_estimate, transform_dataframe


@dataclass(frozen=True)
class EstimateInputs:
    reactor_type: str
    power_mwt: float
    enrichment: float
    interest_rate: float
    discount_rate: float
    construction_duration: float
    debt_to_equity: float
    operation_mode: str
    emergency_shutdowns: float
    startup_duration: float
    startup_duration_refueling: float
    tax_credit_type: str
    tax_credit_value: Optional[float]
    plant_lifetime: float
    n_rings_per_assembly: Optional[int] = None
    active_height: Optional[float] = None
    n_assembly_rings: Optional[int] = None
    n_core_rings: Optional[int] = None
    tax_credit_units: Optional[int] = None
    noak_unit_number: int = 100
    licensing_learning_rate: float = 0.20


@dataclass
class EstimateResult:
    display_df: pd.DataFrame
    enriched_df: pd.DataFrame
    detailed_sorted_df: pd.DataFrame
    params: dict


@dataclass(frozen=True)
class LcoeAtNoakInputs(EstimateInputs):
    noak_unit_number: int = 10
    include_diagnostics: bool = False


@dataclass
class LcoeAtNoakResult:
    mean: float
    std: float
    diag_df: Optional[pd.DataFrame]
    diag_params: dict


def _base_overrides(inputs: EstimateInputs) -> dict:
    overrides = {
        'Interest Rate': inputs.interest_rate,
        'Discount Rate': inputs.discount_rate,
        'Construction Duration': inputs.construction_duration,
        'Debt To Equity Ratio': inputs.debt_to_equity,
        'Escalation Year': ESCALATION_YEAR,
        'Operation Mode': inputs.operation_mode,
        'Emergency Shutdowns Per Year': inputs.emergency_shutdowns,
        'Startup Duration after Emergency Shutdown': inputs.startup_duration,
        'Startup Duration after Refueling': inputs.startup_duration_refueling,
        'Levelization Period': inputs.plant_lifetime,
        'NOAK Unit Number': int(inputs.noak_unit_number),
        'Licensing Learning': float(inputs.licensing_learning_rate),
    }
    if inputs.tax_credit_type == 'PTC':
        overrides['PTC credit value'] = inputs.tax_credit_value
        overrides['PTC credit period'] = 10
        overrides['Tax Rate'] = 0.21
    elif inputs.tax_credit_type == 'ITC':
        overrides['ITC credit level'] = inputs.tax_credit_value
    if inputs.tax_credit_type in ('PTC', 'ITC') and inputs.tax_credit_units is not None:
        overrides['Number of Units Claiming ITC/PTC'] = int(inputs.tax_credit_units)
    return overrides


def _build_app_params(inputs: EstimateInputs, overrides: dict) -> dict:
    return build_params(
        inputs.reactor_type,
        inputs.power_mwt,
        inputs.enrichment,
        overrides,
        n_rings_per_assembly=inputs.n_rings_per_assembly,
        active_height=inputs.active_height,
        n_assembly_rings=inputs.n_assembly_rings,
        n_core_rings=inputs.n_core_rings,
    )


def run_estimate(inputs: EstimateInputs) -> EstimateResult:
    """Run the app's full cost-estimate path for one committed input set."""
    params = _build_app_params(inputs, _base_overrides(inputs))

    # Silence the cost engine's per-account print spam ("For the cost of the
    # Account ..."). Those lines flood Streamlit Cloud's log panel and make
    # actual diagnostics like [mem] hard to find.
    with redirect_stdout(io.StringIO()):
        raw_df = bottom_up_cost_estimate('cost/Cost_Database.xlsx', params)
        enriched_df, detailed_sorted_df = cost_drivers_estimate(raw_df, params)

    return EstimateResult(
        display_df=transform_dataframe(enriched_df),
        enriched_df=enriched_df,
        detailed_sorted_df=detailed_sorted_df,
        params=params,
    )


def _get_mean_std(df: pd.DataFrame, account, which='FOAK') -> Tuple[float, float]:
    mask = df['Account'].astype(str).str.strip() == str(account)
    mean_col = next(
        (c for c in df.columns
         if c.startswith(f'{which} Estimated Cost (') and 'std' not in c),
        None,
    )
    std_col = next(
        (c for c in df.columns if f'{which} Estimated Cost std' in c),
        None,
    )
    if mean_col is None or not mask.any():
        return float('nan'), float('nan')
    mean_vals = pd.to_numeric(df.loc[mask, mean_col], errors='coerce').dropna()
    if mean_vals.empty:
        return float('nan'), float('nan')
    mean = float(mean_vals.iloc[0])
    if std_col is None:
        return mean, 0.0
    std_vals = pd.to_numeric(df.loc[mask, std_col], errors='coerce').fillna(0)
    std = float(std_vals.iloc[0]) if not std_vals.empty else 0.0
    return mean, std


def run_lcoe_at_noak_unit(inputs: LcoeAtNoakInputs) -> LcoeAtNoakResult:
    """Run one NOAK-unit anchor for the Costs-in-Perspective plot."""
    overrides = _base_overrides(inputs)
    overrides['NOAK Unit Number'] = int(inputs.noak_unit_number)

    if inputs.tax_credit_type == 'PTC':
        lcoe_account = 'LCOE with PTC'
    elif inputs.tax_credit_type == 'ITC':
        lcoe_account = 'LCOE (ITC-adjusted)'
    else:
        lcoe_account = 'LCOE'

    params = _build_app_params(inputs, overrides)
    with redirect_stdout(io.StringIO()):
        raw_df = bottom_up_cost_estimate('cost/Cost_Database.xlsx', params)

    # Read mean / std directly from the cost engine's raw output.
    # cost_drivers_estimate enriches with per-account LCOE columns we don't use
    # here, and transform_dataframe only int-truncates floats, so both are
    # skipped to save one full enrichment pass per anchor.
    mean, std = _get_mean_std(raw_df, lcoe_account, 'NOAK')

    diag_params = {
        'NOAK Unit Number': params.get('NOAK Unit Number'),
        'Power MWt': params.get('Power MWt'),
        'Power MWe': params.get('Power MWe'),
        'Thermal Efficiency': params.get('Thermal Efficiency'),
        'Capacity Factor': params.get('Capacity Factor'),
        'Construction Duration': params.get('Construction Duration'),
        'Interest Rate': params.get('Interest Rate'),
        'Discount Rate': params.get('Discount Rate'),
        'Debt To Equity Ratio': params.get('Debt To Equity Ratio'),
        'Levelization Period': params.get('Levelization Period'),
        'Annual Electricity Production': params.get('Annual Electricity Production'),
    }

    diag_df = None
    if inputs.include_diagnostics:
        _foak_cols = [c for c in raw_df.columns if c.startswith('FOAK Estimated Cost (')]
        _noak_cols = [c for c in raw_df.columns if c.startswith('NOAK Estimated Cost (')]
        _keep = [c for c in (['Account', 'Account Title'] + _foak_cols + _noak_cols)
                 if c in raw_df.columns]
        diag_df = raw_df[_keep].copy() if _keep else None

    return LcoeAtNoakResult(
        mean=float(mean) if mean == mean else float('nan'),
        std=float(std) if std == std else 0.0,
        diag_df=diag_df,
        diag_params=diag_params,
    )
