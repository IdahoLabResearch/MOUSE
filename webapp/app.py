# Copyright 2025 Battelle Energy Alliance, LLC
# Released under the MIT License.
"""
MOUSE Streamlit Web App
Microreactor Online Unified Simulation Engine cost estimation without OpenMC.

Run from the MOUSE repo root:
    streamlit run webapp/app.py
"""

# ---------------------------------------------------------------------------
# Ensure the MOUSE repo root is in sys.path so all MOUSE modules resolve.
# ---------------------------------------------------------------------------
import sys
import os

_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

# ---------------------------------------------------------------------------
# Memory tracing (disabled). tracemalloc was enabled briefly to identify
# the leak source; the diagnostic dump showed Python-side allocations are
# minor (~60 MB of ~688 MB total) — the real memory lives in C-extension
# buffers (numpy/pandas/matplotlib/openmc) that tracemalloc can't see.
# Even at depth=1 it adds measurable per-allocation overhead on the
# pandas-heavy cost engine, so it stays off in production. The diagnostic
# helper still works (it checks _tracemalloc.is_tracing() and falls
# back to gc-based counts if tracing is off).
# ---------------------------------------------------------------------------
import tracemalloc as _tracemalloc

# ---------------------------------------------------------------------------
# IMPORTANT: Stub openmc and watts BEFORE any MOUSE import.
# ---------------------------------------------------------------------------
from unittest.mock import MagicMock


class _MaterialStub:
    def __init__(self, name=None, temperature=None):
        self.name = name
        self.temperature = temperature
        self.density = 0.0

    def set_density(self, units, value):
        self.density = value

    def add_nuclide(self, *args, **kwargs):
        pass

    def add_element(self, *args, **kwargs):
        pass

    def add_s_alpha_beta(self, *args, **kwargs):
        pass

    @staticmethod
    def mix_materials(materials, fractions, method, name=None):
        result = _MaterialStub(name=name)
        try:
            result.density = sum(m.density * f for m, f in zip(materials, fractions))
        except Exception:
            result.density = 0.0
        return result


class _MaterialsStub:
    def append(self, mat):
        pass

    def extend(self, mats):
        pass


class _OpenMCStub(MagicMock):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Material = _MaterialStub
        self.Materials = _MaterialsStub


_openmc_stub = _OpenMCStub()
_openmc_stub.Material = _MaterialStub
_openmc_stub.Materials = _MaterialsStub

for _mod in ['openmc', 'openmc.deplete', 'openmc.mgxs']:
    sys.modules[_mod] = _openmc_stub

sys.modules['watts'] = MagicMock()

# matplotlib is imported at module level by core_design/utils.py and
# core_design/correction_factor.py for plot helpers (create_universe_plot,
# keff_comparison_vs_Time, etc.) that the webapp never calls. Stub it so
# those imports succeed without installing the heavy matplotlib package.
_mpl_stub = MagicMock()
for _mod in ['matplotlib', 'matplotlib.pyplot', 'matplotlib.patches',
             'matplotlib.colors']:
    sys.modules[_mod] = _mpl_stub

# ---------------------------------------------------------------------------
# Cache the OpenMC materials database build. collect_materials_data is
# called dozens of times per Run (drums.py alone calls it 7 times per
# cost-engine run, and we have ~8 cost-engine runs).
#
# We replace the openmc-based collect_materials_data with a flat JSON
# lookup. Audit showed only drums.py reads materials_database[X].density;
# every other materials_database[X] access in core_design/openmc_template_*
# is `cell.fill = mat`, a no-op because openmc is stubbed. So a thin
# object with just .name and .density satisfies all callers.
#
# Patch must run BEFORE reactor_config / core_design templates are
# imported, so their `from ... import collect_materials_data` binds to
# the JSON version.
# ---------------------------------------------------------------------------
import json as _json
import core_design.openmc_materials_database as _materials_db_mod


class _ThinMaterial:
    """Minimal stand-in for openmc.Material. drums.py reads .density;
    openmc_template_* assigns these to .fill on stubbed cells (no-op)."""
    __slots__ = ('name', 'density')
    def __init__(self, name, density):
        self.name = name
        self.density = density


_MATERIALS_JSON_PATH = os.path.join(os.path.dirname(__file__),
                                    'materials_densities.json')
with open(_MATERIALS_JSON_PATH) as _f:
    _MATERIALS_RAW = _json.load(_f)

# Pre-build per-reactor lookup dicts at module load. Total objects ~80
# (25-26 per reactor x 3 reactors), each __slots__ wrapper is ~80 bytes,
# so ~6 KB of permanent footprint. Replaces the ~MB-scale per-Run cost
# of repeatedly building openmc.Material objects (each carrying nuclides,
# add_element call records, etc.).
_THIN_MATERIALS_BY_REACTOR = {
    rtype: {name: _ThinMaterial(name, float(density))
            for name, density in mats.items()}
    for rtype, mats in _MATERIALS_RAW.items()
}


def _cached_collect_materials_data(params):
    rtype = params.get('reactor type', 'LTMR')
    return _THIN_MATERIALS_BY_REACTOR.get(rtype,
                                          _THIN_MATERIALS_BY_REACTOR['LTMR'])


_materials_db_mod.collect_materials_data = _cached_collect_materials_data

# ---------------------------------------------------------------------------
# Standard imports (after stubs are in place)
# ---------------------------------------------------------------------------
import base64
import html
import io
import math
import sqlite3
import uuid
import warnings
from datetime import datetime, timezone
from pathlib import Path

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import altair as alt
import streamlit as st
# streamlit_analytics2 removed (memory: ~10-25 MB freed). Replace the
# library with a no-op context manager so the existing
# `with streamlit_analytics.track():` block keeps working without
# requiring a deep re-indent of the entire script body. The custom
# SQLite-backed visit logger below (`_log_visit_once_per_session`) is
# kept; it's lightweight and gives us the analytics we actually use.
from contextlib import nullcontext as _nullcontext

class _NoOpAnalytics:
    @staticmethod
    def track():
        return _nullcontext()

streamlit_analytics = _NoOpAnalytics()
from st_cookies_manager import EncryptedCookieManager

from reactor_config import SubcriticalError, ShortLifetimeError, ESCALATION_YEAR
from webapp.fuel_lifetime_estimator import (
    get_ltmr_keff_curve,
    get_ltmr_peaking_factor,
    get_ltmr_leakage,
)
from webapp.gcmr_fuel_lifetime_estimator import (
    get_gcmr_peaking_factor,
    get_gcmr_leakage,
    get_gcmr_keff_curve,
)
from webapp.hpmr_fuel_lifetime_estimator import (
    get_hpmr_keff_curve,
    get_hpmr_peaking_factor,
    get_hpmr_leakage,
)
from webapp.estimate_service import (
    EstimateInputs,
    LcoeAtNoakInputs,
    run_estimate,
    run_lcoe_at_noak_unit,
)

# ---------------------------------------------------------------------------
# Performance patches: cache Excel reads that would otherwise repeat on every run.
# ---------------------------------------------------------------------------
import functools
import cost.cost_escalation as _ce

_orig_inflation = _ce.calculate_inflation_multiplier

@functools.lru_cache(maxsize=128)
def _cached_inflation_multiplier(file_path, base_dollar_year, cost_type, escalation_year):
    return _orig_inflation(file_path, base_dollar_year, cost_type, escalation_year)

_ce.calculate_inflation_multiplier = _cached_inflation_multiplier

import pandas as _pd_orig
_orig_read_excel = _pd_orig.read_excel

@functools.lru_cache(maxsize=8)
def _cached_read_excel(file_path, sheet_name):
    return _orig_read_excel(file_path, sheet_name=sheet_name)

def _patched_read_excel(file_path, sheet_name=0, **kwargs):
    if isinstance(file_path, str) and file_path.endswith('Cost_Database.xlsx') and not kwargs:
        return _cached_read_excel(file_path, sheet_name).copy()
    return _orig_read_excel(file_path, sheet_name=sheet_name, **kwargs)

_pd_orig.read_excel = _patched_read_excel
_ce.pd.read_excel = _patched_read_excel

# ---------------------------------------------------------------------------
# Memory monitor
# ---------------------------------------------------------------------------
# Called once per rerun, right after the sidebar badge renders, BEFORE any
# st.stop() upstream can short-circuit the script. Two tiers:
#
#   Every rerun (cheap, no UX cost):
#     - gc.collect() reclaims the rerun's transient DataFrames
#     - plt.close('all') drops matplotlib figure references after they've
#       been serialized to PNG by st.pyplot (the single biggest non-cache
#       source of RAM growth)
#     - malloc_trim asks glibc to release freed pages back to the OS
#       (no-op on macOS/Windows; load-bearing on Streamlit Cloud Linux)
#
#   At 600 MB (sweep — clears every cache we own + dumps diagnostics):
#     - st.cache_data (cost engine + anchors)
#     - functools.lru_cache instances (_cached_inflation_multiplier,
#       _cached_read_excel) — st.cache_data.clear() does NOT touch these
#     (the former bounded OrderedDict _materials_cache is gone now that
#     collect_materials_data is a flat JSON lookup with no runtime build)
#     - Logs tracemalloc top 10 allocation sites, gc top object types,
#       and session_state size so we can identify what's leaking past
#       the cache sweep. Threshold is intentionally low (600 vs the
#       previous 800) to give ~350 MB of headroom before the friendly
#       stop, and to surface diagnostics earlier in the leak cycle.
#     Next Run rebuilds cold; this is the cost-warmth-for-survival trade.
#
#   At 950 MB (friendly stop):
#     - Render an actionable banner and st.stop() the rerun cleanly.
#       Avoids Streamlit Cloud's OOM "Oh no" page that requires manual
#       reboot. Do NOT use os._exit() — Streamlit Cloud treats it as a
#       crash, not a restart.
def _log_memory_diagnostics(rss_mb):
    """Dumps targeted diagnostics to find the leak source:
      - top gc object types by count (general orientation)
      - top 5 functools.partial .func qualnames (what kind of partials)
      - top 5 function qualnames (what kind of closures)
      - session_state DEEP size with per-key breakdown for DataFrames
        and ndarrays (shallow getsizeof misses C-allocated buffers)
      - Streamlit cache_data storage size if discoverable
    Output goes to stdout (Streamlit Cloud logs)."""
    print(f"[mem-diag] post-sweep RSS={rss_mb:.0f} MB — diagnostics follow",
          flush=True)

    try:
        import gc as _gc
        from collections import Counter as _Counter
        objs = _gc.get_objects()
        type_counts = _Counter(type(o).__name__ for o in objs)
        print("[mem-diag] top 10 gc object types by count:", flush=True)
        for name, count in type_counts.most_common(10):
            print(f"[mem-diag]   {count:>8}  {name}", flush=True)
    except Exception as e:
        print(f"[mem-diag] gc count failed: {e}", flush=True)
        objs = []

    try:
        import functools as _ft
        partials = [o for o in objs if isinstance(o, _ft.partial)]
        partial_funcs = _Counter(
            getattr(p.func, '__qualname__', None)
            or getattr(p.func, '__name__', None)
            or type(p.func).__name__
            for p in partials
        )
        print(f"[mem-diag] {len(partials)} functools.partial — "
              f"top 5 by .func:", flush=True)
        for name, count in partial_funcs.most_common(5):
            print(f"[mem-diag]   {count:>5}  {name}", flush=True)
    except Exception as e:
        print(f"[mem-diag] partial breakdown failed: {e}", flush=True)

    try:
        import types as _types
        funcs = [o for o in objs if isinstance(o, _types.FunctionType)]
        func_quals = _Counter(getattr(f, '__qualname__', '?') for f in funcs)
        print(f"[mem-diag] {len(funcs)} functions — top 5 by qualname:",
              flush=True)
        for name, count in func_quals.most_common(5):
            print(f"[mem-diag]   {count:>5}  {name}", flush=True)
    except Exception as e:
        print(f"[mem-diag] function breakdown failed: {e}", flush=True)

    try:
        import sys as _sys
        import pandas as _pd_diag
        import numpy as _np_diag
        keys = list(st.session_state.keys())
        breakdown = []
        for k in keys:
            try:
                v = st.session_state[k]
                if isinstance(v, _pd_diag.DataFrame):
                    sz = int(v.memory_usage(deep=True).sum())
                elif isinstance(v, _np_diag.ndarray):
                    sz = int(v.nbytes)
                elif isinstance(v, dict):
                    sz = _sys.getsizeof(v) + sum(
                        _sys.getsizeof(x) for x in v.values()
                    )
                else:
                    sz = _sys.getsizeof(v)
                breakdown.append((str(k), sz, type(v).__name__))
            except Exception:
                continue
        breakdown.sort(key=lambda x: x[1], reverse=True)
        total = sum(sz for _, sz, _ in breakdown)
        print(f"[mem-diag] session_state deep: {len(keys)} keys, "
              f"~{total/1024/1024:.2f} MB total", flush=True)
        for k, sz, tname in breakdown[:8]:
            print(f"[mem-diag]   {sz/1024/1024:6.2f} MB  "
                  f"{tname:<14}  {k}", flush=True)
    except Exception as e:
        print(f"[mem-diag] session_state deep size failed: {e}", flush=True)

    try:
        from streamlit.runtime.caching import cache_data_api as _cd
        store = getattr(_cd, 'CACHE_DATA_MESSAGE_REPLAY_CTX', None)
        registry = getattr(_cd, '_data_caches', None)
        if registry is not None:
            caches = getattr(registry, '_caches', {})
            print(f"[mem-diag] streamlit cache_data: {len(caches)} caches",
                  flush=True)
            for cname, c in list(caches.items())[:5]:
                inner = getattr(c, '_mem_cache', None) or getattr(c, 'cache', None)
                if inner is not None:
                    try:
                        size = len(inner)
                    except Exception:
                        size = '?'
                    print(f"[mem-diag]   {size} entries  {cname}", flush=True)
    except Exception as e:
        print(f"[mem-diag] streamlit cache inspection failed: {e}", flush=True)


def _run_memory_monitor():
    import gc as _gc
    _gc.collect()

    try:
        import ctypes as _ctypes
        _ctypes.CDLL("libc.so.6").malloc_trim(0)
    except (OSError, AttributeError):
        pass

    try:
        import psutil as _psutil
        _rss_mb = _psutil.Process().memory_info().rss / (1024 * 1024)
    except ImportError:
        return

    if _rss_mb > 600:
        st.cache_data.clear()
        try:
            _cached_inflation_multiplier.cache_clear()
        except Exception:
            pass
        try:
            _cached_read_excel.cache_clear()
        except Exception:
            pass
        # Triple gc pass breaks reference cycles that a single collect
        # leaves untouched. (Was originally added for matplotlib Figure
        # <-> Axes <-> Artist mutual refs; now matplotlib is fully gone
        # but the practice of multiple passes is cheap and still helps
        # with pandas/numpy intermediate cycles.)
        _gc.collect()
        _gc.collect()
        _gc.collect()
        try:
            import ctypes as _ctypes
            _ctypes.CDLL("libc.so.6").malloc_trim(0)
        except (OSError, AttributeError):
            pass
        _new_rss = _psutil.Process().memory_info().rss / (1024 * 1024)
        print(f"[mem] {_rss_mb:.0f} MB > 600 MB -> cleared caches, "
              f"now {_new_rss:.0f} MB", flush=True)
        _log_memory_diagnostics(_new_rss)

    if _rss_mb > 950:
        print(f"[mem] {_rss_mb:.0f} MB > 950 MB -> friendly stop",
              flush=True)
        st.error(
            "**Memory limit reached.** The app is approaching its "
            "1 GB memory ceiling. Please reload the page to start a "
            "fresh session. If this keeps happening, the app may "
            "have many concurrent users  try again in a few minutes."
        )
        st.stop()


# ---------------------------------------------------------------------------
# Reactor metadata: full names and design images
# ---------------------------------------------------------------------------
_REACTOR_LABELS = {
    'LTMR': 'Liquid Metal Microreactor (LTMR)',
    'GCMR': 'Gas Cooled Microreactor (GCMR)',
    'HPMR': 'Heat Pipe Microreactor (HPMR)',
}
_LABEL_TO_KEY = {v: k for k, v in _REACTOR_LABELS.items()}

_ASSETS = os.path.join(_repo_root, 'assets', 'Ref_openmc_2d_designs')

_REACTOR_IMAGES = {
    'LTMR': {
        'main': (
            os.path.join(_ASSETS, 'LTMR_core.png'),
            'Simplified scoping-level design used for cost estimation, not a '
            'detailed engineering design. LTMR core cross-section: hexagonal '
            'arrangement of UZrH alloy fuel pins and ZrH moderator pins cooled '
            'by NaK liquid metal, surrounded by a graphite radial reflector '
            'with control drums.',
        ),
        'details': [
            (
                os.path.join(_ASSETS, 'LTMR_fuel_pin_universe.png'),
                'LTMR fuel pin cross-section from center outward: zirconium cladding, '
                'gap, U-ZrH fuel meat, gap, and SS304 outer cladding.',
            ),
            (
                os.path.join(_ASSETS, 'LTMR_moderator_pin_universe.png'),
                'LTMR moderator pin cross-section ZrH hydrogen moderator encased in '
                'SS304 cladding, interspersed between fuel pins to thermalize neutrons.',
            ),
        ],
    },
    'GCMR': {
        # The main image is overridden at runtime with the per-(N_A, N_C)
        # core PNG; this static entry is the fallback when the per-pair
        # image isn't available.
        'main': (
            os.path.join(_ASSETS, 'GCMR_Core.png'),
            'GCMR core cross-section hexagonal fuel assemblies containing TRISO fuel '
            'compacts arranged in a honeycomb pattern, cooled by helium gas flowing '
            'through dedicated coolant channels, with graphite reflector and control drums.',
        ),
        # The fuel assembly entry below is overridden at runtime with the
        # per-N_A image. The TRISO particle and zoomed fuel-assembly
        # images don't depend on geometry kept static.
        'details': [
            (
                os.path.join(_ASSETS, 'GCMR_Fuel Assembly.png'),
                'GCMR fuel assembly cross-section TRISO fuel compacts, helium coolant '
                'channels, and ZrH moderator booster pins embedded in a graphite matrix.',
            ),
            (
                os.path.join(_ASSETS, 'GCMR_fuel_assembly_zoomed.png'),
                'GCMR fuel assembly zoomed individual TRISO particles visible within '
                'the graphite fuel compact at the target packing fraction.',
            ),
            (
                os.path.join(_ASSETS, 'GCMR_TRISO_particle.png'),
                'TRISO fuel particle multi-layer design with UN fuel kernel, buffer '
                'graphite, inner PyC, SiC pressure vessel, and outer PyC coating that '
                'retains fission products up to ~1600 °C.',
            ),
        ],
    },
    'HPMR': {
        'main': (
            os.path.join(_ASSETS, 'HPMR_core.png'),
            'HPMR core cross-section monolithic graphite/metal core with hexagonal '
            'fuel assemblies and embedded alkali-metal heat pipes that passively transfer '
            'heat to the secondary side, with graphite reflector and control drums.',
        ),
        'details': [
            (
                os.path.join(_ASSETS, 'HPMR_fuel_assembly.png'),
                'HPMR fuel assembly cross-section TRISO fuel pins and heat pipes '
                'arranged in a hexagonal pattern within the graphite monolith block.',
            ),
            (
                os.path.join(_ASSETS, 'HPMR_fuel_pin_universe.png'),
                'HPMR fuel pin cross-section homogenized TRISO fuel region surrounded '
                'by a thin helium gap within the monolith.',
            ),
            (
                os.path.join(_ASSETS, 'HPMR_heatpipe_universe.png'),
                'HPMR heat pipe cross-section working fluid region and outer cladding '
                'that passively carries heat from the core to the power conversion system '
                'with no moving parts.',
            ),
        ],
    },
}

# ---------------------------------------------------------------------------
# Cached cost estimate (module-level so cache persists across reruns)
# ---------------------------------------------------------------------------
# Aspect-ratio bounds for the Active Height slider:
# 0.5 ≤ H / D ≤ 2.0
# where D is the active core diameter (= 2 × active_radius, no reflector).
# H is the active fuel height. Reflectors are excluded from both.
ASPECT_RATIO_MIN = 0.5
ASPECT_RATIO_MAX = 2.0

# LTMR per-N geometry:
# - ACTIVE_R = Core Radius − Reflector Thickness (≈ 2.836·N − 1.136 cm)
# - DIAMETER_CM = 2 × ACTIVE_R (active core diameter, NO reflector)
# Trained N values come directly from the parametric study Excel; intermediate
# N values are linearly interpolated from the trained table adequate for
# slider display, but flagged in the UI as interpolated.
LTMR_TRAINED_N = {10, 12, 14, 18, 24}
LTMR_N_TO_ACTIVE_RADIUS_CM = {
    n: round(2.836 * n - 1.136, 1)
    for n in [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
}
LTMR_N_TO_DIAMETER_CM = {
    n: int(round(2 * ar))
    for n, ar in LTMR_N_TO_ACTIVE_RADIUS_CM.items()
}


def _ltmr_diameter_label(n, d):
    star = '' if n in LTMR_TRAINED_N else ' *'
    return f"{d} cm{star}"


LTMR_DIAMETER_LABELS = [_ltmr_diameter_label(n, d)
                        for n, d in LTMR_N_TO_DIAMETER_CM.items()]
LTMR_DIAMETER_LABEL_TO_N = {label: n
                            for n, label in zip(LTMR_N_TO_DIAMETER_CM.keys(),
                                                LTMR_DIAMETER_LABELS)}

# GCMR per-(N_A, N_C) geometry. Formulas (verified against parametric study):
# Assembly_FTF(N_A) = 2.25 × (N_A − 1) × √3
# Reflector(N_A) = Assembly_FTF / 2
# Active_Radius = Assembly_FTF × N_C
# Core_Radius = Active_Radius + Reflector
# Diameter = 2 × Core_Radius
import math as _math
GCMR_TRAINED_PAIRS = {(4, 3), (5, 3), (6, 4), (6, 5), (7, 5)}
GCMR_NA_VALUES = [4, 5, 6, 7]
GCMR_NC_VALUES = [3, 4, 5]


def _gcmr_assembly_ftf(n_a):
    return 2.25 * (n_a - 1) * _math.sqrt(3.0)


def _gcmr_active_radius(n_a, n_c):
    return _gcmr_assembly_ftf(n_a) * n_c


def _gcmr_diameter_cm(n_a, n_c):
    ftf = _gcmr_assembly_ftf(n_a)
    return 2.0 * (ftf * n_c + ftf / 2.0)


# Build the 12-pair lookup of ACTIVE core diameter (no reflector),
# sorted ascending so the slider goes small → large left-to-right.
GCMR_PAIR_TO_DIAMETER_CM = {}
for _na in GCMR_NA_VALUES:
    for _nc in GCMR_NC_VALUES:
        GCMR_PAIR_TO_DIAMETER_CM[(_na, _nc)] = int(round(2 * _gcmr_active_radius(_na, _nc)))
GCMR_PAIR_TO_DIAMETER_CM = dict(sorted(GCMR_PAIR_TO_DIAMETER_CM.items(),
                                       key=lambda kv: kv[1]))


def _gcmr_diameter_label(na, nc, d):
    star = '' if (na, nc) in GCMR_TRAINED_PAIRS else ' *'
    return f"{d} cm{star}"


GCMR_DIAMETER_LABELS = [_gcmr_diameter_label(na, nc, d)
                        for (na, nc), d in GCMR_PAIR_TO_DIAMETER_CM.items()]
GCMR_DIAMETER_LABEL_TO_PAIR = {
    label: pair
    for pair, label in zip(GCMR_PAIR_TO_DIAMETER_CM.keys(),
                           GCMR_DIAMETER_LABELS)
}

# Mass-scaling constant derived from the GCMR parametric study reference
# row (N_A=6, N_C=5, H=215 cm, E=0.1975, P=1 MWt). Total uranium mass per
# unit of [F_A × F_C × H] volume index.
GCMR_G_PER_VOLUME_INDEX = 0.5776


# HPMR Active core diameter per N_C, with N_A locked to 6 (the only
# value present in the HPMR parametric study).
# Active diameter = 2 × Active radius (excludes reflector)
# Active radius = (sqrt(3)/2) × hex_edge
# hex_edge = Assembly_FTF × (N_C - 1) + Assembly_FTF/2 + 6.6
# Assembly_FTF = (Lattice_Pitch × (N_A - 1) + 1.4 × R_pin) × sqrt(3)
# At N_A=6, Lattice_Pitch=3.4, R_pin=1.05 → Assembly_FTF ≈ 31.99 cm.
HPMR_NA_FIXED = 6
_HPMR_LATTICE_PITCH = 3.4
_HPMR_FUEL_PIN_OUTER_R = 1.05


def _hpmr_assembly_ftf():
    return (_HPMR_LATTICE_PITCH * (HPMR_NA_FIXED - 1)
            + 1.4 * _HPMR_FUEL_PIN_OUTER_R) * _math.sqrt(3.0)


def _hpmr_active_radius(n_c):
    ftf = _hpmr_assembly_ftf()
    hex_edge = ftf * (n_c - 1) + ftf / 2.0 + 6.6
    return 0.5 * _math.sqrt(3.0) * hex_edge


HPMR_NC_VALUES = [3, 4, 5, 6, 7]
HPMR_NC_TO_DIAMETER_CM = {
    nc: int(round(2 * _hpmr_active_radius(nc))) for nc in HPMR_NC_VALUES
}
HPMR_DIAMETER_LABELS = [
    f"{HPMR_NC_TO_DIAMETER_CM[nc]} cm"
    for nc in HPMR_NC_VALUES
]
HPMR_DIAMETER_LABEL_TO_NC = {
    label: nc for label, nc in zip(HPMR_DIAMETER_LABELS, HPMR_NC_VALUES)
}


@st.cache_data(show_spinner=False, max_entries=2)
def _run_estimate(reactor_type, power_mwt, enrichment, interest_rate, discount_rate,
                  construction_duration, debt_to_equity, operation_mode,
                  emergency_shutdowns, startup_duration, startup_duration_refueling,
                  tax_credit_type, tax_credit_value, plant_lifetime,
                  n_rings_per_assembly=None, active_height=None,
                  n_assembly_rings=None, n_core_rings=None,
                  tax_credit_units=None):
    result = run_estimate(EstimateInputs(
        reactor_type=reactor_type,
        power_mwt=power_mwt,
        enrichment=enrichment,
        interest_rate=interest_rate,
        discount_rate=discount_rate,
        construction_duration=construction_duration,
        debt_to_equity=debt_to_equity,
        operation_mode=operation_mode,
        emergency_shutdowns=emergency_shutdowns,
        startup_duration=startup_duration,
        startup_duration_refueling=startup_duration_refueling,
        tax_credit_type=tax_credit_type,
        tax_credit_value=tax_credit_value,
        plant_lifetime=plant_lifetime,
        n_rings_per_assembly=n_rings_per_assembly,
        active_height=active_height,
        n_assembly_rings=n_assembly_rings,
        n_core_rings=n_core_rings,
        tax_credit_units=tax_credit_units,
    ))
    return (result.display_df, result.enriched_df,
            result.detailed_sorted_df, result.params)


# Single-point cost-engine call for the Costs-in-Perspective plot.
# Used to add an intermediate anchor (e.g. NOAK Unit Number = 10)
# between the headline FOAK (N=1) and NOAK (N=user_setting). Delegates
# to estimate_service so the app keeps one cost-engine execution path;
# the service still skips full enrichment for this anchor to save work.
@st.cache_data(show_spinner=False, max_entries=2)
def _lcoe_at_noak_unit(reactor_type, power_mwt, enrichment, interest_rate, discount_rate,
                       construction_duration, debt_to_equity, operation_mode,
                       emergency_shutdowns, startup_duration, startup_duration_refueling,
                       tax_credit_type, tax_credit_value, plant_lifetime,
                       n_rings_per_assembly=None, active_height=None,
                       n_assembly_rings=None, n_core_rings=None,
                       noak_unit_number=10,
                       tax_credit_units=None):
    """Returns (mean, std, per_account_df) for one NOAK Unit Number.

    The per_account_df includes Account, Account Title, FOAK and NOAK
    columns for all rows useful for diagnosing which account drives
    a wild LCOE value.
    """
    result = run_lcoe_at_noak_unit(LcoeAtNoakInputs(
        reactor_type=reactor_type,
        power_mwt=power_mwt,
        enrichment=enrichment,
        interest_rate=interest_rate,
        discount_rate=discount_rate,
        construction_duration=construction_duration,
        debt_to_equity=debt_to_equity,
        operation_mode=operation_mode,
        emergency_shutdowns=emergency_shutdowns,
        startup_duration=startup_duration,
        startup_duration_refueling=startup_duration_refueling,
        tax_credit_type=tax_credit_type,
        tax_credit_value=tax_credit_value,
        plant_lifetime=plant_lifetime,
        n_rings_per_assembly=n_rings_per_assembly,
        active_height=active_height,
        n_assembly_rings=n_assembly_rings,
        n_core_rings=n_core_rings,
        tax_credit_units=tax_credit_units,
        noak_unit_number=noak_unit_number,
    ))
    return result.mean, result.std, result.diag_df, result.diag_params


# ---------------------------------------------------------------------------
# Anonymous visitor analytics
# ---------------------------------------------------------------------------
_ANALYTICS_DB_PATH = Path(_repo_root) / 'webapp' / 'analytics.db'
_ANALYTICS_COOKIE_NAME = 'anonymous_id'


def _get_cookie_manager():
    return EncryptedCookieManager(
        prefix='mouse_',
        password=st.secrets.get('COOKIE_PASSWORD', 'change-this-in-streamlit-secrets'),
    )


def _get_analytics_conn():
    conn = sqlite3.connect(_ANALYTICS_DB_PATH, check_same_thread=False)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS visits (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            anonymous_id TEXT NOT NULL,
            session_id TEXT NOT NULL,
            first_seen_utc TEXT NOT NULL,
            visit_time_utc TEXT NOT NULL,
            page_name TEXT,
            reactor_type TEXT,
            user_agent TEXT,
            language TEXT
        )
    """)
    conn.commit()
    return conn


def _get_or_create_anonymous_id(cookies):
    anonymous_id = cookies.get(_ANALYTICS_COOKIE_NAME)
    if not anonymous_id:
        anonymous_id = str(uuid.uuid4())
        cookies[_ANALYTICS_COOKIE_NAME] = anonymous_id
        cookies.save()
    return anonymous_id


def _get_first_seen(conn, anonymous_id, now_utc):
    row = conn.execute(
        "SELECT MIN(first_seen_utc) FROM visits WHERE anonymous_id = ?",
        (anonymous_id,),
    ).fetchone()
    return row[0] if row and row[0] else now_utc


def _log_visit_once_per_session(conn, anonymous_id, reactor_type=None, page_name='main'):
    if 'analytics_session_id' not in st.session_state:
        st.session_state.analytics_session_id = str(uuid.uuid4())

    if 'analytics_visit_logged' not in st.session_state:
        st.session_state.analytics_visit_logged = False

    if st.session_state.analytics_visit_logged:
        return

    headers = dict(st.context.headers) if hasattr(st.context, 'headers') else {}
    user_agent = headers.get('User-Agent', '')
    language = headers.get('Accept-Language', '')
    now_utc = datetime.now(timezone.utc).isoformat()
    first_seen_utc = _get_first_seen(conn, anonymous_id, now_utc)

    conn.execute(
        """
        INSERT INTO visits (
            anonymous_id,
            session_id,
            first_seen_utc,
            visit_time_utc,
            page_name,
            reactor_type,
            user_agent,
            language
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            anonymous_id,
            st.session_state.analytics_session_id,
            first_seen_utc,
            now_utc,
            page_name,
            reactor_type,
            user_agent,
            language,
        ),
    )
    conn.commit()
    st.session_state.analytics_visit_logged = True


def _get_analytics_summary(conn):
    total_visits = conn.execute("SELECT COUNT(*) FROM visits").fetchone()[0]
    unique_visitors = conn.execute(
        "SELECT COUNT(DISTINCT anonymous_id) FROM visits"
    ).fetchone()[0]
    returning_visitors = conn.execute(
        """
        SELECT COUNT(*)
        FROM (
            SELECT anonymous_id
            FROM visits
            GROUP BY anonymous_id
            HAVING COUNT(*) > 1
        )
        """
    ).fetchone()[0]

    visits_per_person = pd.read_sql_query(
        """
        SELECT
            anonymous_id,
            COUNT(*) AS visit_count,
            MIN(first_seen_utc) AS first_seen_utc,
            MAX(visit_time_utc) AS last_visit_utc,
            MAX(COALESCE(reactor_type, '')) AS last_reactor_type
        FROM visits
        GROUP BY anonymous_id
        ORDER BY last_visit_utc DESC
        """,
        conn,
    )

    recent_visits = pd.read_sql_query(
        """
        SELECT
            anonymous_id,
            visit_time_utc,
            reactor_type,
            page_name
        FROM visits
        ORDER BY visit_time_utc DESC
        LIMIT 25
        """,
        conn,
    )

    return total_visits, unique_visitors, returning_visitors, visits_per_person, recent_visits


def _render_analytics_sidebar(conn):
    total_visits, unique_visitors, returning_visitors, visits_per_person, recent_visits = _get_analytics_summary(conn)

    with st.sidebar.expander('📈 Anonymous Visitor Analytics'):
        c1, c2 = st.columns(2)
        c1.metric('Total visits', total_visits)
        c2.metric('Unique visitors', unique_visitors)

        c3, c4 = st.columns(2)
        c3.metric('Returning', returning_visitors)
        c4.metric(
            'Avg visits/user',
            f"{(total_visits / unique_visitors):.2f}" if unique_visitors else '0.00'
        )

        st.caption('Visits per anonymous visitor')
        st.dataframe(
            visits_per_person,
            width='stretch',
            hide_index=True,
            height=220,
        )

        st.caption('Most recent visits')
        st.dataframe(
            recent_visits,
            width='stretch',
            hide_index=True,
            height=220,
        )

# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------
def _get_mean_std(df, account, which='FOAK'):
    row = df[df['Account'] == account]
    if row.empty:
        return float('nan'), float('nan')
    mean_prefix = 'FOAK Estimated Cost (' if which == 'FOAK' else 'NOAK Estimated Cost ('
    std_prefix = 'FOAK Estimated Cost std (' if which == 'FOAK' else 'NOAK Estimated Cost std ('
    mean_cols = [c for c in df.columns if c.startswith(mean_prefix)]
    std_cols = [c for c in df.columns if c.startswith(std_prefix)]
    try:
        mean = float(row[mean_cols[0]].iloc[0]) if mean_cols else float('nan')
    except (TypeError, ValueError):
        mean = float('nan')
    try:
        std = float(row[std_cols[0]].iloc[0]) if std_cols else float('nan')
    except (TypeError, ValueError):
        std = float('nan')
    return mean, std


def _fmt_cost(mean, std):
    if math.isnan(mean):
        return 'N/A'
    # Sub-$1M values (e.g. training costs) round to 1 decimal so
    # they read as "$0.3M" instead of "$0M". Values >= $1M stay as
    # whole millions to keep OCC/TCI/Direct cards uncluttered.
    use_decimal = abs(mean) < 1e6
    fmt = (lambda v: f'${v / 1e6:.1f}M') if use_decimal else (lambda v: f'${round(v / 1e6)}M')
    if math.isnan(std) or std == 0:
        return fmt(mean)
    # "Mean [low – high]" the mean is the headline, the bracketed
    # range (mean ± 1σ) is rendered in lighter slate so the eye anchors
    # on the central number first.
    return (
        f'{fmt(mean)} '
        f'<span style="color:#64748b;font-weight:400;">'
        f'[{fmt(mean - std)} &ndash; {fmt(mean + std)}]'
        f'</span>'
    )


def _fmt_lcoe(mean, std):
    if math.isnan(mean):
        return 'N/A'
    m = int(round(mean))
    if math.isnan(std) or std == 0:
        return f'${m}/MW<sub>e</sub>h'
    lo = int(round(mean - std))
    hi = int(round(mean + std))
    return (
        f'${m}/MW<sub>e</sub>h '
        f'<span style="color:#64748b;font-weight:400;">'
        f'[${lo} &ndash; ${hi}]'
        f'</span>'
    )


def _fmt_lcoh(mean, std):
    if math.isnan(mean):
        return 'N/A'
    m = int(round(mean))
    if math.isnan(std) or std == 0:
        return f'${m}/MW<sub>t</sub>h'
    lo = int(round(mean - std))
    hi = int(round(mean + std))
    return (
        f'${m}/MW<sub>t</sub>h '
        f'<span style="color:#64748b;font-weight:400;">'
        f'[${lo} &ndash; ${hi}]'
        f'</span>'
    )


def _get_lcof(df, which='FOAK'):
    mean_col = 'FOAK LCOE' if which == 'FOAK' else 'NOAK LCOE'
    std_col = 'FOAK LCOE_std' if which == 'FOAK' else 'NOAK LCOE_std'
    if mean_col not in df.columns:
        return float('nan'), float('nan')
    mask = df['Account'].isin([25, 80])
    mean_vals = pd.to_numeric(df.loc[mask, mean_col], errors='coerce').dropna()
    mean = float(mean_vals.sum()) if len(mean_vals) > 0 else float('nan')
    if std_col not in df.columns:
        return mean, float('nan')
    std_vals = pd.to_numeric(df.loc[mask, std_col], errors='coerce').fillna(0)
    std = float(np.sqrt((std_vals ** 2).sum()))
    return mean, std


def _fmt_metric(v):
    """Format a scoping-metric number per Section 3 convention.

    - If value is NaN -> 'N/A'.
    - If exactly 0 -> '0'.
    - If |value| >= 1 -> rounded integer with thousands commas.
    - If 0 < |value| < 1 -> rounded to ONE significant figure
      (e.g. 0.123 -> 0.1, 0.045 -> 0.05, 0.0023 -> 0.002). If rounding
      pushes the value across an order-of-magnitude boundary
      (e.g. 0.97 -> 1), the integer branch is used instead.
    """
    if isinstance(v, float) and math.isnan(v):
        return 'N/A'
    if v == 0:
        return '0'
    av = abs(v)
    if av >= 1:
        return f'{round(v):,}'
    exp = math.floor(math.log10(av))
    decimals = -exp
    rounded = round(v, decimals)
    if abs(rounded) >= 1:
        return f'{round(rounded):,}'
    return f'{rounded:.{decimals}f}'


def _fmt_table_val(x):
    if x == '-' or x is None or x == '':
        return x
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v == 0:
        return '0'
    sign = '-' if v < 0 else ''
    a = abs(v)
    if a >= 1e7:
        return f'{sign}{round(a / 1e6):.0f}M'
    elif a >= 1e6:
        return f'{sign}{round(a / 1e6):.0f}M'
    elif a >= 1e5:
        return f'{sign}{round(a / 1e4) * 10:.0f}K'
    elif a >= 1e4:
        return f'{sign}{round(a / 1e3):.0f}K'
    elif a >= 1e3:
        return f'{sign}{round(a / 1e3):.0f}K'
    else:
        return f'{sign}{int(round(a))}'


# ---------------------------------------------------------------------------
# HTML/CSS card helpers
# ---------------------------------------------------------------------------

def _help_icon(text):
    """Return HTML for a `?` help icon with a hover tooltip.

    `text` is plain text (no HTML); it's HTML-escaped and exposed via
    the `data-help` attribute (used by the CSS `::after` tooltip) and
    `aria-label` (for screen readers). We deliberately do NOT set the
    native `title=""` attribute because that triggers the browser's
    built-in tooltip *on top of* our styled tooltip the result is
    the same text rendering twice on hover.

    Newlines (`\\n` or `\\n\\n`) inside the text are collapsed to a
    single space before being put into the attribute value. Embedded
    newlines in attribute values cause Streamlit's CommonMark parser
    to treat the surrounding HTML as a "broken block" and dump the
    attribute string as visible text into the card. Strip them defensively.
    """
    flat = ' '.join(text.split())
    safe = html.escape(flat)
    return (
        f'<span class="mouse-help-icon" data-help="{safe}" '
        f'aria-label="{safe}" tabindex="0">?</span>'
    )


def _section_header(num, title_html, subtitle_html, top_margin='2.5rem'):
    """Render a Section 0X header as a tinted card.

    The header is one rounded box with:
    - a 4-px INL-blue left stripe (signature motif, same as Run Summary
      and the subsection headers),
    - a light-grey background (#f1f3f5 same as the Run Summary card)
      so the section break is visually distinct from inline content,
    - the section number prefix (`01 —`) in slate so the title pops,
    - the title in INL-blue uppercase, 1.15 rem,
    - the question-style subtitle in slate underneath.

    `num` is a 2-char string ('01'..'05'). `title_html` and
    `subtitle_html` may contain HTML entities; they're rendered as-is.
    `top_margin` defaults to 2.5 rem (between sections); Section 1
    overrides it to 0.5 rem (because the Run Summary card above has
    its own margin-bottom and the gap would otherwise stack).
    """
    st.markdown(
        f'<div style="margin-top:{top_margin};background:#f1f3f5;'
        f'border-radius:8px;'
        f'padding:0.85rem 1.1rem;margin-bottom:1.25rem;">'
        f'<div style="font-size:1.15rem;font-weight:700;color:#1B4F8C;'
        f'letter-spacing:0.04em;text-transform:uppercase;">'
        f'<span style="color:#94a3b8;">{num} &mdash;</span> '
        f'{title_html}</div>'
        f'<p style="color:#64748b;font-size:0.85rem;margin:0.2rem 0 0 0;">'
        f'{subtitle_html}'
        f'</p>'
        f'</div>',
        unsafe_allow_html=True,
    )


def _kpi_card(col, title, foak_val, noak_val, help_text=None):
    help_html = _help_icon(help_text) if help_text else ''
    col.markdown(
        f'''<div style="background:white;border-radius:8px;padding:1.1rem 1.25rem;
                        border:1px solid #bfdbfe;
                        min-height:110px;overflow-wrap:break-word;">
              <div style="font-size:0.85rem;font-weight:600;color:#64748b;
                          text-transform:uppercase;letter-spacing:0.09em;margin-bottom:0.65rem;">{title}{help_html}</div>
              <div style="display:flex;flex-wrap:wrap;align-items:baseline;gap:0.45rem;margin-bottom:0.45rem;">
                <span style="background:#fff3ed;color:#c84b1e;font-size:0.85rem;font-weight:700;
                             padding:0.12rem 0.4rem;border-radius:4px;letter-spacing:0.05em;
                             flex-shrink:0;">FOAK</span>
                <span style="font-size:1rem;font-weight:600;color:#0a2540;overflow-wrap:anywhere;">{foak_val}</span>
              </div>
              <div style="height:1px;background:#f1f3f5;margin:0.3rem 0;"></div>
              <div style="display:flex;flex-wrap:wrap;align-items:baseline;gap:0.45rem;margin-top:0.45rem;">
                <span style="background:#eff6ff;color:#1B4F8C;font-size:0.85rem;font-weight:700;
                             padding:0.12rem 0.4rem;border-radius:4px;letter-spacing:0.05em;
                             flex-shrink:0;">NOAK</span>
                <span style="font-size:1rem;font-weight:600;color:#0a2540;overflow-wrap:anywhere;">{noak_val}</span>
              </div>
            </div>''',
        unsafe_allow_html=True,
    )


def _keff_altair_chart(times, keffs):
    """Altair line+markers chart of k_eff vs Time with a horizontal
    reference line at k_eff=1. Returns the chart; caller renders it.

    Replaces a matplotlib version that contributed to the per-Run leak
    (TransformNode lambdas pinned by Streamlit's PNG cache). Vega-Lite
    renders client-side so the server holds only the data points."""
    _df = pd.DataFrame({'time': times, 'keff': keffs})
    _ymax = max(1.005, float(np.max(keffs)) + 0.005)
    _x_axis = alt.Axis(
        labelFontSize=13, labelColor='#000000', labelFontWeight=500,
        titleFontSize=14, titleColor='#000000', titleFontWeight='bold',
        domain=True, domainColor='#000000', domainWidth=1.5,
        tickColor='#000000', tickWidth=1.2,
    )
    _y_axis = alt.Axis(
        labelFontSize=13, labelColor='#000000', labelFontWeight=500,
        titleFontSize=14, titleColor='#000000', titleFontWeight='bold',
        domain=True, domainColor='#000000', domainWidth=1.5,
        tickColor='#000000', tickWidth=1.2,
    )
    _line = alt.Chart(_df).mark_line(color='#1B4F8C', strokeWidth=2).encode(
        x=alt.X('time:Q', title='Time (days)', axis=_x_axis),
        y=alt.Y('keff:Q', title='k_eff',
                scale=alt.Scale(domain=[0.98, _ymax]),
                axis=_y_axis),
    )
    _pts = alt.Chart(_df).mark_point(
        color='#1B4F8C', size=60, filled=True, stroke='white', strokeWidth=1
    ).encode(x='time:Q', y='keff:Q')
    _ref = alt.Chart(pd.DataFrame({'y': [1.0]})).mark_rule(
        color='#0a2540', strokeDash=[4, 4], strokeWidth=1.5
    ).encode(y='y:Q')
    return (_line + _pts + _ref).properties(height=260)


def _wrap_label_two_lines(text, max_chars=22):
    """Wrap a label string to 2 lines using '|' as the line-break
    delimiter (Altair's labelExpr split on this delimiter renders the
    pieces on separate lines). Returns the single-line input unchanged
    if it already fits in max_chars or has fewer than 2 words. Splits
    at the word boundary that minimizes the longer of the two lines."""
    text = str(text).strip()
    if len(text) <= max_chars:
        return text
    words = text.split()
    if len(words) < 2:
        return text
    best_split, best_max = 1, float('inf')
    for i in range(1, len(words)):
        line1 = ' '.join(words[:i])
        line2 = ' '.join(words[i:])
        m = max(len(line1), len(line2))
        if m < best_max:
            best_max, best_split = m, i
    return ' '.join(words[:best_split]) + '|' + ' '.join(words[best_split:])


def _grouped_lcoe_bars_chart(df, label_col='Account Title', height=400,
                             label_angle=-35, category_col=None):
    """Altair grouped bar chart (FOAK / NOAK side by side) with error
    bars. df must already be in the desired display order (typically
    sorted descending by FOAK LCOE); columns required:
    <label_col>, 'FOAK LCOE', 'NOAK LCOE'; optional: 'FOAK LCOE_std',
    'NOAK LCOE_std', <category_col>.

    When category_col is supplied, the x-axis tick labels are colored
    by COST_DRIVER_CATEGORY_PALETTE so the cost category is conveyed
    without spending the bar color channel (which stays reserved for
    FOAK vs NOAK).

    Replaces a matplotlib version that, multiplied across the
    per-section calls per Run, was a major leak source."""
    _has_foak_std = 'FOAK LCOE_std' in df.columns
    _has_noak_std = 'NOAK LCOE_std' in df.columns
    _rows = []
    for _, _row in df.iterrows():
        _fv = float(_row['FOAK LCOE'])
        _fs = float(_row['FOAK LCOE_std']) if _has_foak_std else 0.0
        _nv = float(_row['NOAK LCOE'])
        _ns = float(_row['NOAK LCOE_std']) if _has_noak_std else 0.0
        _rows.append({label_col: str(_row[label_col]),
                      'Type': 'FOAK', 'LCOE': _fv,
                      'lower': _fv - _fs, 'upper': _fv + _fs})
        _rows.append({label_col: str(_row[label_col]),
                      'Type': 'NOAK', 'LCOE': _nv,
                      'lower': _nv - _ns, 'upper': _nv + _ns})
    _plot_df = pd.DataFrame(_rows)
    # Explicit sort list preserves df's existing row order. Without
    # this, Vega-Lite falls back to alphabetical, which scrambles the
    # "largest driver first" ordering we want. labelOverlap=False
    # forces every x-axis label to render (the default 'greedy' policy
    # hides ~half the labels when they collide).
    _label_order = [str(v) for v in df[label_col]]

    # Per-label color via Vega expression: each x-axis tick gets the
    # color of its cost category. Falls back to the default text color
    # when no category column is provided or a label has no mapping.
    _default_label_color = '#0a2540'
    _label_color = _default_label_color
    if category_col is not None and category_col in df.columns:
        _label_to_color = {
            str(_lbl): COST_DRIVER_CATEGORY_PALETTE.get(
                str(_cat), _default_label_color)
            for _lbl, _cat in zip(df[label_col], df[category_col])
        }
        if _label_to_color:
            _ternaries = ''.join(
                f'datum.value === {_json.dumps(_lbl)} ? '
                f'{_json.dumps(_col)} : '
                for _lbl, _col in _label_to_color.items()
            )
            _label_color = alt.ExprRef(
                f'{_ternaries}{_json.dumps(_default_label_color)}'
            )

    _color_scale = alt.Scale(domain=['FOAK', 'NOAK'],
                             range=['#c84b1e', '#1B4F8C'])
    _x_enc = alt.X(f'{label_col}:N',
                   sort=_label_order,
                   title=None,
                   axis=alt.Axis(labelAngle=label_angle,
                                 labelColor=_label_color,
                                 labelFontSize=12, labelFontWeight=600,
                                 labelLimit=200, labelOverlap=False,
                                 # Multi-line wrap: labels passed through
                                 # _wrap_label_two_lines insert '|' at the
                                 # break point; split() turns that into a
                                 # 2-line label. Labels without '|' render
                                 # on a single line (split returns a
                                 # 1-element array).
                                 labelExpr="split(datum.label, '|')"))
    _bars = alt.Chart(_plot_df).mark_bar().encode(
        x=_x_enc,
        xOffset=alt.XOffset('Type:N', sort=['FOAK', 'NOAK']),
        y=alt.Y('LCOE:Q', title='LCOE Contribution ($/MWh)',
                axis=alt.Axis(titleColor='#0a2540', labelColor='#0a2540',
                              gridDash=[3, 3], gridColor='#cbd5e1')),
        color=alt.Color('Type:N', scale=_color_scale,
                        legend=alt.Legend(title=None, orient='top-right')),
    )
    _errs = alt.Chart(_plot_df).mark_errorbar(thickness=1.8, ticks=True,
                                              color='#0a2540').encode(
        x=_x_enc,
        xOffset=alt.XOffset('Type:N', sort=['FOAK', 'NOAK']),
        y=alt.Y('lower:Q', title='LCOE Contribution ($/MWh)'),
        y2='upper:Q',
    )
    return (_bars + _errs).properties(height=height).configure_view(
        stroke=None, fill='#f8fafc'
    )


CURATED_COST_DRIVER_ACCOUNTS = [
    # 10: pre-construction expanded one level.
    '11', '12', '13', '14', '15',
    # 20: direct cost curated to familiar user-facing groupings.
    '21',
    '221.1', '221.2', '221.31', '221.32', '221.33', '221.34',
    '222', '223.2', '226', '227', '228',
    '23', '24', '25', '26',
    # 30, 40, 60: expanded to available children.
    '31', '32', '33', '34', '35', '36',
    '41',
    '61', '62', '63', '64', '65', '66', '67', '68', '69',
    # 70: O&M curated, with account 75 expanded to replacement drivers.
    '71', '72', '73',
    '751', '752', '753', '754', '755', '756', '757', '758', '759',
    '76', '77', '78', '79',
    # 80: fuel cycle expanded one level.
    '81', '82', '83',
]


COST_DRIVER_TITLE_OVERRIDES = {
    '759': 'Other Direct Plant Maintenance & Replacements',
}


COST_DRIVER_DEFINITION_SOURCE = (
    'Definitions are condensed from the EPRI GN-COA appendix on final account '
    'descriptions. Replacement rows use the same sustaining-capital category.'
)


COST_DRIVER_DEFINITIONS = {
    '11': 'Land needed for the reactor site and related facilities; examples include land purchases and land rights.',
    '12': 'Permits that allow permanent plant construction to proceed; examples include federal, state, and local site permits.',
    '13': 'Licensing work needed to construct and operate the plant; examples include preapplication work, application preparation, and regulator review.',
    '14': 'Plant construction and operating permits beyond the main license; examples include environmental and energy-market approvals.',
    '15': 'Studies supporting plant construction and operation; examples include traffic, archaeological, historical, and site studies.',
    '21': 'Civil work and permanent structures for the plant; examples include site preparation, reactor island structures, buildings, roads, and foundations.',
    '221.1': 'Reactor vessel boundary and supporting accessories; examples include vessel shells, supports, and structural accessories.',
    '221.2': 'Equipment that changes or shuts down reactor reactivity; examples include drums, rods, blades, shutdown rods, and drives.',
    '221.31': 'Core materials that reflect neutrons back into the reactor; examples include reflector assemblies or blocks.',
    '221.32': 'Shielding located inside the core geometry; external biological shielding is counted with civil structures.',
    '221.33': 'Non-coolant material that moderates the neutron spectrum; examples include moderator blocks or assemblies.',
    '221.34': 'Supplemental moderator material within non-fuel core internals; examples include booster moderator elements.',
    '222': 'System that carries heat from the core to the application interface; examples include coolant, piping, pumps, circulators, and heat exchangers.',
    '223.2': 'Safety cooling equipment that removes heat through the reactor vessel; examples include external cavity cooling components.',
    '226': 'Supporting reactor plant systems not otherwise itemized; examples include inert gas, makeup coolant, treatment, auxiliary cooling, and sampling equipment.',
    '227': 'Controls and instruments for monitoring and protecting the reactor; examples include panels, computers, sensors, cables, and software.',
    '228': 'Miscellaneous reactor-plant installation items; examples include painting, welder qualification, and insulation.',
    '23': 'Equipment that converts reactor heat into the final product; examples include turbines, generators, process-heat export, hydrogen equipment, and heat rejection.',
    '24': 'Electrical service and distribution equipment; examples include switchgear, station service, switchboards, protective equipment, raceways, and wiring.',
    '25': 'First core fuel purchased before commissioning; examples include mining, conversion, enrichment, fabrication, and other fissionable feed material.',
    '26': 'Plant equipment not captured in the major systems; examples include cranes, utility systems, communications, furnishings, and remote inspection equipment.',
    '31': 'Temporary construction support and field indirects; examples include tools, vehicles, consumables, utilities, laydown areas, field shops, and insurance.',
    '32': 'Supervision of construction craft work; examples include field engineers, superintendents, and contractor supervisory staff.',
    '33': 'Startup and commissioning work before commercial operation; examples include trial runs, fuel loading, coolant loading, and test runs.',
    '34': 'Transport of major project equipment and modules; examples include fuel, reactor modules, energy-conversion modules, and construction modules.',
    '35': 'Engineering and design services for the plant; examples include layout, site engineering, project engineering, design QA, and owner oversight.',
    '36': 'Project and construction management support; examples include QA, procurement, contracts, controls, HR, medical, and safety support.',
    '41': 'Pre-startup staffing preparation; examples include recruiting and training operators before commissioning or demonstration testing.',
    '61': 'Allowance for capital-cost escalation; normally excluded from fixed-year constant-dollar estimates.',
    '62': 'Financing cost before commercial operation; examples include interest during construction on financed upfront costs.',
    '63': 'Accounting depreciation of plant assets.',
    '64': 'Capitalized financial allowance represented in the estimate.',
    '65': 'Capitalized financial allowance represented in the estimate.',
    '66': 'Capitalized financial allowance represented in the estimate.',
    '67': 'Capitalized financial allowance represented in the estimate.',
    '68': 'Capitalized financial allowance represented in the estimate.',
    '69': 'Financial-cost contingency for uncertainty before operation; examples include schedule-related financing uncertainty.',
    '71': 'Annual staffing cost for operations and maintenance; examples include operators, security, maintenance, engineering, management, and benefits.',
    '72': 'Recurring non-fuel consumables used in operation; examples include borated water, control rods, burnable poisons, and process chemicals.',
    '73': 'Recurring regulatory oversight cost; examples include plant inspections and regulator fees.',
    '751': 'Sustaining-capital cost for maintaining long-term plant capability; this row represents reactor pressure vessel replacement.',
    '752': 'Sustaining-capital cost for maintaining long-term plant capability; this row represents core barrel or inner-vessel replacement.',
    '753': 'Sustaining-capital cost for maintaining long-term plant capability; this row represents moderator replacement.',
    '754': 'Sustaining-capital cost for maintaining long-term plant capability; this row represents reflector replacement.',
    '755': 'Sustaining-capital cost for maintaining long-term plant capability; this row represents reactivity-control equipment replacement.',
    '756': 'Sustaining-capital cost for maintaining long-term plant capability; this row represents integrated heat-transfer equipment replacement.',
    '757': 'Sustaining-capital cost for maintaining long-term plant capability; examples include additional major replacement work.',
    '758': 'Sustaining-capital cost for maintaining long-term plant capability; examples include additional major replacement work.',
    '759': ('Sustaining-capital cost for maintaining or improving plant capability after '
            'operation begins. Examples include replacements or upgrades for civil '
            'structures, safety systems, fuel handling, I&C, energy conversion, '
            'electrical equipment, cranes, utilities, communications, and fixtures.'),
    '76': 'Recurring taxes and insurance during operation; salary-related taxes are counted with staffing.',
    '77': 'Maintenance outage cost during operation; examples include scheduled and unscheduled outages, excluding refueling work.',
    '78': 'Recurring provision for plant end-of-life work; examples include decommissioning, decontamination, dismantling, and close-out engineering.',
    '79': 'Contingency for annual O&M uncertainty; management reserve is excluded.',
    '81': 'Recurring cost of refueling operations; examples include fuel management, licensing support, QA, inspection, handling, and storage planning.',
    '82': 'Recurring cost of buying additional fuel; examples include reload uranium supply, conversion, enrichment, fabrication, and tails disposal.',
    '83': 'Recurring cost of spent-fuel management; examples include interim storage, reprocessing, recovered-material credits, and final disposal.',
}


COST_DRIVER_CATEGORY_ORDER = [
    'Preconstruction',
    'Direct',
    'Indirect',
    'Training',
    'Financial',
    'Annualized O&M',
    'Annualized fuel',
]


COST_DRIVER_CATEGORY_BY_LEAD_DIGIT = {
    '1': 'Preconstruction',
    '2': 'Direct',
    '3': 'Indirect',
    '4': 'Training',
    '6': 'Financial',
    '7': 'Annualized O&M',
    '8': 'Annualized fuel',
}


COST_DRIVER_CATEGORY_PALETTE = {
    'Preconstruction': '#6b7280',
    'Direct':          '#1B4F8C',
    'Indirect':        '#0891b2',
    'Training':        '#7c3aed',
    'Financial':       '#c84b1e',
    'Annualized O&M':  '#15803d',
    'Annualized fuel': '#b45309',
}


def _account_code(x):
    try:
        code = f'{float(x):.6f}'.rstrip('0').rstrip('.')
    except (TypeError, ValueError):
        code = str(x).strip()
    return code


def _cost_driver_title(code, title):
    return COST_DRIVER_TITLE_OVERRIDES.get(code, str(title).strip())


def _cost_driver_treatment(code):
    return 'Annualized' if str(code).strip()[0] in ('7', '8') else 'Capitalized'


def _cost_driver_category(code):
    _lead = str(code).strip()[:1]
    return COST_DRIVER_CATEGORY_BY_LEAD_DIGIT.get(_lead, 'Other')


def _curated_cost_driver_candidates(enriched_df):
    _df = enriched_df.copy()
    _df['_Account Code'] = _df['Account'].apply(_account_code)
    _df = _df[_df['_Account Code'].isin(CURATED_COST_DRIVER_ACCOUNTS)].copy()
    if _df.empty:
        return _df
    _df = (_df.sort_values('FOAK LCOE', ascending=False)
             .drop_duplicates('_Account Code', keep='first'))
    _df['Driver Title'] = [
        _cost_driver_title(code, title)
        for code, title in zip(_df['_Account Code'], _df['Account Title'])
    ]
    _df['Category'] = [_cost_driver_category(code) for code in _df['_Account Code']]
    return _df


def _hidden_cost_driver_lcoe_totals(enriched_df, min_foak_lcoe=10.0):
    _df = _curated_cost_driver_candidates(enriched_df)
    if _df.empty:
        return 0.0, 0.0, 0
    _hidden = _df[pd.to_numeric(_df['FOAK LCOE'], errors='coerce') < min_foak_lcoe]
    if _hidden.empty:
        return 0.0, 0.0, 0
    _foak = pd.to_numeric(_hidden['FOAK LCOE'], errors='coerce').fillna(0.0).sum()
    _noak = pd.to_numeric(_hidden['NOAK LCOE'], errors='coerce').fillna(0.0).sum()
    return float(_foak), float(_noak), int(len(_hidden))


def _curated_cost_drivers(enriched_df, min_foak_lcoe=10.0,
                          max_drivers=28):
    """Return the curated cost-driver rows requested for the app plot.

    The account list intentionally mixes levels in the Code of Accounts:
    some parent accounts stay aggregated, while selected parents are expanded
    to the account level most recognizable to users.
    """
    _df = _curated_cost_driver_candidates(enriched_df)
    if _df.empty:
        return _df

    _df = _df[_df['FOAK LCOE'] >= min_foak_lcoe]
    _df = _df.sort_values('FOAK LCOE', ascending=False).head(max_drivers)
    _df['Driver Label'] = [
        _wrap_label_two_lines(title, max_chars=26)
        for title in _df['Driver Title']
    ]
    return _df


def _render_cost_driver_definition_table(driver_df):
    if driver_df.empty:
        return

    _rows = []
    for _, _row in driver_df.sort_values('FOAK LCOE', ascending=False).iterrows():
        _code = str(_row['_Account Code'])
        _name = str(_row.get('Driver Title') or _row['Account Title']).strip()
        _treatment = _cost_driver_treatment(_code)
        _definition = COST_DRIVER_DEFINITIONS.get(
            _code,
            f'Includes cost allowances for {_name.lower()}.'
        )
        _rows.append(
            '<div style="border:1px solid #e2e8f0;border-radius:6px;'
            'padding:0.65rem 0.75rem;background:#ffffff;">'
            '<div style="font-size:0.85rem;font-weight:700;color:#0a2540;'
            'line-height:1.25;margin-bottom:0.25rem;">'
            f'{html.escape(_name)} ({_treatment})</div>'
            '<div style="font-size:0.82rem;line-height:1.35;color:#3c4257;'
            'display:-webkit-box;-webkit-line-clamp:3;-webkit-box-orient:vertical;'
            'overflow:hidden;">'
            f'{html.escape(_definition)}</div>'
            '</div>'
        )

    st.markdown(
        '<div style="font-size:0.9rem;font-weight:600;color:#0a2540;'
        'margin:1.25rem 0 0.45rem 0;">Account definitions</div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        '<p style="color:#64748b;font-size:0.8rem;line-height:1.35;'
        'margin:0 0 0.65rem 0;">'
        f'{html.escape(COST_DRIVER_DEFINITION_SOURCE)}</p>',
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div style="display:grid;grid-template-columns:repeat(2,minmax(0,1fr));'
        'gap:0.65rem;margin-bottom:0.25rem;">'
        f'{"".join(_rows)}</div>',
        unsafe_allow_html=True,
    )


def _side_view_altair_chart(diameter_cm, active_height_cm,
                            axial_reflector_cm, radial_reflector_cm):
    """Altair to-scale side view of the reactor cylinder. Replaces a
    matplotlib version that created one Figure + Axes + 4 annotations
    per Run; with Streamlit's PNG cache pinning, that cost ~2-3 MB/Run
    of TransformNode/CallbackRegistry state. Two nested rectangles
    (outer envelope + active core) with four dimension labels."""
    total_h = active_height_cm + 2 * axial_reflector_cm
    active_d = diameter_cm - 2 * radial_reflector_cm
    half_d = diameter_cm / 2.0
    half_h = total_h / 2.0
    half_active_d = active_d / 2.0
    half_active_h = active_height_cm / 2.0

    # Domain with padding on the right (for H_active/H_total labels)
    # and below (for D_active/D_total labels) so labels never clip.
    x_pad_l = diameter_cm * 0.05
    x_pad_r = diameter_cm * 0.30
    y_pad_t = total_h * 0.05
    y_pad_b = total_h * 0.18
    x_min = -half_d - x_pad_l
    x_max = half_d + x_pad_r
    y_min = -half_h - y_pad_b
    y_max = half_h + y_pad_t

    _x_scale = alt.Scale(domain=[x_min, x_max], nice=False)
    _y_scale = alt.Scale(domain=[y_min, y_max], nice=False)

    _rect_df = pd.DataFrame([
        {'x': -half_d, 'x2': half_d,
         'y': -half_h, 'y2': half_h},
        {'x': -half_active_d, 'x2': half_active_d,
         'y': -half_active_h, 'y2': half_active_h},
    ])
    _rects = alt.Chart(_rect_df).mark_rect(
        fill='white', stroke='#0a2540', strokeWidth=1.5,
    ).encode(
        x=alt.X('x:Q', scale=_x_scale, axis=None),
        x2='x2:Q',
        y=alt.Y('y:Q', scale=_y_scale, axis=None),
        y2='y2:Q',
    )

    _labels_left = pd.DataFrame([
        {'x': half_active_d + diameter_cm * 0.02,
         'y': half_h * 0.30,
         'text': f'H_active = {active_height_cm:.0f} cm'},
        {'x': half_d + diameter_cm * 0.02,
         'y': -half_h * 0.30,
         'text': f'H_total = {total_h:.0f} cm'},
    ])
    _labels_center = pd.DataFrame([
        {'x': 0,
         'y': -half_active_h - total_h * 0.04,
         'text': f'D_active = {active_d:.0f} cm'},
        {'x': 0,
         'y': -half_h - total_h * 0.06,
         'text': f'D_total = {diameter_cm:.0f} cm'},
    ])
    _text_left = alt.Chart(_labels_left).mark_text(
        fontSize=11, fontWeight='bold', color='#0a2540',
        baseline='middle', align='left',
    ).encode(
        x=alt.X('x:Q', scale=_x_scale, axis=None),
        y=alt.Y('y:Q', scale=_y_scale, axis=None),
        text='text:N',
    )
    _text_center = alt.Chart(_labels_center).mark_text(
        fontSize=11, fontWeight='bold', color='#0a2540',
        baseline='top', align='center',
    ).encode(
        x=alt.X('x:Q', scale=_x_scale, axis=None),
        y=alt.Y('y:Q', scale=_y_scale, axis=None),
        text='text:N',
    )

    # Aspect-preserved sizing: 1 cm in x ~ 1 cm in y on screen.
    aspect = (y_max - y_min) / (x_max - x_min)
    base_w = 380
    chart_h = max(160, min(640, int(base_w * aspect)))

    return (_rects + _text_left + _text_center).properties(
        width=base_w, height=chart_h,
        title=alt.TitleParams('Side View (to scale)',
                              fontSize=11, fontWeight='bold',
                              color='#0a2540', anchor='start'),
    ).configure_view(stroke=None)


def _materials_section(reactor_type, params):
    """Render a 'Materials & Components' panel that lists the basic
    materials of the reactor (fuel, moderator, reflector, drums, …)
    plus, for LTMR, the per-assembly fuel and moderator pin counts.
    All values are read directly from `params` no interpolation."""
    def _pretty(name):
        """Replace underscores with spaces so material names like
        'UZrH_alloy' read as 'UZrH alloy' in the displayed table."""
        if not isinstance(name, str):
            return name
        return name.replace('_', ' ')

    rows = []
    # Fuel row. Per reactor type:
    #   LTMR: UZrH alloy fuel pins; show the U weight fraction.
    #   GCMR: TRISO particles (UCO kernel) modeled explicitly in the
    #         compact (particle radii + packing fraction resolved).
    #   HPMR: TRISO fuel modeled as homogenized (particles smeared into
    #         the pin volume rather than resolved individually).
    _fuel_str = _pretty(params.get('Fuel', ''))
    if reactor_type == 'LTMR' and 'U_met_wo' in params:
        _u_wo = float(params['U_met_wo']) * 100.0
        _fuel_str = f"{_fuel_str} ({_u_wo:.0f} wt% U)"
    elif reactor_type == 'GCMR':
        _fuel_str = f"TRISO &mdash; {_fuel_str} kernel (particles resolved)"
    elif reactor_type == 'HPMR':
        _fuel_str = "TRISO (modeled as homogenized)"
    rows.append(('Fuel', _fuel_str))

    rows.append(('Moderator', _pretty(params.get('Moderator', ''))))
    if params.get('Moderator Booster Materials'):
        rows.append(('Moderator booster',
                     ', '.join(_pretty(m) for m in params['Moderator Booster Materials'])))
    if reactor_type == 'HPMR':
        rows.append(('Cooling device', _pretty(params.get('Cooling Device', 'Heat pipes'))))
    else:
        rows.append(('Coolant', _pretty(params.get('Coolant', ''))))
    rows.append(('Radial reflector', _pretty(params.get('Radial Reflector', ''))))
    rows.append(('Axial reflector', _pretty(params.get('Axial Reflector', ''))))
    rows.append(('Control drum absorber', _pretty(params.get('Control Drum Absorber', ''))))
    rows.append(('Control drum reflector', _pretty(params.get('Control Drum Reflector', ''))))

    # LTMR is a single-assembly core no per-assembly distinction needed.
    if reactor_type == 'LTMR':
        if 'Fuel Pin Count' in params:
            rows.append(('Number of fuel pins',
                         f"{int(params['Fuel Pin Count']):,}"))
        if 'Moderator Pin Count' in params:
            rows.append(('Number of moderator pins',
                         f"{int(params['Moderator Pin Count']):,}"))

    # GCMR-specific TRISO details
    if reactor_type == 'GCMR':
        if 'Packing Fraction' in params:
            _pf = float(params['Packing Fraction']) * 100.0
            rows.append(('TRISO packing fraction', f"{_pf:.1f} %"))
        if 'Total Number of TRISO Particles' in params:
            _ntot = int(params['Total Number of TRISO Particles'])
            rows.append(('Total TRISO particles', f"{_ntot:,}"))

    # Render as a 2-column CSS grid (label col + value col), both sized
    # to their max content. Labels share a right edge, values share a
    # left edge no ragged alignment, no stray table borders.
    body = ''.join(
        f'<div style="color:#3c4257;font-weight:600;white-space:nowrap;">{k}</div>'
        f'<div style="color:#0a2540;font-weight:500;">{v}</div>'
        for k, v in rows
    )
    st.markdown(
        '<div style="font-size:0.85rem;font-weight:600;color:#64748b;text-transform:uppercase;'
        'letter-spacing:0.09em;margin-bottom:0.45rem;">Materials &amp; Components</div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        f'''<div style="background:#f7f8fa;border:1px solid #bfdbfe;border-radius:8px;
                        padding:0.85rem 1.1rem;margin-bottom:0.85rem;">
              <div style="display:grid;grid-template-columns:max-content max-content;
                          row-gap:0.35rem;column-gap:1.25rem;font-size:0.85rem;">
                {body}
              </div>
            </div>''',
        unsafe_allow_html=True,
    )


def _ring_count(rings_over_one_edge):
    rings = int(rings_over_one_edge)
    return 2 * rings * (rings - 1) + 2 * sum(range(1, rings - 1)) + 2 * rings - 1


def _density_g_cm3(params, material_name, default=0.0):
    rtype = params.get('reactor type', 'LTMR')
    return float(_MATERIALS_RAW.get(rtype, {}).get(material_name, default))


def _cylinder_mass_kg(radius_cm, height_cm, density_g_cm3, count=1):
    return count * math.pi * radius_cm ** 2 * height_cm * density_g_cm3 / 1000.0


def _annulus_mass_kg(r_outer_cm, r_inner_cm, height_cm, density_g_cm3, count=1):
    area_cm2 = math.pi * (r_outer_cm ** 2 - r_inner_cm ** 2)
    return count * area_cm2 * height_cm * density_g_cm3 / 1000.0


def _uranium_atomic_mass(params):
    enrichment = float(params.get('Enrichment', 0.1975))
    return enrichment * 235.0439299 + (1.0 - enrichment) * 238.0507882


def _uco_uranium_mass_fraction(params):
    # UCO is represented as a UO2/UC mixture. Per uranium atom, the
    # non-uranium atoms are 2*f oxygen atoms and (1-f) carbon atoms.
    f_uo2 = float(params.get('UO2 atom fraction', 0.7))
    m_u = _uranium_atomic_mass(params)
    return m_u / (m_u + 2.0 * f_uo2 * 15.999 + (1.0 - f_uo2) * 12.011)


def _homog_triso_uranium_mass_fraction(params):
    # Atom densities from core_design/openmc_materials_database.py for
    # the HPMR homogenized TRISO material.
    e = float(params.get('Enrichment', 0.1975))
    atom_densities = {
        'U235': 0.00130037929 * e,
        'U238': 0.00130037929 * (1.0 - e),
        'O16': 2.59371545E-03,
        'O17': 1.05004397E-06,
        'O18': 5.99797186E-06,
        'Si28': 2.76954169E-03,
        'Si29': 1.40694868E-04,
        'Si30': 9.28556098E-05,
        'C12': 7.31619752E-02,
        'C13': 7.58819416E-04,
    }
    atomic_masses = {
        'U235': 235.0439299, 'U238': 238.0507882,
        'O16': 15.9949146, 'O17': 16.9991318, 'O18': 17.9991596,
        'Si28': 27.9769265, 'Si29': 28.9764947, 'Si30': 29.9737701,
        'C12': 12.0, 'C13': 13.0033548,
    }
    total = sum(atom_densities[k] * atomic_masses[k] for k in atom_densities)
    uranium = atom_densities['U235'] * atomic_masses['U235'] + atom_densities['U238'] * atomic_masses['U238']
    return uranium / total


def _heatpipe_steel_mass_kg(params):
    # Heat pipe material is modeled as homogenized SS316 + sodium. For
    # transportability, include only the steel contribution and leave the
    # sodium working fluid out with other coolant inventories.
    atom_densities = {
        'Si28': 1.49701E-02, 'Si29': 7.60143E-04, 'Si30': 5.01090E-04,
        'Cr50': 6.46763E-03, 'Cr52': 1.24724E-01, 'Cr53': 1.41423E-02, 'Cr54': 3.52029E-03,
        'Mn55': 1.66133E-02,
        'Fe54': 3.12186E-02, 'Fe56': 4.90061E-01, 'Fe57': 1.13180E-02, 'Fe58': 1.50617E-03,
        'Ni58': 6.33738E-02, 'Ni60': 2.44119E-02, 'Ni61': 1.06115E-03, 'Ni62': 3.38338E-03, 'Ni64': 8.61654E-04,
        'Mo92': 1.75699E-03, 'Mo94': 1.09514E-03, 'Mo95': 1.88484E-03, 'Mo96': 1.97478E-03,
        'Mo97': 1.13066E-03, 'Mo98': 2.85681E-03, 'Mo100': 1.14011E-03,
    }
    atomic_masses = {
        'Si28': 27.9769265, 'Si29': 28.9764947, 'Si30': 29.9737701,
        'Cr50': 49.9460418, 'Cr52': 51.9405062, 'Cr53': 52.9406481, 'Cr54': 53.9388792,
        'Mn55': 54.9380439,
        'Fe54': 53.9396082, 'Fe56': 55.9349363, 'Fe57': 56.9353928, 'Fe58': 57.9332744,
        'Ni58': 57.9353424, 'Ni60': 59.9307859, 'Ni61': 60.9310556, 'Ni62': 61.9283454, 'Ni64': 63.9279660,
        'Mo92': 91.9068072, 'Mo94': 93.9050836, 'Mo95': 94.9058374, 'Mo96': 95.9046748,
        'Mo97': 96.9060169, 'Mo98': 97.9054036, 'Mo100': 99.9074680,
    }
    amu_to_g_per_cm3 = 1.66053906660
    steel_density = amu_to_g_per_cm3 * sum(atom_densities[k] * atomic_masses[k] for k in atom_densities)
    return _cylinder_mass_kg(
        float(params.get('Heat Pipe Radii', [0.0])[0]),
        float(params.get('Active Height', 0.0)),
        steel_density,
        int(params.get('Number of Heatpipes', 0)),
    )


def _transport_fuel_masses(params, reactor_type):
    uranium_mass_kg = float(params.get('Uranium Mass', 0.0))
    fuel_mass_kg = uranium_mass_kg
    matrix_already_in_moderator_kg = 0.0
    heatpipe_steel_kg = 0.0

    if reactor_type == 'LTMR':
        u_weight_fraction = float(params.get('U_met_wo', 0.0))
        if u_weight_fraction > 0:
            fuel_mass_kg = uranium_mass_kg / u_weight_fraction
        radii = params.get('Fuel Pin Radii', [])
        materials = params.get('Fuel Pin Materials', [])
        height = float(params.get('Active Height', 0.0))
        count = int(params.get('Fuel Pin Count', 0))
        if len(radii) >= 5 and len(materials) >= 5:
            fuel_mass_kg += _cylinder_mass_kg(radii[0], height, _density_g_cm3(params, materials[0]), count)
            fuel_mass_kg += _annulus_mass_kg(radii[4], radii[3], height, _density_g_cm3(params, materials[4]), count)

    elif reactor_type == 'GCMR':
        fuel_mass_kg = uranium_mass_kg / _uco_uranium_mass_fraction(params)
        radii = params.get('Fuel Pin Radii', [])
        materials = params.get('Fuel Pin Materials', [])
        particles = int(params.get('Total Number of TRISO Particles', 0))
        if len(radii) >= 5 and len(materials) >= 5 and particles > 0:
            for i in range(1, 5):
                fuel_mass_kg += _annulus_mass_kg(
                    radii[i], radii[i - 1], 1.0,
                    _density_g_cm3(params, materials[i]), particles,
                )
        compact_count = _ring_count(params.get('Assembly Rings', 1) - 1) * _ring_count(params.get('Core Rings', 1))
        compact_volume = math.pi * float(params.get('Compact Fuel Radius', 0.0)) ** 2 * float(params.get('Active Height', 0.0))
        matrix_already_in_moderator_kg = (
            compact_count
            * compact_volume
            * (1.0 - float(params.get('Packing Fraction', 0.0)))
            * _density_g_cm3(params, params.get('Matrix Material', 'Graphite'))
            / 1000.0
        )
        fuel_mass_kg += matrix_already_in_moderator_kg

    elif reactor_type == 'HPMR':
        fuel_mass_kg = uranium_mass_kg / _homog_triso_uranium_mass_fraction(params)
        heatpipe_steel_kg = _heatpipe_steel_mass_kg(params)

    return fuel_mass_kg, matrix_already_in_moderator_kg, heatpipe_steel_kg


def _info_card(col, title, value, subtitle='', accent='#64748b', bg='white', border='#bfdbfe'):
    sub_html = f'<div style="font-size:0.85rem;color:#64748b;margin-top:0.2rem;">{subtitle}</div>' if subtitle else ''
    col.markdown(
        f'''<div style="background:{bg};border:1px solid {border};border-radius:8px;
                        padding:1.1rem 1.25rem;min-height:80px;">
              <div style="font-size:0.85rem;font-weight:600;color:{accent};
                          text-transform:uppercase;letter-spacing:0.09em;margin-bottom:0.35rem;">{title}</div>
              <div style="font-size:1.5rem;font-weight:700;color:#0a2540;line-height:1.2;">{value}</div>
              {sub_html}
            </div>''',
        unsafe_allow_html=True,
    )


# ---------------------------------------------------------------------------
# Global CSS
# ---------------------------------------------------------------------------
_CSS = """
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');

html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
}

/* ── App background ── */
.stApp { background: #f7f8fa !important; }

/* ── Hide Streamlit chrome ── */
#MainMenu { visibility: hidden; }
footer { visibility: hidden; }
[data-testid="stHeader"] {
    background: transparent !important;
    box-shadow: none !important;
}

/* ── Sidebar ── */
section[data-testid="stSidebar"] > div:first-child {
    background: #f1f3f5;
    border-right: 4px solid #1B4F8C;
}
section[data-testid="stSidebar"] .stMarkdown strong {
    color: #1B4F8C !important;
    font-size: 0.85rem;
    text-transform: uppercase;
    letter-spacing: 0.1em;
}
section[data-testid="stSidebar"] .stMarkdown p {
    color: #64748b !important;
}
section[data-testid="stSidebar"] label,
section[data-testid="stSidebar"] [data-testid="stWidgetLabel"] p {
    color: #3c4257 !important;
    font-weight: 500 !important;
    font-size: 0.85rem !important;
}
section[data-testid="stSidebar"] .stCaption p {
    color: #64748b !important;
}
section[data-testid="stSidebar"] hr {
    border-color: #bfdbfe !important;
    margin: 0.6rem 0 !important;
}
/* Help icons (?) in sidebar — filled INL blue circle with white question mark */
section[data-testid="stSidebar"] [data-testid="stTooltipHoverTarget"],
section[data-testid="stSidebar"] .stTooltipHoverTarget {
    background: #1B4F8C !important;
    border-radius: 50% !important;
    width: 1.15rem !important;
    height: 1.15rem !important;
    display: inline-flex !important;
    align-items: center !important;
    justify-content: center !important;
    position: relative !important;
    opacity: 1 !important;
}
section[data-testid="stSidebar"] [data-testid="stTooltipHoverTarget"] svg,
section[data-testid="stSidebar"] .stTooltipHoverTarget svg {
    display: none !important;
}
section[data-testid="stSidebar"] [data-testid="stTooltipHoverTarget"]::after,
section[data-testid="stSidebar"] .stTooltipHoverTarget::after {
    content: '?' !important;
    color: white !important;
    font-weight: 700 !important;
    font-size: 0.85rem !important;
    line-height: 1 !important;
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif !important;
}
section[data-testid="stSidebar"] .stButton > button {
    background: #c84b1e !important;
    color: white !important;
    font-weight: 600 !important;
    border: none !important;
    border-radius: 8px !important;
    padding: 0.65rem 1rem !important;
    font-size: 1rem !important;
    letter-spacing: 0.02em !important;
    width: 100%;
}
section[data-testid="stSidebar"] .stButton > button:hover {
    background: #9a3412 !important;
    transform: translateY(-1px);
}

/* ── Sidebar inputs (force white bg + visible border on the gray sidebar) ── */
section[data-testid="stSidebar"] [data-baseweb="select"] > div,
section[data-testid="stSidebar"] [data-baseweb="input"],
section[data-testid="stSidebar"] [data-baseweb="input"] > div,
section[data-testid="stSidebar"] [data-testid="stNumberInput"] input,
section[data-testid="stSidebar"] [data-testid="stTextInput"] input,
section[data-testid="stSidebar"] [data-testid="stNumberInput"] button,
section[data-testid="stSidebar"] input[type="text"],
section[data-testid="stSidebar"] input[type="number"] {
    background: white !important;
    border: 1px solid #cbd5e1 !important;
}

/* ── Tabs ── */
.stTabs [data-baseweb="tab-list"] {
    background: transparent;
    gap: 0.25rem;
    border-bottom: 1px solid #1B4F8C;
}
.stTabs [data-baseweb="tab"] {
    background: transparent;
    border: none;
    border-radius: 8px 8px 0 0;
    color: #64748b;
    font-weight: 600;
    font-size: 1rem;
    padding: 0.55rem 1.2rem;
    margin-bottom: -1px;
}
.stTabs [aria-selected="true"] {
    background: white !important;
    color: #1B4F8C !important;
    border-top: 1px solid #1B4F8C !important;
    border-left: 1px solid #1B4F8C !important;
    border-right: 1px solid #1B4F8C !important;
    border-bottom: 1px solid white !important;
}
.stTabs [data-baseweb="tab-panel"] {
    background: white;
    border-radius: 0 8px 8px 8px;
    padding: 1.5rem 1.75rem;
    border: 1px solid #1B4F8C;
}

/* ── Expander ── */
/* Box around the expander is kept subtle so the expanded content has  */
/* a visible container, but the summary (click target) is styled to    */
/* match the popover-trigger caption style: 0.85rem slate text, faint   */
/* underline, hover turns it INL-blue. Consistent affordance vocabulary.*/
[data-testid="stExpander"] {
    border: 1px solid #bfdbfe !important;
    border-radius: 8px !important;
    background: white;
}
[data-testid="stExpander"] summary,
[data-testid="stExpander"] details > summary,
[data-testid="stExpander"] details > summary p,
[data-testid="stExpander"] [role="button"],
[data-testid="stExpander"] [role="button"] p,
[data-testid="stExpander"] [data-testid="stExpanderToggleIcon"] {
    color: #64748b !important;
    font-weight: 400 !important;
    font-size: 0.85rem !important;
}
[data-testid="stExpander"] summary,
[data-testid="stExpander"] details > summary,
[data-testid="stExpander"] details > summary p,
[data-testid="stExpander"] [role="button"] p {
    text-decoration: underline !important;
    text-decoration-color: #cbd5e1 !important;
    text-underline-offset: 0.25rem !important;
}
[data-testid="stExpander"] summary:hover,
[data-testid="stExpander"] details > summary:hover,
[data-testid="stExpander"] details > summary:hover p,
[data-testid="stExpander"] [role="button"]:hover,
[data-testid="stExpander"] [role="button"]:hover p,
[data-testid="stExpander"] summary:hover [data-testid="stExpanderToggleIcon"] {
    color: #1B4F8C !important;
    text-decoration-color: #1B4F8C !important;
    cursor: pointer !important;
}

/* ── Inline help icon ── */
/* Small filled INL-blue circle with a white "?". On hover, a navy     */
/* tooltip with the descriptive text appears below. Reuses the         */
/* existing 999-px pill radius (no new shape). Native title="" attr    */
/* is set alongside so touch / accessibility tools still work.         */
.mouse-help-icon {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 0.95rem;
    height: 0.95rem;
    background: #1B4F8C;
    color: white !important;
    font-size: 0.65rem;
    font-weight: 700;
    border-radius: 999px;
    margin-left: 0.35rem;
    cursor: help;
    position: relative;
    vertical-align: middle;
    flex-shrink: 0;
    line-height: 1;
    text-transform: none;
    letter-spacing: 0;
}
.mouse-help-icon:hover::after,
.mouse-help-icon:focus::after {
    content: attr(data-help);
    position: absolute;
    top: calc(100% + 0.4rem);
    left: 0;
    z-index: 100;
    background: #0a2540;
    color: white;
    padding: 0.6rem 0.85rem;
    border-radius: 4px;
    font-size: 0.85rem;
    font-weight: 400;
    line-height: 1.45;
    width: 260px;
    max-width: 260px;
    text-transform: none;
    letter-spacing: normal;
    text-align: left;
    white-space: normal;
    box-shadow: 0 2px 8px rgba(10, 37, 64, 0.15);
}

/* ── Popover trigger button ── */
/* Restyle the popover button so its label reads like a caption        */
/* (st.caption styling: 0.85rem, slate #64748b, normal weight) instead  */
/* of a bold primary button. The popover behavior is unchanged the     */
/* user clicks the label to expand and read the long explanation. This  */
/* makes "How is Fuel Lifetime estimated?" visually consistent with the */
/* k_eff caption underneath it.                                          */
[data-testid="stPopover"] button,
.stPopover > button,
button[kind="secondary"][data-testid="baseButton-secondary"][aria-haspopup] {
    background: transparent !important;
    border: none !important;
    color: #64748b !important;
    font-size: 0.85rem !important;
    font-weight: 400 !important;
    padding: 0.2rem 0 !important;
    box-shadow: none !important;
    text-decoration: underline !important;
    text-decoration-color: #cbd5e1 !important;
    text-underline-offset: 0.25rem !important;
}
[data-testid="stPopover"] button:hover,
.stPopover > button:hover {
    color: #1B4F8C !important;
    text-decoration-color: #1B4F8C !important;
}

/* ── Download button ── */
.stDownloadButton > button {
    background: #1B4F8C !important;
    color: white !important;
    border: none !important;
    border-radius: 8px !important;
    font-weight: 600 !important;
}
.stDownloadButton > button:hover {
    background: #0a2540 !important;
}

/* ── Spinner ── */
[data-testid="stSpinner"] > div {
    border-top-color: #1B4F8C !important;
}

/* ── Dataframe ── */
[data-testid="stDataFrame"] { border-radius: 8px; overflow: hidden; }
"""

# ---------------------------------------------------------------------------
# Page config (must be first Streamlit call)
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title='MOUSE Microreactor Cost Estimator',
    layout='wide',
)


@st.cache_data(show_spinner=False)
def _load_logo_b64(path):
    """Read a PNG and return its base64-encoded string for inline embedding."""
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('ascii')


@st.cache_data(show_spinner=False)
def _load_state_industrial_price_distribution():
    """Per-state industrial retail price distribution (cents/kWh).

    Precomputed from EIA Form 861 (Sales_Ult_Cust_2024.xlsx) by
    cost_drivers_and_markets/build_state_industrial_price_distribution.py.
    One row per state with n_utilities, weighted_mean (the EIA Table 5.6.A
    sales-weighted definition), and the empirical distribution of
    per-utility prices: min / p10 / q1 / median / q3 / p90 / max.
    """
    _path = (Path(__file__).resolve().parent.parent
             / 'assets' / 'Ref_Results'
             / 'state_industrial_price_distribution_2024.csv')
    return pd.read_csv(_path)

cookies = _get_cookie_manager()
if not cookies.ready():
    st.stop()

analytics_conn = _get_analytics_conn()
anonymous_id = _get_or_create_anonymous_id(cookies)

# ---------------------------------------------------------------------------
# Main app wrapped in analytics tracker
# ---------------------------------------------------------------------------
with streamlit_analytics.track():

    st.markdown(f'<style>{_CSS}</style>', unsafe_allow_html=True)

    # ── Sidebar inputs ──────────────────────────────────────────────────────
    with st.sidebar:
        # MOUSE logo in the sidebar header. Uses the black PNG since
        # Streamlit's sidebar background is light grey by default.
        st.markdown(
            f'<div style="text-align:center;padding:0.5rem 0 1rem 0;">'
            f'<img src="data:image/png;base64,{_load_logo_b64("assets/logos/MOUSE-logo_R1_black.png")}" '
            f'alt="MOUSE" style="max-width:100%;height:auto;">'
            f'<div style="color:#64748b;font-size:0.85rem;letter-spacing:0.06em;'
            f'text-transform:uppercase;line-height:1.4;margin-top:0.6rem;">'
            f'Microreactor Optimization<br>Using Simulation &amp; Economics</div>'
            f'</div>',
            unsafe_allow_html=True,
        )

        st.divider()
        st.markdown('**A Reactor Design**')

        reactor_label = st.selectbox(
            'Reactor Type',
            options=list(_REACTOR_LABELS.values()),
            help='Select a microreactor design to estimate costs for.',
        )
        reactor_type = _LABEL_TO_KEY[reactor_label]

        _log_visit_once_per_session(analytics_conn, anonymous_id, reactor_type=reactor_type, page_name='main')

        # Per-reactor enrichment floor. All three reactors expose a 5%
        # floor; values below the data-validated band (HPMR parametric
        # study covers E in [0.10, 0.1975]) are extrapolated and let
        # users explicitly see the subcritical region.
        _enrichment_min = {'LTMR': 0.05, 'GCMR': 0.05, 'HPMR': 0.05}
        _e_min = _enrichment_min.get(reactor_type, 0.05)
        # Clamp the default to the per-reactor min in case a previous
        # session left a lower value cached.
        _e_default = max(_e_min, 0.1975)
        enrichment = st.slider(
            'Enrichment (U-235 fraction)',
            min_value=_e_min,
            max_value=0.1975,
            value=_e_default,
            step=0.0025,
            format='%.4f',
            key=f'enrichment_{reactor_type}',
            help=('U-235 enrichment fraction. Affects uranium masses and fuel '
                  'lifetime. Higher enrichment means more SWU and higher cost; an '
                  'SWU premium multiplier of 1.15 is applied above 10% enrichment. '
                  'Down-blending (diluting higher-enriched U) is not considered.'),
        )

        _power_defaults = {'LTMR': 20, 'GCMR': 15, 'HPMR': 5}
        _power_max = {'LTMR': 64, 'GCMR': 50, 'HPMR': 60}

        power_mwt = st.slider(
            'Thermal Power (MWₜ)',
            min_value=1,
            max_value=_power_max[reactor_type],
            value=_power_defaults[reactor_type],
            step=1,
            key=f'power_{reactor_type}',
            help='Thermal power output.',
        )

        # LTMR / GCMR / HPMR: extra geometry inputs (diameter + active height).
        n_rings_per_assembly = None # LTMR only
        n_assembly_rings = None # GCMR (varies) and HPMR (locked at 6)
        n_core_rings = None # GCMR / HPMR
        active_height = None
        if reactor_type == 'LTMR':
            # Default to N=12 (95 cm), a mid-range trained geometry.
            _default_diameter_label = next(
                lbl for lbl in LTMR_DIAMETER_LABELS
                if LTMR_DIAMETER_LABEL_TO_N[lbl] == 12
            )
            _diameter_label = st.select_slider(
                'Active Core Diameter',
                options=LTMR_DIAMETER_LABELS,
                value=_default_diameter_label,
                key='ltmr_diameter',
                help=('Active core diameter (does NOT include the radial reflector). '
                      'Discrete values mapped to the number of fuel rings per '
                      'assembly. Values marked with * are interpolated between '
                      'trained geometries.'),
            )
            n_rings_per_assembly = LTMR_DIAMETER_LABEL_TO_N[_diameter_label]

            _ar_ltmr = LTMR_N_TO_ACTIVE_RADIUS_CM[n_rings_per_assembly]
            _ad_ltmr = LTMR_N_TO_DIAMETER_CM[n_rings_per_assembly] # active diameter
            _h_min = max(1, int(round(ASPECT_RATIO_MIN * _ad_ltmr)))
            _h_max = int(round(ASPECT_RATIO_MAX * _ad_ltmr))
            _h_default = int(round(_ad_ltmr)) # H/D = 1.0

            active_height = st.slider(
                'Active Height (cm)',
                min_value=_h_min,
                max_value=_h_max,
                value=_h_default,
                step=1,
                key=f'ltmr_active_height_{n_rings_per_assembly}',
                help=(f'Active fuel height in cm. Bounds correspond to aspect ratio '
                      f'(H / D, where D is the active core diameter) between '
                      f'{ASPECT_RATIO_MIN} and {ASPECT_RATIO_MAX}. For this geometry '
                      f'the active core diameter is {_ad_ltmr} cm.'),
            )

        elif reactor_type == 'GCMR':
            # Default to (N_A=6, N_C=5), the reference GCMR design
            _default_label = next(
                lbl for lbl in GCMR_DIAMETER_LABELS
                if GCMR_DIAMETER_LABEL_TO_PAIR[lbl] == (6, 5)
            )
            _diameter_label = st.select_slider(
                'Active Core Diameter',
                options=GCMR_DIAMETER_LABELS,
                value=_default_label,
                key='gcmr_diameter',
                help=('Active core diameter (does NOT include the radial reflector). '
                      'Discrete values mapped to (Assembly Rings, Core Rings). Values '
                      'marked with * are interpolated between trained geometries.'),
            )
            n_assembly_rings, n_core_rings = GCMR_DIAMETER_LABEL_TO_PAIR[_diameter_label]

            _ar_gcmr = _gcmr_active_radius(n_assembly_rings, n_core_rings)
            _ad_gcmr = 2.0 * _ar_gcmr # active diameter
            _h_min = max(1, int(round(ASPECT_RATIO_MIN * _ad_gcmr)))
            _h_max = int(round(ASPECT_RATIO_MAX * _ad_gcmr))
            _h_default = int(round(_ad_gcmr)) # H/D = 1.0

            active_height = st.slider(
                'Active Height (cm)',
                min_value=_h_min,
                max_value=_h_max,
                value=_h_default,
                step=1,
                key=f'gcmr_active_height_{n_assembly_rings}_{n_core_rings}',
                help=(f'Active fuel height in cm. Bounds correspond to aspect ratio '
                      f'(H / D, where D is the active core diameter) between '
                      f'{ASPECT_RATIO_MIN} and {ASPECT_RATIO_MAX}. For this geometry '
                      f'the active core diameter is {_ad_gcmr:.0f} cm.'),
            )

        elif reactor_type == 'HPMR':
            # HPMR's parametric study has N_A locked at 6, so the
            # active diameter slider varies only N_C. H is selected
            # independently with H/D in [ASPECT_RATIO_MIN, ASPECT_RATIO_MAX],
            # matching LTMR / GCMR.
            _default_label = next(
                lbl for lbl in HPMR_DIAMETER_LABELS
                if HPMR_DIAMETER_LABEL_TO_NC[lbl] == 5
            )
            _diameter_label = st.select_slider(
                'Active Core Diameter',
                options=HPMR_DIAMETER_LABELS,
                value=_default_label,
                key='hpmr_diameter',
                help=('Active core diameter (does NOT include the radial '
                      'reflector). Discrete values mapped to Core Rings '
                      '(N_C); Assembly Rings (N_A) is locked at 6 in the '
                      'HPMR parametric study.'),
            )
            n_assembly_rings = HPMR_NA_FIXED # always 6
            n_core_rings = HPMR_DIAMETER_LABEL_TO_NC[_diameter_label]

            _ar_hpmr = _hpmr_active_radius(n_core_rings)
            _ad_hpmr = 2.0 * _ar_hpmr # active diameter
            _h_min = max(1, int(round(ASPECT_RATIO_MIN * _ad_hpmr)))
            _h_max = int(round(ASPECT_RATIO_MAX * _ad_hpmr))
            _h_default = int(round(_ad_hpmr)) # H/D = 1.0

            active_height = st.slider(
                'Active Height (cm)',
                min_value=_h_min,
                max_value=_h_max,
                value=_h_default,
                step=1,
                key=f'hpmr_active_height_{n_core_rings}',
                help=(f'Active fuel height in cm. Bounds correspond to aspect '
                      f'ratio (H / D, where D is the active core diameter) '
                      f'between {ASPECT_RATIO_MIN} and {ASPECT_RATIO_MAX}. For '
                      f'this geometry the active core diameter is {_ad_hpmr:.0f} cm.'),
            )

        st.divider()
        st.markdown('**B Operation Parameters**')

        operation_mode = st.selectbox(
            'Operation Mode',
            options=['Remotely Monitored', 'On-Site Staffed'],
            help=(
                '**Remotely Monitored:** Operators monitor the reactor remotely and are '
                'required on-site only for emergencies or shutdown.\n\n'
                '**On-Site Staffed:** Operators must be physically present in the control room 24/7.'
            ),
        )
        emergency_shutdowns = st.slider(
            'Number of emergency shutdowns per year',
            min_value=0.1, max_value=10.0, value=2.0, step=0.1, format='%.1f',
            help=('Average number of unplanned shutdowns per year. Often higher '
                  'early in operation and decreases as operating experience accumulates.'),
        )
        startup_duration = st.slider(
            'Startup Duration after Emergency Shutdown (days)',
            min_value=1, max_value=180, value=21, step=1,
            help=('Days the reactor is offline after an unplanned emergency shutdown. '
                  'Varies by event; this is a rough average.'),
        )
        startup_duration_refueling = st.slider(
            'Startup Duration after Refueling (days)',
            min_value=1, max_value=30, value=14, step=1,
            help=('Planned days to bring the reactor back to full power after a '
                  'refueling shutdown at end of fuel lifetime.'),
        )

        st.divider()
        st.markdown('**C Economic Parameters**')

        debt_to_equity = st.slider(
            'Debt-to-Equity Ratio',
            min_value=0.0, max_value=5.0, value=1.0, step=0.1, format='%.1f',
            help='Ratio of debt to equity financing (e.g. 1.0 = equal debt and equity).',
        )
        interest_rate = st.slider(
            'Interest Rate (%)',
            min_value=0.0, max_value=20.0, value=7.0, step=0.5, format='%.1f',
            help=(
                "Annual interest rate on the loan used to build the reactor (the project's "
                "cost of debt). Used only to compute Interest During Construction (IDC, "
                "Account 62): only the debt-financed portion of the OCC (set via "
                "Debt-to-Equity Ratio) accrues this interest, accumulating over the "
                "construction period. Higher rate means larger IDC and higher Total "
                "Capital Investment. Does NOT appear in LCOE. Typical range: 5% to 12% "
                "for nuclear projects."
            ),
        )
        discount_rate = st.slider(
            'Discount Rate (%)',
            min_value=3.0, max_value=15.0, value=7.0, step=0.5, format='%.1f',
            help=(
                "Annual rate used to convert future cash flows into present value. This "
                "is the project's Weighted Average Cost of Capital (WACC), the blended "
                "return required by debt and equity investors. Applied to the LCOE/LCOH "
                "calculation (discounting future revenues and costs) and to Capital "
                "Recovery Factors (annualizing replacements and fuel costs). A higher "
                "rate makes near-term costs weigh more and far-future revenues weigh "
                "less, raising LCOE. Should be greater than or equal to the Interest "
                "Rate. Typical: 3% to 7% (public or government) or 8% to 15% (private "
                "nuclear)."
            ),
        )
        construction_duration = st.slider(
            'Construction Duration (months)',
            min_value=1, max_value=120, value=12, step=1,
            help='Number of months from ground-break to commercial operation.',
        )
        plant_lifetime = st.slider(
            'Plant Lifetime (years)',
            min_value=10, max_value=100, value=60, step=1,
            help=(
                'Economic lifetime of the plant used to levelize costs over time. '
                'A longer lifetime spreads capital costs over more years, reducing LCOE. '
                'Minimum: 10 years. Maximum: 100 years.'
            ),
        )

        st.markdown('**Government Subsidy (IRA Tax Credits)**')
        tax_credit_type = st.selectbox(
            'Tax Credit',
            options=['None', 'PTC', 'ITC'],
            index=2,
            help=(
                "**None:** no tax credit applied.\n\n"
                "**PTC (Production Tax Credit):** a per-MWh credit earned for "
                "each MWh of electricity sold during the credit period. Reduces LCOE.\n\n"
                "**ITC (Investment Tax Credit):** a one-time credit applied to "
                "the Overnight Capital Cost (OCC) of the plant. Reduces OCC.\n\n"
                "Note: ITC and PTC are mutually exclusive; only one can be selected per project."
            ),
        )
        tax_credit_value = None
        if tax_credit_type == 'PTC':
            tax_credit_value = st.select_slider(
                'PTC Credit Value ($/MWh)',
                options=[3.0, 3.3, 3.6, 15.0, 16.5, 18.0],
                value=15.0,
                format_func=lambda x: f'${x:.1f}/MWh',
                help=(
                    "The Production Tax Credit (PTC) is a per-MWh credit earned for every "
                    "MWh of electricity produced and sold during the credit period "
                    "(typically 10 years). Under the Inflation Reduction Act (IRA, "
                    "Section 45Y), advanced nuclear facilities placed in service after "
                    "December 31, 2024 may qualify for the Clean Electricity PTC.\n\n"
                    "**Base credit rate:**\n"
                    "* **$3/MWh:** when prevailing wage requirements are NOT met.\n"
                    "* **$15/MWh:** when prevailing wage and apprenticeship requirements "
                    "ARE met (5x multiplier).\n\n"
                    "**Optional stackable bonuses (added on top of the base rate):**\n"
                    "* **+10%** for domestic content (US-made iron, steel, manufactured products).\n"
                    "* **+10%** for siting in an 'energy community' (areas affected by coal "
                    "plant closures or fossil fuel employment decline).\n\n"
                    "**Typical values to choose from:**\n"
                    "* **$3.0/MWh:** base, no bonuses, no prevailing wage.\n"
                    "* **$3.3/MWh:** base + one 10% bonus.\n"
                    "* **$3.6/MWh:** base + both bonuses.\n"
                    "* **$15.0/MWh:** prevailing wage met, no bonuses.\n"
                    "* **$16.5/MWh:** prevailing wage + one 10% bonus.\n"
                    "* **$18.0/MWh:** prevailing wage + both bonuses."
                ),
            )
        elif tax_credit_type == 'ITC':
            tax_credit_value = st.select_slider(
                'ITC Credit Level',
                options=[0.06, 0.30, 0.40, 0.50],
                value=0.30,
                format_func=lambda x: f'{x*100:.0f}%',
                help=(
                    "The Investment Tax Credit (ITC) is a one-time credit that reduces "
                    "the Overnight Capital Cost (OCC) of the plant. Under the Inflation "
                    "Reduction Act (IRA, Section 48E), advanced nuclear facilities placed "
                    "in service after December 31, 2024 may qualify for the Clean "
                    "Electricity ITC.\n\n"
                    "**The level you can claim depends on which requirements your project meets:**\n"
                    "* **6%:** base rate, when prevailing wage requirements are NOT met.\n"
                    "* **30%:** when prevailing wage and apprenticeship requirements ARE met.\n"
                    "* **40%:** 30% + 10% bonus for domestic content (US-made iron, steel, "
                    "manufactured products).\n"
                    "* **50%:** 30% + 10% domestic content + 10% bonus for siting in an "
                    "'energy community' (areas affected by coal plant closures or fossil "
                    "fuel employment decline)."
                ),
            )

        # IRA sunset cutoff — only relevant when a credit is selected.
        # FOAK = unit 1, NOAK column = unit 'NOAK Unit Number'. A unit
        # qualifies only if its position in the deployment sequence is
        # <= this cutoff. Past the cutoff, ITC/PTC-adjusted outputs
        # fall back to the un-subsidized values.
        tax_credit_units = None
        if tax_credit_type in ('PTC', 'ITC'):
            tax_credit_units = st.slider(
                f'Number of Units Claiming {tax_credit_type}',
                min_value=1,
                max_value=50,
                value=10,
                step=1,
                help=(
                    f"Number of units in the order book that may claim the selected "
                    f"{tax_credit_type} before the IRA sunset year is reached. A unit is "
                    f"eligible only if its position in the deployment sequence is "
                    f"≤ this cutoff.\n\n"
                    f"When a unit is past the cutoff, the {tax_credit_type}-adjusted metrics "
                    f"fall back to the un-subsidized values, producing a step in "
                    f"the LCOE-vs-deployment-scale curve at the sunset point."
                ),
            )

        st.divider()
        run_button = st.button('⚡ Run Analysis', type='primary', width='stretch')
        if run_button:
            # Run a generational GC pass before kicking off a new
            # cost-engine run. Frees DataFrames from the previous Run
            # that Python's incremental GC hasn't reclaimed yet. Keeps
            # the Streamlit Cloud memory footprint from drifting up
            # over a long session of repeated Runs.
            import gc as _gc
            _gc.collect()
            st.session_state.has_run = True
            # Snapshot every input the downstream compute / render reads.
            # Sidebar widgets keep updating as the user tinkers, but results
            # are frozen to whatever was committed on the last click — so
            # tweaking a slider no longer auto-triggers a recompute.
            st.session_state.committed_inputs = {
                'reactor_type': reactor_type,
                'enrichment': enrichment,
                'power_mwt': power_mwt,
                'n_rings_per_assembly': n_rings_per_assembly,
                'active_height': active_height,
                'n_assembly_rings': n_assembly_rings,
                'n_core_rings': n_core_rings,
                'operation_mode': operation_mode,
                'emergency_shutdowns': emergency_shutdowns,
                'startup_duration': startup_duration,
                'startup_duration_refueling': startup_duration_refueling,
                'debt_to_equity': debt_to_equity,
                'interest_rate': interest_rate,
                'discount_rate': discount_rate,
                'construction_duration': construction_duration,
                'plant_lifetime': plant_lifetime,
                'tax_credit_type': tax_credit_type,
                'tax_credit_value': tax_credit_value,
                'tax_credit_units': tax_credit_units,
            }

        st.divider()
        st.markdown('**💬 Feedback**')
        st.caption('Help us improve MOUSE by sharing your thoughts.')
        st.link_button(
            '📝 Give Feedback',
            'https://qualtricsxm69xy9s7vm.qualtrics.com/jfe/form/SV_4Pb0vub9xCcsVV4',
            width='stretch',
        )

        if st.secrets.get("SHOW_ANALYTICS_PANEL", False):
            st.divider()
            _render_analytics_sidebar(analytics_conn)

        # ── Subtle memory badge ──────────────────────────────────
        # Last thing in the sidebar. Intentionally very small and
        # faint slate text on the slate-grey sidebar; a casual user
        # won't notice it, but a heavy user can glance to see how
        # close the process is to the 1 GB Streamlit Cloud ceiling.
        try:
            import psutil as _psutil_sidebar
            _rss_mb = _psutil_sidebar.Process().memory_info().rss / (1024 * 1024)
            st.markdown(
                f'<div style="font-size:0.7rem;color:#94a3b8;'
                f'text-align:right;margin-top:1rem;">'
                f'RAM {int(_rss_mb)} / 1024 MB</div>',
                unsafe_allow_html=True,
            )
        except Exception:
            pass

    # ── Memory monitor ──────────────────────────────────────────────────────
    # Must run BEFORE any st.stop() upstream can short-circuit the rerun.
    # Always cheap (gc + plt.close + malloc_trim); above 800 MB also sweeps
    # caches; above 950 MB renders a friendly stop banner.
    _run_memory_monitor()

    # ── Welcome banner ──────────────────────────────────────────────────────
    # Wrap all welcome content in a single slot so we can explicitly
    # clear it when the user clicks Run; otherwise Streamlit can leave
    # the previous rerun's welcome content visible during the long
    # computation, alongside the new "Computing all results" banner.
    welcome_slot = st.empty()

    if not st.session_state.get('has_run', False):
      with welcome_slot.container():
        # ── INL / repository credit (rendered ABOVE the welcome banner
        # so the repo link is the first thing the user sees) ──────────
        st.markdown(
            f'''<div style="background:#f1f3f5;border:1px solid #94a3b8;
                           border-left:4px solid #1B4F8C;border-radius:8px;
                           padding:1.8rem 2.2rem;color:#3c4257;
                           margin-bottom:1.5rem;">
                 <p style="font-size:1rem;line-height:1.55;margin:0 0 1.1rem;color:#3c4257;
                           padding-bottom:1.1rem;border-bottom:1px solid #93c5fd;">
                   This webapp is built on the <strong style="color:#0a2540;">MOUSE</strong> tool by
                   <strong style="color:#0a2540;">Idaho National Laboratory</strong> at
                   <a href="https://github.com/IdahoLabResearch/MOUSE" target="_blank"
                      style="color:#1B4F8C;font-weight:600;text-decoration:none;">
                     github.com/IdahoLabResearch/MOUSE</a>, designed for
                   <strong style="color:#0a2540;">accessibility and ease of use</strong>. Use the repository
                   directly for full control over custom geometries, materials, and advanced analyses.
                   <strong style="color:#0a2540;">Contributions via pull requests are welcome:</strong>
                   new reactor designs, improved cost data, or feature enhancements.
                 </p>
                 <img src="data:image/png;base64,{_load_logo_b64("assets/logos/MOUSE-logo_R1_black.png")}"
                      alt="MOUSE"
                      style="height:64px;width:auto;display:block;margin:0 0 0.4rem;">
                 <div style="font-size:1rem;font-weight:600;color:#64748b;margin-bottom:0.7rem;">
                   Microreactor Optimization Using Simulation and Economics
                 </div>
                 <p style="font-size:1rem;margin:0 0 1.1rem;max-width:780px;color:#3c4257;line-height:1.55;">
                   MOUSE is a first-of-its-kind framework that connects nuclear microreactor design choices
                   to bottom-up technoeconomic outcomes, helping engineers, investors, and policy makers see
                   how geometry, fuel, enrichment, and operating assumptions drive capital cost and LCOE before
                   major decisions are locked in. The web app requires no simulation environment or code setup:
                   it serves precomputed OpenMC-based parametric studies through fast nearest-neighbor interpolation,
                   updating estimates within seconds as users adjust LTMR, GCMR, and HPMR inputs. Its five sections cover
                   design and cost summaries, transportability, scoping metrics, cost decomposition, and market
                   competitiveness.
                 </p>
                 <div style="display:flex;gap:1rem;flex-wrap:wrap;">
                   <div style="background:white;border:1px solid #bfdbfe;border-radius:8px;padding:0.85rem 1.1rem;">
                     <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;
                                 font-weight:700;color:#1B4F8C;margin-bottom:0.3rem;">Reactor Types</div>
                     <div style="font-weight:600;font-size:1rem;color:#0a2540;">LTMR · GCMR · HPMR</div>
                     <div style="font-size:0.85rem;color:#64748b;margin-top:0.25rem;">Liquid Metal · Gas Cooled · Heat Pipe</div>
                   </div>
                   <div style="background:white;border:1px solid #bfdbfe;border-radius:8px;padding:0.85rem 1.1rem;">
                     <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;
                                 font-weight:700;color:#1B4F8C;margin-bottom:0.3rem;">Costs</div>
                     <div style="font-weight:600;font-size:1rem;color:#0a2540;">OCC · TCI · LCOE · LCOH · LCOF</div>
                     <div style="font-size:0.85rem;color:#64748b;margin-top:0.25rem;">Bottom-up estimation · Cost drivers · IRA credits</div>
                   </div>
                   <div style="background:white;border:1px solid #bfdbfe;border-radius:8px;padding:0.85rem 1.1rem;">
                     <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;
                                 font-weight:700;color:#1B4F8C;margin-bottom:0.3rem;">Neutronics &amp; Thermal Hydraulics at a glance</div>
                     <div style="font-weight:600;font-size:1rem;color:#0a2540;">Peaking Factor · Leakage · Power Density · Coolant Inventory</div>
                     <div style="font-size:0.85rem;color:#64748b;margin-top:0.25rem;">first-order scoping</div>
                   </div>
                   <div style="background:white;border:1px solid #bfdbfe;border-radius:8px;padding:0.85rem 1.1rem;">
                     <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;
                                 font-weight:700;color:#1B4F8C;margin-bottom:0.3rem;">Fuel Cycle</div>
                     <div style="font-weight:600;font-size:1rem;color:#0a2540;">U-235 / U-238 Mass · Lifetime · Discharge Burnup</div>
                   </div>
                   <div style="background:white;border:1px solid #bfdbfe;border-radius:8px;padding:0.85rem 1.1rem;">
                     <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;
                                 font-weight:700;color:#1B4F8C;margin-bottom:0.3rem;">Transportability</div>
                     <div style="font-weight:600;font-size:1rem;color:#0a2540;">Component Dims · Mass · Truck · Rail · Sea</div>
                     <div style="font-size:0.85rem;color:#64748b;margin-top:0.25rem;">first-order scoping</div>
                   </div>
                   <div style="background:white;border:1px solid #bfdbfe;border-radius:8px;padding:0.85rem 1.1rem;">
                     <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;
                                 font-weight:700;color:#1B4F8C;margin-bottom:0.3rem;">Costs in Perspective</div>
                     <div style="font-weight:600;font-size:1rem;color:#0a2540;">NOAK LCOE vs. market benchmarks</div>
                     <div style="font-size:0.85rem;color:#64748b;margin-top:0.25rem;">vs. wholesale and retail electricity prices</div>
                   </div>
                 </div>
               </div>''',
            unsafe_allow_html=True,
        )

        st.markdown(
            '''<div style="background:#fffbeb;border:1px solid #fcd34d;
                           border-left:4px solid #f59e0b;border-radius:8px;
                           padding:1.1rem 1.25rem;margin-bottom:1rem;">
                 <div style="display:flex;align-items:center;gap:0.5rem;margin-bottom:0.75rem;">
                   <span style="font-size:1.15rem;">⚠️</span>
                   <span style="font-size:1rem;font-weight:700;color:#92400e;
                                text-transform:uppercase;letter-spacing:0.07em;">
                     Important Caveats Please Read Before Use
                   </span>
                 </div>
                 <ul style="margin:0;padding-left:1.2rem;color:#92400e;font-size:1rem;line-height:1.75;">
                   <li><strong>Pre-conceptual designs:</strong> The three reactor designs (LTMR, GCMR, HPMR)
                       are pre-conceptual. They have not been fully optimized and do not represent
                       final or licensed configurations.</li>
                   <li><strong>Incomplete information:</strong> Cost estimates were developed with incomplete
                       engineering and procurement data. Significant uncertainties exist across all
                       accounts, particularly for novel or first-of-a-kind components.</li>
                   <li><strong>Not for investment decisions:</strong> These estimates must <em>not</em> be
                       used as the basis for financial, investment, or procurement decisions. They are
                       intended solely for research, screening, and comparative analysis.</li>
                 </ul>
               </div>''',
            unsafe_allow_html=True,
        )

        st.markdown(
            '''<div style="background:#eff6ff;border:1px solid #cbd5e1;
                           border-left:4px solid #1B4F8C;border-radius:8px;
                           padding:1.1rem 1.25rem;margin-top:1rem;margin-bottom:1rem;
                           font-size:1.15rem;font-weight:600;color:#1B4F8C;">
                 👈 How to use: configure your reactor in the sidebar on the left,
                 then click <strong>⚡ Run Analysis</strong> at the bottom of the sidebar.
               </div>''',
            unsafe_allow_html=True,
        )

        st.markdown(
            """
            <div style='text-align: center; font-size: 1rem; color: #64748b; padding-top: 2rem; padding-bottom: 1rem;'>
                © 2025 Battelle Energy Alliance, LLC. MOUSE is released under the MIT License.
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.stop()

    # When has_run=True, force-clear any leftover welcome content
    # from the previous rerun so only the computing banner is visible
    # during the long computation.
    welcome_slot.empty()

    # Compare live sidebar widgets against the last committed snapshot.
    # If anything has changed since the last Run click, blank the results
    # area and prompt the user to click Run again. This is what makes the
    # results disappear the moment a slider moves.
    _committed = st.session_state.committed_inputs
    _current_inputs = {
        'reactor_type': reactor_type,
        'enrichment': enrichment,
        'power_mwt': power_mwt,
        'n_rings_per_assembly': n_rings_per_assembly,
        'active_height': active_height,
        'n_assembly_rings': n_assembly_rings,
        'n_core_rings': n_core_rings,
        'operation_mode': operation_mode,
        'emergency_shutdowns': emergency_shutdowns,
        'startup_duration': startup_duration,
        'startup_duration_refueling': startup_duration_refueling,
        'debt_to_equity': debt_to_equity,
        'interest_rate': interest_rate,
        'discount_rate': discount_rate,
        'construction_duration': construction_duration,
        'plant_lifetime': plant_lifetime,
        'tax_credit_type': tax_credit_type,
        'tax_credit_value': tax_credit_value,
        'tax_credit_units': tax_credit_units,
    }
    if _current_inputs != _committed:
        st.markdown(
            '<div style="background:#eff6ff;border:1px solid #cbd5e1;'
            'border-left:4px solid #1B4F8C;border-radius:8px;'
            'padding:1.1rem 1.25rem;margin-top:1rem;margin-bottom:1rem;'
            'font-size:1.15rem;font-weight:600;color:#1B4F8C;">'
            '👈 Inputs changed. Click <strong>⚡ Run Analysis</strong> '
            'in the sidebar to update results.'
            '</div>',
            unsafe_allow_html=True,
        )
        st.stop()

    # Inputs match the snapshot — restore from committed to be explicit
    # about which values drive the compute / render below.
    reactor_type = _committed['reactor_type']
    enrichment = _committed['enrichment']
    power_mwt = _committed['power_mwt']
    n_rings_per_assembly = _committed['n_rings_per_assembly']
    active_height = _committed['active_height']
    n_assembly_rings = _committed['n_assembly_rings']
    n_core_rings = _committed['n_core_rings']
    operation_mode = _committed['operation_mode']
    emergency_shutdowns = _committed['emergency_shutdowns']
    startup_duration = _committed['startup_duration']
    startup_duration_refueling = _committed['startup_duration_refueling']
    debt_to_equity = _committed['debt_to_equity']
    interest_rate = _committed['interest_rate']
    discount_rate = _committed['discount_rate']
    construction_duration = _committed['construction_duration']
    plant_lifetime = _committed['plant_lifetime']
    tax_credit_type = _committed['tax_credit_type']
    tax_credit_value = _committed['tax_credit_value']
    tax_credit_units = _committed.get('tax_credit_units')

    # ── Show single progress banner covering BOTH the basic estimate
    # and the NOAK deployment-scale sweep that follows. ─────────────────────
    _precompute_slot = st.empty()
    _precompute_slot.markdown(
        '<div style="background:#fffbeb;'
        'border:2px solid #f59e0b;border-radius:8px;padding:1.8rem 2.2rem;'
        'margin-bottom:1rem;color:#92400e;text-align:center;">'
        '<div style="font-size:1.15rem;font-weight:700;margin-bottom:0.4rem;'
        'letter-spacing:0.02em;">⏳ Computing all results...</div>'
        '<div style="font-size:1rem;font-weight:500;line-height:1.5;">'
        'Roughly <strong>10-15 seconds</strong>. All results will appear at once when finished.'
        '</div></div>',
        unsafe_allow_html=True,
    )

    # ── Run cost estimate ───────────────────────────────────────────────────
    try:
        display_df, enriched_df, detailed_sorted_df, params = _run_estimate(
            reactor_type, power_mwt, enrichment,
            interest_rate / 100.0, discount_rate / 100.0, construction_duration, debt_to_equity,
            operation_mode, emergency_shutdowns, startup_duration, startup_duration_refueling,
            tax_credit_type, tax_credit_value, plant_lifetime,
            n_rings_per_assembly=n_rings_per_assembly,
            active_height=active_height,
            n_assembly_rings=n_assembly_rings,
            n_core_rings=n_core_rings,
            tax_credit_units=tax_credit_units,
        )
    except SubcriticalError as exc:
        _precompute_slot.empty()
        st.error('### ⚠ Reactor is Subcritical')
        st.warning(str(exc))
        ca, cb, cc = st.columns(3)
        _info_card(ca, 'Fuel Lifetime', '0 days', accent='#dc2626', bg='#fef2f2', border='#fecaca')
        _info_card(cb, 'Thermal Power', f'{power_mwt} MW<sub>t</sub>', accent='#9a3412', bg='#fffbeb', border='#fcd34d')
        _info_card(cc, 'Enrichment', f'{enrichment*100:.2f}%', accent='#9a3412', bg='#fffbeb', border='#fcd34d')
        st.info('No cost estimate is available for a subcritical operating point. '
                'Try reducing the power or increasing the enrichment.')

        st.markdown(
            """
            <div style='text-align: center; font-size: 1rem; color: #64748b; padding-top: 2rem; padding-bottom: 1rem;'>
                © 2025 Battelle Energy Alliance, LLC. MOUSE is released under the MIT License.
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.stop()
    except ShortLifetimeError as exc:
        _precompute_slot.empty()
        # Extract the estimated lifetime in days from the message ("only N days")
        import re
        m = re.search(r'only (\d+) days', str(exc))
        est_lt_days = int(m.group(1)) if m else 0
        st.error('### ⚠ Fuel Lifetime Too Short')
        st.warning(str(exc))
        ca, cb, cc = st.columns(3)
        _info_card(ca, 'Fuel Lifetime', f'{est_lt_days} days',
                   accent='#dc2626', bg='#fef2f2', border='#fecaca')
        _info_card(cb, 'Thermal Power', f'{power_mwt} MW<sub>t</sub>',
                   accent='#9a3412', bg='#fffbeb', border='#fcd34d')
        _info_card(cc, 'Enrichment', f'{enrichment*100:.2f}%',
                   accent='#9a3412', bg='#fffbeb', border='#fcd34d')
        st.info('No cost estimate is performed when the fuel lifetime is below 90 days. '
                'Try increasing the diameter, the height, or the enrichment.')

        st.markdown(
            """
            <div style='text-align: center; font-size: 1rem; color: #64748b; padding-top: 2rem; padding-bottom: 1rem;'>
                © 2025 Battelle Energy Alliance, LLC. MOUSE is released under the MIT License.
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.stop()
    except Exception as exc:
        _precompute_slot.empty()
        st.error(f'Cost estimation failed: {exc}')
        import traceback
        st.code(traceback.format_exc())

        st.markdown(
            """
            <div style='text-align: center; font-size: 1rem; color: #64748b; padding-top: 2rem; padding-bottom: 1rem;'>
                © 2025 Battelle Energy Alliance, LLC. MOUSE is released under the MIT License.
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.stop()

    # ── Extract key params ──────────────────────────────────────────────────
    fuel_lifetime = params.get('Fuel Lifetime', float('nan'))
    power_mwe = params.get('Power MWe', float('nan'))
    capacity_factor = params.get('Capacity Factor', float('nan'))

    try:
        _capacity_factor_value = float(capacity_factor)
    except (TypeError, ValueError):
        _capacity_factor_value = float('nan')

    if not math.isnan(_capacity_factor_value) and _capacity_factor_value < 0.10:
        _precompute_slot.empty()
        _cf_percent = max(0.0, _capacity_factor_value * 100.0)
        _cf_value_text = 'zero' if _capacity_factor_value < 0 else f'{_cf_percent:.1f}%'
        st.error('### ⚠ Capacity Factor Too Low')
        st.warning(
            f'The calculated capacity factor is {_cf_value_text}, which is too low for meaningful results.'
        )
        st.info(
            'No results are shown because cost and production outputs are not meaningful '
            'with this operating profile. Reduce outage assumptions such as emergency '
            'shutdown frequency or startup duration and run the analysis again.'
        )

        st.markdown(
            """
            <div style='text-align: center; font-size: 1rem; color: #64748b; padding-top: 2rem; padding-bottom: 1rem;'>
                © 2025 Battelle Energy Alliance, LLC. MOUSE is released under the MIT License.
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.stop()

    _fl_str = f'{fuel_lifetime / 365:.1f} yrs' if not math.isnan(float(fuel_lifetime)) else 'N/A'
    _fl_days = f'{int(fuel_lifetime):,} days' if not math.isnan(float(fuel_lifetime)) else ''
    _cf_str = f'{capacity_factor * 100:.1f}%' if not math.isnan(float(capacity_factor)) else 'N/A'

    # ── Collect all summary values ──────────────────────────────────────────
    if tax_credit_type == 'ITC':
        occ_account = 'OCC (ITC-adjusted)'
        tci_account = 'TCI (ITC-adjusted)'
        lcoe_account = 'LCOE (ITC-adjusted)'
    elif tax_credit_type == 'PTC':
        occ_account = 'OCC'
        tci_account = 'TCI'
        lcoe_account = 'LCOE with PTC'
    else:
        occ_account = 'OCC'
        tci_account = 'TCI'
        lcoe_account = 'LCOE'

    occ_f, occ_f_std = _get_mean_std(display_df, occ_account, 'FOAK')
    tci_f, tci_f_std = _get_mean_std(display_df, tci_account, 'FOAK')
    lcoe_f, lcoe_f_std = _get_mean_std(display_df, lcoe_account, 'FOAK')
    lcoh_f, lcoh_f_std = _get_mean_std(display_df, 'LCOH', 'FOAK')
    lcof_f, lcof_f_std = _get_lcof(enriched_df, 'FOAK')

    occ_n, occ_n_std = _get_mean_std(display_df, occ_account, 'NOAK')
    tci_n, tci_n_std = _get_mean_std(display_df, tci_account, 'NOAK')
    lcoe_n, lcoe_n_std = _get_mean_std(display_df, lcoe_account, 'NOAK')
    lcoh_n, lcoh_n_std = _get_mean_std(display_df, 'LCOH', 'NOAK')
    lcof_n, lcof_n_std = _get_lcof(enriched_df, 'NOAK')

    # ── Pre-compute the slow cost-engine sweep BEFORE rendering any
    # tab, so all results show up at once when the spinner clears.
    # Without this, the user sees half the tab populate, then a long
    # gap, then the rest as the LCOE sweep finishes (15-25s on a
    # cold cache).  Re-using _mid_results below in the Costs in
    # Perspective block keeps the cost-engine call count the same.
    _N_user_pre = int(round(float(params.get('NOAK Unit Number', 100))))
    if _N_user_pre < 2:
        _N_user_pre = 100
    _N_mids_pre = [2, 10]
    _mid_results = {}
    for _N in _N_mids_pre:
        if _N <= 1 or _N >= _N_user_pre:
            continue
        try:
            _m_i, _s_i, _, _ = _lcoe_at_noak_unit(
                reactor_type, power_mwt, enrichment,
                interest_rate / 100.0, discount_rate / 100.0,
                construction_duration,
                debt_to_equity, operation_mode, emergency_shutdowns,
                startup_duration, startup_duration_refueling,
                tax_credit_type, tax_credit_value, plant_lifetime,
                n_rings_per_assembly=n_rings_per_assembly,
                active_height=active_height,
                n_assembly_rings=n_assembly_rings,
                n_core_rings=n_core_rings,
                noak_unit_number=_N,
                tax_credit_units=tax_credit_units,
            )
            _mid_results[_N] = (_m_i, _s_i)
        except Exception as _e:
            st.warning(f'Could not compute N={_N} anchor: {_e}')
    _precompute_slot.empty()

    # ── Result hero banner (rendered AFTER all computation finishes) ────────
    _badge_style = ('background:#eff6ff;border:1px solid #bfdbfe;color:#1B4F8C;'
                    'border-radius:999px;font-size:0.85rem;font-weight:600;'
                    'padding:0.2rem 0.6rem;letter-spacing:0.06em;margin-left:0.6rem;')
    if tax_credit_type == 'PTC':
        _units_text = f' (for the first {int(tax_credit_units)} units)' if tax_credit_units is not None else ''
        _credit_badge = f'<span style="{_badge_style}">PTC ${tax_credit_value}/MWh{_units_text}</span>'
    elif tax_credit_type == 'ITC':
        _units_text = f' (for the first {int(tax_credit_units)} units)' if tax_credit_units is not None else ''
        _credit_badge = f'<span style="{_badge_style}">ITC {int(tax_credit_value*100)}%{_units_text}</span>'
    else:
        _credit_badge = ('<span style="background:#f1f3f5;border:1px solid #cbd5e1;'
                         'color:#64748b;border-radius:999px;font-size:0.85rem;'
                         'font-weight:600;padding:0.2rem 0.6rem;letter-spacing:0.06em;'
                         'margin-left:0.6rem;">No Tax Credit</span>')

    st.markdown(
        f'''<div style="background:#f1f3f5;border:1px solid #94a3b8;
                       border-left:4px solid #1B4F8C;border-radius:8px;
                       padding:1.8rem 2.2rem;color:#3c4257;
                       margin-bottom:1.5rem;">
             <div style="font-size:0.85rem;font-weight:600;letter-spacing:0.12em;
                         text-transform:uppercase;color:#64748b;margin-bottom:0.5rem;">
               Run Summary
             </div>
             <h2 style="font-size:1.5rem;font-weight:700;margin:0 0 1.1rem;color:#0a2540;line-height:1.2;">
               {reactor_label} {_credit_badge}
             </h2>
             <div style="font-size:0.85rem;font-weight:600;text-transform:uppercase;letter-spacing:0.09em;color:#94a3b8;margin:0 0 0.4rem;">Design</div>
             <div style="display:flex;gap:2rem;flex-wrap:wrap;margin-bottom:1rem;">
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Thermal Power</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{power_mwt} MW<sub>t</sub></div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Enrichment</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{enrichment*100:.2f} wt%</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Active Diameter</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{_diameter_label.replace(' *', '').strip()}</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Active Height</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{active_height} cm</div>
               </div>
             </div>
             <div style="font-size:0.85rem;font-weight:600;text-transform:uppercase;letter-spacing:0.09em;color:#94a3b8;margin:0 0 0.4rem;">Operations</div>
             <div style="display:flex;gap:2rem;flex-wrap:wrap;margin-bottom:1rem;">
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Operation Mode</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{operation_mode}</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Emergency Shutdowns / Yr</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{emergency_shutdowns:.1f}</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Emergency Outage</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{startup_duration} days</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Refuel Outage</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{startup_duration_refueling} days</div>
               </div>
             </div>
             <div style="font-size:0.85rem;font-weight:600;text-transform:uppercase;letter-spacing:0.09em;color:#94a3b8;margin:0 0 0.4rem;">Economics</div>
             <div style="display:flex;gap:2rem;flex-wrap:wrap;">
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Debt-to-Equity</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{debt_to_equity:.1f}</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Interest Rate</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{interest_rate:.1f}%</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Discount Rate</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{discount_rate:.1f}%</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Construction</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{construction_duration} months</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">Plant Lifetime</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{plant_lifetime} yrs</div>
               </div>
               <div>
                 <div style="font-size:0.85rem;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;">All Costs In</div>
                 <div style="font-weight:600;font-size:1rem;color:#0a2540;">{ESCALATION_YEAR} USD</div>
               </div>
             </div>
           </div>''',
        unsafe_allow_html=True,
    )

    # Tabs were removed; bands render in sequential document order.

    # Card queue shared across containers. Cards are appended during
    # the compute block in tab_design (via _fuel_card with a section=
    # tag), then dispatched into the right band when rendered.
    # 3 = Fuel cycle metrics (Band 3, third column)
    # 4 = Core neutronics (Band 3, first column)
    # 5 = Core thermal-hydraulics (Band 3, second column)
    _section_cards = {3: [], 4: [], 5: []}

    # ═══════════════════════════════════════════════════════════════
    # BAND 1 Design & Cost Headlines
    # Sections 1-3: What does it cost? · What does it look like? ·
    # How much fuel and how long does it last?
    # ═══════════════════════════════════════════════════════════════
    # Margin-top is smaller than the other bands because the Run
    # Summary card immediately above already carries its own
    # margin-bottom; otherwise the gap stacks visibly.
    _section_header(
        '01',
        'Design &amp; Cost Summary',
        'What it looks like, what it costs, how it performs',
        top_margin='0.5rem',
    )

    # ── Subsection 1: Geometry & Layout ─────────────────────────
    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:0 0 0.85rem 0;">Geometry &amp; Layout</div>',
        unsafe_allow_html=True,
    )
    geo_left, geo_right = st.columns([1, 1], gap='large')
    with geo_left:
        main_img, main_caption = _REACTOR_IMAGES[reactor_type]['main']

        # For LTMR, pick the per-N core cross-section that matches
        # the user's chosen diameter (Number of Rings per Assembly).
        if reactor_type == 'LTMR':
            _n = int(params.get('Number of Rings per Assembly', 12))
            _per_n_path = os.path.join(_ASSETS, f'LTMR_core_N{_n}.png')
            if os.path.exists(_per_n_path):
                main_img = _per_n_path
                main_caption = (
                    f"Simplified scoping-level design used for cost estimation, "
                    f"not a detailed engineering design. LTMR core cross-section "
                    f"for N = {_n} rings per assembly "
                    f"(active core radius {LTMR_N_TO_ACTIVE_RADIUS_CM[_n]:.1f} cm, "
                    f"active core diameter {LTMR_N_TO_DIAMETER_CM[_n]} cm). "
                    f"Hexagonal arrangement of UZrH alloy fuel pins and ZrH "
                    f"moderator pins cooled by NaK liquid metal, surrounded by a "
                    f"graphite radial reflector with 18 control drums."
                )

        # For GCMR, pick the per-(N_A, N_C) core cross-section.
        if reactor_type == 'GCMR':
            _na = int(params.get('Assembly Rings', 6))
            _nc = int(params.get('Core Rings', 5))
            _per_pair_path = os.path.join(
                _ASSETS, f'GCMR_core_NA{_na}_NC{_nc}.png'
            )
            if os.path.exists(_per_pair_path):
                main_img = _per_pair_path
                main_caption = (
                    f"Simplified scoping-level design used for cost estimation, "
                    f"not a detailed engineering design. GCMR core cross-section "
                    f"for (N_A={_na} assembly rings, N_C={_nc} core rings). "
                    f"Hexagonal fuel assemblies containing TRISO fuel compacts "
                    f"arranged in a honeycomb pattern, cooled by helium gas, "
                    f"with graphite reflector and control drums."
                )

        # For HPMR, pick the per-(N_A, N_C) core cross-section
        # (N_A is locked at 6 in the parametric study, but the
        # filename pattern matches the GCMR convention so future
        # N_A variation is supported).
        if reactor_type == 'HPMR':
            _na = int(params.get('Number of Rings per Assembly', 6))
            _nc = int(params.get('Number of Rings per Core', 5))
            _per_pair_path = os.path.join(
                _ASSETS, f'HPMR_core_NA{_na}_NC{_nc}.png'
            )
            if os.path.exists(_per_pair_path):
                main_img = _per_pair_path
                main_caption = (
                    f"Simplified scoping-level design used for cost estimation, "
                    f"not a detailed engineering design. HPMR core cross-section "
                    f"for (N_A={_na} assembly rings, N_C={_nc} core rings). "
                    f"TRISO fuel modeled as a homogenized fuel pin (TRISO particles "
                    f"smeared into the pin volume rather than resolved individually) "
                    f"interleaved with sodium heat pipes inside monolith-graphite "
                    f"moderator blocks, surrounded by a graphite reflector with "
                    f"12 control drums."
                )

        # Cross-section on top, side-view directly below (LTMR/GCMR
        # only both have geometry-driven height inputs). Both render
        # at the column width, so the side view's vertical extent
        # relative to the cross-section's diameter makes H/D visually
        # obvious.
        st.image(main_img, width='stretch')
        st.caption(main_caption)

        if reactor_type in ('LTMR', 'GCMR', 'HPMR'):
            # The figure has built-in right padding to mirror the
            # cross-section's legend area, so the side-view rectangle
            # aligns horizontally with the cross-section circle above.
            # Save with bbox_inches=None (st.pyplot's default 'tight'
            # would crop the right-side empty space, defeating the
            # alignment).
            _diam = float(2.0 * params['Core Radius'])
            _h = float(params['Active Height'])
            _ax_r = float(params.get('Axial Reflector Thickness', 0.0))
            _rad_r = float(params.get('Radial Reflector Thickness', 0.0))
            st.altair_chart(
                _side_view_altair_chart(_diam, _h, _ax_r, _rad_r),
                use_container_width=False,
            )

    with geo_right:
        # Materials & components table at the top of the right column,
        # followed by the aspect-ratio readout and the detailed cross-
        # section gallery (both moved here from under the side view).
        _materials_section(reactor_type, params)

        _active_radius_cm = (float(2.0 * params['Core Radius'])
                             - 2.0 * float(params.get('Radial Reflector Thickness', 0.0))
                             ) / 2.0
        _active_diam_cm = 2.0 * _active_radius_cm
        if _active_diam_cm > 0:
            _aspect = float(params['Active Height']) / _active_diam_cm
            st.markdown(
                f'<div style="margin-top:0.5rem;font-size:0.85rem;color:#0a2540;">'
                f'<strong>Aspect ratio (H/D):</strong> {_aspect:,.2f}'
                f'</div>',
                unsafe_allow_html=True,
            )
            with st.popover('Tell me more', width='content'):
                st.markdown(
                    'Active height &divide; active diameter, using the active fissioning '
                    'region only (no reflector, no vessel). A value of 1.0 means a '
                    '"cube-like" cylinder where height equals diameter, which '
                    'minimises neutron leakage per unit core volume. Tall, slender '
                    'cores (H/D > 1.5) lose more axial neutrons; short, fat cores '
                    '(H/D < 0.5) lose more radial neutrons. Most microreactor '
                    'designs sit between 0.7 and 1.5. The webapp constrains H/D '
                    'to [0.5, 2.0] via the active-height slider.'
                )

        with st.expander('View detailed cross-sections'):
            for img_path, img_caption in _REACTOR_IMAGES[reactor_type]['details']:
                # GCMR: replace the static fuel-assembly image with the
                # per-N_A version that matches the user's selection.
                if reactor_type == 'GCMR' and 'Fuel Assembly.png' in img_path:
                    _na = int(params.get('Assembly Rings', 6))
                    _per_na_path = os.path.join(
                        _ASSETS, f'GCMR_fuel_assembly_NA{_na}.png'
                    )
                    if os.path.exists(_per_na_path):
                        img_path = _per_na_path
                        img_caption = (
                            f"GCMR fuel assembly cross-section for N_A={_na} "
                            f"(assembly rings). TRISO fuel compacts, helium coolant "
                            f"channels, and ZrH moderator booster pins embedded in "
                            f"a graphite matrix."
                        )
                # HPMR: same pattern, swap to per-N_A assembly image.
                elif reactor_type == 'HPMR' and 'fuel_assembly' in img_path.lower() \
                                          and 'zoomed' not in img_path.lower():
                    _na = int(params.get('Number of Rings per Assembly', 6))
                    _per_na_path = os.path.join(
                        _ASSETS, f'HPMR_fuel_assembly_NA{_na}.png'
                    )
                    if os.path.exists(_per_na_path):
                        img_path = _per_na_path
                        img_caption = (
                            f"HPMR fuel assembly cross-section for N_A={_na} "
                            f"(assembly rings). Homogenized-TRISO fuel pins "
                            f"and sodium heat pipes embedded in a monolith-graphite "
                            f"moderator block."
                        )
                st.image(img_path, width='stretch')
                st.caption(img_caption)

    # ── Subsection 2: Performance ───────────────────────────────
    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:1.25rem 0 0.85rem 0;">Performance</div>',
        unsafe_allow_html=True,
    )
    perf_left, perf_right = st.columns([1, 1], gap='large')
    with perf_left:
        ic1, ic2, ic3 = st.columns(3)
        _info_card(ic1, 'Electric Power Output', f'{power_mwe:.1f} MW<sub>e</sub>',
                   subtitle=f'Thermal input: {power_mwt} MW<sub>t</sub>')
        _info_card(ic2, 'Fuel Lifetime', _fl_str, subtitle=_fl_days)
        _info_card(ic3, 'Capacity Factor', _cf_str,
                   subtitle='Accounts for refueling & shutdowns')

        with st.popover('How is Fuel Lifetime estimated?', width='content'):
            if reactor_type == 'LTMR':
                st.markdown(
                    'We pre-ran a parametric set of **OpenMC depletion simulations** '
                    'of the LTMR that swept active heights, active diameters '
                    '(ring counts), enrichments, and thermal powers. From each '
                    'simulation we extracted the time at which the core falls to '
                    'k_eff = 1 (end of cycle). For your inputs, the Fuel Lifetime '
                    'shown here is a **KNN local regression** over those simulations '
                    '(4 nearest cases, distance-weighted). Typical error is '
                    '**5-15%** for cases close to the trained design space, and '
                    'may be larger for ring counts marked * in the diameter slider '
                    '(outside the trained grid) or near the criticality boundary.'
                )
            elif reactor_type == 'GCMR':
                st.markdown(
                    'We pre-ran a parametric set of **OpenMC depletion simulations** '
                    'of the GCMR that swept active heights, (Assembly Rings, Core '
                    'Rings) geometry pairs, enrichments, and thermal powers. From '
                    'each simulation we extracted the time at which the core falls '
                    'to k_eff = 1 (end of cycle). For your inputs, the Fuel '
                    'Lifetime shown here is a **KNN local regression** over those '
                    'simulations (4 nearest cases, distance-weighted). Typical '
                    'error is **5-15%** for cases close to the trained design '
                    'space, and may be larger for (N_A, N_C) pairs marked * in '
                    'the diameter slider (outside the trained grid) or near the '
                    'criticality boundary.'
                )
            elif reactor_type == 'HPMR':
                st.markdown(
                    'We pre-ran a parametric set of **OpenMC depletion simulations** '
                    'of the HPMR (104 cases) that swept active heights, core '
                    'rings (N_C), enrichments, and thermal powers, with N_A '
                    'locked at 6. From each simulation we extracted the time at '
                    'which the core falls to k_eff = 1 (end of cycle). For your '
                    'inputs, the Fuel Lifetime shown here is a **KNN local '
                    'regression** over those simulations (4 nearest cases, '
                    'distance-weighted). Typical error is **5-15%** inside the '
                    'trained envelope (N_C ∈ [3..7], H = 136-1056 cm, E = 10-19.75%, '
                    'P = 1-60 MWₜ). Mass U-235 / U-238 are derived from the exact '
                    'HPMR mass formula, so off-grid geometry queries are mass-'
                    'consistent. The 8-10% enrichment band lets you see the '
                    'subcritical region explicitly &mdash; the model returns 0 '
                    'there because the nearest training point at E = 10% is itself '
                    'subcritical.'
                )

    with perf_right:
        # ── Reactivity vs Time (LTMR only) ──────────────────────────────
        _reactivity_swing_pct = None # filled in for LTMR if curve available
        if reactor_type == 'LTMR':
            _times, _keffs = get_ltmr_keff_curve(
                n_rings_per_assembly = params['Number of Rings per Assembly'],
                active_height = params['Active Height'],
                enrichment = params['Enrichment'],
                power_mwt = params['Power MWt'],
                anchor_lifetime_days = params.get('Fuel Lifetime', None),
            )
            if _times.size >= 2:
                _k_bol = float(_keffs[0])
                if _k_bol > 1.0:
                    _reactivity_swing_pct = (_k_bol - 1.0) / _k_bol * 100.0
            if _times.size >= 2:
                st.markdown(
                    '<div style="font-size:0.85rem;font-weight:600;color:#64748b;'
                    'text-transform:uppercase;letter-spacing:0.09em;'
                    'margin-bottom:0.6rem;">k_eff vs Time</div>',
                    unsafe_allow_html=True,
                )
                st.altair_chart(_keff_altair_chart(_times, _keffs),
                                use_container_width=True)
                st.caption(
                    '**k_eff** (the neutron multiplication factor) measures '
                    'whether the chain reaction sustains itself: > 1 it grows, '
                    '= 1 critical (steady), < 1 it dies down. The curve shown '
                    'comes from the parametric set of **OpenMC depletion '
                    'simulations** of the LTMR (sweeping active heights, '
                    'diameters, enrichments, and thermal powers). For your '
                    'inputs, k_eff at each timestep is the distance-weighted '
                    'average of the 4 nearest training cases. The time axis '
                    'is anchored so the k_eff = 1 crossing matches the Fuel '
                    'Lifetime above; the subcritical tail is omitted.'
                )
                st.markdown('<div style="height:1rem"></div>', unsafe_allow_html=True)

        # ── Reactivity vs Time (GCMR same pattern as LTMR) ────────────
        if reactor_type == 'GCMR':
            _times, _keffs = get_gcmr_keff_curve(
                assembly_rings = params['Assembly Rings'],
                core_rings = params['Core Rings'],
                active_height = params['Active Height'],
                enrichment = params['Enrichment'],
                power_mwt = params['Power MWt'],
                anchor_lifetime_days = params.get('Fuel Lifetime', None),
            )
            if _times.size >= 2:
                _k_bol = float(_keffs[0])
                if _k_bol > 1.0:
                    _reactivity_swing_pct = (_k_bol - 1.0) / _k_bol * 100.0
            if _times.size >= 2:
                st.markdown(
                    '<div style="font-size:0.85rem;font-weight:600;color:#64748b;'
                    'text-transform:uppercase;letter-spacing:0.09em;'
                    'margin-bottom:0.6rem;">k_eff vs Time</div>',
                    unsafe_allow_html=True,
                )
                st.altair_chart(_keff_altair_chart(_times, _keffs),
                                use_container_width=True)
                st.caption(
                    '**k_eff** (the neutron multiplication factor) measures '
                    'whether the chain reaction sustains itself: > 1 it grows, '
                    '= 1 critical (steady), < 1 it dies down. The curve shown '
                    'comes from the parametric set of **OpenMC depletion '
                    'simulations** of the GCMR (sweeping enrichment, assembly '
                    'and core rings, active height, and thermal power). For '
                    'your inputs, k_eff at each timestep is the distance-'
                    'weighted average of the 4 nearest training cases in the '
                    '(E, N_A, N_C, H, P) feature space. The time axis is '
                    'anchored so the k_eff = 1 crossing matches the Fuel '
                    'Lifetime above; the subcritical tail is omitted.'
                )
                st.markdown('<div style="height:1rem"></div>', unsafe_allow_html=True)

        # ── Reactivity vs Time (HPMR same pattern as LTMR / GCMR) ─────
        if reactor_type == 'HPMR':
            _times, _keffs = get_hpmr_keff_curve(
                n_rings_per_assembly = params['Number of Rings per Assembly'],
                n_rings_per_core = params['Number of Rings per Core'],
                active_height = params['Active Height'],
                enrichment = params['Enrichment'],
                power_mwt = params['Power MWt'],
                anchor_lifetime_days = params.get('Fuel Lifetime', None),
            )
            if _times.size >= 2:
                _k_bol = float(_keffs[0])
                if _k_bol > 1.0:
                    _reactivity_swing_pct = (_k_bol - 1.0) / _k_bol * 100.0
            if _times.size >= 2:
                st.markdown(
                    '<div style="font-size:0.85rem;font-weight:600;color:#64748b;'
                    'text-transform:uppercase;letter-spacing:0.09em;'
                    'margin-bottom:0.6rem;">k_eff vs Time</div>',
                    unsafe_allow_html=True,
                )
                st.altair_chart(_keff_altair_chart(_times, _keffs),
                                use_container_width=True)
                st.caption(
                    '**k_eff** (the neutron multiplication factor) measures '
                    'whether the chain reaction sustains itself: > 1 it grows, '
                    '= 1 critical (steady), < 1 it dies down. The curve shown '
                    'comes from the parametric set of **OpenMC depletion '
                    'simulations** of the HPMR (sweeping enrichment, core '
                    'rings N_C, active height, and thermal power; N_A is '
                    'locked at 6 in the training set). For your inputs, '
                    'k_eff at each timestep is the distance-weighted average '
                    'of the 4 nearest training cases in the (E, N_C, H, P) '
                    'feature space. The time axis is anchored so the k_eff = '
                    '1 crossing matches the Fuel Lifetime above; the '
                    'subcritical tail is omitted.'
                )
                st.markdown('<div style="height:1rem"></div>', unsafe_allow_html=True)

    # ── Subsection 3: Electricity Costs ───────────────────────
    # Capital story on the left (OCC + TCI + OCC decomposition);
    # operational + market story on the right (Annualized Cost +
    # AC decomposition + Levelized Costs).
    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:1.25rem 0 0.4rem 0;">Electricity Costs</div>'
        '<p style="color:#64748b;font-size:0.85rem;margin:0 0 0.85rem 0;">'
        'Ranges shown are <strong>mean &minus; 1&sigma;</strong> to '
        '<strong>mean + 1&sigma;</strong>.'
        '</p>',
        unsafe_allow_html=True,
    )
    econ_left, econ_right = st.columns([1, 1], gap='large')
    with econ_left:
        st.markdown(
            '<div style="font-size:0.85rem;font-weight:600;color:#64748b;text-transform:uppercase;'
            'letter-spacing:0.09em;margin-bottom:0.6rem;">Capital Costs</div>',
            unsafe_allow_html=True,
        )
        cc1, cc2 = st.columns(2)
        _kpi_card(
            cc1, 'Overnight Capital Cost (OCC)',
            _fmt_cost(occ_f, occ_f_std), _fmt_cost(occ_n, occ_n_std),
            help_text=(
                'Overnight Capital Cost: the total cost to build the plant '
                'in reference-year dollars as if construction could happen '
                'instantly. Sum of pre-construction, direct, indirect, and '
                'training costs (Code of Accounts 10-40). Excludes financing '
                'interest and escalation, so reactor designs can be compared '
                'on a like-for-like basis.'
            ),
        )
        _kpi_card(
            cc2, 'Total Capital Investment (TCI)',
            _fmt_cost(tci_f, tci_f_std), _fmt_cost(tci_n, tci_n_std),
            help_text=(
                'Total Capital Investment: OCC plus capitalized financial '
                'costs accrued during construction (interest on debt-financed '
                'portion, escalation). This is the real, time-value-of-money '
                'amount needed to bring the plant to commercial operation.'
            ),
        )

        # ── OCC decomposition (accounts 10/20/30/40) ──────────────
        # Layout: a single 4-column × 6-row CSS grid so every element
        # at the same grid-row (titles, FOAK badges, FOAK values,
        # separator, NOAK badges, NOAK values) lines up horizontally
        # across all 4 components no matter how tall any single cell
        # gets (e.g. when "Pre-construction" wraps).
        _pc_f, _pc_f_std = _get_mean_std(display_df, 10, 'FOAK')
        _dc_f, _dc_f_std = _get_mean_std(display_df, 20, 'FOAK')
        _ic_f, _ic_f_std = _get_mean_std(display_df, 30, 'FOAK')
        _tc_f, _tc_f_std = _get_mean_std(display_df, 40, 'FOAK')
        _pc_n, _pc_n_std = _get_mean_std(display_df, 10, 'NOAK')
        _dc_n, _dc_n_std = _get_mean_std(display_df, 20, 'NOAK')
        _ic_n, _ic_n_std = _get_mean_std(display_df, 30, 'NOAK')
        _tc_n, _tc_n_std = _get_mean_std(display_df, 40, 'NOAK')

        _cats = [
            ('Preconstruction Cost',
             _fmt_cost(_pc_f, _pc_f_std), _fmt_cost(_pc_n, _pc_n_std),
             'Site selection, characterization, licensing and permitting, '
             'owner’s costs, and early engineering studies incurred '
             'before physical construction begins (Code of Accounts 10).'),
            ('Direct Cost',
             _fmt_cost(_dc_f, _dc_f_std), _fmt_cost(_dc_n, _dc_n_std),
             'Materials, equipment, and labor for the physical plant: '
             'reactor and containment, balance-of-plant systems, '
             'structures, electrical and I&C. The largest bucket in most '
             'cost estimates (Code of Accounts 20).'),
            ('Indirect Cost',
             _fmt_cost(_ic_f, _ic_f_std), _fmt_cost(_ic_n, _ic_n_std),
             'Construction-support services that are not part of the '
             'physical plant: engineering and design services, construction '
             'management, temporary facilities, insurance and taxes during '
             'construction (Code of Accounts 30).'),
            ('Training Cost',
             _fmt_cost(_tc_f, _tc_f_std), _fmt_cost(_tc_n, _tc_n_std),
             'Initial training program for operations and maintenance '
             'staff before commercial operation: operator training, '
             'simulator costs, qualification programs, and training '
             'materials (Code of Accounts 40).'),
        ]
        _title_row = ''.join(
            f'<div style="font-size:0.85rem;font-weight:600;color:#64748b;'
            f'line-height:1.25;">{c[0]}</div>'
            for c in _cats
        )
        # Help icons live on their own grid row so they align across
        # columns regardless of whether the title above wrapped.
        _help_row = ''.join(
            f'<div>{_help_icon(c[3])}</div>'
            for c in _cats
        )
        _foak_badge_row = ''.join(
            '<div><span style="background:#fff3ed;color:#c84b1e;font-size:0.85rem;'
            'font-weight:700;padding:0.12rem 0.4rem;border-radius:4px;'
            'letter-spacing:0.05em;">FOAK</span></div>'
            for _ in _cats
        )
        _foak_val_row = ''.join(
            f'<div style="font-size:0.85rem;font-weight:600;color:#0a2540;'
            f'overflow-wrap:break-word;">{c[1]}</div>'
            for c in _cats
        )
        _sep_row = ''.join(
            '<div style="height:1px;background:#f1f3f5;"></div>'
            for _ in _cats
        )
        _noak_badge_row = ''.join(
            '<div><span style="background:#eff6ff;color:#1B4F8C;font-size:0.85rem;'
            'font-weight:700;padding:0.12rem 0.4rem;border-radius:4px;'
            'letter-spacing:0.05em;">NOAK</span></div>'
            for _ in _cats
        )
        _noak_val_row = ''.join(
            f'<div style="font-size:0.85rem;font-weight:600;color:#0a2540;'
            f'overflow-wrap:break-word;">{c[2]}</div>'
            for c in _cats
        )
        _decomp_html = (
            f'<div style="background:white;border-radius:8px;padding:1.1rem 1.25rem;'
            f'border:1px solid #bfdbfe;margin-top:0.75rem;">'
            f'<div style="font-size:0.85rem;font-weight:600;color:#64748b;'
            f'text-transform:uppercase;letter-spacing:0.09em;margin-bottom:1rem;">'
            f'Overnight Capital Cost Decomposition</div>'
            f'<div style="display:grid;grid-template-columns:repeat(4, minmax(0, 1fr));'
            f'column-gap:1.25rem;row-gap:0.4rem;align-items:start;">'
            f'{_title_row}'
            f'{_help_row}'
            f'{_foak_badge_row}'
            f'{_foak_val_row}'
            f'{_sep_row}'
            f'{_noak_badge_row}'
            f'{_noak_val_row}'
            f'</div>'
            f'</div>'
        )
        st.markdown(_decomp_html, unsafe_allow_html=True)

    with econ_right:
        # ── Annualized costs (total + O&M / Fuel decomposition) ────
        st.markdown(
            '<div style="font-size:0.85rem;font-weight:600;color:#64748b;text-transform:uppercase;'
            'letter-spacing:0.09em;margin-bottom:0.6rem;">Annualized Costs</div>',
            unsafe_allow_html=True,
        )

        _ac_f, _ac_f_std = _get_mean_std(display_df, 'AC', 'FOAK')
        _ac_n, _ac_n_std = _get_mean_std(display_df, 'AC', 'NOAK')
        _om_f, _om_f_std = _get_mean_std(display_df, 70, 'FOAK')
        _om_n, _om_n_std = _get_mean_std(display_df, 70, 'NOAK')
        _fuel_f, _fuel_f_std = _get_mean_std(display_df, 80, 'FOAK')
        _fuel_n, _fuel_n_std = _get_mean_std(display_df, 80, 'NOAK')

        _kpi_card(
            st, 'Annualized Cost',
            _fmt_cost(_ac_f, _ac_f_std), _fmt_cost(_ac_n, _ac_n_std),
            help_text=(
                'Total recurring yearly cost to operate the plant: O&M, fuel, '
                'and capital plant expenditures, summed and annualized. The '
                'numerator of LCOE before discounting.'
            ),
        )

        _ac_cats = [
            ('O&M Cost',
             _fmt_cost(_om_f, _om_f_std), _fmt_cost(_om_n, _om_n_std),
             'Annual operations and maintenance: staff salaries, routine '
             'maintenance, consumables, insurance, regulatory fees, and '
             'central-facility costs (Code of Accounts 70).'),
            ('Fuel Cost',
             _fmt_cost(_fuel_f, _fuel_f_std), _fmt_cost(_fuel_n, _fuel_n_std),
             'Annual fuel cycle costs: fuel purchase, fabrication, '
             'enrichment, and back-end disposal / storage charges '
             '(Code of Accounts 80).'),
        ]
        _ac_title_row = ''.join(
            f'<div style="font-size:0.85rem;font-weight:600;color:#64748b;'
            f'line-height:1.25;">{c[0]}</div>'
            for c in _ac_cats
        )
        _ac_help_row = ''.join(
            f'<div>{_help_icon(c[3])}</div>'
            for c in _ac_cats
        )
        _ac_foak_badge_row = ''.join(
            '<div><span style="background:#fff3ed;color:#c84b1e;font-size:0.85rem;'
            'font-weight:700;padding:0.12rem 0.4rem;border-radius:4px;'
            'letter-spacing:0.05em;">FOAK</span></div>'
            for _ in _ac_cats
        )
        _ac_foak_val_row = ''.join(
            f'<div style="font-size:0.85rem;font-weight:600;color:#0a2540;'
            f'overflow-wrap:break-word;">{c[1]}</div>'
            for c in _ac_cats
        )
        _ac_sep_row = ''.join(
            '<div style="height:1px;background:#f1f3f5;"></div>'
            for _ in _ac_cats
        )
        _ac_noak_badge_row = ''.join(
            '<div><span style="background:#eff6ff;color:#1B4F8C;font-size:0.85rem;'
            'font-weight:700;padding:0.12rem 0.4rem;border-radius:4px;'
            'letter-spacing:0.05em;">NOAK</span></div>'
            for _ in _ac_cats
        )
        _ac_noak_val_row = ''.join(
            f'<div style="font-size:0.85rem;font-weight:600;color:#0a2540;'
            f'overflow-wrap:break-word;">{c[2]}</div>'
            for c in _ac_cats
        )
        _ac_decomp_html = (
            f'<div style="background:white;border-radius:8px;padding:1.1rem 1.25rem;'
            f'border:1px solid #bfdbfe;margin-top:0.75rem;">'
            f'<div style="font-size:0.85rem;font-weight:600;color:#64748b;'
            f'text-transform:uppercase;letter-spacing:0.09em;margin-bottom:1rem;">'
            f'Annualized Cost Decomposition</div>'
            f'<div style="display:grid;grid-template-columns:repeat(2, minmax(0, 1fr));'
            f'column-gap:1.25rem;row-gap:0.4rem;align-items:start;">'
            f'{_ac_title_row}'
            f'{_ac_help_row}'
            f'{_ac_foak_badge_row}'
            f'{_ac_foak_val_row}'
            f'{_ac_sep_row}'
            f'{_ac_noak_badge_row}'
            f'{_ac_noak_val_row}'
            f'</div>'
            f'</div>'
        )
        st.markdown(_ac_decomp_html, unsafe_allow_html=True)

        st.markdown('<div style="height:0.75rem"></div>', unsafe_allow_html=True)

        st.markdown(
            '<div style="font-size:0.85rem;font-weight:600;color:#64748b;text-transform:uppercase;'
            'letter-spacing:0.09em;margin-bottom:0.6rem;">Levelized Costs</div>',
            unsafe_allow_html=True,
        )
        # LCOE + LCOF only; LCOH moved to the "Heat Costs" subsection
        # below so the electricity / heat tracks read cleanly apart.
        lc1, lc2 = st.columns(2)
        _kpi_card(lc1, 'LCOE ($/MW<sub>e</sub>h)',
                  _fmt_lcoe(lcoe_f, lcoe_f_std), _fmt_lcoe(lcoe_n, lcoe_n_std))
        _kpi_card(lc2, 'LCOF Fuel Only ($/MW<sub>e</sub>h)',
                  _fmt_lcoe(lcof_f, lcof_f_std), _fmt_lcoe(lcof_n, lcof_n_std))

        st.markdown('<div style="height:0.5rem"></div>', unsafe_allow_html=True)
        with st.popover('What do FOAK and NOAK mean?', width='content'):
            st.markdown(
                '**FOAK** — *First-of-a-Kind*. Includes learning-curve and '
                'contingency premiums for an initial deployment.\n\n'
                '**NOAK** — *Nth-of-a-Kind*. Reflects mature, serial-production '
                'cost reductions at scale.'
            )

        # SECTION 3 fuel-cycle header was here; rendering moved
        # to Band 3 (Fuel cycle column). Computation continues
        # below so the _fuel_card queue is still populated.

        # ── Fuel Inventory (uses st.metric so we get the built-in help icon) ──
        _u235_g = float(params.get('Mass U235', 0.0))
        _u238_g = float(params.get('Mass U238', 0.0))
        _hm_g = _u235_g + _u238_g
        _mwe = float(params.get('Power MWe', 0.0)) or float('nan')
        _hm_kg = _hm_g / 1.0e3 # g → kg
        _hm_kg_per_mwe = _hm_kg / _mwe # kg / MWe
        _fis_kg = _u235_g / 1.0e3 # g → kg
        _fis_kg_per_mwe = _fis_kg / _mwe # kg / MWe

        # _section_cards is initialized above the tab declarations
        # so it's visible from both tab_design and tab_physics.
        def _fuel_card(title, value_str, help_text,
                       accent=None, bg=None, border=None,
                       status=None, section=3):
            # QUEUE-MODE: appends to _section_cards[section] instead
            # of rendering immediately. The actual render happens
            # later inside the appropriate tab. Status/style choice
            # is computed at queue time so callers don't need to
            # think about it.
            if status == 'warning':
                _bg, _br, _tc = '#fffbeb', '#fcd34d', '#92400e'
            elif status == 'error':
                _bg, _br, _tc = '#fef2f2', '#fecaca', '#b91c1c'
            else:
                _bg, _br, _tc = '#ffffff', '#bfdbfe', '#1B4F8C'
            _help_html = html.escape(help_text) if help_text else ''
            _html = (
                f'<div style="background:{_bg};border:1px solid {_br};'
                f'border-radius:8px;padding:1.1rem 1.25rem;'
                f'margin-bottom:0.85rem;">'
                f'<div style="font-size:0.85rem;font-weight:600;color:{_tc};'
                f'text-transform:uppercase;letter-spacing:0.09em;margin-bottom:0.35rem;">'
                f'{title}</div>'
                f'<div style="font-size:1rem;font-weight:700;color:#0a2540;'
                f'line-height:1.3;margin-bottom:0.45rem;">{value_str}</div>'
                f'<div style="font-size:0.85rem;font-weight:400;color:#64748b;'
                f'line-height:1.45;">{_help_html}</div>'
                f'</div>'
            )
            _section_cards.setdefault(section, []).append(_html)

        def _scoping_callout(text):
            st.markdown(
                f'<div style="background:#fffbeb;border:1.5px solid #f59e0b;'
                f'border-radius:8px;padding:0.85rem 1.1rem;margin-bottom:1rem;'
                f'font-size:0.85rem;line-height:1.55;color:#92400e;">'
                f'<strong style="color:#92400e;">Scoping-only note. </strong>'
                f'{text}</div>',
                unsafe_allow_html=True,
            )

        def _render_section_cards(section_num):
            for _h in _section_cards.get(section_num, []):
                st.markdown(_h, unsafe_allow_html=True)

        _fuel_card(
            'Fuel loading',
            f'{_fmt_metric(_hm_kg)} kgHM | {_fmt_metric(_hm_kg_per_mwe)} kgHM/MWe',
            ('Total uranium mass in the core (HM = Heavy Metal = '
             'U235 + U238). kgHM/MWe normalises by net electric '
             'output, useful for comparing fuel inventory across '
             'reactor types. Microreactors typically run 100 to '
             '1,000 kgHM/MWe.'),
            section=3,
        )
        _fuel_card(
            'Fissile loading',
            f'{_fmt_metric(_fis_kg)} kg | {_fmt_metric(_fis_kg_per_mwe)} kg/MWe',
            ('Mass of fissile material (U235 for LEU/HALEU fuel). '
             'Drives enrichment cost (HALEU runs $3k to $15k/kg) and '
             'safeguards requirements. kg/MWe is much higher for '
             'microreactors than for commercial LWRs (~5 kg/MWe).'),
            section=3,
        )

        # Peaking factor + discharge burnup only meaningful for
        # critical cases. Skip if HM mass or lifetime aren't available.
        _lifetime_days = float(params.get('Fuel Lifetime', 0.0))
        if _hm_kg > 0 and _lifetime_days > 0:
            _pf = 0.0
            if reactor_type == 'LTMR':
                _pf = get_ltmr_peaking_factor(
                    n_rings_per_assembly = params['Number of Rings per Assembly'],
                    active_height = params['Active Height'],
                    enrichment = params['Enrichment'],
                    power_mwt = params['Power MWt'],
                )
            elif reactor_type == 'GCMR':
                _pf = get_gcmr_peaking_factor(
                    assembly_rings = params['Assembly Rings'],
                    core_rings = params['Core Rings'],
                    active_height = params['Active Height'],
                    enrichment = params['Enrichment'],
                    power_mwt = params['Power MWt'],
                )
            elif reactor_type == 'HPMR':
                _pf = get_hpmr_peaking_factor(
                    n_rings_per_assembly = params['Number of Rings per Assembly'],
                    n_rings_per_core = params['Number of Rings per Core'],
                    active_height = params['Active Height'],
                    enrichment = params['Enrichment'],
                    power_mwt = params['Power MWt'],
                )

            # Average discharge burnup: total energy / total HM mass.
            # MWt × days = MW·d, divided by kg → MWd/kg.
            _bu_avg = (float(params['Power MWt']) * _lifetime_days) / _hm_kg
            _bu_max = _bu_avg * _pf if _pf > 0 else 0.0

            _fuel_card(
                'Peaking factor',
                f'{_pf:.1f}' if _pf > 0 else 'N/A',
                ('Ratio of the maximum local fission rate to the core '
                 'average. Drives the gap between average and peak fuel '
                 'pin temperatures and burnup. Lower is better; typical '
                 'microreactor values are 1.5 to 3.0. Interpolated from '
                 'the parametric study.'),
                section=4,
            )

            # Axial + total leakage (BOL, %). Inside the trained
            # H range we use the KNN-interpolated value; outside we
            # fall back to a one-group migration-area physics
            # formula (which actually responds to the user's H).
            _ax_lk, _tot_lk, _lk_src = 0.0, 0.0, None
            _active_radius_cm = float(2.0 * params['Core Radius']) / 2.0
            # active_radius excludes the reflector recompute from active diameter
            _active_radius_cm = (float(2.0 * params['Core Radius'])
                                 - 2.0 * float(params.get('Radial Reflector Thickness', 0.0))
                                 ) / 2.0
            _r_refl = float(params.get('Radial Reflector Thickness', 0.0))
            _z_refl = float(params.get('Axial Reflector Thickness', 0.0))
            if reactor_type == 'LTMR':
                _ax_lk, _tot_lk, _lk_src = get_ltmr_leakage(
                    n_rings_per_assembly = params['Number of Rings per Assembly'],
                    active_height = params['Active Height'],
                    enrichment = params['Enrichment'],
                    power_mwt = params['Power MWt'],
                    active_radius_cm = _active_radius_cm,
                    radial_reflector_cm = _r_refl,
                    axial_reflector_cm = _z_refl,
                )
            elif reactor_type == 'GCMR':
                _ax_lk, _tot_lk, _lk_src = get_gcmr_leakage(
                    assembly_rings = params['Assembly Rings'],
                    core_rings = params['Core Rings'],
                    active_height = params['Active Height'],
                    enrichment = params['Enrichment'],
                    power_mwt = params['Power MWt'],
                    active_radius_cm = _active_radius_cm,
                    radial_reflector_cm = _r_refl,
                    axial_reflector_cm = _z_refl,
                )
            elif reactor_type == 'HPMR':
                _ax_lk, _tot_lk, _lk_src = get_hpmr_leakage(
                    n_rings_per_assembly = params['Number of Rings per Assembly'],
                    n_rings_per_core = params['Number of Rings per Core'],
                    active_height = params['Active Height'],
                    enrichment = params['Enrichment'],
                    power_mwt = params['Power MWt'],
                    active_radius_cm = _active_radius_cm,
                    radial_reflector_cm = _r_refl,
                    axial_reflector_cm = _z_refl,
                )

            _src_note = (
                ' Source: interpolated from the parametric study '
                '(inside the trained design space).'
                if _lk_src == 'interpolated' else
                ' Source: computed from a one group physics formula '
                '(outside the trained range).'
            )

            if _ax_lk > 0:
                _fuel_card(
                    'Axial leakage (BOL)',
                    f'{_fmt_metric(_ax_lk)} %',
                    ('Fraction of neutrons that escape through the top '
                     'or bottom of the active core at beginning of life. '
                     'Driven by Active Height (shorter cores leak more) '
                     'and axial reflector thickness. Reflector and '
                     'control drums are accounted for.'
                     + _src_note),
                    section=4,
                )
            if _tot_lk > 0:
                _fuel_card(
                    'Total leakage (BOL)',
                    f'{_fmt_metric(_tot_lk)} %',
                    ('Total fraction of neutrons that escape the active '
                     'core (axial + radial) at beginning of life. Driven '
                     'by core dimensions (Active Height and active '
                     'radius). Microreactors typically have higher total '
                     'leakage (5 to 35 percent) than commercial LWRs '
                     '(~3 percent) due to their small size.'
                     + _src_note),
                    section=4,
                )
            _fuel_card(
                'Discharge burnup (avg)',
                f'{_fmt_metric(_bu_avg)} MWd/kgHM',
                ('Average burnup of the fuel at end of life (when the '
                 'reactor first becomes subcritical). Computed as '
                 'Power [MWt] × Fuel Lifetime [days] / Heavy Metal '
                 'mass [kg]. Headline economic metric: higher means '
                 'more energy per kg of fuel, lowering fuel cost per '
                 'MWh.'),
                section=3,
            )
            if _bu_max > 0:
                _fuel_card(
                    'Discharge burnup (max)',
                    f'{_fmt_metric(_bu_max)} MWd/kgHM',
                    ('Peak burnup of the most heavily depleted region '
                     '(= average × peaking factor). The design limiting '
                     'value: cladding integrity, fission gas release, '
                     'and dimensional change all depend on the peak. '
                     'Commercial LWRs are licensed to ~62 MWd/kgU peak; '
                     'TRISO fueled designs allow much higher.'),
                    section=3,
                )

            # Mining intensity uses MOUSE's existing 'Natural Uranium Mass'
            # (computed in fuel_calculations using T=0.25% tails and natural
            # U feed F=0.71%, the standard mass-balance formula).
            _nat_u_kg = float(params.get('Natural Uranium Mass', 0.0))
            _mwh_total = _mwe * _lifetime_days * 24.0
            if _nat_u_kg > 0 and _mwh_total > 0:
                _mining = (_nat_u_kg * 1000.0) / _mwh_total # kg→g, /MWh
                _fuel_card(
                    'Mining intensity',
                    f'{_fmt_metric(_mining)} gU/MWh',
                    ('Mass of natural uranium mined and milled per MWh '
                     'of electricity produced. Computed using tails '
                     'enrichment 0.25 percent and feed 0.71 percent. '
                     'Typical: commercial LWRs ~17 to 25 gU/MWh, HALEU '
                     'microreactors ~30 to 80, natural U reactors '
                     '(CANDU) ~150 to 200. Lower is better.'),
                    section=3,
                )

            if _reactivity_swing_pct is not None:
                _fuel_card(
                    'Reactivity swing',
                    f'{_fmt_metric(_reactivity_swing_pct)} %Δk/k',
                    ('Total reactivity consumed by burnup over the fuel '
                     'cycle, in %Δk/k. Drives control drum sizing: drum '
                     'worth must exceed the swing to keep cold clean '
                     'k_eff below 1 with all drums in. Typical: '
                     'commercial LWRs ~10 to 15 percent, HALEU '
                     'microreactors ~15 to 40 percent.'),
                    section=4,
                )

            # Heat flux at the fuel-pin surface (MW/m²) already in
            # params from calculate_heat_flux: Power / total pin
            # cylindrical surface area.
            _hflux = float(params.get('Heat Flux', 0.0))
            if _hflux > 0:
                _fuel_card(
                    'Heat flux (avg)',
                    f'{_fmt_metric(_hflux * 100.0)} W/cm² ({_fmt_metric(_hflux)} MW/m²)',
                    ('Average heat flux at the outer surface of the '
                     'fuel pins = Power / (π × pin diameter × H × pin '
                     'count). Sets the convective heat removal '
                     'requirement on the coolant. Typical microreactor '
                     'values: 0.1 to 1 MW/m² (10 to 100 W/cm²).'),
                    section=5,
                )

            # Separative Work Units (SWU): kg-SWU total and per MWh.
            _swu = float(params.get('SWU', 0.0))
            if _swu > 0 and _mwh_total > 0:
                _swu_per_mwh = _swu * 1000.0 / _mwh_total # kg → g, /MWh
                _fuel_card(
                    'Enrichment SWU',
                    f'{_fmt_metric(_swu)} kg-SWU | {_fmt_metric(_swu_per_mwh)} g-SWU/MWh',
                    ('Separative Work Units consumed to produce the '
                     'fuel for this core, the standard metric for '
                     'enrichment effort. Computed using the value '
                     'function method (tails 0.25 percent, feed 0.71 '
                     'percent). Drives enrichment cost ($50 to $200/kg '
                     'SWU); higher enrichment products need more SWU '
                     'per kg.'),
                    section=3,
                )

            # Onsite coolant inventory power-scaled, reactor-specific.
            # LTMR: 1833 kg/MWₜ of NaK (Creys-Malville scaling).
            # GCMR: 3.3 kg/MWₜ of helium (UNT 919556 tables 17 & 18).
            # HPMR: 0 (heat pipes individually sealed; no bulk inventory).
            #
            # Display in tons with deliberately coarse rounding to signal
            # this is a rough estimate, not a precision figure:
            # >= 1 ton -> integer tons
            # < 1 ton -> 1 significant figure
            if reactor_type in ('LTMR', 'GCMR'):
                _inv_kg = float(params.get('Onsite Coolant Inventory', 0.0))
                if _inv_kg > 0:
                    _coolant_label = {'LTMR': 'NaK', 'GCMR': 'Helium'}[reactor_type]
                    _tons = _inv_kg / 1000.0
                    if _tons >= 1:
                        _val_str = f'{round(_tons):,d} ton'
                    else:
                        _val_str = f'{_tons:.1g} ton'
                    _inv_help = {
                        'LTMR': (
                            'Rough estimate of the on site primary NaK '
                            'inventory (filled core plus storage). '
                            'Scales linearly with thermal power at '
                            '~1833 kg/MWt, derived from the Creys '
                            'Malville sodium plant. NaK is sealed and '
                            'not consumed, so no periodic replacement. '
                            'Drives the coolant procurement line in OCC.'
                        ),
                        'GCMR': (
                            'Rough estimate of the on site helium '
                            'inventory. Scales linearly with thermal '
                            'power at ~3.3 kg/MWt. He has a steady ~10 '
                            'percent/year leakage rate, so one tenth '
                            'is replaced annually. Drives both the OCC '
                            'coolant line and the OPEX make up term.'
                        ),
                    }[reactor_type]
                    _fuel_card(
                        'Coolant inventory',
                        f'{_val_str} ({_coolant_label})',
                        _inv_help,
                        section=5,
                    )

            # Coolant mass flow rate + primary-loop power fraction.
            # LTMR: NaK pump. GCMR: He compressor. HPMR: heat
            # pipes (no flow, no card).
            if reactor_type == 'LTMR':
                _mdot = float(params.get('Coolant Mass Flow Rate', 0.0))
                if _mdot > 0:
                    _fuel_card(
                        'Coolant mass flow rate',
                        f'{_fmt_metric(_mdot)} kg/s',
                        ('Primary coolant mass flow rate from '
                         'm_dot = Power_MWt × 10⁶ / (ΔT × c_p) with '
                         'the LTMR\'s fixed ΔT = 90 °C (Tin 430, Tout '
                         '520 °C) and the heat capacity of NaK. Sets '
                         'the primary pump and heat exchanger sizing. '
                         'Typical: 5 to 100 kg/s for 1 to 20 MWt.'),
                        section=5,
                    )

                _pump_kW = float(params.get('Primary Pump Mechanical Power', 0.0))
                if _pump_kW > 0 and _mwe > 0:
                    _pump_pct = 100.0 * _pump_kW / (_mwe * 1000.0)
                    _fuel_card(
                        'Primary pump fraction',
                        f'{_fmt_metric(_pump_pct)} % ({_fmt_metric(_pump_kW)} kW)',
                        ('Primary loop pump mechanical power as a '
                         'fraction of gross electric output. Computed '
                         'from P = m_dot × Δp / (ρ × η) with Δp = 250 '
                         'kPa, ρ = 750 kg/m³ (NaK at ~500 °C), η = '
                         '0.75. Typical liquid metal values are 1 to '
                         '3 percent.'),
                        section=5,
                    )

            elif reactor_type == 'GCMR':
                # Coolant Mass Flow Rate is set by mass_flow_rate()
                # in tools.py and represents the TOTAL helium flow
                # across all primary loops (loop_factor adjustment
                # already applied). ΔT is 250 °C for the GCMR
                # configuration (300 → 550 °C).
                _mdot_gcmr = float(params.get('Coolant Mass Flow Rate', 0.0))
                if _mdot_gcmr > 0:
                    _fuel_card(
                        'Coolant mass flow rate',
                        f'{_fmt_metric(_mdot_gcmr)} kg/s',
                        ('Primary helium mass flow rate, total across '
                         'all primary loops. Computed from m_dot = '
                         'Power_MWt × 10⁶ / (ΔT × c_p) per loop with '
                         'GCMR\'s ΔT = 250 °C (Tin 300, Tout 550 °C) '
                         'and the c_p of He. Helium\'s low density '
                         'keeps these small: 1 to 10 kg/s for 1 to 20 '
                         'MWt.'),
                        section=5,
                    )

                # Primary Loop Compressor Power is the per-loop
                # power in W (compressor_power in tools.py uses
                # Primary Loop Mass Flow Rate, which is per-loop).
                # Multiply by Primary Loop Count (2 for GCMR) to
                # get the total compressor power, then express as
                # a fraction of gross electric output.
                _comp_w_per_loop = float(params.get('Primary Loop Compressor Power', 0.0))
                _n_loops = float(params.get('Primary Loop Count', 1))
                _comp_kW_total = (_comp_w_per_loop * _n_loops) / 1000.0
                if _comp_kW_total > 0 and _mwe > 0:
                    _comp_pct = 100.0 * _comp_kW_total / (_mwe * 1000.0)
                    _fuel_card(
                        'Primary compressor fraction',
                        f'{_fmt_metric(_comp_pct)} % ({_fmt_metric(_comp_kW_total)} kW total)',
                        ('Primary loop helium compressor power as a '
                         'fraction of gross electric output. Computed '
                         'from Δp · m_dot / (η · ρ_He) with ρ_He = '
                         '3.33 kg/m³, Δp = 50 kPa, η = 0.8, times '
                         'Primary Loop Count (2 for GCMR). Typical: '
                         '2 to 6 percent (higher than liquid metals '
                         'due to He\'s low density).'),
                        section=5,
                    )

            elif reactor_type == 'HPMR':
                # HPMR uses sealed Na heat pipes (no flowing primary
                # loop), so the relevant design metrics are per-pipe
                # thermal duty and linear heat rate rather than mass
                # flow / pump power.
                _n_hp = int(params.get('Number of Heatpipes', 0))
                if _n_hp > 0:
                    _fuel_card(
                        'Number of heat pipes',
                        f'{_n_hp:,}',
                        ('Total Na heat pipe count across the core = '
                         '(heat pipes per assembly) × (fuel assemblies '
                         'per core). Typical HPMR designs have 500 to '
                         '5,000 heat pipes (eVinci class ~1,000 to '
                         '2,000). More pipes lower per pipe duty but '
                         'raise fabrication cost.'),
                        section=5,
                    )

                    _kw_per_hp = float(params['Power MWt']) * 1000.0 / _n_hp
                    _fuel_card(
                        'Power per heat pipe',
                        f'{_fmt_metric(_kw_per_hp)} kW/HP',
                        ('Thermal duty per Na heat pipe = Power_MWt × '
                         '1000 / N_HP. Capacity is bounded by four '
                         'physical limits (capillary, sonic, '
                         'entrainment, boiling). At 700 to 900 °C the '
                         'practical envelope is 5 to 15 kW/HP for 1.5 '
                         'to 3 cm pipes; above ~20 kW/HP needs exotic '
                         'wick geometries.'),
                        section=5,
                    )

                    _h_cm = float(params.get('Active Height', 0.0))
                    if _h_cm > 0:
                        _kw_per_cm = _kw_per_hp / _h_cm
                        _fuel_card(
                            'Linear heat rate',
                            f'{_fmt_metric(_kw_per_cm)} kW/cm',
                            ('Heat pipe linear heat rate = (kW per HP) '
                             '/ Active Height. THE key sizing metric '
                             'for sodium heat pipes. Published Na HP '
                             'test data (LANL, INL) shows steady state '
                             'up to 0.20 to 0.30 kW/cm for 1.5 to 2.5 '
                             'cm pipes at 700 to 900 °C. Typical '
                             'designs target 0.05 to 0.20 kW/cm for '
                             'margin.'),
                            section=5,
                        )

                    _n_pins = int(params.get('Fuel Pin Count', 0))
                    if _n_pins > 0:
                        _pin_hp_ratio = _n_pins / _n_hp
                        _fuel_card(
                            'Pin to heat pipe ratio',
                            f'{_fmt_metric(_pin_hp_ratio)} (pins / HP)',
                            ('Fuel pins served per heat pipe = Fuel '
                             'Pin Count / N_HP. HPMR designs typically '
                             'run 1 to 6 pins per HP depending on '
                             'lattice. Above ~6, each HP must drain '
                             'heat from many pins, raising local '
                             'cladding to HP flux and forcing larger '
                             'HP diameter or smaller pin pitch.'),
                            section=5,
                        )

            # Power density = Power_MWₜ / Active Core Volume.
            # Active core is a hex prism with apothem = active_radius
            # and height = Active Height. _active_radius_cm was already
            # computed earlier in this block (Core Radius − reflector).
            _vol_cm3 = (2.0 * _math.sqrt(3.0)
                        * (_active_radius_cm ** 2)
                        * float(params['Active Height']))
            _vol_m3 = _vol_cm3 * 1.0e-6
            if _vol_m3 > 0:
                _pd = float(params['Power MWt']) / _vol_m3
                _fuel_card(
                    'Power density',
                    f'{_fmt_metric(_pd)} MW/m³',
                    ('Thermal power per unit active core volume = '
                     'Power_MWt / (2√3 · R² · H), where R is the hex '
                     'apothem (no reflector). Tells you how '
                     'aggressively the fuel is loaded. Typical: '
                     'eVinci/HPMR class ~10 to 15 MW/m³, LTMR ~10 to '
                     '30, metal fuel microreactors ~50 to 100, '
                     'commercial LWRs ~100.'),
                    section=5,
                )

            # Aspect ratio moved out of the TH section into the
            # "What does it look like" panel in img_col (Band 1).

        # Section 3 cards now rendered in Band 3 (Fuel cycle
        # column) instead of inline here.

        st.markdown('<div style="height:0.75rem"></div>', unsafe_allow_html=True)

    # ── Subsection 4: Heat Costs ────────────────────────────────
    # Heat-application costs derived from the electricity numbers
    # via factors from INL/RPT-23/72866 §8.3 Cost-Adjustment
    # Methodology (data-set averages across 30+ designs). The same
    # constants are used in cost/non_direct_cost.py to compute LCOH,
    # so the values shown here are consistent with the LCOE→LCOH
    # pipeline. Decomposition is not shown because the factors are
    # applied as totals (Capital and O&M+fuel as aggregates), not
    # per-account.
    _HEAT_OCC_FACTOR = 0.795   # OCC_heat = OCC × 0.795
    _HEAT_ANN_FACTOR = 0.966   # AC_heat = AC × 0.966
    _ho_f,  _ho_f_std  = occ_f * _HEAT_OCC_FACTOR, occ_f_std * _HEAT_OCC_FACTOR
    _ho_n,  _ho_n_std  = occ_n * _HEAT_OCC_FACTOR, occ_n_std * _HEAT_OCC_FACTOR
    _hac_f, _hac_f_std = _ac_f * _HEAT_ANN_FACTOR, _ac_f_std * _HEAT_ANN_FACTOR
    _hac_n, _hac_n_std = _ac_n * _HEAT_ANN_FACTOR, _ac_n_std * _HEAT_ANN_FACTOR

    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:1.25rem 0 0.4rem 0;">Heat Costs</div>'
        '<p style="color:#64748b;font-size:0.85rem;margin:0 0 0.85rem 0;">'
        'Ranges shown are <strong>mean &minus; 1&sigma;</strong> to '
        '<strong>mean + 1&sigma;</strong>.'
        '</p>',
        unsafe_allow_html=True,
    )
    # Same LCOH value, rendered in two unit systems: $/MWhth and
    # $/MMBtu. Conversion: 1 MWh thermal = 3.4121 MMBtu (exact).
    _MMBTU_PER_MWHT = 3.4121

    def _fmt_lcoh_mmbtu(mean, std):
        if math.isnan(mean):
            return 'N/A'
        m = mean / _MMBTU_PER_MWHT
        if math.isnan(std) or std == 0:
            return f'${m:.1f}/MMBtu'
        lo = (mean - std) / _MMBTU_PER_MWHT
        hi = (mean + std) / _MMBTU_PER_MWHT
        return (
            f'${m:.1f}/MMBtu '
            f'<span style="color:#64748b;font-weight:400;">'
            f'[${lo:.1f} &ndash; ${hi:.1f}]'
            f'</span>'
        )

    hc1, hc2, hc3, hc4 = st.columns(4)
    _kpi_card(
        hc1, 'Capital Cost',
        _fmt_cost(_ho_f, _ho_f_std), _fmt_cost(_ho_n, _ho_n_std),
        help_text=(
            'Capital cost for heat-only operation: electricity OCC × 0.795. '
            'The 0.795 factor is the data-set average from INL/RPT-23/72866 '
            '§8.3, accounting for removing the power-conversion system and '
            'the related reduction in direct and indirect costs.'
        ),
    )
    _kpi_card(
        hc2, 'Annualized Cost',
        _fmt_cost(_hac_f, _hac_f_std), _fmt_cost(_hac_n, _hac_n_std),
        help_text=(
            'Annual O&M + fuel for heat-only operation: electricity '
            'Annualized Cost × 0.966. The 0.966 factor is the data-set '
            'average from INL/RPT-23/72866 §8.3.'
        ),
    )
    _kpi_card(
        hc3, 'LCOH ($/MW<sub>t</sub>h)',
        _fmt_lcoh(lcoh_f, lcoh_f_std), _fmt_lcoh(lcoh_n, lcoh_n_std),
        help_text=(
            'Levelized Cost of Heat in $/MWh thermal. Computed via the '
            'full PV formula using OCC_heat and AC_heat, then multiplied '
            'by the reactor\'s thermal efficiency to convert from '
            '$/MWh-electric to $/MWh-thermal.'
        ),
    )
    _kpi_card(
        hc4, 'LCOH ($/MMBtu)',
        _fmt_lcoh_mmbtu(lcoh_f, lcoh_f_std), _fmt_lcoh_mmbtu(lcoh_n, lcoh_n_std),
        help_text=(
            'Same LCOH expressed in $ per million BTU of thermal energy. '
            'Conversion: 1 MWh-thermal = 3.4121 MMBtu, so $/MMBtu = '
            '($/MWh-thermal) ÷ 3.4121. $/MMBtu is the conventional unit '
            'for industrial process-heat purchasing.'
        ),
    )

    # ═══════════════════════════════════════════════════════════
    # BAND 2 Can we ship it?
    # ═══════════════════════════════════════════════════════════
    _section_header(
        '02',
        'Transportability Check',
        'Can we ship it?',
    )

    # ─────────────────────────────────────────────────────────────
    # Transportability Considerations
    # ─────────────────────────────────────────────────────────────
    # Per-component dimensions (height, diameter) and dry mass for
    # the four nested envelopes that make up a microreactor module.
    # Per-mode badges compare the outermost (RVACS) envelope to
    # three transport-mode dimensional limits. No total mass is
    # computed (per design choice each component is shown alone).
    #
    # Rounding rules (match the coolant-inventory card):
    # weights ≥ 1 ton -> integer tons
    # weights < 1 ton -> 1 significant figure
    # dimensions -> meters with 1 decimal
    # Inner "Can we ship it?" sub-header removed; the band-level
    # header above already announces the section.
    st.markdown(
        '<div style="background:#f7f8fa;border:1px solid #bfdbfe;border-radius:8px;'
        'padding:0.85rem 1.1rem;margin-bottom:0.9rem;font-size:0.85rem;line-height:1.45;color:#3c4257;">'
        'Transportability is one of the headline features of microreactors '
        'they can move by truck, rail, or sea container, which differentiates '
        'them from large NPPs. The numbers below report the component geometry '
        'and then screen the assembled reactor module against standard road, '
        'rail, and sea shipping limits.'
        '</div>',
        unsafe_allow_html=True,
    )

    # ── Module Geometry subsection header ──────────────────────
    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:0 0 0.4rem 0;">Module Geometry</div>'
        '<p style="color:#64748b;font-size:0.85rem;margin:0 0 0.85rem 0;">'
        'The reactor is built up from up to four nested layers, listed '
        'below from innermost to outermost. Their masses are shown separately '
        'for transparency, but the transportability check treats them as one '
        'assembled reactor module.'
        '</p>',
        unsafe_allow_html=True,
    )

    # ── Build per-component rows (height_m, diameter_m, mass_kg) ──
    def _lb_str(kg):
        """Convert modeled component mass from kilograms to whole pounds."""
        if kg is None:
            return ''
        try:
            kg = float(kg)
        except (TypeError, ValueError):
            return ''
        if kg <= 0:
            return ''
        return f'{kg * 2.2046226218:,.0f} lb'

    def _m1(cm):
        if cm is None or cm <= 0:
            return ''
        return f'{cm/100.0:.1f} m'

    # MOUSE label disambiguation:
    # - LTMR: 'Vessel' = reactor vessel, 'Guard Vessel' = guard vessel
    # - GCMR: 'Vessel' = core barrel (internal),
    # 'Guard Vessel' = the actual RPV (He pressure boundary)
    # - HPMR: 'Vessel' = reactor vessel, no guard vessel
    _vessel_height_cm = float(params.get('Vessel Height', 0.0))
    _bottom_depth_cm = float(params.get('Vessel Bottom Depth', 0.0))
    _vessel_thk_cm = float(params.get('Vessel Thickness', 0.0))
    _guard_thk_cm = float(params.get('Guard Vessel Thickness', 0.0))
    _gap_v_g_cm = float(params.get('Gap Between Vessel And Guard Vessel', 0.0))

    # Reactor (core + reflectors + drums). All component masses
    # are in kg. Use full modeled fuel mass instead of uranium
    # heavy-metal mass alone. For GCMR, graphite matrix inside fuel
    # compacts is already counted by the moderator-mass routine, so
    # subtract it from the moderator term to avoid double counting.
    _fuel_mass_kg, _fuel_matrix_overlap_kg, _heatpipe_steel_kg = _transport_fuel_masses(params, reactor_type)
    _reactor_dia_cm = 2.0 * float(params.get('Core Radius', 0.0))
    _reactor_h_cm = (float(params.get('Active Height', 0.0))
                        + 2.0 * float(params.get('Axial Reflector Thickness', 0.0)))
    _reactor_mass_kg = (
        _fuel_mass_kg
        + max(0.0, float(params.get('Moderator Mass', 0.0)) - _fuel_matrix_overlap_kg)
        + float(params.get('Moderator Booster Mass', 0.0))
        + float(params.get('Radial Reflector Mass', 0.0))
        + float(params.get('Axial Reflector Mass', 0.0))
        + float(params.get('Control Drums Mass', 0.0))
        + _heatpipe_steel_kg
    )

    # Reactor vessel (the pressure boundary)
    if reactor_type == 'GCMR':
        # MOUSE 'Guard Vessel' is the RPV for GCMR
        _rv_outer_r_cm = (float(params.get('Guard Vessel Radius', 0.0))
                          + _guard_thk_cm)
        _rv_height_cm = (_vessel_height_cm
                          + _bottom_depth_cm + _vessel_thk_cm + _gap_v_g_cm
                          + _guard_thk_cm)
        _rv_mass_kg = float(params.get('Guard Vessel Mass', 0.0))
    else:
        _rv_outer_r_cm = (float(params.get('Vessel Radius', 0.0)) + _vessel_thk_cm)
        _rv_height_cm = _vessel_height_cm + _bottom_depth_cm
        _rv_mass_kg = float(params.get('Vessel Mass', 0.0))
    _rv_dia_cm = 2.0 * _rv_outer_r_cm

    # Guard vessel only for LTMR
    _has_guard = (reactor_type == 'LTMR'
                  and float(params.get('Guard Vessel Thickness', 0.0)) > 0)
    if _has_guard:
        _gv_outer_r_cm = (float(params.get('Guard Vessel Radius', 0.0)) + _guard_thk_cm)
        _gv_dia_cm = 2.0 * _gv_outer_r_cm
        _gv_height_cm = (_vessel_height_cm
                          + _bottom_depth_cm + _vessel_thk_cm + _gap_v_g_cm
                          + _guard_thk_cm)
        _gv_mass_kg = float(params.get('Guard Vessel Mass', 0.0))
    else:
        _gv_dia_cm = _gv_height_cm = _gv_mass_kg = 0.0

    # RVACS (cooling vessel + intake vessel combined)
    _rvacs_outer_r_cm = float(params.get('Vessels Total Radius', 0.0))
    _rvacs_dia_cm = 2.0 * _rvacs_outer_r_cm
    _rvacs_height_cm = float(params.get('Vessels Total Height', 0.0))
    _rvacs_mass_kg = (float(params.get('Cooling Vessel Mass', 0.0))
                         + float(params.get('Intake Vessel Mass', 0.0)))

    # ── Render component table ──
    # Each row carries an always-visible small-font description
    # under the component name explaining what is included and
    # what is excluded. The hover-only title= tooltip approach
    # was unreliable across browsers/themes, so the notes are
    # shown inline in the table itself. Explicit colors
    # throughout so the table is readable regardless of the
    # Streamlit theme.
    _CELL = ('padding:0.55rem 0.8rem;color:#0a2540;'
             'border-bottom:1px solid #bfdbfe;'
             'vertical-align:top;')
    _CELL_C = _CELL + 'text-align:center;'
    _CELL_NAME = _CELL + 'font-weight:600;'
    _DESC = ('font-size:0.85rem;font-weight:400;color:#64748b;'
             'line-height:1.4;margin-top:0.2rem;')

    _moderator_for_type = {
        'LTMR': 'ZrH',
        'GCMR': 'graphite (with ZrH booster pins)',
        'HPMR': 'monolith graphite',
    }.get(reactor_type, 'moderator')
    _reactor_desc = (
        f'Includes: modeled fuel mass, moderator ({_moderator_for_type}), '
        'radial + axial reflector, and control drums.'
        + (' Fuel mass includes UZrH fuel meat plus modeled Zr and SS304 fuel-pin metal regions.'
           if reactor_type == 'LTMR' else '')
        + (' Fuel mass includes UCO kernels, TRISO coating layers, and graphite fuel-compact matrix.'
           if reactor_type == 'GCMR' else '')
        + (' Fuel mass includes homogenized TRISO fuel; heat-pipe steel is included, Na working fluid is excluded.'
           if reactor_type == 'HPMR' else '')
    )
    _rv_desc_extra = (
        ' <span style="color:#1B4F8C;">For GCMR this maps to '
        'MOUSE\'s internal "Guard Vessel" field (the RPV).</span>'
        if reactor_type == 'GCMR' else ''
    )
    _rv_desc = (
        'Diameter = 2 × (vessel radius + thickness). '
        'Height = active core + axial reflector + lower '
        'plenum + upper plenum + bottom dish. '
        'Top closure dome not modeled. Mass = vessel wall only.'
        + _rv_desc_extra
    )
    _gv_desc = (
        'Secondary containment shell around the reactor vessel '
        'for primary coolant leak containment. Mass = '
        'guard vessel wall only (no internals).'
    )
    _gv_na_desc = (
        'Intentionally omitted for this reactor type He is '
        'inert (GCMR) and heat pipes are individually sealed '
        '(HPMR), so neither has a bulk primary coolant '
        'requiring secondary containment.'
    )
    _rvacs_desc = (
        'The two outer vessels that comprise the Reactor Vessel '
        'Auxiliary Cooling System: the cooling vessel + the intake '
        'vessel, treated here as one shipping envelope. Diameter = '
        '2 × intake vessel outer radius; height is the full external '
        'envelope. Mass = cooling vessel wall + intake vessel wall '
        'only (no air, no insulation, no support structure).'
    )

    _rows_html = []
    _rows_html.append(
        '<tr style="background:#ffffff;">'
        f'<td style="{_CELL_NAME}">Reactor (core + reflectors + drums)'
        f'<div style="{_DESC}">{_reactor_desc}</div></td>'
        f'<td style="{_CELL_C}">{_m1(_reactor_h_cm)}</td>'
        f'<td style="{_CELL_C}">{_m1(_reactor_dia_cm)}</td>'
        f'<td style="{_CELL_C}">{_lb_str(_reactor_mass_kg)}</td>'
        '</tr>'
    )
    _rows_html.append(
        '<tr style="background:#f7f8fa;">'
        f'<td style="{_CELL_NAME}">Reactor vessel'
        f'<div style="{_DESC}">{_rv_desc}</div></td>'
        f'<td style="{_CELL_C}">{_m1(_rv_height_cm)}</td>'
        f'<td style="{_CELL_C}">{_m1(_rv_dia_cm)}</td>'
        f'<td style="{_CELL_C}">{_lb_str(_rv_mass_kg)}</td>'
        '</tr>'
    )
    if _has_guard:
        _rows_html.append(
            '<tr style="background:#ffffff;">'
            f'<td style="{_CELL_NAME}">Guard vessel'
            f'<div style="{_DESC}">{_gv_desc}</div></td>'
            f'<td style="{_CELL_C}">{_m1(_gv_height_cm)}</td>'
            f'<td style="{_CELL_C}">{_m1(_gv_dia_cm)}</td>'
            f'<td style="{_CELL_C}">{_lb_str(_gv_mass_kg)}</td>'
            '</tr>'
        )
    else:
        _rows_html.append(
            '<tr style="background:#ffffff;">'
            f'<td style="{_CELL_NAME};color:#64748b;">Guard vessel'
            f'<div style="{_DESC}">{_gv_na_desc}</div></td>'
            f'<td colspan="3" style="{_CELL_C};color:#64748b;">N/A not used for this reactor type</td>'
            '</tr>'
        )
    _rows_html.append(
        '<tr style="background:#f7f8fa;">'
        f'<td style="{_CELL_NAME}">Reactor Vessel Auxiliary Cooling System (cooling vessel + intake vessel)'
        f'<div style="{_DESC}">{_rvacs_desc}</div></td>'
        f'<td style="{_CELL_C}">{_m1(_rvacs_height_cm)}</td>'
        f'<td style="{_CELL_C}">{_m1(_rvacs_dia_cm)}</td>'
        f'<td style="{_CELL_C}">{_lb_str(_rvacs_mass_kg)}</td>'
        '</tr>'
    )

    _TH = ('padding:0.55rem 0.8rem;font-size:0.85rem;'
           'text-transform:uppercase;letter-spacing:0.06em;'
           'color:#3c4257;font-weight:600;')
    st.markdown(
        '<div style="margin-bottom:0.9rem;">'
        '<table style="width:100%;border-collapse:collapse;'
        'font-size:0.85rem;background:#ffffff;color:#0a2540;'
        'border:1px solid #bfdbfe;border-radius:8px;overflow:hidden;">'
        '<thead style="background:#f1f3f5;">'
        '<tr>'
        f'<th style="{_TH};text-align:left;">Component</th>'
        f'<th style="{_TH};text-align:center;">Height</th>'
        f'<th style="{_TH};text-align:center;">Diameter</th>'
        f'<th style="{_TH};text-align:center;">Mass (lb)</th>'
        '</tr></thead><tbody>'
        + ''.join(_rows_html) +
        '</tbody></table>'
        '</div>',
        unsafe_allow_html=True,
    )

    # ── Notes panel ──
    # Only items NOT already covered by the per-row descriptions in
    # the table above: global exclusions (shielding, coolant, support
    # gear), reactor-specific caveats (GCMR labeling, HPMR heat-pipe
    # treatment), and the "Guard vessel N/A" note for reactor
    # types without a bulk primary coolant.
    _gcmr_note = (
        '<li><strong>GCMR labeling:</strong> what is shown as '
        '"Reactor vessel" here maps to MOUSE\'s internal '
        '<em>Guard Vessel</em> field, because for the GCMR the '
        'outer pressure shell is the RPV, not the inner core '
        'barrel.</li>'
        if reactor_type == 'GCMR' else ''
    )
    _hpmr_note = (
        '<li><strong>HPMR heat pipes:</strong> steel mass is included '
        'using the modeled heat-pipe material and geometry; the Na '
        'working fluid is excluded with coolant inventories.</li>'
        if reactor_type == 'HPMR' else ''
    )
    _gv_na_note = (
        ''
        if _has_guard else
        '<li><strong>Guard vessel intentionally omitted:</strong> '
        'helium is inert (GCMR) and each heat pipe is individually '
        'sealed (HPMR), so neither has a bulk primary coolant '
        'requiring secondary containment.</li>'
    )
    st.markdown(
        '<div style="background:#eff6ff;border:1px solid #bfdbfe;'
        'border-radius:8px;padding:0.85rem 1.1rem;margin-bottom:0.9rem;'
        'font-size:0.85rem;line-height:1.55;color:#1B4F8C;">'
        '<div style="font-weight:600;font-size:0.85rem;'
        'text-transform:uppercase;letter-spacing:0.06em;'
        'color:#1B4F8C;margin-bottom:0.45rem;">Notes &amp; assumptions</div>'
        '<ul style="margin:0;padding-left:1.2rem;color:#1B4F8C;">'
        '<li><strong>Shielding excluded:</strong> in vessel '
        'shielding (B<sub>4</sub>C) and out of vessel shielding '
        '(WEP / concrete biological shield) are not included. '
        'Shielding adds significant mass and outer dimension to '
        'the module as shipped or as installed, depending on '
        'whether it ships with the reactor or is built on site.</li>'
        '<li><strong>Coolant excluded:</strong> primary coolant '
        'inventory (NaK for LTMR, He for GCMR, heat pipe Na for '
        'HPMR) is not included in the mass column.</li>'
        + _gcmr_note
        + _hpmr_note
        + _gv_na_note +
        '</ul>'
        '</div>',
        unsafe_allow_html=True,
    )

    # ── Transport-mode subsection header ──────────────────────
    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:1.25rem 0 0.85rem 0;">Transport Mode Compatibility</div>',
        unsafe_allow_html=True,
    )

    # Caveat panel: the fit-check is a scoping geometry/mass test
    # against generic shipping envelopes; it does NOT model the
    # real constraints that govern actual nuclear shipments.
    st.markdown(
        '<div style="background:#fffbeb;border:1px solid #fcd34d;border-radius:8px;'
        'padding:0.85rem 1.1rem;margin-bottom:0.9rem;'
        'font-size:0.85rem;line-height:1.45;color:#92400e;">'
        '<strong>Caveat:</strong> geometry and mass check only. '
        'Fueled modules, or modules removed after operation, '
        'require shielded casks that add tens of tonnes and can '
        'exceed the ISO container envelope.'
        '</div>',
        unsafe_allow_html=True,
    )

    # ── "How to read this" panel ──
    # Plain-English primer on what the cards below actually check, so
    # non-specialist readers don't have to reverse-engineer the
    # comparison from the badge labels.
    st.markdown(
        '<div style="background:#f7f8fa;border:1px solid #bfdbfe;border-radius:8px;'
        'padding:0.85rem 1.1rem;margin-bottom:0.9rem;'
        'font-size:0.85rem;line-height:1.45;color:#3c4257;">'
        '<strong>How to read this:</strong> the three columns below '
        '(Road, Rail, Sea) each list standard shipping limits used in '
        'that mode. Each card compares your reactor module\'s outer '
        'dimensions and mass to one limit. The reactor module is evaluated '
        'in a horizontal shipping orientation. The road summary uses green for no '
        'permit, yellow for a permitted overweight/oversize shipment, and '
        'orange for a superload above the MOUSE 150,000 lb screening '
        'threshold. The rail summary uses green for an ISO-container '
        'intermodal shipment, yellow for a dimensional shipment requiring '
        'railroad clearance, and orange for specialized railcar or '
        'special-train planning. The sea summary uses green for a standard '
        'ISO-container shipment and yellow when breakbulk or heavy-lift '
        'service is required.'
        '</div>',
        unsafe_allow_html=True,
    )

    # ── Transport-mode envelope limits ──
    # Organised by mode (Road / Rail / Sea), each column listing
    # the envelopes that apply to that mode. ISO containers appear
    # under all three modes (they're multimodal); AAR Plate F is
    # rail-only; US road no-permit is truck-only. Project-specific
    # heavy options (Schnabel rail, breakbulk sea) get a closing
    # note per mode rather than a fit-check.
    # Each envelope carries:
    # - `gross_t`: the number SHOWN to the user (container structural
    #   rating for ISOs, Gross Vehicle Weight for US road).
    # - `payload_t`: what the fit-check actually compares against
    #   (container gross − tare for ISOs; GVW − tractor+trailer tare
    #   for US road). This is the realistic cargo mass available.
    # - `height_note`: optional small qualifier shown after the
    #   height number (e.g. "(route-dep.)" for state-set limits).
    _envelopes = {
        'iso20': {
            'name': 'ISO 20 ft container',
            'width_m': 2.35,
            'height_m': 2.39,
            'length_m': 5.90,
            'payload_t': 21.7,  # 24 t total minus ~2.3 t container weight
            'help_text': (
                'A standard 20 ft shipping container, used worldwide on '
                'trucks, trains, and ships. Maximum cargo weight is '
                '~21.7 t (the 24 t total limit minus ~2.3 t for the '
                'container itself). Key advantage: the same container '
                'moves between truck, train, and ship without ever '
                'being unpacked.'
            ),
            'cite_html': (
                'Source: <a href="https://www.iso.org/standard/76912.html" '
                'target="_blank" style="color:#1B4F8C;">ISO 668:2020</a>'
            ),
            'road_note_html': (
                'Typical cargo limit on US trucks: 21 to 22 t '
                '(24 t total minus ~2.3 t container weight).'
            ),
        },
        'iso40': {
            'name': 'ISO 40 ft container',
            'width_m': 2.35,
            'height_m': 2.39,
            'length_m': 12.03,
            'payload_t': 26.78,  # 30.48 t total minus ~3.7 t container weight
            'help_text': (
                'A standard 40 ft shipping container with the same '
                'width and height as the 20 ft, roughly double the '
                'length. Maximum cargo weight is ~26.8 t (30.5 t total '
                'minus ~3.7 t container). On US roads, federal truck '
                'weight rules often reduce what you can actually load '
                'below this number.'
            ),
            'cite_html': (
                'Source: <a href="https://www.iso.org/standard/76912.html" '
                'target="_blank" style="color:#1B4F8C;">ISO 668:2020</a>'
            ),
            'road_note_html': (
                'Typical cargo limit on US trucks: 26 to 28 t. '
                'US federal 80,000 lb truck weight law often reduces '
                'this further in practice.'
            ),
        },
        'iso40hc': {
            'name': 'ISO 40 ft High Cube container',
            'width_m': 2.35,
            'height_m': 2.70,
            'length_m': 12.03,
            'payload_t': 26.58,  # 30.48 t total minus ~3.9 t container weight
            'help_text': (
                'Like a standard 40 ft container but ~30 cm taller '
                'inside. The usual choice when the module is too tall '
                'for a regular container. Maximum cargo weight ~26.6 t. '
                'On US trucks, loaded height (~4.1 m) sits right at the '
                'no permit clearance limit, so the truck chassis '
                'matters.'
            ),
            'cite_html': (
                'Source: <a href="https://www.iso.org/standard/76912.html" '
                'target="_blank" style="color:#1B4F8C;">ISO 668:2020</a>'
            ),
            'road_note_html': (
                'Typical cargo limit on US trucks: 26 to 27 t. '
                'Loaded height on a standard chassis is ~4.1 m, right '
                'at the US no permit clearance limit.'
            ),
        },
        'us_no_permit': {
            'name': 'US road transport category',
            'width_m': 2.59,
            'height_m': 4.11,
            'height_note': '(loaded vehicle planning limit)',
            'length_m': None,
            'payload_t': None,
            'road_category': True,
            # MOUSE screening assumptions. The legal 80,000 lb threshold
            # applies to gross vehicle weight, not reactor-module mass.
            'tractor_trailer_tare_lb': 32000.0,
            'deck_height_m': 1.52,
            'no_permit_gvw_lb': 80000.0,
            'superload_gvw_lb': 150000.0,
            'help_text': (
                'The reactor module is evaluated horizontally. Estimated loaded '
                'vehicle weight equals the assembled reactor-module mass plus a '
                '32,000 lb tractor/trailer planning allowance; the shipping frame '
                'or cradle is not yet included. Loaded height equals the module '
                'diameter plus a 1.52 m trailer deck. Green means no permit is '
                'needed under the MOUSE screening limits. Yellow means an '
                'overweight or oversize permit is needed. Orange means the '
                'estimated loaded vehicle weight exceeds the MOUSE 150,000 lb '
                'superload screening threshold. Actual requirements vary by '
                'state and route.'
            ),
            'cite_html': (
                'Sources: <a href="https://www.law.cornell.edu/uscode/text/23/127" '
                'target="_blank" style="color:#1B4F8C;">23 USC § 127 '
                '(federal weight limits)</a>; '
                '<a href="https://www.ecfr.gov/current/title-23/chapter-I/'
                'subchapter-G/part-658/section-658.15" target="_blank" '
                'style="color:#1B4F8C;">23 CFR § 658.15 (width)</a>. '
                'The 4.11 m loaded-height value and 150,000 lb superload '
                'threshold are MOUSE planning assumptions.'
            ),
        },
        'aar_plate_f': {
            'name': 'Rail flatcar (oversized cargo)',
            'width_m': 3.25,
            'height_m': 5.18,
            'height_note': '(above rail)',
            'length_m': None,
            'payload_t': None,
            'help_text': (
                'Preliminary dimensional-envelope check for an open rail '
                'flatcar. A geometric fit does not mean the shipment is '
                'approved. Dimensional rail loads require railroad '
                'clearance, and final acceptance depends on the complete '
                'loaded-car profile, axle and wheel loading, center of '
                'gravity, securement, route obstructions, interchange '
                'railroads, and loading/unloading tracks. The simplified '
                'MOUSE check does not yet apply a route-specific mass or '
                'length limit and does not model the narrowing upper '
                'corners of the loading gauge.'
            ),
            'cite_html': (
                'Screening envelope; clearance requirements: '
                '<a href="https://www.up.com/shipping/machinery/seven-steps" '
                'target="_blank" style="color:#1B4F8C;">Union Pacific '
                'dimensional-load process</a>'
            ),
        },
    }

    _mode_groups = [
        {'name': 'Road',
         'envelope_keys': ['iso20', 'iso40', 'iso40hc', 'us_no_permit']},
        {'name': 'Rail',
         'envelope_keys': ['iso20', 'iso40', 'iso40hc', 'aar_plate_f']},
        {'name': 'Sea',
         'envelope_keys': ['iso20', 'iso40', 'iso40hc']},
    ]

    # Reactor envelope + total mass used by every fit-check.
    _rvacs_dia_m = _rvacs_dia_cm / 100.0
    _rvacs_h_m = _rvacs_height_cm / 100.0
    _badge_total_kg = (_reactor_mass_kg + _rv_mass_kg
                       + _gv_mass_kg + _rvacs_mass_kg)
    _badge_total_t = _badge_total_kg / 1000.0

    def _module_fits_envelope(env):
        """Return a preliminary horizontal mass-and-envelope fit."""
        _horizontal = {
            'width_m': _rvacs_dia_m,
            'height_m': _rvacs_dia_m,
            'length_m': _rvacs_h_m,
        }
        _geometry_ok = (
            _horizontal['width_m'] <= env['width_m']
            and _horizontal['height_m'] <= env['height_m']
            and (
                env['length_m'] is None
                or _horizontal['length_m'] <= env['length_m']
            )
        )
        _weight_ok = (
            env.get('payload_t') is None
            or _badge_total_t <= env['payload_t']
        )
        return _geometry_ok and _weight_ok

    def _render_rail_sea_category_card(mode_name):
        """Render mode-level logistics categories for rail and sea.

        These colors describe logistics complexity, not a statutory permit
        ladder. Detailed carrier, route, port, vessel, and securement review
        remains necessary for an actual shipment.
        """
        _iso_keys = ('iso20', 'iso40', 'iso40hc')
        _fitting_isos = [
            _envelopes[key]['name']
            for key in _iso_keys
            if _module_fits_envelope(_envelopes[key])
        ]

        if mode_name == 'Rail':
            _title = 'Rail transport category'
            _help_text = (
                'This is a preliminary logistics-complexity classification, '
                'not a permit determination. Green means at least one ISO '
                'container fit check passes, allowing a standard intermodal '
                'concept. Yellow means no ISO container fits but the '
                'simplified rail-flatcar envelope fits; a dimensional-load '
                'clearance proposal is required. Orange means the simplified '
                'flatcar envelope also fails, so specialized equipment such '
                'as a heavy-duty or Schnabel car and possibly special-train '
                'planning may be needed. Even a green result still requires '
                'carrier acceptance, route and interchange review, loading '
                'and unloading track approval, axle-loading and center-of-'
                'gravity checks, securement design, and a commercial quote.'
            )
            if _fitting_isos:
                _badge_text = '✓ Rail (standard ISO-container/intermodal shipment)'
                _bc = ('#15803d', '#dcfce7', '#bbf7d0')
                _reason = 'Fits: ' + ', '.join(_fitting_isos) + '.'
            elif _module_fits_envelope(_envelopes['aar_plate_f']):
                _badge_text = '⚠ Rail (dimensional shipment — clearance required)'
                _bc = ('#a16207', '#fef9c3', '#fde047')
                _reason = (
                    'ISO-container limits exceeded; rail flatcar fits. '
                    'Railroad clearance required.'
                )
            else:
                _badge_text = '⚠ Rail (specialized railcar or special-train planning required)'
                _bc = ('#c2410c', '#ffedd5', '#fdba74')
                _reason = (
                    'Standard container and flatcar limits exceeded. '
                    'Specialized rail planning required.'
                )
            _source_html = (
                'Clearance process: <a href="https://www.up.com/shipping/'
                'machinery/seven-steps" target="_blank" '
                'style="color:#1B4F8C;">Union Pacific dimensional loads</a>'
            )
        elif mode_name == 'Sea':
            _title = 'Sea transport category'
            _help_text = (
                'This is a preliminary logistics-complexity classification, '
                'not a maritime permit determination. Green means at least '
                'one ISO-container fit check passes. Yellow means no ISO '
                'container fits and the module would require a breakbulk, '
                'multipurpose, or heavy-lift shipping concept. A yellow '
                'result requires vessel and port suitability checks, lift '
                'planning, stowage and securing design, and confirmation of '
                'terminal handling and inland road or rail legs. MOUSE does '
                'not yet model crane capacity, hatch and deck limits, port '
                'access, or shipment-specific handling costs.'
            )
            if _fitting_isos:
                _badge_text = '✓ Sea (standard ISO-container shipment)'
                _bc = ('#15803d', '#dcfce7', '#bbf7d0')
                _reason = 'Fits: ' + ', '.join(_fitting_isos) + '.'
            else:
                _badge_text = '⚠ Sea (breakbulk or heavy-lift shipment required)'
                _bc = ('#a16207', '#fef9c3', '#fde047')
                _reason = (
                    'ISO-container limits exceeded. Breakbulk or heavy-lift '
                    'planning required.'
                )
            _source_html = (
                'Cargo planning and securing: '
                '<a href="https://www.imo.org/en/ourwork/safety/pages/'
                'css-code.aspx" target="_blank" style="color:#1B4F8C;">'
                'IMO CSS Code</a>'
            )
        else:
            return ''

        return (
            '<div style="background:#ffffff;border:1px solid #bfdbfe;'
            'border-radius:8px;padding:0.7rem 0.85rem;margin-bottom:0.6rem;'
            'color:#0a2540;">'
            f'<div style="font-weight:600;font-size:0.85rem;color:#0a2540;'
            f'margin-bottom:0.35rem;">{_title}'
            f'{_help_icon(_help_text)}</div>'
            f'<div style="display:inline-block;background:{_bc[1]};'
            f'border:1px solid {_bc[2]};color:{_bc[0]};font-size:0.85rem;'
            f'font-weight:600;padding:0.15rem 0.5rem;border-radius:8px;'
            f'margin-bottom:0.4rem;max-width:100%;white-space:normal;'
            f'line-height:1.3;">{_badge_text}</div>'
            f'<div style="font-size:0.82rem;color:#64748b;line-height:1.35;'
            f'margin-bottom:0.25rem;">{_reason}</div>'
            f'<div style="font-size:0.85rem;color:#64748b;line-height:1.4;">'
            f'{_source_html}</div>'
            '</div>'
        )

    def _render_road_category_card(env):
        """Render the three-level US road screening classification."""
        _tare_lb = float(env['tractor_trailer_tare_lb'])
        _deck_m = float(env['deck_height_m'])
        _loaded_vehicle_weight_lb = _badge_total_kg * 2.2046226218 + _tare_lb

        # Horizontal reactor shipment: diameter controls width and package
        # height; modeled module height controls transport length.
        _package_width_m = _rvacs_dia_m
        _package_height_m = _rvacs_dia_m
        _loaded_height_m = _package_height_m + _deck_m
        _width_ok = _package_width_m <= env['width_m']
        _height_ok = _loaded_height_m <= env['height_m']
        _geometry_ok = _width_ok and _height_ok

        if _loaded_vehicle_weight_lb > env['superload_gvw_lb']:
            _badge_text = '⚠ US road (superload)'
            _bc = ('#c2410c', '#ffedd5', '#fdba74')
            _reason = (
                f'Estimated loaded vehicle weight: '
                f'{_loaded_vehicle_weight_lb:,.0f} lb — superload.'
            )
        elif _geometry_ok and _loaded_vehicle_weight_lb <= env['no_permit_gvw_lb']:
            _badge_text = '✓ US road (no permits needed)'
            _bc = ('#15803d', '#dcfce7', '#bbf7d0')
            _reason = (
                f'Estimated loaded vehicle weight: '
                f'{_loaded_vehicle_weight_lb:,.0f} lb. Width and loaded '
                'height are within no-permit limits.'
            )
        else:
            _badge_text = '⚠ US road (permit required)'
            _bc = ('#a16207', '#fef9c3', '#fde047')
            _triggers = []
            if _loaded_vehicle_weight_lb > env['no_permit_gvw_lb']:
                _triggers.append('weight')
            if not _width_ok:
                _triggers.append('width')
            if not _height_ok:
                _triggers.append('loaded height')
            _reason = (
                f'Estimated loaded vehicle weight: '
                f'{_loaded_vehicle_weight_lb:,.0f} lb. Permit required due '
                'to ' + ', '.join(_triggers) + '.'
            )

        return (
            '<div style="background:#ffffff;border:1px solid #bfdbfe;'
            'border-radius:8px;padding:0.7rem 0.85rem;margin-bottom:0.6rem;'
            'color:#0a2540;">'
            f'<div style="font-weight:600;font-size:0.85rem;color:#0a2540;'
            f'margin-bottom:0.35rem;">{env["name"]}'
            f'{_help_icon(env["help_text"])}</div>'
            f'<div style="display:inline-block;background:{_bc[1]};'
            f'border:1px solid {_bc[2]};color:{_bc[0]};font-size:0.85rem;'
            f'font-weight:600;padding:0.15rem 0.5rem;border-radius:8px;'
            f'margin-bottom:0.4rem;max-width:100%;white-space:normal;'
            f'line-height:1.3;">{_badge_text}</div>'
            f'<div style="font-size:0.85rem;color:#3c4257;line-height:1.4;'
            f'margin-bottom:0.25rem;">w ≤ {env["width_m"]:.2f} m '
            f'&nbsp;|&nbsp; loaded h ≤ {env["height_m"]:.2f} m '
            f'&nbsp;|&nbsp; no permit ≤ {env["no_permit_gvw_lb"]:,.0f} lb '
            f'&nbsp;|&nbsp; permit required: &gt;{env["no_permit_gvw_lb"]:,.0f}–'
            f'{env["superload_gvw_lb"]:,.0f} lb '
            f'&nbsp;|&nbsp; superload &gt; {env["superload_gvw_lb"]:,.0f} lb</div>'
            f'<div style="font-size:0.82rem;color:#64748b;line-height:1.35;'
            f'margin-bottom:0.25rem;">{_reason}</div>'
            f'<div style="font-size:0.85rem;color:#64748b;line-height:1.4;">'
            f'{env["cite_html"]}</div>'
            '</div>'
        )

    def _render_envelope_card(env, mode_name=''):
        # Horizontal reactor shipment: width = diameter, height = diameter,
        # and length = modeled module height.
        _horizontal = {
            'width_m': _rvacs_dia_m,
            'height_m': _rvacs_dia_m,
            'length_m': _rvacs_h_m,
        }
        _width_ok = _horizontal['width_m'] <= env['width_m']
        _height_ok = _horizontal['height_m'] <= env['height_m']
        _length_ok = (
            env['length_m'] is None
            or _horizontal['length_m'] <= env['length_m']
        )
        _geometry_ok = _width_ok and _height_ok and _length_ok
        _weight_ok = (
            env.get('payload_t') is None
            or _badge_total_t <= env['payload_t']
        )
        _fits = _geometry_ok and _weight_ok
        _bc = (
            ('#15803d', '#dcfce7', '#bbf7d0')
            if _fits else ('#b91c1c', '#fee2e2', '#fecaca')
        )

        if _fits:
            _badge_text = '✓ fits'
            _fit_note = ''
        elif _geometry_ok and not _weight_ok:
            _badge_text = '✗ exceeds weight'
            _fit_note = 'Geometry fits, but the cargo-mass limit is exceeded.'
        else:
            _badge_text = '✗ does not fit'
            _failures = []
            if not _width_ok:
                _failures.append('width')
            if not _height_ok:
                _failures.append('height')
            if not _length_ok:
                _failures.append('length')
            _fit_note = 'Exceeds ' + ', '.join(_failures) + '.'

        _height_note = env.get('height_note', '')
        _height_note_str = (
            f' <span style="color:#64748b;">{_height_note}</span>'
            if _height_note else ''
        )
        _len_str = (
            f' &nbsp;|&nbsp; len ≤ {env["length_m"]:.2f} m'
            if env['length_m'] is not None else ''
        )
        _wt_str = (
            f' &nbsp;|&nbsp; weight ≤ {env["payload_t"]:.1f} t'
            if env.get('payload_t') is not None else ''
        )
        _mode_key = (mode_name or '').lower() + '_note_html'
        _mode_note = env.get(_mode_key, '')
        return (
            '<div style="background:#ffffff;border:1px solid #bfdbfe;'
            'border-radius:8px;padding:0.7rem 0.85rem;margin-bottom:0.6rem;'
            'color:#0a2540;">'
            f'<div style="font-weight:600;font-size:0.85rem;color:#0a2540;'
            f'margin-bottom:0.35rem;">'
            f'{env["name"]}{_help_icon(env["help_text"])}'
            f'</div>'
            f'<div style="display:inline-block;background:{_bc[1]};'
            f'border:1px solid {_bc[2]};color:{_bc[0]};font-size:0.85rem;'
            f'font-weight:600;padding:0.15rem 0.5rem;border-radius:8px;'
            f'margin-bottom:0.4rem;max-width:100%;white-space:normal;'
            f'line-height:1.3;">{_badge_text}</div>'
            f'<div style="font-size:0.85rem;color:#3c4257;margin-bottom:0.25rem;">'
            f'w ≤ {env["width_m"]:.2f} m &nbsp;|&nbsp; '
            f'h ≤ {env["height_m"]:.2f} m{_height_note_str}{_len_str}{_wt_str}'
            f'</div>'
            + (f'<div style="font-size:0.82rem;color:#64748b;line-height:1.35;'
               f'margin-bottom:0.25rem;">{_fit_note}</div>'
               if _fit_note else '')
            + (f'<div style="font-size:0.85rem;color:#3c4257;line-height:1.4;'
               f'margin-bottom:0.25rem;">{env["note_html"]}</div>'
               if env.get('note_html') else '')
            + (f'<div style="font-size:0.85rem;color:#3c4257;line-height:1.4;'
               f'margin-bottom:0.25rem;">{_mode_note}</div>'
               if _mode_note else '')
            + f'<div style="font-size:0.85rem;color:#64748b;line-height:1.4;">'
              f'{env["cite_html"]}</div>'
            + '</div>'
        )

    # Three columns: Road / Rail / Sea. Each column gets a small
    # uppercase mode-group label followed by the envelope cards and
    # a project-specific closing note for that mode.
    _mode_cols = st.columns(3, gap='medium')
    for _col, _group in zip(_mode_cols, _mode_groups):
        with _col:
            st.markdown(
                f'<div style="font-size:0.85rem;font-weight:600;'
                f'color:#64748b;text-transform:uppercase;'
                f'letter-spacing:0.09em;margin-bottom:0.6rem;">'
                f'{_group["name"]}</div>',
                unsafe_allow_html=True,
            )
            _cards_html = ''.join(
                _render_road_category_card(_envelopes[k])
                if _envelopes[k].get('road_category')
                else _render_envelope_card(
                    _envelopes[k], mode_name=_group['name']
                )
                for k in _group['envelope_keys']
            )
            if _group['name'] in ('Rail', 'Sea'):
                _cards_html += _render_rail_sea_category_card(_group['name'])
            st.markdown(_cards_html, unsafe_allow_html=True)

    # Footnote — what each badge does and doesn't check, plus the
    # oversized fallback for every mode. Rendered as two parallel
    # bullets so they read as siblings, not as text + an afterthought.
    st.markdown(
        '<div style="background:#fffbeb;border:1px solid #fcd34d;border-radius:8px;'
        'padding:0.85rem 1.1rem;margin-bottom:1rem;font-size:0.85rem;line-height:1.45;color:#92400e;">'
        '<strong>Notes:</strong>'
        '<ul style="margin:0.4rem 0 0 1.2rem;padding:0;">'
        '<li>Each badge compares the assembled reactor module: the '
        'outermost RVACS envelope (diameter and height) and the sum of all '
        'component masses. The reactor module is evaluated horizontally. '
        'Per-component dimensions and masses are shown above.</li>'
        '<li>The US road classification estimates loaded vehicle weight '
        'using the reactor-module mass plus a 32,000 lb tractor/trailer '
        'tare and calculates loaded height using a 1.52 m trailer deck. '
        'Shipping-frame mass is not yet included.</li>'
        '<li>The rail and sea colors describe logistics complexity, not '
        'a legal permit ladder. Rail uses green for an ISO-container '
        'intermodal concept, yellow for a dimensional shipment requiring '
        'railroad clearance, and orange for specialized railcar or '
        'special-train planning. Sea uses green for ISO-container shipping '
        'and yellow for breakbulk or heavy-lift planning.</li>'
        '<li>For modules larger than these envelopes, road uses permitted '
        'heavy haul, rail uses project-specific heavy-duty or Schnabel cars, '
        'and sea uses breakbulk or heavy-lift vessels. Actual feasibility '
        'depends on the selected carrier, route, railroads, ports, vessel, '
        'lifting plan, cargo securement, and first/last-mile access.</li>'
        '</ul>'
        '</div>',
        unsafe_allow_html=True,
    )

    # Cost hero (Capital + Levelized cards + FOAK/NOAK legend)
    # was here; moved up into info_col under the k_eff plot as
    # part of the Band 1 layout.

    # ═══════════════════════════════════════════════════════════════
    # TAB 2 PHYSICS & THERMAL
    # Sections 4-5: Core neutronics at glance · What are the
    # thermal-hydraulic conditions at a glance (Section 6 Transportability
    # is rendered in the Design & Economics tab.)
    # ═══════════════════════════════════════════════════════════════
    # ═══════════════════════════════════════════════════════════
    # BAND 3 Microreactor Design & Performance Metrics
    # 3 columns: Neutronics · Thermal-Hydraulics · Fuel cycle.
    # ═══════════════════════════════════════════════════════════
    _section_header(
        '03',
        'Microreactor Scoping Metrics',
        'How does the design perform across neutronics, thermal-hydraulics, '
        'and fuel cycle?',
    )

    # Helper: render a stacked subsection with a header, optional
    # scoping-caveat popover, and the queued cards in a 2-column grid.
    def _render_metric_subsection(title, section_num, caveat_md=None,
                                  empty_msg=None, top_margin='0'):
        st.markdown(
            f'<div style="font-size:1rem;font-weight:700;color:#0a2540;'
            f'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
            f'margin:{top_margin} 0 0.6rem 0;">'
            f'{title}</div>',
            unsafe_allow_html=True,
        )
        if caveat_md is not None:
            with st.popover('⚠ Scoping caveat', width='content'):
                st.markdown(caveat_md)
        cards = _section_cards.get(section_num, [])
        if not cards:
            st.info(empty_msg or 'No metrics available for this case.')
            return
        # Cards split alternately between the left and right columns
        # so reading order goes top-left → top-right → next row → …
        _left, _right = st.columns(2, gap='medium')
        for i, _h in enumerate(cards):
            target = _left if i % 2 == 0 else _right
            with target:
                st.markdown(_h, unsafe_allow_html=True)

    _render_metric_subsection(
        'Neutronics',
        section_num=4,
        caveat_md=(
            'These are first-pass scoping designs for initial economic '
            'analysis. They have **NOT** been checked for shutdown '
            'margin, reactivity coefficients, kinetics, control-drum '
            'worth, or transient response. The k_eff, peaking factor, '
            'and leakage values are sufficient to size the fuel cycle '
            'and bound discharge burnup, but a full safety / licensing '
            'analysis would require additional Monte Carlo and '
            'depletion calculations.'
        ),
        empty_msg='No neutronics metrics available — typically means '
                  'the case is subcritical.',
        top_margin='0',
    )

    _render_metric_subsection(
        'Thermal-Hydraulics',
        section_num=5,
        caveat_md=(
            'Heat flux, power density, and per-component thermal duty '
            'are reported for first-order sizing only. The app does '
            '**NOT** verify passive heat-removal capability, peak '
            'fuel / cladding / coolant temperatures, transient cooling '
            'under loss-of-flow or loss-of-heat-sink, or material '
            'temperature limits. Steady-state design margins must be '
            'confirmed with a coupled thermal-hydraulics analysis '
            'before relying on these numbers.'
        ),
        empty_msg='No thermal-hydraulic metrics available for this case.',
        top_margin='1.25rem',
    )

    _render_metric_subsection(
        'Fuel cycle',
        section_num=3,
        caveat_md=None,
        empty_msg='No fuel-cycle metrics available for this case.',
        top_margin='1.25rem',
    )

    # ═══════════════════════════════════════════════════════════════
    # BAND 4 Cost Detail (cost drivers, full breakdown, costs in
    # perspective). Rendered as a single column with successive rows.
    # ═══════════════════════════════════════════════════════════════
    _section_header(
        '04',
        'Cost Decomposition',
        'What&#39;s driving the cost?',
    )

    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:0 0 0.85rem 0;">Cost drivers</div>',
        unsafe_allow_html=True,
    )
    _drv = _curated_cost_drivers(enriched_df, min_foak_lcoe=10.0,
                                 max_drivers=28)

    if _drv.empty:
        st.info('No curated cost drivers with FOAK LCOE >= 10 $/MWh found.')
    else:
        _hidden_foak, _hidden_noak, _hidden_count = (
            _hidden_cost_driver_lcoe_totals(enriched_df, min_foak_lcoe=10.0)
        )
        _hidden_note = (
            f' The {_hidden_count} hidden small drivers sum to '
            f'{_hidden_foak:.1f} $/MWh FOAK and {_hidden_noak:.1f} $/MWh NOAK.'
            if _hidden_count else
            ' No curated cost drivers are hidden by the 10 $/MWh threshold.'
        )
        st.markdown(
            '<p style="color:#64748b;font-size:0.85rem;margin-bottom:0.6rem;">'
            'Curated cost-driver view using account levels familiar to most '
            'users. Drivers are sorted by FOAK LCOE contribution; rows below '
            '10 $/MWh are hidden. Bars show FOAK (orange) and NOAK (navy). '
            'Error bars show +/-1 standard deviation across Monte Carlo '
            f'samples.{_hidden_note}</p>',
            unsafe_allow_html=True,
        )
        _category_inline = '&nbsp; '.join(
            f'<span style="color:{COST_DRIVER_CATEGORY_PALETTE[_c]};'
            f'font-weight:700;">{html.escape(_c)}</span>'
            for _c in COST_DRIVER_CATEGORY_ORDER
        )
        st.markdown(
            '<div style="border:1px solid #e2e8f0;border-radius:6px;'
            'background:#ffffff;padding:0.55rem 0.75rem;'
            'margin:0 0 1rem 0;font-size:0.85rem;color:#3c4257;'
            'line-height:1.45;">'
            'Each x-axis label is colored by cost category: '
            f'{_category_inline}.'
            '</div>',
            unsafe_allow_html=True,
        )

        _chunk_size = 7
        _chunks = [_drv.iloc[i:i + _chunk_size].copy()
                   for i in range(0, len(_drv), _chunk_size)]
        _plot_col, _ = st.columns([7, 3])
        with _plot_col:
            for _i, _chunk in enumerate(_chunks[:4]):
                _start = _i * _chunk_size + 1
                _end = _start + len(_chunk) - 1
                _title = (f'Top {_end} cost drivers'
                          if _i == 0 else f'Cost drivers ranked {_start}-{_end}')
                st.markdown(
                    f'<div style="font-size:0.9rem;font-weight:600;color:#0a2540;'
                    f'margin:{0 if _i == 0 else 1.25}rem 0 0.45rem 0;">'
                    f'{_title}</div>',
                    unsafe_allow_html=True,
                )
                st.altair_chart(
                    _grouped_lcoe_bars_chart(
                        _chunk,
                        label_col='Driver Label',
                        height=390,
                        label_angle=-35,
                        category_col='Category',
                    ),
                    use_container_width=True,
                )
        _render_cost_driver_definition_table(_drv)

    # ─── Band 4 (cont.) Full breakdown ──────────────────────────────
    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:1.25rem 0 0.85rem 0;">Full breakdown</div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        f'<p style="color:#64748b;font-size:0.85rem;margin-bottom:1rem;">'
        f'Full cost breakdown by Code of Accounts. Highlighted rows are parent accounts. '
        f'First-of-a-kind (FOAK) columns represent the first deployed unit; '
        f'nth-of-a-kind (NOAK) columns represent the mature deployment unit selected in the inputs. '
        f'The cost columns include both capitalized costs and annualized operating, maintenance, '
        f'replacement, decommissioning, and fuel costs. '
        f'The levelized cost of energy (LCOE) columns estimate how much each account contributes '
        f'to the total electricity cost. '
        f'Values are shown as <strong>mean</strong> '
        f'<span style="color:#64748b;">[mean &minus; 1&sigma; to mean + 1&sigma;]</span>; '
        f'a single value means no uncertainty was computed for that account. '
        f'All dollar values in {ESCALATION_YEAR} USD.</p>',
        unsafe_allow_html=True,
    )

    _foak_col = next((c for c in display_df.columns
                      if c.startswith('FOAK Estimated Cost (') and 'std' not in c), None)
    _noak_col = next((c for c in display_df.columns
                      if c.startswith('NOAK Estimated Cost (') and 'std' not in c), None)
    _foak_std = next((c for c in display_df.columns if 'FOAK Estimated Cost std' in c), None)
    _noak_std = next((c for c in display_df.columns if 'NOAK Estimated Cost std' in c), None)
    _have_lcoe = 'FOAK LCOE' in enriched_df.columns

    def _fmt_table_money(x):
        v = _fmt_table_val(x)
        return v if v == '-' else f'${v}'

    def _value_cell(main, low=None, high=None):
        if main == '-':
            return '<span class="mouse-table-empty">-</span>'
        if low is None or high is None or low == high:
            return f'<span class="mouse-value-main">{main}</span>'
        return (
            f'<span class="mouse-value-main">{main}</span>'
            f'<span class="mouse-value-range">[{low} - {high}]</span>'
        )

    def _pm(mean_val, std_val):
        m = _fmt_table_money(mean_val)
        if m == '-':
            return _value_cell('-')
        try:
            mn = float(mean_val)
            sd = float(std_val)
        except (TypeError, ValueError):
            return _value_cell(m)
        if sd == 0 or np.isnan(sd):
            return _value_cell(m)
        lo = _fmt_table_money(max(0, mn - sd))
        hi = _fmt_table_money(mn + sd)
        return _value_cell(m, lo, hi)

    def _fmt_lcoe_tab(v):
        try:
            v = float(v)
        except (TypeError, ValueError):
            return '-'
        if np.isnan(v) or v <= 0:
            return '-'
        if v < 1:
            return '&lt; $1'
        return f'${int(round(v))}'

    def _pm_lcoe(mean_val, std_val):
        m = _fmt_lcoe_tab(mean_val)
        if m == '-':
            return _value_cell('-')
        try:
            mn = float(mean_val)
            sd = float(std_val)
        except (TypeError, ValueError):
            return _value_cell(m)
        if np.isnan(sd) or sd <= 0:
            return _value_cell(m)
        lo = _fmt_lcoe_tab(max(0.0, mn - sd))
        hi = _fmt_lcoe_tab(mn + sd)
        return _value_cell(m, lo, hi)

    def _fmt_account(x):
        if isinstance(x, str):
            return x
        try:
            v = float(x)
            return str(int(v)) if v == int(v) else f'{v:g}'
        except (TypeError, ValueError):
            return str(x)

    table_df = display_df[['Account', 'Account Title']].copy()
    table_df['Account'] = table_df['Account'].apply(_fmt_account)

    _sf = display_df[_foak_std] if _foak_std else pd.Series('-', index=display_df.index)
    _sn = display_df[_noak_std] if _noak_std else pd.Series('-', index=display_df.index)
    table_df['FOAK Cost ($)'] = [_pm(m, s) for m, s in zip(display_df[_foak_col], _sf)]
    table_df['NOAK Cost ($)'] = [_pm(m, s) for m, s in zip(display_df[_noak_col], _sn)]

    if _have_lcoe:
        _ei = display_df.index
        _e_fl = enriched_df['FOAK LCOE'].reindex(_ei)
        _e_nl = enriched_df['NOAK LCOE'].reindex(_ei)
        _e_fls = enriched_df['FOAK LCOE_std'].reindex(_ei) \
                 if 'FOAK LCOE_std' in enriched_df.columns else pd.Series(np.nan, index=_ei)
        _e_nls = enriched_df['NOAK LCOE_std'].reindex(_ei) \
                 if 'NOAK LCOE_std' in enriched_df.columns else pd.Series(np.nan, index=_ei)
        table_df['FOAK LCOE ($/MWh)'] = [_pm_lcoe(m, s) for m, s in zip(_e_fl, _e_fls)]
        table_df['NOAK LCOE ($/MWh)'] = [_pm_lcoe(m, s) for m, s in zip(_e_nl, _e_nls)]

    def _account_level(acct_str):
        try:
            v = float(str(acct_str).strip())
        except (ValueError, TypeError):
            return '-'
        ip = int(v)
        n = len(str(ip))
        if n <= 2:
            return 0 if ip % 10 == 0 else 1
        elif n == 3:
            return 3 if v != ip else 2
        return 4

    _acct_levels = [_account_level(a) for a in table_df['Account']]
    _idx_to_pos = {idx: pos for pos, idx in enumerate(table_df.index)}

    _EM = '\u2003'
    _PREFIX = {'-': '', 0: '', 1: _EM, 2: f'{_EM}{_EM}',
               3: f'{_EM}{_EM}{_EM}', 4: f'{_EM}{_EM}{_EM}{_EM}'}
    table_df['Account Title'] = [
        f"{_PREFIX.get(lv, _EM * 4)}{title}"
        for lv, title in zip(_acct_levels, display_df['Account Title'])
    ]

    _LEVEL_STYLE = {
        '-': ('background-color:#fffaf0', 'color:#7c2d12', 'font-weight:600'),
         0: ('background-color:#0a2540', 'color:#ffffff', 'font-weight:600'),
         1: ('background-color:#eaf4fb', 'color:#0a2540', 'font-weight:600'),
         2: ('background-color:#f7fbfe', 'color:#0a2540', 'font-weight:500'),
         3: ('background-color:#ffffff', 'color:#3c4257', 'font-weight:400'),
         4: ('background-color:#ffffff', 'color:#3c4257', 'font-weight:400'),
    }

    # Hand-rolled HTML table instead of pandas Styler. Each call to
    # df.style.apply() allocates ~one functools.partial(_default_formatter)
    # per cell that Streamlit's forward-message cache pins across reruns,
    # accumulating ~700-1000 partials per Run with no mechanism to release
    # them. The HTML approach uses inline styles directly and creates no
    # pandas internals at all. Tradeoff: loses st.dataframe's interactive
    # column resize / sort UI. The Excel download below still gives users
    # interactive analysis.
    import html as _html_lib
    _num_cols = [c for c in table_df.columns if c not in ('Account', 'Account Title')]
    _num_col_set = set(_num_cols)
    _thead_style = (
        'background:#0a2540;color:#fff;font-weight:600;'
        'font-size:0.85rem;text-transform:uppercase;letter-spacing:0.05em;'
        'padding:10px 12px;text-align:left;position:sticky;top:0;z-index:5;'
    )
    _HEADER_LABELS = {
        'FOAK Cost ($)': 'FOAK Cost<br>($)',
        'NOAK Cost ($)': 'NOAK Cost<br>($)',
        'FOAK LCOE ($/MWh)': 'FOAK LCOE<br>($/MW<sub>e</sub>h)',
        'NOAK LCOE ($/MWh)': 'NOAK LCOE<br>($/MW<sub>e</sub>h)',
    }
    _thead_cells = []
    for col in table_df.columns:
        _sticky = ''
        if col == 'Account':
            _sticky = 'left:0;z-index:7;min-width:6.25rem;width:6.25rem;'
        elif col == 'Account Title':
            _sticky = 'left:6.25rem;z-index:7;min-width:31rem;width:31rem;'
        elif col == 'FOAK LCOE ($/MWh)':
            _sticky = 'border-left:2px solid #8fb3d1;'
        _thead_cells.append(
            f'<th style="{_thead_style}{_sticky}">'
            f'{_HEADER_LABELS.get(col, _html_lib.escape(str(col)))}</th>'
        )
    _thead_html = ''.join(_thead_cells)

    _rows_html = []
    _cols_list = list(table_df.columns)
    for pos, (_idx, row) in enumerate(table_df.iterrows()):
        lv = _acct_levels[pos]
        bg, fg, fw = _LEVEL_STYLE.get(
            lv, ('background-color:#ffffff', 'color:#3c4257', 'font-weight:400')
        )
        row_style = f'{bg};{fg};{fw}'
        row_class = f'mouse-level-{lv}' if lv != '-' else 'mouse-level-summary'
        cells = []
        for col in _cols_list:
            align = 'right' if col in _num_col_set else 'left'
            sticky = ''
            extra_class = ' mouse-num-cell' if col in _num_col_set else ''
            if col == 'Account':
                sticky = f'position:sticky;left:0;z-index:3;{bg};'
            elif col == 'Account Title':
                sticky = f'position:sticky;left:6.25rem;z-index:3;{bg};'
            elif col == 'FOAK LCOE ($/MWh)':
                sticky = 'border-left:2px solid #c5d7e6;'
            if col in _num_col_set:
                val = str(row[col])
            elif col == 'Account Title':
                val = str(row[col])
            else:
                val = _html_lib.escape(str(row[col]))
            cells.append(
                f'<td class="{row_class}{extra_class}" style="{row_style};{sticky}'
                f'text-align:{align};padding:7px 12px;border-bottom:1px solid #e8eef5;">'
                f'{val}</td>'
            )
        _rows_html.append('<tr>' + ''.join(cells) + '</tr>')

    _html_table = (
        '<style>'
        '.mouse-breakdown-wrap table{border-collapse:separate;border-spacing:0;'
        'width:100%;font-size:0.88rem;line-height:1.25;}'
        '.mouse-breakdown-wrap th,.mouse-breakdown-wrap td{vertical-align:middle;}'
        '.mouse-breakdown-wrap th:nth-child(n+3),.mouse-breakdown-wrap .mouse-num-cell{text-align:right;}'
        '.mouse-breakdown-wrap .mouse-value-main{display:block;font-weight:700;color:inherit;white-space:nowrap;}'
        '.mouse-breakdown-wrap .mouse-value-range{display:block;margin-top:0.08rem;'
        'font-size:0.78rem;font-weight:500;color:#64748b;white-space:nowrap;}'
        '.mouse-breakdown-wrap .mouse-level-0 .mouse-value-range{color:#dbeafe;}'
        '.mouse-breakdown-wrap .mouse-level-summary .mouse-value-range{color:#b45309;}'
        '.mouse-breakdown-wrap .mouse-table-empty{color:#94a3b8;}'
        '.mouse-breakdown-wrap .mouse-level-0 td{border-top:2px solid #0a2540;}'
        '.mouse-breakdown-wrap .mouse-level-1 td{border-top:1px solid #bfdbfe;}'
        '.mouse-breakdown-wrap .mouse-level-3 td,.mouse-breakdown-wrap .mouse-level-4 td{'
        'font-size:0.84rem;}'
        '</style>'
        '<div style="max-height:580px;overflow-y:auto;'
        'border:1px solid #d8e3ee;border-radius:8px;margin-bottom:0.5rem;'
        'box-shadow:0 1px 2px rgba(15,23,42,0.06);" class="mouse-breakdown-wrap">'
        '<table>'
        f'<thead><tr>{_thead_html}</tr></thead>'
        f'<tbody>{"".join(_rows_html)}</tbody>'
        '</table></div>'
    )
    st.markdown(_html_table, unsafe_allow_html=True)
    del _rows_html, _html_table, _thead_html

    st.markdown('<div style="height:1rem"></div>', unsafe_allow_html=True)

    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        display_df.to_excel(writer, index=False, sheet_name='Cost Estimate')
    buffer.seek(0)

    dl_col, _ = st.columns([1, 3])
    dl_col.download_button(
        label='⬇ Download Full Excel',
        data=buffer.getvalue(),
        file_name=f'MOUSE_cost_estimate_{reactor_type}.xlsx',
        mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
        width='stretch',
    )
    # ═══════════════════════════════════════════════════════════════
    # SECTION 05 Costs in Perspective
    # ═══════════════════════════════════════════════════════════════
    _section_header(
        '05',
        'Market &amp; Geographic Competitiveness',
        'Is it competitive?',
    )
    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:0 0 0.85rem 0;">LCOE vs Units Deployed</div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div style="background:#f7f8fa;border:1px solid #bfdbfe;border-radius:8px;'
        'padding:0.85rem 1.1rem;margin-bottom:0.9rem;font-size:0.85rem;line-height:1.45;color:#3c4257;">'
        'The curve is anchored at four deployment scales: N=1 (FOAK from the headline '
        'card), N=2, N=10 (two extra cost-engine calls, cached), and N=user-set NOAK '
        'Unit Number (default 100, NOAK from the headline card). The shaded band '
        'reflects the one-sigma uncertainty around each anchor and is connected '
        'piecewise-linearly between them. Overlaid against seven US-relevant '
        'electricity-market price ranges to indicate where the reactor would be '
        'cost-competitive at each scale.'
        '</div>',
        unsafe_allow_html=True,
    )

    # Use the headline FOAK (N=1) and NOAK (N=user_setting) values
    # as anchors, and add ONE extra cost-engine call at N=10 to
    # give an intermediate point. Same extraction path as the
    # headline cards (bottom_up_cost_estimate ->
    # cost_drivers_estimate -> transform_dataframe ->
    # _get_mean_std), so values are consistent by construction.
    # _mid_results and _N_user_pre were pre-computed above the
    # tabs (under a single spinner) so the whole page renders at
    # once. Re-bind the local names this block expects.
    _N_user = _N_user_pre
    _N_mids = _N_mids_pre

    try:
        _foak_m = float(lcoe_f)
        _noak_m = float(lcoe_n)
    except (TypeError, ValueError):
        _foak_m = _noak_m = float('nan')
    try:
        _foak_s = float(lcoe_f_std) if lcoe_f_std == lcoe_f_std else 0.0
        _noak_s = float(lcoe_n_std) if lcoe_n_std == lcoe_n_std else 0.0
    except (TypeError, ValueError):
        _foak_s = _noak_s = 0.0

    # Build the anchor list: FOAK + intermediates + NOAK, sorted
    # by N ascending. Skip any intermediate with NaN.
    _units = [1] + [_N for _N in _N_mids
                    if _N in _mid_results
                    and _mid_results[_N][0] == _mid_results[_N][0]] + [_N_user]
    _means = [_foak_m]
    _stds = [_foak_s]
    for _N in _N_mids:
        if _N in _mid_results and _mid_results[_N][0] == _mid_results[_N][0]:
            _means.append(_mid_results[_N][0])
            _stds.append(_mid_results[_N][1])
    _means.append(_noak_m)
    _stds.append(_noak_s)

    # NaN-tolerant gating: drop points where the LCOE came back
    # NaN, plot whatever valid points remain. Scatter markers +
    # mean line are always drawn (they don't need the spline).
    _u_arr = np.array(_units, dtype=float) if _units else np.array([])
    _m_arr = np.array(_means, dtype=float) if _means else np.array([])
    _s_arr = np.array(_stds, dtype=float) if _stds else np.array([])
    if _m_arr.size:
        _valid = ~np.isnan(_m_arr)
        _u_arr = _u_arr[_valid]
        _m_arr = _m_arr[_valid]
        _s_arr = _s_arr[_valid]
        _s_arr = np.nan_to_num(_s_arr, nan=0.0)

    if _u_arr.size >= 2:
        _fill_color, _edge_color = '#1B4F8C', '#0a2540'

        # Plain linear interpolation between anchor points  no spline,
        # no overshoot; np.interp gives straight line segments which is
        # the honest representation when we only have 2-3 anchors.
        try:
            _x_smooth = np.linspace(_u_arr.min(), _u_arr.max(), 300)
            _m_smooth = np.interp(_x_smooth, _u_arr, _m_arr)
            _s_smooth = np.interp(_x_smooth, _u_arr, _s_arr)
        except Exception:
            _x_smooth = _u_arr
            _m_smooth = _m_arr
            _s_smooth = _s_arr

        # Y-axis sizing  90th-percentile of (mean + std) rather than
        # max, so a single outlier point doesn't blow up the y-range
        # and squash the rest of the data. Hard cap at $800 so even
        # outliers stay readable.
        if _m_arr.size:
            _band_vals = _m_arr + _s_arr
            _band_p90 = float(np.nanpercentile(_band_vals, 90))
        else:
            _band_p90 = 0.0
        _ymax = max(410.0, _band_p90 * 1.15)
        _ymax = min(_ymax, 800.0)

        # Market benchmark ranges. Each market has a top (high) and
        # bottom (low) y-value bracketing an indicative cost range.
        # Within a category all entries share one color.
        _PREMIUM = '#92400e'   # Remote, Defense, Island & Mining
        _ALASKA = '#64748b'    # Alaska railbelt
        _US_GRID = '#15803d'   # U.S. grid
        _markets = [
            ('Remote communities', 2, 4, 400, 290, _PREMIUM, True),
            ('Defense', 8, 4, 316, 296, _PREMIUM, False),
            ('Island & Mining', 30, 12, 380, 190, _PREMIUM, False),
            ('Alaska railbelt electricity', 48, 8, 313, 182, _ALASKA, False),
            ('Alaska railbelt generation', 60, 8, 166, 62, _ALASKA, False),
            ('U.S. grid electricity', 75, 12, 270, 79, _US_GRID, False),
            ('U.S. grid generation', 88, 10, 55, 29, _US_GRID, False),
        ]
        _label_offsets = {
            'Remote communities': 8, 'Defense': 25,
            'Alaska railbelt generation': -15, 'U.S. grid generation': 20,
        }
        _market_rows = []
        for _name, _x, _w, _ys, _ye, _c, _down in _markets:
            _ys_capped = min(_ys, 400)
            if _name == 'Alaska railbelt generation':
                _ly = _ye + _label_offsets[_name]
            else:
                _ly = _ys_capped + _label_offsets.get(_name, 15)
            _market_rows.append({
                'name': _name, 'x_left': _x, 'x_right': _x + _w,
                'x_mid': _x + _w / 2.0, 'color': _c,
                'y_start': _ys_capped, 'y_end': _ye, 'label_y': _ly,
                'down_only': _down,
            })
        _market_df = pd.DataFrame(_market_rows)

        # ----- LCOE band + mean line + anchor points -----
        _band_df = pd.DataFrame({
            'units': _x_smooth,
            'lower': _m_smooth - _s_smooth,
            'upper': _m_smooth + _s_smooth,
            'mean': _m_smooth,
        })
        _anchor_df = pd.DataFrame({'units': _u_arr, 'lcoe': _m_arr})

        _x_axis = alt.Axis(
            values=[1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
            title='Number of Units Deployed',
            titleFontSize=12, titleFontWeight='bold', titleColor='#000000',
            labelFontSize=11, labelColor='#000000', labelFontWeight=500,
            domain=True, domainColor='#000000', domainWidth=1.5,
            tickColor='#000000',
            gridDash=[3, 3], gridColor='#cbd5e1',
        )
        _y_axis = alt.Axis(
            title='LCOE ($/MWh)',
            titleFontSize=12, titleFontWeight='bold', titleColor='#000000',
            labelFontSize=11, labelColor='#000000', labelFontWeight=500,
            domain=True, domainColor='#000000', domainWidth=1.5,
            tickColor='#000000',
            gridDash=[3, 3], gridColor='#cbd5e1',
        )

        _band = alt.Chart(_band_df).mark_area(
            color=_fill_color, opacity=0.45,
            stroke=_edge_color, strokeWidth=1.5,
        ).encode(
            x=alt.X('units:Q', scale=alt.Scale(domain=[1, 100]), axis=_x_axis),
            y=alt.Y('lower:Q', scale=alt.Scale(domain=[0, _ymax]),
                    axis=_y_axis),
            y2='upper:Q',
        )
        _mean_line = alt.Chart(_band_df).mark_line(
            color=_edge_color, strokeWidth=2.0,
        ).encode(
            x=alt.X('units:Q', scale=alt.Scale(domain=[1, 100])),
            y=alt.Y('mean:Q', scale=alt.Scale(domain=[0, _ymax])),
        )
        _anchors = alt.Chart(_anchor_df).mark_point(
            color=_edge_color, size=80, filled=True,
            stroke='white', strokeWidth=1.2,
        ).encode(
            x=alt.X('units:Q', scale=alt.Scale(domain=[1, 100])),
            y=alt.Y('lcoe:Q', scale=alt.Scale(domain=[0, _ymax])),
        )

        # ----- Market benchmark bars -----
        # Top bar (skip rows where down_only=True)
        _top_df = _market_df[~_market_df['down_only']]
        _top = alt.Chart(_top_df).mark_rule(
            strokeWidth=4, strokeCap='round',
        ).encode(
            x=alt.X('x_left:Q', scale=alt.Scale(domain=[1, 100])),
            x2='x_right:Q',
            y=alt.Y('y_start:Q', scale=alt.Scale(domain=[0, _ymax])),
            color=alt.Color('color:N', scale=None, legend=None),
        )
        _bottom = alt.Chart(_market_df).mark_rule(
            strokeWidth=4, strokeCap='round',
        ).encode(
            x=alt.X('x_left:Q', scale=alt.Scale(domain=[1, 100])),
            x2='x_right:Q',
            y=alt.Y('y_end:Q', scale=alt.Scale(domain=[0, _ymax])),
            color=alt.Color('color:N', scale=None, legend=None),
        )
        _vert = alt.Chart(_market_df).mark_rule(
            strokeWidth=2.4, strokeCap='round',
        ).encode(
            x=alt.X('x_mid:Q', scale=alt.Scale(domain=[1, 100])),
            y=alt.Y('y_start:Q', scale=alt.Scale(domain=[0, _ymax])),
            y2='y_end:Q',
            color=alt.Color('color:N', scale=None, legend=None),
        )
        _market_labels = alt.Chart(_market_df).mark_text(
            fontSize=9, fontWeight='bold', baseline='middle', align='center',
        ).encode(
            x=alt.X('x_mid:Q', scale=alt.Scale(domain=[1, 100])),
            y=alt.Y('label_y:Q', scale=alt.Scale(domain=[0, _ymax])),
            text='name:N',
            color=alt.Color('color:N', scale=None, legend=None),
        )

        # If the LCOE band sits entirely above the y-axis cap, warn
        # the user before rendering the (invisible) curve.
        _band_min_visible = (_m_arr - _s_arr).min() if _m_arr.size else 0.0
        if _band_min_visible > _ymax:
            _band_lo = (_m_arr - _s_arr).min()
            _band_hi = (_m_arr + _s_arr).max()
            st.markdown(
                f'<div style="background:#fffbeb;border:1px solid #f59e0b;'
                f'border-radius:8px;padding:0.85rem 1.1rem;margin-bottom:0.6rem;'
                f'color:#92400e;font-size:0.85rem;line-height:1.55;">'
                f'<strong style="color:#92400e;">⚠️ Reactor LCOE off the chart.</strong> '
                f'The reactor LCOE band ranges roughly '
                f'<strong>${_band_lo:,.0f}-${_band_hi:,.0f}/MWh</strong>, which is above '
                f'the chart\'s <strong>${int(_ymax)}/MWh</strong> ceiling, so the curve is '
                f'not visible on the plot below. The market benchmarks remain visible for '
                f'reference. To bring the curve into the chart, try a higher reactor power, '
                f'higher enrichment, longer plant lifetime, or a larger NOAK Unit Number '
                f'any of those reduces the LCOE.'
                f'</div>',
                unsafe_allow_html=True,
            )

        _lcoe_chart = (
            _band + _mean_line + _anchors + _top + _bottom + _vert + _market_labels
        ).properties(height=420).configure_view(stroke=None)
        st.altair_chart(_lcoe_chart, use_container_width=True)

        # Market definitions panel (matches the user-provided spec)
        st.markdown(
            '<div style="background:#eff6ff;border:1px solid #bfdbfe;'
            'border-radius:8px;padding:0.85rem 1.1rem;margin-bottom:0.9rem;'
            'font-size:0.85rem;line-height:1.55;color:#1B4F8C;">'
            '<div style="font-weight:600;font-size:0.85rem;'
            'text-transform:uppercase;letter-spacing:0.06em;'
            'color:#1B4F8C;margin-bottom:0.45rem;">Market definitions</div>'
            '<ul style="margin:0;padding-left:1.2rem;">'
            '<li><strong>U.S. Grid Generation:</strong> regional average wholesale price; '
            'excludes transmission, distribution, and customer charges.</li>'
            '<li><strong>U.S. Grid Electricity:</strong> state-level average retail '
            'electricity price, all sectors; excludes Alaska and Hawaii.</li>'
            '<li><strong>Alaska Railbelt Generation:</strong> wholesale generation cost; '
            'excludes transmission, distribution, and customer charges.</li>'
            '<li><strong>Alaska Railbelt Electricity:</strong> retail price including '
            'generation, transmission, distribution, and adjustments.</li>'
            '<li><strong>Island &amp; Mining:</strong> all-in diesel- or LNG-based '
            'electricity cost, including fuel delivery.</li>'
            '<li><strong>Defense:</strong> remote base electricity cost plus a premium '
            'for reliability and security.</li>'
            '<li><strong>Remote Communities:</strong> off-road community diesel '
            'generation cost.</li>'
            '</ul>'
            '</div>',
            unsafe_allow_html=True,
        )

    # ─────────────────────────────────────────────────────────
    # Behind-the-meter / Distributed Generation comparison:
    # within-state industrial retail price distribution vs LCOE
    # ─────────────────────────────────────────────────────────
    # Data source: EIA Form 861 "Sales to Ultimate Customers",
    # 2024 final data. Per-utility industrial revenue / sales
    # gives an effective $/kWh for every utility in a state;
    # aggregating yields a true within-state distribution
    # rather than a single state-average number. Filters:
    # drop Utility Number 99999 (EIA imputation), keep
    # Service Type == 'Bundled' so generation + T&D are
    # bundled (apples-to-apples with reactor LCOE), drop rows
    # with non-positive industrial revenue or sales.
    # Precomputed by:
    # cost_drivers_and_markets/build_state_industrial_price_distribution.py
    _STATE_FULL_NAMES = {
        'AL': 'Alabama', 'AK': 'Alaska', 'AZ': 'Arizona',
        'AR': 'Arkansas', 'CA': 'California', 'CO': 'Colorado',
        'CT': 'Connecticut', 'DE': 'Delaware', 'DC': 'D.C.',
        'FL': 'Florida', 'GA': 'Georgia', 'HI': 'Hawaii',
        'ID': 'Idaho', 'IL': 'Illinois', 'IN': 'Indiana',
        'IA': 'Iowa', 'KS': 'Kansas', 'KY': 'Kentucky',
        'LA': 'Louisiana', 'ME': 'Maine', 'MD': 'Maryland',
        'MA': 'Massachusetts', 'MI': 'Michigan', 'MN': 'Minnesota',
        'MS': 'Mississippi', 'MO': 'Missouri', 'MT': 'Montana',
        'NE': 'Nebraska', 'NV': 'Nevada', 'NH': 'New Hampshire',
        'NJ': 'New Jersey', 'NM': 'New Mexico', 'NY': 'New York',
        'NC': 'North Carolina', 'ND': 'North Dakota', 'OH': 'Ohio',
        'OK': 'Oklahoma', 'OR': 'Oregon', 'PA': 'Pennsylvania',
        'RI': 'Rhode Island', 'SC': 'South Carolina', 'SD': 'South Dakota',
        'TN': 'Tennessee', 'TX': 'Texas', 'UT': 'Utah',
        'VT': 'Vermont', 'VA': 'Virginia', 'WA': 'Washington',
        'WV': 'West Virginia', 'WI': 'Wisconsin', 'WY': 'Wyoming',
    }

    st.markdown('<div style="height:1.2rem"></div>', unsafe_allow_html=True)
    st.markdown(
        '<div style="font-size:1rem;font-weight:700;color:#0a2540;'
        'border-left:4px solid #0a2540;padding:0.4rem 0 0.4rem 0.75rem;'
        'margin:1.25rem 0 0.85rem 0;">Behind-the-Meter Comparison: '
        'Reactor LCOE vs State Retail Electricity Prices</div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div style="background:#f7f8fa;border:1px solid #bfdbfe;border-radius:8px;'
        'padding:0.85rem 1.1rem;margin-bottom:0.9rem;font-size:0.85rem;line-height:1.55;color:#3c4257;">'

        '<p style="margin:0 0 0.55rem 0;">'
        '<strong>What this shows.</strong> '
        'For behind-the-meter (BTM) and bring-your-own-generator (BYOG) deployments '
        'data centers, industrial loads, mining, defense a useful '
        '<em>first-order</em> cost benchmark is the retail electricity price the '
        'customer would otherwise pay.'
        '</p>'

        '<p style="margin:0 0 0.4rem 0;"><strong>How to read it.</strong></p>'
        '<ul style="margin:0 0 0.55rem 1.2rem;padding:0;">'
        '<li><strong>Bar height</strong> = sales-weighted mean industrial retail '
        'price in that state (the EIA Table 5.6.A state-average definition).</li>'
        '<li><strong>Whisker</strong> = 10th to 90th percentile of individual '
        'utility prices in that state, showing how much variation exists within '
        'the state.</li>'
        '<li><strong>Horizontal bands</strong> = reactor FOAK and NOAK LCOE '
        '±1&sigma; ranges.</li>'
        '</ul>'

        '<p style="margin:0 0 0.4rem 0;"><strong>Bar color (competitiveness of '
        'the sales-weighted mean vs the LCOE bands).</strong></p>'
        '<ul style="margin:0 0 0.55rem 1.2rem;padding:0;">'
        '<li><span style="color:#15803d;font-weight:600;">Green:</span> state-average '
        'retail price exceeds the FOAK lower band reactor wins even at FOAK.</li>'
        '<li><span style="color:#92400e;font-weight:600;">Yellow:</span> state-average '
        'price between NOAK and FOAK bands reactor wins only at NOAK scale.</li>'
        '<li><span style="color:#64748b;font-weight:600;">Gray:</span> state-average '
        'price below the NOAK lower band not competitive even at scale.</li>'
        '</ul>'

        '<p style="margin:0 0 0.55rem 0;border-left:2px solid #cbd5e1;'
        'padding:0 0 0 0.6rem;font-size:0.85rem;color:#3c4257;">'
        '<strong>Note.</strong> The bar top is usually not in the middle of the '
        'whisker. The bar (sales-weighted mean) and the whisker (10th-90th '
        'percentile of utility prices) answer different questions. Each utility '
        'counts equally in the percentiles, but the weighted mean is dominated '
        'by whichever utilities sell the most MWh. In most states a few big '
        'investor-owned utilities sell the bulk of industrial power at '
        'lower-than-typical rates, so the weighted mean sits in the lower part '
        'of the whisker.'
        '</p>'

        '<p style="margin:0 0 0.55rem 0;font-size:0.8rem;color:#64748b;line-height:1.4;">'
        '<strong>Caveat.</strong> First-order benchmark only; actual avoided cost '
        'depends on demand charges, time-of-use structure, standby tariffs, and '
        'PPAs. D.C. is omitted (no bundled industrial utilities in 2024).'
        '</p>'

        '<ul style="margin:0 0 0 1.2rem;padding:0;font-size:0.85rem;color:#64748b;">'
        '<li>Source: U.S. EIA Form 861, 2024 final data, '
        '<code>Sales_Ult_Cust_2024.xlsx</code> (sheet <code>States</code>). '
        '<a href="https://www.eia.gov/electricity/data/eia861/" target="_blank" '
        'style="color:#1B4F8C;">eia.gov/electricity/data/eia861</a>.</li>'
        '<li>Per-utility price = INDUSTRIAL Revenues / Sales. '
        'Per state: weighted mean = &Sigma;(revenue) / &Sigma;(sales); '
        'p10 / p90 are unweighted percentiles of the per-utility prices.</li>'
        '</ul>'

        '</div>',
        unsafe_allow_html=True,
    )

    try:
        _foak_m_btm = float(lcoe_f)
        _noak_m_btm = float(lcoe_n)
        _foak_s_btm = float(lcoe_f_std) if lcoe_f_std == lcoe_f_std else 0.0
        _noak_s_btm = float(lcoe_n_std) if lcoe_n_std == lcoe_n_std else 0.0
    except (TypeError, ValueError):
        _foak_m_btm = _noak_m_btm = float('nan')
        _foak_s_btm = _noak_s_btm = 0.0

    if (_foak_m_btm == _foak_m_btm and _noak_m_btm == _noak_m_btm):
        # Load per-state distribution and convert cents/kWh -> $/MWh
        # (×10) for every percentile / mean column so chart units match
        # the reactor LCOE.
        _dist_df = _load_state_industrial_price_distribution().copy()
        _PRICE_COLS = ['weighted_mean', 'min', 'p10', 'q1', 'median',
                       'q3', 'p90', 'max']
        for _c in _PRICE_COLS:
            _dist_df[_c] = _dist_df[_c] * 10.0
        _dist_df['state_name'] = _dist_df['state'].map(_STATE_FULL_NAMES)
        _dist_df = _dist_df.sort_values('weighted_mean', ascending=False).reset_index(drop=True)

        # Competitiveness thresholds: use the LOWER edge of each
        # uncertainty band so any state whose sales-weighted mean
        # falls inside the FOAK or NOAK band counts as competitive.
        # Sales-weighted mean is used (not median) because it is the
        # direct analog of the EIA Table 5.6.A state-average number
        # the chart used to display as a single bar.
        _foak_low = _foak_m_btm - _foak_s_btm
        _noak_low = _noak_m_btm - _noak_s_btm

        _GREEN = '#16a34a'
        _YELLOW = '#d97706'
        _GRAY = '#94a3b8'

        def _tier(_v):
            if _v >= _foak_low:
                return 'Beats FOAK'
            if _v >= _noak_low:
                return 'Beats NOAK'
            return 'Loses'

        _dist_df['tier'] = _dist_df['weighted_mean'].apply(_tier)

        # Composed box plot per state:
        #   1) ±1σ FOAK/NOAK background bands (mark_rect)
        #   2) Whisker rule from p10 to p90 (tier color, slim)
        #   3) Box from q1 to q3 (tier color)
        #   4) Median line (mark_tick, dark, slightly thicker)
        #   5) Sales-weighted mean marker (diamond)
        #   6) Dashed FOAK / NOAK mean reference lines
        #   7) State labels rotated at the base of the chart
        _band_top_btm = max(
            _foak_m_btm + _foak_s_btm,
            _noak_m_btm + _noak_s_btm,
            float(_dist_df['p90'].max()),
            float(_dist_df['weighted_mean'].max()),
        )
        _ymax_chart = _band_top_btm * 1.18

        _tier_scale = alt.Scale(
            domain=['Beats FOAK', 'Beats NOAK', 'Loses'],
            range=[_GREEN, _YELLOW, _GRAY],
        )

        _bands_df = pd.DataFrame([
            {'low': _foak_m_btm - _foak_s_btm,
             'high': _foak_m_btm + _foak_s_btm, 'band': 'FOAK'},
            {'low': _noak_m_btm - _noak_s_btm,
             'high': _noak_m_btm + _noak_s_btm, 'band': 'NOAK'},
        ])
        _bands = alt.Chart(_bands_df).mark_rect(opacity=0.16).encode(
            y='low:Q', y2='high:Q',
            color=alt.Color(
                'band:N',
                scale=alt.Scale(domain=['FOAK', 'NOAK'],
                                range=['#c84b1e', '#1B4F8C']),
                legend=None,
            ),
        )

        _lines_df = pd.DataFrame([
            {'mean': _foak_m_btm,
             'label': f'FOAK LCOE ${_foak_m_btm:.0f}/MWh'},
            {'mean': _noak_m_btm,
             'label': f'NOAK LCOE ${_noak_m_btm:.0f}/MWh'},
        ])
        _lines = alt.Chart(_lines_df).mark_rule(
            strokeDash=[6, 4], strokeWidth=1.8,
        ).encode(
            y='mean:Q',
            color=alt.Color(
                'label:N',
                scale=alt.Scale(
                    domain=[f'FOAK LCOE ${_foak_m_btm:.0f}/MWh',
                            f'NOAK LCOE ${_noak_m_btm:.0f}/MWh'],
                    range=['#c84b1e', '#1B4F8C'],
                ),
                legend=alt.Legend(title=None, orient='top-right'),
            ),
        )

        _state_order = list(_dist_df['state'])
        _x_state = alt.X(
            'state:N',
            sort=_state_order,
            title='U.S. states (sorted by sales-weighted mean, high to low)',
            axis=alt.Axis(
                labels=False, ticks=False,
                titleFontSize=13, titleFontWeight='bold',
                titleColor='#000000',
                domain=True, domainColor='#000000', domainWidth=1.5,
            ),
        )
        _y_axis = alt.Axis(
            labelFontSize=12, labelColor='#000000', labelFontWeight=500,
            titleFontSize=13, titleFontWeight='bold', titleColor='#000000',
            domain=True, domainColor='#000000', domainWidth=1.5,
            tickColor='#000000', tickWidth=1.2,
            gridDash=[3, 3], gridColor='#cbd5e1',
        )
        _y_scale = alt.Scale(domain=[0, _ymax_chart])

        _tooltip = [
            alt.Tooltip('state_name:N', title='State'),
            alt.Tooltip('n_utilities:Q', title='Utilities (n)'),
            alt.Tooltip('weighted_mean:Q', title='Weighted mean $/MWh',
                        format='.1f'),
            alt.Tooltip('median:Q', title='Median utility $/MWh', format='.1f'),
            alt.Tooltip('p10:Q', title='p10 $/MWh', format='.1f'),
            alt.Tooltip('p90:Q', title='p90 $/MWh', format='.1f'),
        ]

        _bars = alt.Chart(_dist_df).mark_bar().encode(
            x=_x_state,
            y=alt.Y('weighted_mean:Q',
                    title='Industrial retail price ($/MWh, EIA 861, 2024)',
                    scale=_y_scale, axis=_y_axis),
            color=alt.Color('tier:N', scale=_tier_scale, legend=None),
            tooltip=_tooltip,
        )
        _whiskers = alt.Chart(_dist_df).mark_rule(
            color='#0a2540', strokeWidth=1.2,
        ).encode(
            x=_x_state,
            y=alt.Y('p10:Q', scale=_y_scale, axis=_y_axis),
            y2='p90:Q',
            tooltip=_tooltip,
        )
        _whisker_caps = alt.Chart(_dist_df).transform_fold(
            ['p10', 'p90'], as_=['percentile', 'value'],
        ).mark_tick(color='#0a2540', thickness=1.2, size=8).encode(
            x=_x_state,
            y=alt.Y('value:Q', scale=_y_scale, axis=_y_axis),
            tooltip=_tooltip,
        )

        _label_df = pd.DataFrame({
            'state': _state_order,
            'state_name': [_STATE_FULL_NAMES.get(c, c) for c in _state_order],
            'y_anchor': [_ymax_chart * 0.012] * len(_state_order),
        })
        _labels = alt.Chart(_label_df).mark_text(
            angle=270, baseline='middle', align='left',
            fontSize=11, color='#000000',
        ).encode(
            x=alt.X('state:N', sort=_state_order),
            y='y_anchor:Q',
            text='state_name:N',
        )

        _state_chart = (
            _bands + _bars + _whiskers + _whisker_caps + _lines + _labels
        ).properties(height=460).resolve_scale(color='independent')
        st.altair_chart(_state_chart, use_container_width=True)

        # Distribution-aware competitiveness counts. For each LCOE tier
        # (FOAK, NOAK) we ask three questions of the within-state utility
        # price distribution:
        #   beats every utility    = state p10    >= LCOE_low
        #     (even the cheapest utility's price is at or above the
        #      reactor's lower band: reactor wins regardless of which
        #      utility serves the site)
        #   beats the median       = state median >= LCOE_low
        #     (reactor wins against more than half the utilities in
        #      the state)
        #   beats some utilities   = state p90    >= LCOE_low
        #     (reactor wins only against the priciest utilities the
        #      expensive end of the state's distribution: competitive
        #      in pockets of the state, not state-wide)
        # NOAK lists are inclusive of FOAK (any FOAK-competitive state
        # is automatically NOAK-competitive since NOAK_low < FOAK_low).
        _foak_beats_p10    = _dist_df.loc[_dist_df['p10']    >= _foak_low, 'state'].tolist()
        _foak_beats_median = _dist_df.loc[_dist_df['median'] >= _foak_low, 'state'].tolist()
        _foak_beats_p90    = _dist_df.loc[_dist_df['p90']    >= _foak_low, 'state'].tolist()
        _noak_beats_p10    = _dist_df.loc[_dist_df['p10']    >= _noak_low, 'state'].tolist()
        _noak_beats_median = _dist_df.loc[_dist_df['median'] >= _noak_low, 'state'].tolist()
        _noak_beats_p90    = _dist_df.loc[_dist_df['p90']    >= _noak_low, 'state'].tolist()

        def _state_list_html(states):
            return ', '.join(states) if states else '<em>none</em>'

        st.markdown(
            f'<div style="background:#f7f8fa;border:1px solid #bfdbfe;'
            f'border-radius:8px;padding:0.85rem 1.1rem;margin-bottom:0.9rem;'
            f'font-size:0.85rem;line-height:1.55;color:#3c4257;">'
            f'Out of <strong>{len(_dist_df)}</strong> states, the reactor\'s '
            f'competitiveness depends on which utility serves the site:'
            f'<div style="margin-top:0.55rem;">'
            f'<span style="display:inline-block;width:0.85rem;height:0.85rem;'
            f'background:{_GREEN};border-radius:3px;margin-right:0.35rem;'
            f'vertical-align:middle;"></span>'
            f'<strong>FOAK competitive</strong>'
            f'<div style="margin-left:1.2rem;margin-top:0.2rem;">'
            f'More cost-competitive than <em>every</em> utility in '
            f'<strong>{len(_foak_beats_p10)}</strong> states: '
            f'{_state_list_html(_foak_beats_p10)}'
            f'</div>'
            f'<div style="margin-left:1.2rem;margin-top:0.2rem;">'
            f'More cost-competitive than the <em>median</em> utility in '
            f'<strong>{len(_foak_beats_median)}</strong> states: '
            f'{_state_list_html(_foak_beats_median)}'
            f'</div>'
            f'<div style="margin-left:1.2rem;margin-top:0.2rem;">'
            f'More cost-competitive than <em>at least the priciest</em> utilities in '
            f'<strong>{len(_foak_beats_p90)}</strong> states (competitive in '
            f'pockets only): {_state_list_html(_foak_beats_p90)}'
            f'</div>'
            f'</div>'
            f'<div style="margin-top:0.55rem;">'
            f'<span style="display:inline-block;width:0.85rem;height:0.85rem;'
            f'background:{_YELLOW};border-radius:3px;margin-right:0.35rem;'
            f'vertical-align:middle;"></span>'
            f'<strong>NOAK competitive</strong>'
            f'<div style="margin-left:1.2rem;margin-top:0.2rem;">'
            f'More cost-competitive than <em>every</em> utility in '
            f'<strong>{len(_noak_beats_p10)}</strong> states: '
            f'{_state_list_html(_noak_beats_p10)}'
            f'</div>'
            f'<div style="margin-left:1.2rem;margin-top:0.2rem;">'
            f'More cost-competitive than the <em>median</em> utility in '
            f'<strong>{len(_noak_beats_median)}</strong> states: '
            f'{_state_list_html(_noak_beats_median)}'
            f'</div>'
            f'<div style="margin-left:1.2rem;margin-top:0.2rem;">'
            f'More cost-competitive than <em>at least the priciest</em> utilities in '
            f'<strong>{len(_noak_beats_p90)}</strong> states (competitive in '
            f'pockets only): {_state_list_html(_noak_beats_p90)}'
            f'</div>'
            f'</div>'
            f'</div>',
            unsafe_allow_html=True,
        )


    # ── End-of-report footer (gives users a clear "I'm done" signal
    #    now that everything is on a single page instead of in tabs) ──
    st.markdown(
        '<div style="margin-top:2rem;padding:0.8rem 0;border-top:2px solid #cbd5e1;'
        'border-bottom:2px solid #cbd5e1;text-align:center;color:#64748b;'
        'font-size:0.85rem;letter-spacing:0.16em;text-transform:uppercase;">'
        '— End of report —</div>',
        unsafe_allow_html=True,
    )

    # ── Caveat reminder (rendered at the very end of the results page) ──────
    st.markdown(
        '''<div style="background:#fffbeb;border:1px solid #fcd34d;
                       border-left:4px solid #f59e0b;border-radius:8px;
                       padding:0.85rem 1.1rem;margin-top:1.5rem;margin-bottom:1.2rem;
                       display:flex;align-items:flex-start;gap:0.6rem;">
             <span style="font-size:1rem;flex-shrink:0;margin-top:0.05rem;">⚠️</span>
             <span style="font-size:1rem;color:#92400e;line-height:1.55;">
               <strong>Caveats:</strong> Reactor designs are pre-conceptual and not fully
               optimized. Cost estimates were produced with incomplete information and
               <strong>must not be used for investment or procurement decisions</strong>.
               Results are intended for research and comparative screening only.
             </span>
           </div>''',
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div style='text-align: center; font-size: 1rem; color: #64748b; padding-top: 2rem; padding-bottom: 1rem;'>
            © 2025 Battelle Energy Alliance, LLC. MOUSE is released under the MIT License.
        </div>
        """,
        unsafe_allow_html=True,
    )
