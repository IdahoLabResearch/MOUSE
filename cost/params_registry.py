# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
Parameter Registry for MOUSE output reporting.

Each entry maps a parameter name to its metadata:
  - group       : Section it belongs to in the Parameters sheet
  - units       : Physical units (empty string if unitless)
  - description : Human-readable description of the parameter
  - source      : 'User Input' or 'Calculated'
  - hidden      : If True, the parameter is excluded from the Parameters sheet
  - array_mode  : How to display list/array values:
                    'summary'  → show BOL, EOL, min, max as separate rows
                    'steps'    → show first step, last step, number of steps
                    'as_is'    → show the value directly (short lists)
                    None       → not an array

Parameters not in this registry are placed in the 'Uncategorized' group
with no units or description and a warning comment.
"""

# Order in which groups appear in the output sheet
GROUP_ORDER = [
    'Settings',
    'Materials',
    'Geometry',
    'Control Drums',
    'Core Design',
    'Overall System',
    'OpenMC Settings',
    'Physics Results',
    'Primary Loop & Balance of Plant',
    'Shielding',
    'Vessels',
    'Operation',
    'Economic Parameters',
    'Central Facility',
    'Learning Rates',
    'Tax Credits',
    'Debug / Intermediate Values',
    'Uncategorized',
]

PARAMS_REGISTRY = {

    # =========================================================
    # Settings
    # =========================================================
    'plotting': {
        'group': 'Settings', 'units': '',
        'description': 'Whether to generate and save geometry plots (Y = yes, N = no)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'cross_sections_xml_location': {
        'group': 'Settings', 'units': '',
        'description': 'File path to the OpenMC cross sections XML file on HPC',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'simplified_chain_thermal_xml': {
        'group': 'Settings', 'units': '',
        'description': 'File path to the simplified depletion chain XML file for thermal reactors',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Materials
    # =========================================================
    'reactor type': {
        'group': 'Materials', 'units': '',
        'description': 'Reactor type identifier: LTMR (liquid metal thermal), GCMR (gas-cooled), or HPMR (heat pipe)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'TRISO Fueled': {
        'group': 'Materials', 'units': '',
        'description': 'Whether the reactor uses TRISO fuel particles (Yes or No)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Fuel': {
        'group': 'Materials', 'units': '',
        'description': 'Primary fuel material type (e.g. TRIGA_fuel, UO2, UN, UC, homog_TRISO)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Enrichment': {
        'group': 'Materials', 'units': 'fraction',
        'description': 'Uranium-235 enrichment fraction (0 to 1)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'H_Zr_ratio': {
        'group': 'Materials', 'units': 'atomic ratio',
        'description': 'Ratio of hydrogen atoms to zirconium atoms in TRIGA fuel (ZrH hydride)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'U_met_wo': {
        'group': 'Materials', 'units': 'weight fraction',
        'description': 'Weight fraction of uranium metal in TRIGA fuel (remainder is ZrH matrix)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'er_wo': {
        'group': 'Materials', 'units': 'weight fraction',
        'description': 'Weight fraction of erbium burnable poison in TRIGA fuel (0 means no erbium)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'UO2 atom fraction': {
        'group': 'Materials', 'units': 'fraction',
        'description': 'Atom fraction of UO2 in the UCO mixed fuel (remainder is UC)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Coolant': {
        'group': 'Materials', 'units': '',
        'description': 'Primary coolant material type (e.g. NaK, Helium)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Radial Reflector': {
        'group': 'Materials', 'units': '',
        'description': 'Radial reflector material type (e.g. Graphite, BeO)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Axial Reflector': {
        'group': 'Materials', 'units': '',
        'description': 'Axial reflector material type — may differ from radial reflector',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Moderator': {
        'group': 'Materials', 'units': '',
        'description': 'Neutron moderator material type (e.g. ZrH, Graphite, monolith_graphite)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Moderator Booster Materials': {
        'group': 'Materials', 'units': '',
        'description': 'List of moderator booster material names, from innermost to outermost region (e.g. [\'ZrH\'] or [\'ZrH\', \'Graphite\'])',
        'source': 'User Input', 'hidden': False, 'array_mode': 'list'},

    'Matrix Material': {
        'group': 'Materials', 'units': '',
        'description': 'Background matrix material surrounding TRISO particles in the fuel compact',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Control Drum Absorber': {
        'group': 'Materials', 'units': '',
        'description': 'Neutron absorber material in the control drums (e.g. B4C_enriched, B4C_natural)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Control Drum Reflector': {
        'group': 'Materials', 'units': '',
        'description': 'Reflector material in the control drums (e.g. Graphite)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Common Temperature': {
        'group': 'Materials', 'units': 'K',
        'description': 'Uniform operating temperature applied to all materials in isothermal simulations',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'HX Material': {
        'group': 'Materials', 'units': '',
        'description': 'Heat exchanger structural material (e.g. SS316)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Secondary Coolant': {
        'group': 'Materials', 'units': '',
        'description': 'Secondary coolant filling the gap between fuel and moderator or heat pipe and moderator (HPMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Cooling Device': {
        'group': 'Materials', 'units': '',
        'description': 'Primary heat removal device type (e.g. heatpipe for HPMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Geometry
    # =========================================================
    'Fuel Pin Materials': {
        'group': 'Geometry', 'units': '',
        'description': 'Ordered list of materials for each fuel pin region from innermost to outermost',
        'source': 'User Input', 'hidden': False, 'array_mode': 'as_is'},

    'Fuel Pin Radii': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Ordered list of radii defining fuel pin region boundaries from innermost to outermost',
        'source': 'User Input', 'hidden': False, 'array_mode': 'as_is'},

    'Heat Pipe Materials': {
        'group': 'Geometry', 'units': '',
        'description': 'Ordered list of materials for each heat pipe region (HPMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': 'as_is'},

    'Heat Pipe Radii': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Ordered list of radii defining heat pipe region boundaries (HPMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': 'as_is'},

    'Moderator Pin Materials': {
        'group': 'Geometry', 'units': '',
        'description': 'Ordered list of materials for each moderator pin region from innermost to outermost',
        'source': 'User Input', 'hidden': False, 'array_mode': 'as_is'},

    'Moderator Pin Radii': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Ordered list of radii defining moderator pin region boundaries',
        'source': 'User Input', 'hidden': False, 'array_mode': 'as_is'},

    'Moderator Pin Inner Radius': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Inner radius of the moderator pin cladding',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Pin Gap Distance': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Gap distance between adjacent pins (center-to-center minus two outer radii)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Pins Arrangement': {
        'group': 'Geometry', 'units': '',
        'description': '2D pin arrangement map (hidden for readability — see input file)',
        'source': 'User Input', 'hidden': True, 'array_mode': None},

    'Number of Rings per Assembly': {
        'group': 'Geometry', 'units': '',
        'description': 'Number of pin rings along the side of the hexagonal assembly (LTMR/HPMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Number of Rings per Core': {
        'group': 'Geometry', 'units': '',
        'description': 'Number of assembly rings along the side of the hexagonal core (HPMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Radial Reflector Thickness': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Thickness of the radial neutron reflector surrounding the active core',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Axial Reflector Thickness': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Thickness of each axial neutron reflector (top and bottom)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Active Height': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Axial height of the active fuel region',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Compact Fuel Radius': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Radius of the cylindrical fuel compact region containing TRISO particles',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Packing Fraction': {
        'group': 'Geometry', 'units': 'fraction',
        'description': 'Volume fraction of TRISO particles within the fuel compact',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Coolant Channel Radius': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Radius of the coolant channels in the GCMR assembly',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Moderator Booster Radii': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'List of cumulative outer radii for moderator booster pin regions, innermost first (e.g. [0.55] or [0.40, 0.55])',
        'source': 'User Input', 'hidden': False, 'array_mode': 'list'},

    'Lattice Pitch': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Center-to-center distance between adjacent lattice positions',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Assembly Rings': {
        'group': 'Geometry', 'units': '',
        'description': 'Number of hexagonal lattice rings in one fuel assembly (GCMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Core Rings': {
        'group': 'Geometry', 'units': '',
        'description': 'Number of assembly rings along the side of the hexagonal core (GCMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Lattice Apothem': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Distance from the hex lattice center to the midpoint of a flat face of the outermost pin lattice.',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Lattice Radius': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Legacy LTMR name for the hex lattice apothem: distance from the hex lattice center to the midpoint of a flat face.',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Assembly FTF': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Assembly flat-to-flat distance (distance between opposite parallel faces of the hexagonal assembly).',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Core Radius': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Total core radius including the radial reflector',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'hexagonal Core Edge Length': {
        'group': 'Geometry', 'units': 'cm',
        'description': 'Edge length of the hexagonal core (HPMR)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Fuel Pin Count': {
        'group': 'Geometry', 'units': '',
        'description': 'Total number of fuel pins (or fuel elements) in the entire core',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Moderator Pin Count': {
        'group': 'Geometry', 'units': '',
        'description': 'Total number of moderator pins in the entire core',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Fuel Pin Count per Assembly': {
        'group': 'Geometry', 'units': '',
        'description': 'Number of fuel pins per hexagonal assembly (HPMR)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Fuel Assemblies Count': {
        'group': 'Geometry', 'units': '',
        'description': 'Total number of fuel assemblies in the core (HPMR)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Number of Heatpipes': {
        'group': 'Geometry', 'units': '',
        'description': 'Total number of heat pipes in the entire core (HPMR)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Number of Heatpipes per Assembly': {
        'group': 'Geometry', 'units': '',
        'description': 'Number of heat pipes per hexagonal assembly (HPMR)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Control Drums
    # =========================================================
    'Number of Drums': {
        'group': 'Control Drums', 'units': '',
        'description': 'Number of control drums in the core (e.g. 6, 12)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Drum Absorber Arc Degrees': {
        'group': 'Control Drums', 'units': 'degrees',
        'description': 'Angular span of the absorber arc on each control drum',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Drum Radius': {
        'group': 'Control Drums', 'units': 'cm',
        'description': 'Outer radius of each control drum. If not specified, it is automatically set to the largest allowable value that avoids drum overlap.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Drum Absorber Thickness': {
        'group': 'Control Drums', 'units': 'cm',
        'description': 'Thickness of the absorber arc section within each control drum',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Drum Height': {
        'group': 'Control Drums', 'units': 'cm',
        'description': 'Height of the control drums (typically active height + 2x axial reflector thickness)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Drum Tube Radius': {
        'group': 'Control Drums', 'units': 'cm',
        'description': 'Outer radius of the control drum tube (slightly larger than drum radius)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Drum Count': {
        'group': 'Control Drums', 'units': '',
        'description': 'Total number of control drums in the core',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'All Drums Volume': {
        'group': 'Control Drums', 'units': 'cm³',
        'description': 'Total volume of all control drums combined',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'All Drums Area': {
        'group': 'Control Drums', 'units': 'cm²',
        'description': 'Total cross-sectional area of all control drums',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Control Drum Absorber Mass': {
        'group': 'Control Drums', 'units': 'kg',
        'description': 'Total mass of absorber material in all control drums',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Control Drum Reflector Mass': {
        'group': 'Control Drums', 'units': 'kg',
        'description': 'Total mass of reflector material in all control drums',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Control Drums Mass': {
        'group': 'Control Drums', 'units': 'kg',
        'description': 'Total combined mass of all control drums (absorber + reflector)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Core Design
    # =========================================================
    'Radial Reflector Mass': {
        'group': 'Core Design', 'units': 'kg',
        'description': 'Total mass of the radial neutron reflector',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Axial Reflector Mass': {
        'group': 'Core Design', 'units': 'kg',
        'description': 'Total mass of both axial reflectors (top and bottom combined)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Reflector Mass': {
        'group': 'Core Design', 'units': 'kg',
        'description': 'Total reflector mass (used in some reactor types as combined radial mass)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Moderator Mass': {
        'group': 'Core Design', 'units': 'kg',
        'description': 'Total mass of the neutron moderator material in the core',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Moderator Booster Mass': {
        'group': 'Core Design', 'units': 'kg',
        'description': 'Total mass of moderator booster material in the core (GCMR)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Moderator Total Area': {
        'group': 'Core Design', 'units': 'cm²',
        'description': 'Total cross-sectional area occupied by moderator material in the core',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Number Of TRISO Particles Per Compact Fuel': {
        'group': 'Core Design', 'units': '',
        'description': 'Number of TRISO particles per fuel compact element',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Total Number of TRISO Particles': {
        'group': 'Core Design', 'units': '',
        'description': 'Total number of TRISO particles in the entire core',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Lattice Compact Volume': {
        'group': 'Core Design', 'units': 'cm³',
        'description': 'Volume of the fuel compact cylinder used in the TRISO lattice model',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Overall System
    # =========================================================
    'Power MWt': {
        'group': 'Overall System', 'units': 'MWt',
        'description': 'Reactor thermal power output',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Power MWe': {
        'group': 'Overall System', 'units': 'MWe',
        'description': 'Reactor net electric power output (thermal power × thermal efficiency)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Power kWe': {
        'group': 'Overall System', 'units': 'kWe',
        'description': 'Reactor net electric power output in kilowatts',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Thermal Efficiency': {
        'group': 'Overall System', 'units': 'fraction',
        'description': 'Net thermal-to-electric conversion efficiency of the power conversion system',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Heat Flux Criteria': {
        'group': 'Overall System', 'units': 'MW/m²',
        'description': 'Maximum allowable fuel surface heat flux (design constraint)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Heat Flux': {
        'group': 'Overall System', 'units': 'MW/m²',
        'description': 'Calculated average fuel surface heat flux at full power',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Burnup Steps': {
        'group': 'Overall System', 'units': 'MWd/kg',
        'description': 'Cumulative burnup values at each depletion time step',
        'source': 'User Input', 'hidden': False, 'array_mode': 'as_is'},  # was 'steps'

    'Time Steps': {
        'group': 'Overall System', 'units': 's',
        'description': 'User-specified depletion step durations used by the integrator',
        'source': 'User Input', 'hidden': False, 'array_mode': 'as_is'},

    # =========================================================
    # OpenMC Settings
    # =========================================================
   'Shutdown Margin Calc': {
        'group': 'OpenMC Settings', 'units': '',
        'description': 'Whether shutdown margin is calculated. When True, the shutdown-state keff is evaluated with the absorber facing the core (ARI) at Cold Shutdown Temperature, and shutdown margin metrics are reported from ((1 - k_s) / k_s).',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Cold Shutdown Temperature': {
    'group': 'OpenMC Settings', 'units': 'K',
    'description': 'Temperature used for the shutdown-state (ARI) keff calculation during shutdown margin evaluation.',
    'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Isothermal Temperature Coefficients': {
        'group': 'OpenMC Settings', 'units': '',
        'description': 'Whether isothermal temperature reactivity coefficients were calculated',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Temperature Perturbation': {
        'group': 'OpenMC Settings', 'units': 'K',
        'description': 'Temperature increase applied for isothermal temperature coefficient calculation',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Particles': {
        'group': 'OpenMC Settings', 'units': '',
        'description': 'Number of neutron histories per OpenMC batch',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Physics Results
    # =========================================================
    'Fuel Lifetime': {
        'group': 'Physics Results', 'units': 'days',
        'description': 'Estimated fuel cycle length — time for the 3D-corrected keff to fall to 1.0',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Depletion Time Steps': {
        'group': 'Physics Results',
        'units': 'days',
        'description': 'Calculated cumulative depletion times corresponding to the reported keff values and used as the x-axis in keff_comparison_vs_Time.png.',
        'source': 'Calculated',
        'hidden': False,
        'array_mode': 'as_is'},

    'BOL Axial Non-Leakage Probability': {
        'group': 'Physics Results',
        'units': 'fraction',
        'description': 'Beginning-of-life axial non-leakage probability from the 2D-to-3D buckling correction model.',
        'source': 'Calculated',
        'hidden': False,
        'array_mode': None},

    'Estimated Axial Leakage (%)': {
        'group': 'Physics Results',
        'units': '%',
        'description': 'Estimated beginning-of-life axial neutron leakage percentage from the axial buckling correction model, computed as 100 × (1 - axial non-leakage probability).',
        'source': 'Calculated',
        'hidden': False,
        'array_mode': None},

    'BOL Total Non-Leakage Probability': {
        'group': 'Physics Results',
        'units': 'fraction',
        'description': 'Beginning-of-life total non-leakage probability from the axial-plus-radial buckling model. Returns NaN when Core Radius is not available.',
        'source': 'Calculated',
        'hidden': False,
        'array_mode': None},

    'Estimated Total Leakage (%)': {
        'group': 'Physics Results',
        'units': '%',
        'description': 'Estimated beginning-of-life total neutron leakage percentage from the axial-plus-radial buckling model, computed as 100 × (1 - total non-leakage probability). Returns NaN when Core Radius is not available.',
        'source': 'Calculated',
        'hidden': False,
        'array_mode': None},

    'Mass U235': {
        'group': 'Physics Results', 'units': 'g',
        'description': 'Initial mass of U-235 in the fresh fuel load',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Mass U238': {
        'group': 'Physics Results', 'units': 'g',
        'description': 'Initial mass of U-238 in the fresh fuel load',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Uranium Mass': {
        'group': 'Physics Results', 'units': 'kg',
        'description': 'Total initial uranium mass in the fresh fuel load (U-235 + U-238)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Natural Uranium Mass': {
        'group': 'Physics Results', 'units': 'kg',
        'description': 'Mass of natural uranium feed required to produce the enriched fuel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Fuel Tail Waste Mass': {
        'group': 'Physics Results', 'units': 'kg',
        'description': 'Mass of depleted uranium tails produced during the enrichment process',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'SWU': {
        'group': 'Physics Results', 'units': 'kg-SWU',
        'description': 'Separative work units required to produce the enriched fuel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Temp Coeff 2D': {
        'group': 'Physics Results', 'units': 'pcm/K',
        'description': 'Most limiting (least negative) isothermal temperature coefficient from the 2D simulation. '
                       'Negative values indicate self-stabilizing behavior. '
                       'Note: 2D simulation overestimates keff — use 3D corrected value for safety conclusions.',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Temp Coeff 3D (2D corrected)': {
        'group': 'Physics Results', 'units': 'pcm/K',
        'description': 'Most limiting (least negative) isothermal temperature coefficient corrected from 2D to 3D '
                       'using axial neutron leakage correction. Negative values confirm self-stabilizing behavior. '
                       'Use this value for safety analysis.',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Most Limiting Shutdown Margin 2D': {
    'group': 'Physics Results', 'units': 'pcm',
    'description': 'Most limiting shutdown margin from the 2D simulation, computed as min((1 - k_s) / k_s) over depletion steps, where k_s is the shutdown-state keff with the absorber facing the core (ARI) at Cold Shutdown Temperature.',
    'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Maximum Shutdown Margin 2D': {
        'group': 'Physics Results', 'units': 'pcm',
        'description': 'Maximum shutdown margin from the 2D simulation, computed as max((1 - k_s) / k_s) over depletion steps, where k_s is the shutdown-state keff with the absorber facing the core (ARI) at Cold Shutdown Temperature.',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Most Limiting Shutdown Margin 3D (2D corrected)': {
        'group': 'Physics Results', 'units': 'pcm',
        'description': 'Most limiting shutdown margin corrected from 2D to 3D using axial neutron leakage correction, computed as min((1 - k_s) / k_s) over depletion steps, where k_s is the shutdown-state corrected keff with the absorber facing the core (ARI) at Cold Shutdown Temperature.',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Maximum Shutdown Margin 3D (2D corrected)': {
        'group': 'Physics Results', 'units': 'pcm',
        'description': 'Maximum shutdown margin corrected from 2D to 3D using axial neutron leakage correction, computed as max((1 - k_s) / k_s) over depletion steps, where k_s is the shutdown-state corrected keff with the absorber facing the core (ARI) at Cold Shutdown Temperature.',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},    

    'Max Peaking Factor': {
        'group': 'Physics Results', 'units': 'unitless',
        'description': 'Maximum pin peaking factor across all depletion steps — the most limiting value over the fuel lifetime',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Step with Max Peaking Factor': {
        'group': 'Physics Results', 'units': 'unitless',
        'description': 'Depletion step number at which the maximum peaking factor occurs',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Region ID with Max Peaking Factor': {
        'group': 'Physics Results', 'units': 'unitless',
        'description': 'Region ID with the highest peaking factor across all depletion steps; for LTMR this is typically a pin/distribcell ID, while for GCMR it is a mesh-region index',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Max Peaking Factors per Step': {
        'group': 'Physics Results', 'units': 'unitless',
        'description': 'List of the maximum pin peaking factor at each depletion step, in step order',
        'source': 'Calculated', 'hidden': False, 'array_mode': 'as_is'},

    # PF Summary is a nested dict — hidden from the Parameters sheet.
    # The key scalars above (Max Peaking Factor, Step with Max Peaking Factor,
    # Rod ID with Max Peaking Factor) carry the actionable information.
    'PF Summary': {
        'group': 'Physics Results', 'units': '',
        'description': 'Full peaking factor summary table across all depletion steps — see log output for per-pin results',
        'source': 'Calculated', 'hidden': True, 'array_mode': None},

    'keff 2D': {
        'group': 'Physics Results', 'units': '',
        'description': 'keff at each burnup step for the operating configuration (ARO) from the 2D OpenMC simulation, evaluated at Common Temperature.',
        'source': 'Calculated', 'hidden': False, 'array_mode': 'as_is'},

    'keff 3D (2D corrected)': {
        'group': 'Physics Results', 'units': '',
        'description': 'keff at each burnup step for the operating configuration (ARO), corrected from 2D to 3D using axial leakage correction and evaluated at Common Temperature.',
        'source': 'Calculated', 'hidden': False, 'array_mode': 'as_is'},

    # =========================================================
    # Primary Loop & Balance of Plant
    # =========================================================
    'Secondary HX Mass': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'kg',
        'description': 'Mass of the secondary heat exchanger (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Primary HX Mass': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'kg',
        'description': 'Calculated mass of the primary heat exchanger',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Primary Pump': {
        'group': 'Primary Loop & Balance of Plant', 'units': '',
        'description': 'Whether a primary coolant pump is present (Yes or No)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Secondary Pump': {
        'group': 'Primary Loop & Balance of Plant', 'units': '',
        'description': 'Whether a secondary coolant pump is present (Yes or No)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Pump Isentropic Efficiency': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'fraction',
        'description': 'Isentropic efficiency of the primary coolant pump',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Primary Loop Inlet Temperature': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'K',
        'description': 'Primary coolant temperature at the core inlet',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Primary Loop Outlet Temperature': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'K',
        'description': 'Primary coolant temperature at the core outlet',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Secondary Loop Inlet Temperature': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'K',
        'description': 'Secondary loop coolant temperature at the heat exchanger inlet',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Secondary Loop Outlet Temperature': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'K',
        'description': 'Secondary loop coolant temperature at the heat exchanger outlet',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Primary Loop Purification': {
        'group': 'Primary Loop & Balance of Plant', 'units': '',
        'description': 'Whether primary coolant purification is used (reduces coolant replacement frequency)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Primary Loop Count': {
        'group': 'Primary Loop & Balance of Plant', 'units': '',
        'description': 'Number of primary coolant loops',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Primary Loop per loop load fraction': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'fraction',
        'description': 'Fraction of total thermal load handled by each primary coolant loop',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Primary Loop Pressure Drop': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'Pa',
        'description': 'Pressure drop across the primary coolant loop',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Compressor Pressure Ratio': {
        'group': 'Primary Loop & Balance of Plant', 'units': '',
        'description': 'Compressor pressure ratio for gas-cooled reactor (GCMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Compressor Isentropic Efficiency': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'fraction',
        'description': 'Isentropic efficiency of the gas compressor (GCMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'BoP Count': {
        'group': 'Primary Loop & Balance of Plant', 'units': '',
        'description': 'Number of Balance of Plant (power conversion) systems',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'BoP per loop load fraction': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'fraction',
        'description': 'Fraction of total electric load handled by each BoP system',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'BoP Power kWe': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'kWe',
        'description': 'Electric power handled by each BoP system',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Coolant Mass Flow Rate': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'kg/s',
        'description': 'Total primary coolant mass flow rate',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Primary Loop Mass Flow Rate': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'kg/s',
        'description': 'Coolant mass flow rate per primary loop',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Primary Pump Mechanical Power': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'kW',
        'description': 'Mechanical power required to drive the primary coolant pump',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Primary Loop Compressor Power': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'W',
        'description': 'Mechanical power required to drive the primary loop gas compressor (GCMR)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Integrated Heat Transfer Vessel Thickness': {
        'group': 'Primary Loop & Balance of Plant', 'units': 'cm',
        'description': 'Wall thickness of the integrated heat transfer vessel (GCMR, 0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Integrated Heat Transfer Vessel Material': {
        'group': 'Primary Loop & Balance of Plant', 'units': '',
        'description': 'Material of the integrated heat transfer vessel (GCMR)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Shielding
    # =========================================================
    'In Vessel Shield Thickness': {
        'group': 'Shielding', 'units': 'cm',
        'description': 'Thickness of the in-vessel radiation shield (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'In Vessel Shield Inner Radius': {
        'group': 'Shielding', 'units': 'cm',
        'description': 'Inner radius of the in-vessel radiation shield (equal to core radius)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'In Vessel Shield Outer Radius': {
        'group': 'Shielding', 'units': 'cm',
        'description': 'Outer radius of the in-vessel radiation shield',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'In Vessel Shield Material': {
        'group': 'Shielding', 'units': '',
        'description': 'Material of the in-vessel radiation shield',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'In Vessel Shield Mass': {
        'group': 'Shielding', 'units': 'kg',
        'description': 'Total mass of the in-vessel radiation shield',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Out Of Vessel Shield Thickness': {
        'group': 'Shielding', 'units': 'cm',
        'description': 'Thickness of the out-of-vessel radiation shield',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Out Of Vessel Shield Material': {
        'group': 'Shielding', 'units': '',
        'description': 'Material of the out-of-vessel radiation shield',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Out Of Vessel Shield Effective Density Factor': {
        'group': 'Shielding', 'units': 'fraction',
        'description': 'Effective density factor for the out-of-vessel shield — accounts for partial fill (e.g. 0.5 = 50% fill)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Out Of Vessel Shield Mass': {
        'group': 'Shielding', 'units': 'kg',
        'description': 'Total mass of the out-of-vessel radiation shield',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Outer Shield Outer Radius': {
        'group': 'Shielding', 'units': 'cm',
        'description': 'Outer radius of the out-of-vessel shield (outermost boundary of the shielded assembly)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Vessels
    # =========================================================
    'Vessel Radius': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Inner radius of the reactor pressure vessel (or core barrel)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Vessel Thickness': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Wall thickness of the reactor pressure vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Vessel Lower Plenum Height': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Height of the lower plenum region below the active core',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Vessel Upper Plenum Height': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Height of the upper plenum region above the active core',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Vessel Upper Gas Gap': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Height of the upper gas gap above the upper plenum (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Vessel Bottom Depth': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Depth of the vessel bottom head (ellipsoidal)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Vessel Material': {
        'group': 'Vessels', 'units': '',
        'description': 'Structural material of the reactor pressure vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Vessel Height': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Total height of the reactor pressure vessel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Vessel Mass': {
        'group': 'Vessels', 'units': 'kg',
        'description': 'Total mass of the reactor pressure vessel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'RPV Outer Radius': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Outer radius of the reactor pressure vessel (inner radius + wall thickness)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'RPV Outer Height': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Total outer height of the reactor pressure vessel including heads',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Integrated Heat Transfer Vessel Mass': {
        'group': 'Vessels', 'units': 'kg',
        'description': 'Total mass of the integrated heat transfer vessel (GCMR, 0 if none)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Integrated Heat Transfer Vessel Outer Volume': {
        'group': 'Vessels', 'units': 'cm³',
        'description': 'Outer volume of the integrated heat transfer vessel (GCMR, 0 if none)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Gap Between Vessel And Guard Vessel': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Radial gap between the reactor vessel and guard vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Guard Vessel Thickness': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Wall thickness of the guard vessel (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Guard Vessel Material': {
        'group': 'Vessels', 'units': '',
        'description': 'Structural material of the guard vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Guard Vessel Radius': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Outer radius of the guard vessel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Guard Vessel Mass': {
        'group': 'Vessels', 'units': 'kg',
        'description': 'Total mass of the guard vessel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Gap Between Guard Vessel And Cooling Vessel': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Radial gap between the guard vessel and cooling vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Cooling Vessel Thickness': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Wall thickness of the cooling vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Cooling Vessel Material': {
        'group': 'Vessels', 'units': '',
        'description': 'Structural material of the cooling vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Cooling Vessel Radius': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Outer radius of the cooling vessel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Cooling Vessel Mass': {
        'group': 'Vessels', 'units': 'kg',
        'description': 'Total mass of the cooling vessel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Gap Between Cooling Vessel And Intake Vessel': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Radial gap between the cooling vessel and intake vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Intake Vessel Thickness': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Wall thickness of the intake vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Intake Vessel Material': {
        'group': 'Vessels', 'units': '',
        'description': 'Structural material of the intake vessel',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Intake Vessel Radius': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Outer radius of the intake vessel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Intake Vessel Mass': {
        'group': 'Vessels', 'units': 'kg',
        'description': 'Total mass of the intake vessel',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Vessels Total Radius': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Total outer radius of the complete vessel assembly',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Vessels Total Height': {
        'group': 'Vessels', 'units': 'cm',
        'description': 'Total height of the complete vessel assembly',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Total Vessels Mass': {
        'group': 'Vessels', 'units': 'kg',
        'description': 'Combined mass of all vessel components (pressure vessel, guard, cooling, intake)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Operation
    # =========================================================
    'Operation Mode': {
        'group': 'Operation', 'units': '',
        'description': 'Reactor operation mode: Remotely Monitored or On-Site Staffed',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Number of Operators': {
        'group': 'Operation', 'units': '',
        'description': 'Number of operators assigned per shift',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Levelization Period': {
        'group': 'Operation', 'units': 'years',
        'description': 'Plant lifetime used for LCOE levelization (economic analysis period)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Refueling Period': {
        'group': 'Operation', 'units': 'days',
        'description': 'Duration of the scheduled refueling outage',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Emergency Shutdowns Per Year': {
        'group': 'Operation', 'units': 'events/year',
        'description': 'Expected number of unplanned emergency shutdowns per year',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Startup Duration after Refueling': {
        'group': 'Operation', 'units': 'days',
        'description': 'Duration of the startup sequence following a refueling outage',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Startup Duration after Emergency Shutdown': {
        'group': 'Operation', 'units': 'days',
        'description': 'Duration of the startup sequence following an emergency shutdown',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Reactors Monitored Per Operator': {
        'group': 'Operation', 'units': '',
        'description': 'Number of reactors that can be monitored by a single operator',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Security Staff Per Shift': {
        'group': 'Operation', 'units': '',
        'description': 'Number of security personnel required per shift',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Onsite Coolant Inventory': {
        'group': 'Operation', 'units': 'kg',
        'description': 'Total mass of coolant stored on site (fresh + used tanks)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Replacement Coolant Inventory': {
        'group': 'Operation', 'units': 'kg',
        'description': 'Mass of coolant replaced per refueling cycle',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Annual Coolant Supply Frequency': {
        'group': 'Operation', 'units': 'deliveries/year',
        'description': 'Number of coolant resupply deliveries per year (1 if purified, 6 if not)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Capacity Factor': {
        'group': 'Operation', 'units': 'fraction',
        'description': 'Annual plant capacity factor — fraction of the year the plant operates at full power',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Annual Electricity Production': {
        'group': 'Operation', 'units': 'MWh/year',
        'description': 'Total annual electricity production at the calculated capacity factor',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'A75: Vessel Replacement Period (cycles)': {
        'group': 'Operation', 'units': 'cycles',
        'description': 'Number of fuel cycles between scheduled reactor pressure vessel replacements',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'A75: Core Barrel Replacement Period (cycles)': {
        'group': 'Operation', 'units': 'cycles',
        'description': 'Number of fuel cycles between scheduled core barrel replacements',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'A75: Reflector Replacement Period (cycles)': {
        'group': 'Operation', 'units': 'cycles',
        'description': 'Number of fuel cycles between scheduled reflector replacements',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'A75: Drum Replacement Period (cycles)': {
        'group': 'Operation', 'units': 'cycles',
        'description': 'Number of fuel cycles between scheduled control drum replacements',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'A75: Integrated HX Replacement Period (cycles)': {
        'group': 'Operation', 'units': 'cycles',
        'description': 'Number of fuel cycles between scheduled integrated heat exchanger replacements',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Maintenance to Direct Cost Ratio': {
        'group': 'Operation', 'units': 'fraction',
        'description': 'Ratio of annual non-replacement maintenance cost to total direct capital cost',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'A78: CAPEX to Decommissioning Cost Ratio': {
        'group': 'Operation', 'units': 'fraction',
        'description': 'Ratio of total decommissioning cost to total capital cost (typical range: 9%–15%)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Economic Parameters
    # =========================================================
    'Interest Rate': {
        'group': 'Economic Parameters', 'units': 'fraction',
        'description': 'Annual cost of debt used to calculate interest during construction (Account 62)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Discount Rate': {
        'group': 'Economic Parameters', 'units': 'fraction',
        'description': 'Annual discount rate (WACC) used for LCOE levelization and cost annualization',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Construction Duration': {
        'group': 'Economic Parameters', 'units': 'months',
        'description': 'Duration of plant construction used for financing cost calculation',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Debt To Equity Ratio': {
        'group': 'Economic Parameters', 'units': 'fraction',
        'description': 'Ratio of debt financing to equity financing for the project',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Annual Return': {
        'group': 'Economic Parameters', 'units': 'fraction',
        'description': 'Annual return rate assumed for the decommissioning reserve fund',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'NOAK Unit Number': {
        'group': 'Economic Parameters', 'units': '',
        'description': 'Nth-of-a-kind unit number used to calculate the learning curve cost multiplier',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Assumed Number Of Units For Onsite Learning': {
        'group': 'Economic Parameters', 'units': '',
        'description': 'Number of units assumed for onsite learning (typically 2 × NOAK Unit Number)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Land Area': {
        'group': 'Economic Parameters', 'units': 'acres',
        'description': 'Total land area required for the plant site',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Escalation Year': {
        'group': 'Economic Parameters', 'units': 'year',
        'description': 'Reference year to which all costs are escalated for consistent comparison',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Number of Samples': {
        'group': 'Economic Parameters', 'units': '',
        'description': 'Number of Monte Carlo samples used for cost uncertainty quantification',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'indirect to direct field-related cost': {
        'group': 'Economic Parameters', 'units': 'fraction',
        'description': 'Ratio of indirect field (site) costs to total direct field costs — covers site supervision, '
                       'temporary facilities, and other indirect site activities',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'FTEs Per Onsite Operator (24/7)': {
        'group': 'Economic Parameters', 'units': 'FTE/operator',
        'description': 'Full-time equivalents needed per onsite operator position to maintain 24/7 coverage',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'FTEs Per Offsite Operator (24/7)': {
        'group': 'Economic Parameters', 'units': 'FTE/operator',
        'description': 'Full-time equivalents needed per offsite operator position to maintain 24/7 coverage',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'FTEs Per Security Staff (24/7)': {
        'group': 'Economic Parameters', 'units': 'FTE/staff',
        'description': 'Full-time equivalents needed per security staff position to maintain 24/7 coverage',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Hours Per FTE': {
        'group': 'Economic Parameters', 'units': 'hours/year',
        'description': 'Annual working hours per full-time equivalent employee',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Work Hours Per Shift': {
        'group': 'Economic Parameters', 'units': 'hours',
        'description': 'Duration of each work shift',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'FTEs Per Operator Per Year Per Refueling': {
        'group': 'Economic Parameters', 'units': 'FTE/operator/year',
        'description': 'Fraction of an FTE spent per operator per year on refueling-related activities',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'FTEs Per Onsite Operator Per Year': {
        'group': 'Economic Parameters', 'units': 'FTE/operator/year',
        'description': 'Total FTE-years per operator per year including all activities',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Excavation Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of earth to be excavated for plant construction',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Reactor Building
    'Reactor Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the reactor building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Reactor Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the reactor building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Reactor Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the reactor building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Reactor Building Superstructure Area': {
        'group': 'Economic Parameters', 'units': 'm²',
        'description': 'Floor area of the reactor building superstructure',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Integrated Heat Exchanger Building
    'Integrated Heat Exchanger Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the integrated heat exchanger building slab roof (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Integrated Heat Exchanger Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the integrated heat exchanger building basement (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Integrated Heat Exchanger Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the integrated heat exchanger building exterior walls (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Integrated Heat Exchanger Building Superstructure Area': {
        'group': 'Economic Parameters', 'units': 'm²',
        'description': 'Floor area of the integrated heat exchanger building superstructure (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Turbine Building
    'Turbine Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the turbine building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Turbine Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the turbine building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Turbine Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the turbine building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Control Building
    'Control Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the control building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Control Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the control building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Control Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the control building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Manipulator Building
    'Manipulator Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the manipulator building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manipulator Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the manipulator building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manipulator Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the manipulator building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Refueling Building
    'Refueling Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the refueling building slab roof (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Refueling Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the refueling building basement (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Refueling Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the refueling building exterior walls (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Spent Fuel Building
    'Spent Fuel Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the spent fuel building slab roof (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Spent Fuel Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the spent fuel building basement (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Spent Fuel Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the spent fuel building exterior walls (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Emergency Building
    'Emergency Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the emergency building slab roof (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Emergency Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the emergency building basement (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Emergency Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the emergency building exterior walls (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Storage Building
    'Storage Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the storage building slab roof (houses operational spares)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Storage Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the storage building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Storage Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the storage building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Building volumes — Radwaste Building
    'Radwaste Building Slab Roof Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the radioactive waste building slab roof (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radwaste Building Basement Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the radioactive waste building basement (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radwaste Building Exterior Walls Volume': {
        'group': 'Economic Parameters', 'units': 'm³',
        'description': 'Volume of concrete in the radioactive waste building exterior walls (0 if none)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # =========================================================================
    # CENTRAL FACILITY PARAMETERS
    # =========================================================================

    # --- Overall Central Facility ---
    'Estimate Central Facility': {
        'group': 'Central Facility', 'units': '',
        'description': 'Enable central facility cost estimation (True/False). '
                       'When True, reads from Central Facility Database sheet.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Maximum Number of Operating Reactors': {
        'group': 'Central Facility', 'units': '',
        'description': 'Maximum number of reactor units the central facility supports. '
                       'Used to normalize costs per kW of total fleet capacity.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Central Facility Construction Duration': {
        'group': 'Central Facility', 'units': 'months',
        'description': 'Construction duration for central facility. '
                       'Used for financing cost (interest during construction).',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Central Facility Power MWe': {
        'group': 'Central Facility', 'units': 'MWe',
        'description': 'Electrical capacity of the central facility for its own operations.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Site Perimeter': {
        'group': 'Central Facility', 'units': 'm',
        'description': 'Perimeter of entire central facility site for security fencing.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Maintenance Staff Per Shift': {
        'group': 'Central Facility', 'units': 'FTEs',
        'description': 'Number of maintenance staff per shift at the central facility.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Power Mwt of Operating Fleet': {
        'group': 'Central Facility', 'units': 'MWt',
        'description': 'Total thermal power of operating reactor fleet '
                       '(Power MWt × Maximum Number of Operating Reactors).',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    # --- Servicing Facility ---
    'Servicing Facility Perimeter': {
        'group': 'Central Facility', 'units': 'm',
        'description': 'Perimeter of servicing facility for security fencing.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Total Servicing Rate': {
        'group': 'Central Facility', 'units': 'reactors/year',
        'description': 'Number of reactors serviced (refueled/maintained) per year.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Servicing Hot Cell Count': {
        'group': 'Central Facility', 'units': '',
        'description': 'Number of hot cells for reactor servicing operations.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Servicing Hot Cell Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of each servicing hot cell.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Total Servicing Hot Cell Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total volume of all servicing hot cells '
                       '(Servicing Hot Cell Count × Servicing Hot Cell Volume).',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Defueling/Refueling Line Count': {
        'group': 'Central Facility', 'units': '',
        'description': 'Number of defueling/refueling lines in servicing facility.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Power Mwt Processed by Servicing': {
        'group': 'Central Facility', 'units': 'MWt',
        'description': 'Thermal power processed by servicing facility '
                       '(assumes 5% power for low-power testing per hot cell).',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    # Servicing Facility Building Volumes
    'Servicing Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility main building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility main building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility main building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the servicing facility main building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Helium Purification and Storage Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the helium purification and storage building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Helium Purification and Storage Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the helium purification and storage building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Helium Purification and Storage Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the helium purification and storage building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Helium Purification and Storage Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the helium purification and storage building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Servicing Facility Integrated Control Room Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility integrated control room slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Facility Integrated Control Room Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility integrated control room basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Facility Integrated Control Room Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility integrated control room exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Facility Integrated Control Room Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the servicing facility integrated control room',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Servicing Facility Admin Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility admin building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Facility Admin Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility admin building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Facility Admin Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility admin building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Facility Admin Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the servicing facility admin building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Servicing Facility Security Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility security building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Facility Security Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility security building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Facility Security Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the servicing facility security building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Servicing Facility Security Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the servicing facility security building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # --- Manufacturing/Factory Facility ---
    'Factory Perimeter': {
        'group': 'Central Facility', 'units': 'm',
        'description': 'Perimeter of manufacturing/factory facility.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'New Reactor Production Rate': {
        'group': 'Central Facility', 'units': 'reactors/year',
        'description': 'Number of new reactors produced per year.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # Manufacturing Facility Building Volumes
    'Fabrication Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the fabrication building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Fabrication Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the fabrication building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Fabrication Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the fabrication building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Fabrication Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the fabrication building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Feed and Product Warehouse Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the feed and product warehouse building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Feed and Product Warehouse Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the feed and product warehouse building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Feed and Product Warehouse Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the feed and product warehouse building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Feed and Product Warehouse Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the feed and product warehouse building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Manufacturing Facility Integrated Control Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the manufacturing facility integrated control building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manufacturing Facility Integrated Control Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the manufacturing facility integrated control building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manufacturing Facility Integrated Control Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the manufacturing facility integrated control building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manufacturing Facility Integrated Control Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the manufacturing facility integrated control building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Manufacturing Facility Admin Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the manufacturing facility admin building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manufacturing Facility Admin Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the manufacturing facility admin building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manufacturing Facility Admin Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the manufacturing facility admin building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manufacturing Facility Admin Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the manufacturing facility admin building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Manufacturing Facility Security Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the manufacturing facility security building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manufacturing Facility Security Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the manufacturing facility security building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manufacturing Facility Security Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the manufacturing facility security building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Manufacturing Facility Security Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the manufacturing facility security building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # --- New Reactor Facility ---
    'New Reactor Facility Perimeter': {
        'group': 'Central Facility', 'units': 'm',
        'description': 'Perimeter of new reactor fueling and testing facility.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'New Reactor Fueling Line Count': {
        'group': 'Central Facility', 'units': '',
        'description': 'Number of fueling lines for new reactors.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'New Reactor Testing Line Count': {
        'group': 'Central Facility', 'units': '',
        'description': 'Number of testing lines for new reactors.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'New Reactor Testing Hot Cell Count': {
        'group': 'Central Facility', 'units': '',
        'description': 'Number of hot cells for new reactor testing.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'New Reactor Testing Hot Cell Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total volume of all new reactor testing hot cells.',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Power Mwe Processed by New Reactor Facility': {
        'group': 'Central Facility', 'units': 'MWe',
        'description': 'Electrical capacity processed by new reactor facility '
                       '(Power MWe × New Reactor Production Rate).',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    # New Reactor Facility Building Volumes
    'Fresh Fuel Storage and Inspection Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the fresh fuel storage and inspection building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Fresh Fuel Storage and Inspection Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the fresh fuel storage and inspection building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Fresh Fuel Storage and Inspection Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the fresh fuel storage and inspection building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Fresh Fuel Storage and Inspection Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the fresh fuel storage and inspection building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Reactor Fueling Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the reactor fueling building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Reactor Fueling Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the reactor fueling building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Reactor Fueling Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the reactor fueling building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Reactor Fueling Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the reactor fueling building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Reactor Testing Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the reactor testing building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Reactor Testing Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the reactor testing building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Reactor Testing Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the reactor testing building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Reactor Testing Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the reactor testing building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'New Reactor Fuel and Testing Facility Admin Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the new reactor facility admin building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'New Reactor Fuel and Testing Facility Admin Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the new reactor facility admin building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'New Reactor Fuel and Testing Facility Admin Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the new reactor facility admin building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'New Reactor Fuel and Testing Facility Admin Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the new reactor facility admin building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'New Reactor Fuel and Testing Facility Security Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the new reactor facility security building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'New Reactor Fuel and Testing Facility Security Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the new reactor facility security building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'New Reactor Fuel and Testing Facility Security Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the new reactor facility security building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'New Reactor Fuel and Testing Facility Security Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the new reactor facility security building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # --- Radioactive Waste Management Facility ---
    'Radioactive Waste Management Facility Perimeter': {
        'group': 'Central Facility', 'units': 'm',
        'description': 'Perimeter of radioactive waste management facility.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Radioactive Waste Processing Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total volume of radioactive waste processing building.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Processing Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the radioactive waste processing building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Processing Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the radioactive waste processing building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Processing Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the radioactive waste processing building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Radioactive Waste Storage Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total volume of radioactive waste storage building.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Storage Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the radioactive waste storage building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Storage Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the radioactive waste storage building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Storage Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the radioactive waste storage building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Radioactive Waste Management Facility Integrated Control Room Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the waste management facility integrated control room slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Management Facility Integrated Control Room Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the waste management facility integrated control room basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Management Facility Integrated Control Room Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the waste management facility integrated control room exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Management Facility Integrated Control Room Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the waste management facility integrated control room',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Radioactive Waste Management Facility Admin Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the waste management facility admin building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Management Facility Admin Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the waste management facility admin building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Management Facility Admin Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the waste management facility admin building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Management Facility Admin Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the waste management facility admin building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Radioactive Waste Management Facility Security Building Roof Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the waste management facility security building slab roof',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Management Facility Security Building Basement Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the waste management facility security building basement',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Management Facility Security Building Walls Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Volume of concrete in the waste management facility security building exterior walls',
        'source': 'User Input', 'hidden': False, 'array_mode': None},
    'Radioactive Waste Management Facility Security Building Volume': {
        'group': 'Central Facility', 'units': 'm³',
        'description': 'Total enclosed volume of the waste management facility security building',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # --- Transportation ---
    'Local Transport Vehicle Count': {
        'group': 'Central Facility', 'units': '',
        'description': 'Number of local transport vehicles within central facility.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Reactor Transport Vehicle Count': {
        'group': 'Central Facility', 'units': '',
        'description': 'Number of vehicles for transporting reactor units to/from field sites.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Spares/Waste Transport Vehicle Count': {
        'group': 'Central Facility', 'units': '',
        'description': 'Number of vehicles for transporting spare parts and radioactive waste.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'General Transport Vehicle Count': {
        'group': 'Central Facility', 'units': '',
        'description': 'Number of general purpose transport vehicles.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # --- Casks ---
    'Annual Spent Fuel Cask Replacement': {
        'group': 'Central Facility', 'units': 'casks/year',
        'description': 'Annual replacement rate for spent fuel transport casks.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Annual Reactor Cask Replacement': {
        'group': 'Central Facility', 'units': 'casks/year',
        'description': 'Annual replacement rate for reactor transport casks.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Annual Rad Waste Cask Replacement': {
        'group': 'Central Facility', 'units': 'casks/year',
        'description': 'Annual replacement rate for radioactive waste transport casks.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Learning Rates
    # =========================================================
    'No Learning': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'Learning rate for components with no learning benefit (fixed costs, rate = 0)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Licensing Learning': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'Learning rate for licensing costs per doubling of cumulative units produced',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Factory Primary Structure': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'Learning rate for factory-manufactured primary structural components',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Factory Drums': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'Learning rate for factory-manufactured control drum components',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Factory Other': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'Learning rate for other factory-manufactured components',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Factory Be': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'Learning rate for factory-manufactured beryllium components',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Factory BeO': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'Learning rate for factory-manufactured beryllium oxide components',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Non-nuclear off-the-shelf': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'Learning rate for non-nuclear commercial off-the-shelf components',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Onsite Learning': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'Learning rate for onsite construction and installation activities',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'No Learning Cost Multiplier': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'NOAK cost multiplier for components with no learning (always 1.0)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Licensing Learning Cost Multiplier': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'NOAK cost multiplier for licensing costs based on Wright\'s Law',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Factory Primary Structure Cost Multiplier': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'NOAK cost multiplier for factory-built primary structure based on Wright\'s Law',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Factory Drums Cost Multiplier': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'NOAK cost multiplier for factory-built control drums based on Wright\'s Law',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Factory Other Cost Multiplier': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'NOAK cost multiplier for other factory-built components based on Wright\'s Law',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Factory Be Cost Multiplier': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'NOAK cost multiplier for beryllium components based on Wright\'s Law',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Factory BeO Cost Multiplier': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'NOAK cost multiplier for beryllium oxide components based on Wright\'s Law',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Non-nuclear off-the-shelf Cost Multiplier': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'NOAK cost multiplier for non-nuclear commercial components based on Wright\'s Law',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'Onsite Learning Cost Multiplier': {
        'group': 'Learning Rates', 'units': 'fraction',
        'description': 'NOAK cost multiplier for onsite construction activities based on Wright\'s Law',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Tax Credits
    # =========================================================
    'ITC credit level': {
        'group': 'Tax Credits', 'units': 'fraction',
        'description': 'Investment Tax Credit level under IRA Section 48E. '
                       'Values: 0.06 (base), 0.30 (prevailing wage), 0.40 (+domestic content), 0.50 (+energy community). '
                       'Applied as a one-time reduction to the Overnight Capital Cost.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'PTC credit value': {
        'group': 'Tax Credits', 'units': '$/MWh',
        'description': 'Production Tax Credit value per MWh under IRA Section 45Y. '
                       'Values: $3/MWh (base) or $15/MWh (with prevailing wage requirements met).',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'PTC credit period': {
        'group': 'Tax Credits', 'units': 'years',
        'description': 'Duration of the Production Tax Credit eligibility period (typically 10 years under IRA)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Tax Rate': {
        'group': 'Tax Credits', 'units': 'fraction',
        'description': 'Corporate tax rate used to gross up the PTC to its before-tax equivalent in LCOE calculation. '
                       'US federal rate is 0.21 (21%). Use 0.0 for tax-exempt entities.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'domestic_content_bonus': {
        'group': 'Tax Credits', 'units': 'fraction',
        'description': 'IRA stackable bonus for meeting domestic content requirements (+10% on base PTC)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'energy_community_bonus': {
        'group': 'Tax Credits', 'units': 'fraction',
        'description': 'IRA stackable bonus for siting the plant in a designated energy community (+10% on base PTC)',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    'Number of Units Claiming ITC/PTC': {
        'group': 'Tax Credits', 'units': '',
        'description': 'Number of units in the order book that may claim ITC or PTC before the IRA sunset. '
                       'A unit is eligible only if its position in the deployment sequence is <= this cutoff. '
                       'FOAK is unit 1; the NOAK column represents unit "NOAK Unit Number". '
                       'When a unit is past the cutoff, ITC/PTC-adjusted outputs fall back to the un-subsidized values.',
        'source': 'User Input', 'hidden': False, 'array_mode': None},

    # =========================================================
    # Debug / Intermediate Values
    # =========================================================
    'keff 2D high temp': {
        'group': 'Debug / Intermediate Values', 'units': '',
        'description': 'keff at each burnup step from 2D simulation at elevated temperature (Common Temperature + Perturbation) — '
                       'used to calculate isothermal temperature coefficient',
        'source': 'Calculated', 'hidden': False, 'array_mode': 'summary'},

    'keff 3D (2D corrected) high temp': {
        'group': 'Debug / Intermediate Values', 'units': '',
        'description': 'keff at elevated temperature corrected to 3D — used to calculate temperature coefficient',
        'source': 'Calculated', 'hidden': False, 'array_mode': 'summary'},

    'keff 2D ARI': {
        'group': 'Debug / Intermediate Values', 'units': '',
        'description': 'keff at each burnup step for the shutdown configuration with all control drums inserted (ARI) from the 2D simulation, evaluated at Cold Shutdown Temperature for shutdown margin calculation.',
        'source': 'Calculated', 'hidden': False, 'array_mode': 'summary'},

    'keff 3D (2D corrected) ARI': {
        'group': 'Debug / Intermediate Values', 'units': '',
        'description': 'keff at each burnup step for the shutdown configuration with all control drums inserted (ARI), corrected from 2D to 3D and evaluated at Cold Shutdown Temperature for shutdown margin calculation.',
        'source': 'Calculated', 'hidden': False, 'array_mode': 'summary'},

    'number of drums': {
        'group': 'Debug / Intermediate Values', 'units': '',
        'description': 'Total number of control drum positions in the core lattice (GCMR internal variable)',
        'source': 'Calculated', 'hidden': False, 'array_mode': None},

    'keff 2D ARO': {
        'group': 'Debug / Intermediate Values', 'units': '',
        'description': 'keff at each burnup step for the operating configuration with all control drums withdrawn (ARO) from the 2D simulation, evaluated at Common Temperature.',
        'source': 'Calculated', 'hidden': False, 'array_mode': 'summary'},

    'keff 3D (2D corrected) ARO': {
        'group': 'Debug / Intermediate Values', 'units': '',
        'description': 'keff at each burnup step for the operating configuration with all control drums withdrawn (ARO), corrected from 2D to 3D and evaluated at Common Temperature.',
        'source': 'Calculated', 'hidden': False, 'array_mode': 'summary'},

    'Constant': {
        'group': 'Debug / Intermediate Values', 'units': '',
        'description': 'Internal scaling constant set to 1, used as a no-op scaling variable '
                       'for fixed-cost central facility accounts',
        'source': 'Calculated', 'hidden': True, 'array_mode': None},
}
