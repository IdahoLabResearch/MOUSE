# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
# Importing libraries
import openmc

def collect_materials_data(params):
    
    # **************************************************************************************************************************
    #                                               Sec. 1 : MATERIALS
    # **************************************************************************************************************************
    materials = openmc.Materials()
    materials_database = {}
    print("Reading the Materials Database")

    # """""""""""""""""""""
    # Sec. 1.1 : Fuels: TRIGA Fuel and UO2, Uranium Carbide and Nitride
    # """""""""""""""""""""

    # The fuel is declared following a mechanism close to the logic of actual TRIGA fuel
    # fabrication. This type of fuel is specified, for example, as "45/20" fuel,
    # meaning 45% by weight is uranium metal (20% enriched) in a ZrH matrix.
    # TRIGA fuel can also contain 3% by weight of erbium as a burnable absorber in the fuel meat.

    # Declare the individual components of the fuel
    try:
        U_met = openmc.Material(name="U_met")
        U_met.set_density("g/cm3", 19.05)
        U_met.add_nuclide("U235", params['Enrichment'])
        U_met.add_nuclide("U238", 1 - params['Enrichment'])

        ZrH_fuel = openmc.Material(name="ZrH_fuel")
        ZrH_fuel.set_density("g/cm3", 5.63)
        ZrH_fuel.add_element("zirconium", 1.0)
        ZrH_fuel.add_nuclide("H1", params["H_Zr_ratio"])

        Er_bp = openmc.Material(name="Er_bp")
        Er_bp.set_density("g/cm3", 9.07)
        Er_bp.add_element("erbium", 1.0)

        er_wo = params['er_wo'] # burnable poison

        UZrH_alloy = openmc.Material.mix_materials(
            [U_met, ZrH_fuel, Er_bp],
            [params['U_met_wo'], 1 - params['U_met_wo'] - er_wo, er_wo],
            "wo", name="UZrH")
        UZrH_alloy.temperature = params['Common Temperature']
        UZrH_alloy.add_s_alpha_beta("c_H_in_ZrH")
        materials.append(UZrH_alloy)
        materials_database.update({'UZrH_alloy': UZrH_alloy})
    
    except KeyError as e:
        print(f"Skipping UZrH_alloy due to missing parameter: {e}")    

    # UO2
    try:
        UO2 = openmc.Material(name='UO2')
        UO2.set_density('g/cm3', 10.41)
        UO2.add_element('U', 1.0, enrichment=100 * params['Enrichment'])
        UO2.add_nuclide('O16', 2.0)
        UO2.add_s_alpha_beta("c_U_in_UO2")
        UO2.add_s_alpha_beta("c_O_in_UO2")
        materials.append(UO2)
        materials_database.update({'UO2': UO2})
    except KeyError as e:
        print(f"Skipping UO2 due to missing parameter: {e}")     

    # Uranium Carbide
    try:
        UC = openmc.Material(name='UC')
        UC.set_density('g/cm3', 13.0)
        UC.add_element('U', 1.0, enrichment=100 * params['Enrichment'])
        UC.add_element('C', 1.0)
        materials.append(UC)
        materials_database.update({'UC': UC})
    except KeyError as e:
        print(f"Skipping UC due to missing parameter: {e}")    

    # UCO: Mixed uranium dioxide (UO2) and uranium carbide (UC)
    # OpenMC cannot mix materials that already have S(α,β) tables attached.
    # To work around this, a separate UO2_for_mix material is created without
    # S(α,β) tables for mixing purposes only. The S(α,β) tables are then added
    # to the resulting UCO material after mixing.
    # The standalone UO2 material (used directly as fuel) is unaffected.
    # UCO: Mixed uranium dioxide (UO2) and uranium carbide (UC)
    try:
        UO2_for_mix = openmc.Material(name='UO2_for_mix')
        UO2_for_mix.set_density('g/cm3', 10.41)
        UO2_for_mix.add_element('U', 1.0, enrichment=100 * params['Enrichment'])
        UO2_for_mix.add_nuclide('O16', 2.0)

        UCO = openmc.Material.mix_materials(
            [UO2_for_mix, UC],
            [params['UO2 atom fraction'], 1 - params['UO2 atom fraction']],
            'ao',
            name='UCO'
        )
        UCO.temperature = params['Common Temperature']

        # Optional approximation:
        UCO.add_s_alpha_beta("c_U_in_UO2")
        UCO.add_s_alpha_beta("c_O_in_UO2")

        materials.append(UCO)
        materials_database.update({'UCO': UCO})

    except (KeyError, NameError) as e:
        print(f"Skipping UCO due to missing parameter/material: {e}")
    
    # Uranium Nitride
    try:
        UN = openmc.Material(name='UN')  # creates a new material named 'UN'
        UN.set_density('g/cm3', 14.0)
        UN.add_element('U', 1.0, enrichment=100 * params['Enrichment'])
        UN.add_element('N', 1.0)  # adds nitrogen (N) to the material
        UN.add_s_alpha_beta("c_U_in_UN")
        UN.add_s_alpha_beta("c_N_in_UN")
        materials.append(UN)
        materials_database.update({'UN': UN})
    except KeyError as e:
        print(f"Skipping UN due to missing parameter: {e}") 
        
    # U-10Zr
    try:
        UZr = openmc.Material(name='UZr') 
        UZr.set_density('g/cm3', 16.0)
        UZr.add_element('U', 10, 'wo', enrichment=100 * params['Enrichment'])
        UZr.add_element('Zr', 90, 'wo')
        materials.append(UZr)
        materials_database.update({'UZr': UZr})
    except KeyError as e:
        print(f"Skipping U-10Zr due to missing parameter: {e}")

    # Homogenized TRISO fuel
    try:
        U_total = 0.00130037929          # Total U atom density (U235+U238)
        density = 8.08250295E-02  # Total density (atom/b-cm)
        U235_frac = params['Enrichment'] * U_total
        U238_frac = (1 - params['Enrichment']) * U_total
        homog_TRISO = openmc.Material(name='homog_TRISO')
        homog_TRISO.set_density('atom/b-cm', density)
        homog_TRISO.temperature = params['Common Temperature']
        homog_TRISO.add_nuclide('U235', U235_frac, 'ao')
        homog_TRISO.add_nuclide('U238', U238_frac, 'ao')
        homog_TRISO.add_nuclide('O16', 2.59371545E-03, 'ao')
        homog_TRISO.add_nuclide('O17', 1.05004397E-06, 'ao')
        homog_TRISO.add_nuclide('O18', 5.99797186E-06, 'ao')
        homog_TRISO.add_nuclide('Si28', 2.76954169E-03, 'ao')
        homog_TRISO.add_nuclide('Si29', 1.40694868E-04, 'ao')
        homog_TRISO.add_nuclide('Si30', 9.28556098E-05, 'ao')
        homog_TRISO.add_nuclide('C12', 7.31619752E-02, 'ao')
        homog_TRISO.add_nuclide('C13', 7.58819416E-04, 'ao')
        homog_TRISO.add_s_alpha_beta('c_Graphite')
        materials.append(homog_TRISO)
        materials_database.update({'homog_TRISO': homog_TRISO})
    except KeyError as e:
        print(f"Skipping homog_TRISO due to missing parameter: {e}")  

    # """""""""""""""""""""
    # Sec. 1.2 : Hydrides: Zirconium Hydride and yttrium hydride (YHx)
    # """""""""""""""""""""
       
    ZrH = openmc.Material(name="ZrH", temperature=params['Common Temperature'])
    # Internal override used by density-aware temperature-coefficient cases.
    # Other MOUSE calculations retain the original 5.6 g/cm3 value.
    zrh_density = float(params.get('_ZrH Density Override', 5.6))
    ZrH.set_density("g/cm3", zrh_density)
    ZrH.add_nuclide("H1", 1.85)
    ZrH.add_element("zirconium", 1.0)
    ZrH.add_s_alpha_beta("c_H_in_ZrH")

    # Yttrium hydride (YHx) 
    YHx = openmc.Material(name="YHx")
    YHx.set_density("g/cm3", 4.28)
    YHx.add_nuclide("H1", 1.5)  # adds the hydrogen-1 (H-1) nuclide with an atomic ratio of 1.5
    YHx.add_element("yttrium", 1.0)  # adds yttrium with an atomic ratio of 1.0
    # Add thermal scattering data for hydrogen in yttrium hydride (YH2).
    # The add_s_alpha_beta method specifies the S(α,β) thermal scattering treatment for specific materials.
    YHx.add_s_alpha_beta("c_H_in_YH2")
    YHx.temperature = params['Common Temperature']

    materials.extend([ZrH, YHx])
    materials_database.update({'ZrH': ZrH, 'YHx': YHx})

    # """""""""""""""""""""
    # Sec. 1.3 : Coolants: NaK and Helium
    # """""""""""""""""""""
    
    NaK = openmc.Material(name="NaK", temperature=params['Common Temperature'])
    # Internal override used by density-aware temperature-coefficient cases.
    # Other MOUSE calculations retain the original 0.85 g/cm3 value.
    nak_density = float(params.get('_NaK Density Override', 0.85))
    NaK.set_density("g/cm3", nak_density)
    NaK.add_nuclide("Na23", 2.20000e-01)
    NaK.add_nuclide("K39", 7.27413e-01)
    NaK.add_nuclide("K41", 5.24956e-02)

    Helium = openmc.Material(name='Helium')
    # Internal override used by density-aware GCMR temperature-coefficient
    # cases. Other calculations retain the original 0.000166 g/cm3 value.
    helium_density = float(
        params.get('_Helium Density Override', 0.000166)
    )
    Helium.set_density('g/cm3', helium_density)
    Helium.temperature = params['Common Temperature']
    Helium.add_element('He', 1.0)
    
    materials.extend([NaK, Helium])
    materials_database.update({'NaK': NaK, 'Helium': Helium})

    # """""""""""""""""""""
    # Sec. 1.4 : Beryllium and Beryllium Oxide
    # """""""""""""""""""""
    Be = openmc.Material(name="Be")
    Be.add_element("beryllium", 1.0)
    Be.add_s_alpha_beta("c_Be")
    Be.set_density("g/cm3", 1.84)
    Be.temperature = params['Common Temperature']

    BeO = openmc.Material(name="BeO", temperature=params['Common Temperature'])
    BeO.set_density("g/cm3", 3.01)
    BeO.add_element("beryllium", 1.0)
    BeO.add_element("oxygen", 1.0)
    BeO.add_s_alpha_beta("c_Be_in_BeO")
    BeO.add_s_alpha_beta("c_O_in_BeO")
    

    materials.extend([Be, BeO])
    materials_database.update({'Be': Be, 'BeO': BeO})

    # """""""""""""""""""""
    # Sec. 1.5 : Zirconium
    # """""""""""""""""""""
    
    Zr = openmc.Material(name="Zr", temperature=params['Common Temperature'])
    Zr.set_density("g/cm3", 6.49)
    Zr.add_element("zirconium", 1.0)

    materials.append(Zr)
    materials_database.update({'Zr': Zr})
    
    # """""""""""""""""""""
    # Sec. 1.6 : SS304
    # """""""""""""""""""""
    
    SS304 = openmc.Material(name="SS304", temperature=params['Common Temperature'])
    SS304.set_density("g/cm3", 7.98)
    SS304.add_element("carbon", 0.04, "wo")
    SS304.add_element("silicon", 0.50, "wo")
    SS304.add_element("phosphorus", 0.023, "wo")
    SS304.add_element("sulfur", 0.015, "wo")
    SS304.add_element("chromium", 19.00, "wo")
    SS304.add_element("manganese", 1.00, "wo")
    SS304.add_element("iron", 70.173, "wo")
    SS304.add_element("nickel", 9.25, "wo")

    materials.append(SS304)
    materials_database.update({'SS304': SS304})
  
    # """""""""""""""""""""
    # Sec. 1.7 : Carbides: Boron Carbide and Silicon Carbide
    # """""""""""""""""""""  

    # Natural B4C
    B4C_natural = openmc.Material(name="B4C_natural", temperature=params['Common Temperature'])
    B4C_natural.add_element("boron", 4)
    B4C_natural.add_element("carbon", 1)
    B4C_natural.set_density("g/cm3", 2.52)

    # Enriched B4C
    B4C_enriched = openmc.Material(name="B4C_enriched", temperature=params['Common Temperature'])
    B4C_enriched.add_element("boron", 4, enrichment=95.0, enrichment_target='B10', enrichment_type='ao')
    B4C_enriched.add_element("carbon", 1)
    B4C_enriched.set_density("g/cm3", 2.52)

    SiC = openmc.Material(name='SiC')
    SiC.set_density('g/cm3', 3.18)
    SiC.add_element('Si', 0.5)
    SiC.add_element('C', 0.5)

    ZrC = openmc.Material(name='ZrC')
    ZrC.set_density('g/cm3', 6.73)
    ZrC.add_element('Zr', 1.0)
    ZrC.add_element('C', 1.0)

    materials.extend([B4C_natural, B4C_enriched, SiC])
    materials_database.update({'B4C_natural':  B4C_natural, 
                               'B4C_enriched': B4C_enriched, 
                               'SiC': SiC,
                               'ZrC': ZrC})

    # """""""""""""""""""""
    # Sec. 1.8 : Carbon Based Materials : Graphite (Buffer) & pyrolytic carbon (PyC) 
    # """""""""""""""""""""
   
    # Graphite
    Graphite = openmc.Material(
        name='Graphite',
        temperature=params['Common Temperature']
    )
    # Internal override used by density-aware GCMR temperature-coefficient
    # cases. Buffer graphite and PyC retain their fabrication densities.
    graphite_density = float(
        params.get('_Graphite Density Override', 1.60)
    )
    Graphite.set_density('g/cm3', graphite_density)
    Graphite.add_element('C', 1.0)
    # This adds thermal scattering data for graphite.
    Graphite.add_s_alpha_beta('c_Graphite')

    # Graphite of lower density (buffer graphite)
    buffer_graphite = openmc.Material(name='Buffer')
    buffer_graphite.set_density('g/cm3', 0.95)
    buffer_graphite.add_element('C', 1.0)
    buffer_graphite.add_s_alpha_beta('c_Graphite') 

    # Pyrolytic carbon (PyC)
    PyC = openmc.Material(name='PyC')
    PyC.set_density('g/cm3', 1.9)
    PyC.add_element('C', 1.0)
    PyC.add_s_alpha_beta('c_Graphite') 

    materials.extend([Graphite, buffer_graphite, PyC])
    materials_database.update({'Graphite': Graphite, 'buffer_graphite': buffer_graphite, 'PyC': PyC})

    # """""""""""""""""""""
    # Sec. 1.9 : Magnesium Oxide
    # """""""""""""""""""""
    MgO = openmc.Material(name='MgO')
    MgO.set_density('g/cm3', 3.58)
    MgO.add_element('Mg', 1.0)
    MgO.add_element('O', 1.0)
    materials_database.update({'MgO': MgO})

    # """""""""""""""""""""
    # Sec. 1.10 : Tungsten Based Materials: WB, W2B, WB4, WC
    # """""""""""""""""""""
    WB = openmc.Material(name='WB')
    WB.set_density('g/cm3', 15.43)
    WB.add_element('W', 1.0)
    WB.add_element('B', 1.0)

    W2B = openmc.Material(name='W2B')
    W2B.set_density('g/cm3', 16.75)  # doi.org/10.1016/j.jnucmat.2020.152062
    W2B.add_element('W', 2.0)
    W2B.add_element('B', 1.0)

    WB4 = openmc.Material(name='WB4')
    WB4.set_density('g/cm3', 8.23)
    WB4.add_element('W', 1.0)
    WB4.add_element('B', 4.0)

    WC = openmc.Material(name='WC')
    WC.set_density('g/cm3', 15.32)
    WC.add_element('W', 1.0)
    WC.add_element('C', 1.0)

    materials_database.update({'WB': WB, 'W2B': W2B, 'WB4': WB4, 'WC': WC})

    # """""""""""""""""""""
    # Sec. 1.11 : Heat Pipe Microreactor
    # """""""""""""""""""""
    
    # Homogenized heat pipe (SS316 + sodium mixture)
    heatpipe = openmc.Material(name='heatpipe')
    heatpipe.set_density('atom/b-cm', 2.74917E-02)
    heatpipe.temperature = params['Common Temperature']
    heatpipe.add_nuclide('Si28',  1.49701E-02, 'ao')
    heatpipe.add_nuclide('Si29',  7.60143E-04, 'ao')
    heatpipe.add_nuclide('Si30',  5.01090E-04, 'ao')
    heatpipe.add_nuclide('Cr50',  6.46763E-03, 'ao')
    heatpipe.add_nuclide('Cr52',  1.24724E-01, 'ao')
    heatpipe.add_nuclide('Cr53',  1.41423E-02, 'ao')
    heatpipe.add_nuclide('Cr54',  3.52029E-03, 'ao')
    heatpipe.add_nuclide('Mn55',  1.66133E-02, 'ao')
    heatpipe.add_nuclide('Fe54',  3.12186E-02, 'ao')
    heatpipe.add_nuclide('Fe56',  4.90061E-01, 'ao')
    heatpipe.add_nuclide('Fe57',  1.13180E-02, 'ao')
    heatpipe.add_nuclide('Fe58',  1.50617E-03, 'ao')
    heatpipe.add_nuclide('Ni58',  6.33738E-02, 'ao')
    heatpipe.add_nuclide('Ni60',  2.44119E-02, 'ao')
    heatpipe.add_nuclide('Ni61',  1.06115E-03, 'ao')
    heatpipe.add_nuclide('Ni62',  3.38338E-03, 'ao')
    heatpipe.add_nuclide('Ni64',  8.61654E-04, 'ao')
    heatpipe.add_nuclide('Mo92',  1.75699E-03, 'ao')
    heatpipe.add_nuclide('Mo94',  1.09514E-03, 'ao')
    heatpipe.add_nuclide('Mo95',  1.88484E-03, 'ao')
    heatpipe.add_nuclide('Mo96',  1.97478E-03, 'ao')
    heatpipe.add_nuclide('Mo97',  1.13066E-03, 'ao')
    heatpipe.add_nuclide('Mo98',  2.85681E-03, 'ao')
    heatpipe.add_nuclide('Mo100', 1.14011E-03, 'ao')
    heatpipe.add_nuclide('Na23',  1.79266E-01, 'ao')
   
    materials.append(heatpipe)
    materials_database.update({'heatpipe': heatpipe})

    # Monolith graphite
    monolith_graphite = openmc.Material(name='monolith_graphite')
    monolith_graphite.set_density('g/cm3', 1.63)
    monolith_graphite.temperature = params['Common Temperature']
    monolith_graphite.add_nuclide('C12', 0.9893, 'ao')
    monolith_graphite.add_nuclide('C13', 0.0107, 'ao')
    monolith_graphite.add_s_alpha_beta('c_Graphite')
    materials.append(monolith_graphite)
    materials_database.update({'monolith_graphite': monolith_graphite})

    # """""""""""""""""""""
    # Sec. 1.12 : HPMR VTB-specific materials
    # """""""""""""""""""""
    def add_nuclides(material, nuclides):
        for name, amount in nuclides:
            material.add_nuclide(name, amount, percent_type='ao')

    vtb_temperature = params['Common Temperature']

    vtb_air = openmc.Material(name='vtb_air')
    vtb_air.add_nuclide('He4', 1.0, percent_type='ao')
    vtb_air.set_density('g/cm3', 0.18e-3)
    vtb_air.temperature = vtb_temperature

    vtb_shell_ss = openmc.Material(name='vtb_shell_ss')
    add_nuclides(vtb_shell_ss, [
        ('C12', 1.9010E-03), ('Si28', 9.2693E-03), ('Si29', 4.7251E-04),
        ('Si30', 3.1166E-04), ('P31', 4.1322E-04), ('S32', 2.4710E-04),
        ('S33', 1.9511E-06), ('S34', 1.1056E-05), ('S36', 3.0016E-08),
        ('Cr50', 7.9116E-03), ('Cr52', 1.5257E-01), ('Cr53', 1.7300E-02),
        ('Cr54', 4.3063E-03), ('Mn55', 1.0280E-02), ('Fe54', 3.9029E-02),
        ('Fe56', 6.1213E-01), ('Fe57', 1.4144E-02), ('Fe58', 1.3343E-03),
        ('Ni58', 7.7516E-02), ('Ni60', 2.9859E-02), ('Ni61', 1.2981E-03),
        ('Ni62', 4.1389E-03), ('Ni64', 1.0544E-03), ('Mo92', 2.1259E-03),
        ('Mo94', 1.3299E-03), ('Mo95', 2.3030E-03), ('Mo96', 2.4191E-03),
        ('Mo97', 1.3903E-03), ('Mo98', 3.5249E-03), ('Mo100', 1.4135E-03),
    ])
    vtb_shell_ss.set_density('g/cm3', 7.90)
    vtb_shell_ss.temperature = vtb_temperature

    vtb_shell_air_mod = openmc.Material(name='vtb_shell_air_mod')
    add_nuclides(vtb_shell_air_mod, [
        ('He4', 1.983E-05), ('C12', 4.349E-05), ('Si28', 2.121E-04),
        ('Si29', 1.081E-05), ('Si30', 7.130E-06), ('P31', 9.454E-06),
        ('S32', 5.653E-06), ('S33', 4.464E-08), ('S34', 2.530E-07),
        ('S36', 6.867E-10), ('Cr50', 1.810E-04), ('Cr52', 3.491E-03),
        ('Cr53', 3.958E-04), ('Cr54', 9.852E-05), ('Mn55', 2.352E-04),
        ('Fe54', 8.929E-04), ('Fe56', 1.400E-02), ('Fe57', 3.236E-04),
        ('Fe58', 3.053E-05), ('Ni58', 1.773E-03), ('Ni60', 6.831E-04),
        ('Ni61', 2.970E-05), ('Ni62', 9.469E-05), ('Ni64', 2.412E-05),
        ('Mo92', 4.864E-05), ('Mo94', 3.043E-05), ('Mo95', 5.269E-05),
        ('Mo96', 5.535E-05), ('Mo97', 3.181E-05), ('Mo98', 8.065E-05),
        ('Mo100', 3.234E-05),
    ])
    vtb_shell_air_mod.set_density('atom/b-cm', 2.290E-02)
    vtb_shell_air_mod.temperature = vtb_temperature

    vtb_shell_air_hp = openmc.Material(name='vtb_shell_air_hp')
    add_nuclides(vtb_shell_air_hp, [
        ('He4', 5.629E-06), ('C12', 1.287E-04), ('Si28', 6.276E-04),
        ('Si29', 3.199E-05), ('Si30', 2.110E-05), ('P31', 2.798E-05),
        ('S32', 1.673E-05), ('S33', 1.321E-07), ('S34', 7.486E-07),
        ('S36', 2.032E-09), ('Cr50', 5.357E-04), ('Cr52', 1.033E-02),
        ('Cr53', 1.171E-03), ('Cr54', 2.916E-04), ('Mn55', 6.960E-04),
        ('Fe54', 2.643E-03), ('Fe56', 4.145E-02), ('Fe57', 9.576E-04),
        ('Fe58', 9.034E-05), ('Ni58', 5.248E-03), ('Ni60', 2.022E-03),
        ('Ni61', 8.789E-05), ('Ni62', 2.802E-04), ('Ni64', 7.139E-05),
        ('Mo92', 1.439E-04), ('Mo94', 9.004E-05), ('Mo95', 1.559E-04),
        ('Mo96', 1.638E-04), ('Mo97', 9.413E-05), ('Mo98', 2.387E-04),
        ('Mo100', 9.570E-05),
    ])
    vtb_shell_air_hp.set_density('atom/b-cm', 6.771E-02)
    vtb_shell_air_hp.temperature = vtb_temperature

    vtb_shell_air_center = openmc.Material(name='vtb_shell_air_center')
    add_nuclides(vtb_shell_air_center, [
        ('He4', 2.68E-05), ('C12', 1.62E-06), ('Si28', 7.92E-06),
        ('Si29', 4.04E-07), ('Si30', 2.66E-07), ('P31', 3.53E-07),
        ('S32', 2.11E-07), ('S33', 1.67E-09), ('S34', 9.45E-09),
        ('S36', 2.57E-11), ('Cr50', 6.76E-06), ('Cr52', 1.30E-04),
        ('Cr53', 1.48E-05), ('Cr54', 3.68E-06), ('Mn55', 8.79E-06),
        ('Fe54', 3.34E-05), ('Fe56', 5.23E-04), ('Fe57', 1.21E-05),
        ('Fe58', 1.14E-06), ('Ni58', 6.63E-05), ('Ni60', 2.55E-05),
        ('Ni61', 1.11E-06), ('Ni62', 3.54E-06), ('Ni64', 9.01E-07),
        ('Mo92', 1.82E-06), ('Mo94', 1.14E-06), ('Mo95', 1.97E-06),
        ('Mo96', 2.07E-06), ('Mo97', 1.19E-06), ('Mo98', 3.01E-06),
        ('Mo100', 1.21E-06),
    ])
    vtb_shell_air_center.set_density('atom/b-cm', 8.815E-04)
    vtb_shell_air_center.temperature = vtb_temperature

    vtb_potassium_vapor = openmc.Material(name='vtb_potassium_vapor')
    add_nuclides(vtb_potassium_vapor, [('K39', 0.93258), ('K40', 0.00012), ('K41', 0.06730)])
    vtb_potassium_vapor.set_density('g/cm3', 1.11e-4)
    vtb_potassium_vapor.temperature = vtb_temperature

    vtb_potassium_liquid = openmc.Material(name='vtb_potassium_liquid')
    add_nuclides(vtb_potassium_liquid, [('K39', 0.93258), ('K40', 0.00012), ('K41', 0.06730)])
    vtb_potassium_liquid.set_density('g/cm3', 0.705)
    vtb_potassium_liquid.temperature = vtb_temperature

    vtb_wick = openmc.Material(name='vtb_wick')
    add_nuclides(vtb_wick, [
        ('C12', 5.589E-04), ('Si28', 2.725E-03), ('Si29', 1.389E-04),
        ('Si30', 9.163E-05), ('P31', 1.215E-04), ('S32', 7.265E-05),
        ('S33', 5.736E-07), ('S34', 3.250E-06), ('S36', 8.825E-09),
        ('Cr50', 2.326E-03), ('Cr52', 4.485E-02), ('Cr53', 5.086E-03),
        ('Cr54', 1.266E-03), ('Mn55', 3.022E-03), ('Fe54', 1.147E-02),
        ('Fe56', 1.800E-01), ('Fe57', 4.158E-03), ('Fe58', 3.923E-04),
        ('Ni58', 2.279E-02), ('Ni60', 8.778E-03), ('Ni61', 3.816E-04),
        ('Ni62', 1.217E-03), ('Ni64', 3.100E-04), ('Mo92', 6.250E-04),
        ('Mo94', 3.910E-04), ('Mo95', 6.771E-04), ('Mo96', 7.112E-04),
        ('Mo97', 4.087E-04), ('Mo98', 1.036E-03), ('Mo100', 4.156E-04),
        ('K39', 6.584E-01), ('K40', 8.472E-05), ('K41', 4.751E-02),
    ])
    vtb_wick.set_density('g/cm3', 2.753)
    vtb_wick.temperature = vtb_temperature

    vtb_hp_vapor_liquid_wick = openmc.Material(name='vtb_hp_vapor_liquid_wick')
    add_nuclides(vtb_hp_vapor_liquid_wick, [
        ('C12', 3.808E-06), ('Si28', 1.856E-05), ('Si29', 9.463E-07),
        ('Si30', 6.242E-07), ('P31', 8.277E-07), ('S32', 4.949E-07),
        ('S33', 3.908E-09), ('S34', 2.214E-08), ('S36', 6.012E-11),
        ('Cr50', 1.585E-05), ('Cr52', 3.055E-04), ('Cr53', 3.465E-05),
        ('Cr54', 8.625E-06), ('Mn55', 2.059E-05), ('Fe54', 7.814E-05),
        ('Fe56', 1.226E-03), ('Fe57', 2.833E-05), ('Fe58', 2.673E-06),
        ('Ni58', 1.553E-04), ('Ni60', 5.980E-05), ('Ni61', 2.600E-06),
        ('Ni62', 8.291E-06), ('Ni64', 2.112E-06), ('Mo92', 4.258E-06),
        ('Mo94', 2.664E-06), ('Mo95', 4.613E-06), ('Mo96', 4.845E-06),
        ('Mo97', 2.784E-06), ('Mo98', 7.058E-06), ('Mo100', 2.831E-06),
        ('K39', 5.895E-03), ('K40', 7.586E-07), ('K41', 4.254E-04),
    ])
    vtb_hp_vapor_liquid_wick.set_density('atom/b-cm', 8.324E-03)
    vtb_hp_vapor_liquid_wick.temperature = vtb_temperature

    vtb_yh_moderator = openmc.Material(name='vtb_yh_moderator')
    add_nuclides(vtb_yh_moderator, [('Y89', 0.357142857), ('H1', 0.642857143)])
    vtb_yh_moderator.add_s_alpha_beta('c_H_in_YH2')
    vtb_yh_moderator.add_s_alpha_beta('c_Y_in_YH2')
    vtb_yh_moderator.set_density('g/cm3', 4.0850)
    vtb_yh_moderator.temperature = vtb_temperature

    vtb_buffer = openmc.Material(name='vtb_buffer')
    vtb_buffer.add_nuclide('C12', 1.0, percent_type='ao')
    vtb_buffer.add_s_alpha_beta('c_Graphite')
    vtb_buffer.set_density('g/cm3', 1.0400)
    vtb_buffer.temperature = vtb_temperature

    vtb_PyC = openmc.Material(name='vtb_PyC')
    vtb_PyC.add_nuclide('C12', 1.0, percent_type='ao')
    vtb_PyC.add_s_alpha_beta('c_Graphite')
    vtb_PyC.set_density('g/cm3', 1.8820)
    vtb_PyC.temperature = vtb_temperature

    vtb_SiC = openmc.Material(name='vtb_SiC')
    add_nuclides(vtb_SiC, [('Si28', 0.4611), ('Si29', 0.0234), ('Si30', 0.0154), ('C12', 0.5)])
    vtb_SiC.set_density('g/cm3', 3.1710)
    vtb_SiC.temperature = vtb_temperature

    vtb_matrix_graphite = openmc.Material(name='vtb_matrix_graphite')
    add_nuclides(vtb_matrix_graphite, [('C12', 0.9999997), ('B10', 0.0000003)])
    vtb_matrix_graphite.add_s_alpha_beta('c_Graphite')
    vtb_matrix_graphite.set_density('g/cm3', 1.8060)
    vtb_matrix_graphite.temperature = vtb_temperature

    vtb_beryllium = openmc.Material(name='vtb_beryllium')
    vtb_beryllium.add_nuclide('Be9', 1.0, percent_type='ao')
    vtb_beryllium.add_s_alpha_beta('c_Be')
    vtb_beryllium.set_density('g/cm3', 1.848)
    vtb_beryllium.temperature = vtb_temperature

    b10_at_frac = params.get('B10_at_frac_B', 0.95)
    vtb_B4C_drum = openmc.Material(name='vtb_B4C_drum')
    vtb_B4C_drum.add_nuclide('B10', b10_at_frac * 0.8, percent_type='ao')
    vtb_B4C_drum.add_nuclide('B11', (1 - b10_at_frac) * 0.8, percent_type='ao')
    vtb_B4C_drum.add_nuclide('C12', 0.2, percent_type='ao')
    vtb_B4C_drum.set_density('g/cm3', 2.510)
    vtb_B4C_drum.temperature = vtb_temperature

    vtb_B4C_central = openmc.Material(name='vtb_B4C_central')
    add_nuclides(vtb_B4C_central, [('B10', 0.76), ('B11', 0.04), ('C12', 0.2)])
    vtb_B4C_central.set_density('g/cm3', 1.25)
    vtb_B4C_central.temperature = vtb_temperature

    vtb_materials = [
        vtb_air,
        vtb_shell_ss,
        vtb_shell_air_mod,
        vtb_shell_air_hp,
        vtb_shell_air_center,
        vtb_potassium_vapor,
        vtb_potassium_liquid,
        vtb_wick,
        vtb_hp_vapor_liquid_wick,
        vtb_yh_moderator,
        vtb_buffer,
        vtb_PyC,
        vtb_SiC,
        vtb_matrix_graphite,
        vtb_beryllium,
        vtb_B4C_drum,
        vtb_B4C_central,
    ]
    materials.extend(vtb_materials)
    materials_database.update({
        'vtb_air': vtb_air,
        'vtb_shell_ss': vtb_shell_ss,
        'vtb_shell_air_mod': vtb_shell_air_mod,
        'vtb_shell_air_hp': vtb_shell_air_hp,
        'vtb_shell_air_center': vtb_shell_air_center,
        'vtb_potassium_vapor': vtb_potassium_vapor,
        'vtb_potassium_liquid': vtb_potassium_liquid,
        'vtb_wick': vtb_wick,
        'vtb_hp_vapor_liquid_wick': vtb_hp_vapor_liquid_wick,
        'vtb_yh_moderator': vtb_yh_moderator,
        'vtb_buffer': vtb_buffer,
        'vtb_PyC': vtb_PyC,
        'vtb_SiC': vtb_SiC,
        'vtb_matrix_graphite': vtb_matrix_graphite,
        'vtb_beryllium': vtb_beryllium,
        'vtb_B4C_drum': vtb_B4C_drum,
        'vtb_B4C_central': vtb_B4C_central,
    })

    return materials_database
