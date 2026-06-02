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
    ZrH.set_density("g/cm3", 5.6)
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
    NaK.set_density("g/cm3", 0.85)
    NaK.add_nuclide("Na23", 2.20000e-01)
    NaK.add_nuclide("K39", 7.27413e-01)
    NaK.add_nuclide("K41", 5.24956e-02)

    Helium = openmc.Material(name='Helium')
    Helium.set_density('g/cm3', 0.000166)
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
    Graphite = openmc.Material(name='Graphite')
    Graphite.set_density('g/cm3', 1.60)
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

    return materials_database