# Importing libraries
import openmc
from core_design.utils import list_to_dict

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

    # The fuel is declared following a mechanism that is close to the logic of the
    # fabrication of the actual TRIGA fuel. This type of fuel is specified for example
    # as "45/20" fuel, which implies that 45% weight is composed of Uranium (metal)
    # that is 20% enriched in a matrix of ZrH. TRIGA fuel can also contain 3% weight of
    # Erbium as a burnable absorber in the fuel meat.

    # First let's declare the individual components of the fuel
    try:
        U_met = openmc.Material(name="U_met")
        U_met.set_density("g/cm3", 19.05)

        U_met.add_nuclide("U235", params['enrichment'])
        U_met.add_nuclide("U238", 1 - params['enrichment'])

        ZrH_fuel = openmc.Material(name="ZrH_fuel")
        ZrH_fuel.set_density("g/cm3", 5.63)

        ZrH_fuel.add_element("zirconium", 1.0)
        ZrH_fuel.add_nuclide("H1", params["H_Zr_ratio"]) # The proportion of hydrogen atoms relative to zirconium atoms


        # Now we mix the components in the right weight %
        # The NRC seems to license up to TRIGA fuel with up to 45% weight Uranium

        TRIGA_fuel = openmc.Material.mix_materials(
            [U_met, ZrH_fuel], [params['U_met_wo'], 1 - params['U_met_wo']], "wo", name="UZrH"
        )

        TRIGA_fuel.temperature = params['common_temperature']
        TRIGA_fuel.add_s_alpha_beta("c_H_in_ZrH")
        materials.append(TRIGA_fuel)
        materials_database.update({ 'TRIGA_fuel': TRIGA_fuel})
    
    except KeyError as e:
        print(f"Skipping TRIGA_fuel due to missing parameter: {e}")    
    
    # # Let's also make a version with 3% Erbium in the meat
    # Er = openmc.Material(name="Er", temperature= params['common_temperature'])
    # Er.set_density("g/cm3", 9.2)
    # Er.add_element("erbium", 1.0)


    #UO2
    try:
        uo2 = openmc.Material(name='UO2')
        uo2.set_density('g/cm3', 10.41)
        uo2.add_element('U', 1.0, enrichment= 100 * params['enrichment'])
        uo2.add_nuclide('O16', 2.0)
        materials.append(uo2)
        materials_database.update({ 'uo2': uo2})
    except KeyError as e:
        print(f"Skipping UO2 due to missing parameter: {e}")     

    
    # Uranium Carbide
    try:
        uc = openmc.Material(name='UC')
        uc.set_density('g/cm3', 13.0)
        uc.add_element('U', 1.0,  enrichment= 100 * params['enrichment'])
        uc.add_element('N', 1.0)
        materials.append(uc)
        materials_database.update({ 'uc': uc})
    except KeyError as e:
        print(f"Skipping UC due to missing parameter: {e}")    


    #Mixed uranium dioxide (UO2) and uranium carbide (UC)
    try:
        uco = openmc.Material.mix_materials([uo2, uc], [params['UO2 atom fraction'], 1- params['UO2 atom fraction']], 'ao') # mixing UO2 and UC by atom fraction
        materials.append(uco)
        materials_database.update({ 'uco': uco})
    except KeyError as e:
        print(f"Skipping UCO due to missing parameter: {e}") 
    
    # Uranium Nitride
    try:
        un = openmc.Material(name='UN') # This creates a new material named 'UN'.
        un.set_density('g/cm3', 14.0)
        un.add_element('U', 1.0, enrichment=100 * params['enrichment'])
        un.add_element('N', 1.0) # This adds nitrogen (N) to the material.
        materials.append(un)
        materials_database.update({ 'un': un})
    except KeyError as e:
        print(f"Skipping UN due to missing parameter: {e}") 


    # """""""""""""""""""""
    # Sec. 1.2 : Hydrides: Zirconium Hydride and yttrium hydride (YHx)
    # """""""""""""""""""""
       
    ZrH = openmc.Material(name="ZrH", temperature= params['common_temperature'])
    ZrH.set_density("g/cm3", 5.6)

    ZrH.add_nuclide("H1", 1.85)
    ZrH.add_element("zirconium", 1.0)
    ZrH.add_s_alpha_beta("c_H_in_ZrH")

    #Yttrium hydride (YHx) 
    YHx = openmc.Material(name="YHx")
    YHx.set_density("g/cm3", 4.28)
    YHx.add_nuclide("H1", 1.5) # This adds hydrogen-1 (H-1) nuclide to the material with an atomic ratio of 1.5.
    YHx.add_element("yttrium", 1.0) #  This adds yttrium to the material with an atomic ratio of 1.0.
    # Adding thermal scattering data for hydrogen in yttrium hydride (YH2). 
    # The add_s_alpha_beta method is used to specify the S(α,β) thermal scattering treatment for specific materials. 
    YHx.add_s_alpha_beta("c_H_in_YH2")
    YHx.temperature = params['common_temperature']

    materials.extend([ZrH, YHx])
    materials_database.update({ 'ZrH': ZrH, 'YHx' :YHx})

    # """""""""""""""""""""
    # Sec. 1.3 : Coolants: NaK and Helium
    # """""""""""""""""""""
    
    NaK = openmc.Material(name="NaK", temperature= params['common_temperature'])
    NaK.set_density("g/cm3", 0.75)
    NaK.add_nuclide("Na23", 2.20000e-01)
    NaK.add_nuclide("K39", 7.27413e-01)
    NaK.add_nuclide("K41", 5.24956e-02)

    helium = openmc.Material(name='Helium')
    helium.set_density('g/cm3', 0.000166)
    helium.add_element('He', 1.0)
    
    materials.extend([NaK, helium])
    materials_database.update({ 'NaK': NaK, 'helium' :helium})

    # """""""""""""""""""""
    # Sec. 1.4 : Beryllium and Beryllium Oxide
    # """""""""""""""""""""
    Be = openmc.Material(name="Be")
    Be.add_element("beryllium", 1.0)
    Be.add_s_alpha_beta("c_Be")
    Be.set_density("g/cm3", 1.84)
    Be.temperature =  params['common_temperature']

    BeO = openmc.Material(name="BeO", temperature= params['common_temperature'])
    BeO.set_density("g/cm3", 3.01)
    BeO.add_element("beryllium", 1.0)
    BeO.add_element("oxygen", 1.0)
    BeO.add_s_alpha_beta("c_Be_in_BeO")

    materials.extend([Be, BeO])
    materials_database.update({ 'Be': Be, 'BeO': BeO})

    # """""""""""""""""""""
    # Sec. 1.5 : Zirconium
    # """""""""""""""""""""
    
    Zr = openmc.Material(name="Zr", temperature= params['common_temperature'])
    Zr.set_density("g/cm3", 6.49)
    Zr.add_element("zirconium", 1.0)

    materials.append(Zr)
    materials_database.update({'Zr': Zr})
    
    # """""""""""""""""""""
    # Sec. 1.6 : SS304
    # """""""""""""""""""""
    
    SS304 = openmc.Material(name="SS304", temperature= params['common_temperature'])
    SS304.set_density("g/cm3", 7.98)
    SS304.add_element("carbon", 0.04)
    SS304.add_element("silicon", 0.50)
    SS304.add_element("phosphorus", 0.023)
    SS304.add_element("sulfur", 0.015)
    SS304.add_element("chromium", 19.00)
    SS304.add_element("manganese", 1.00)
    SS304.add_element("iron", 70.173)
    SS304.add_element("nickel", 9.25)

    materials.append(SS304)
    materials_database.update({'SS304': SS304})
  
    # """""""""""""""""""""
    # Sec. 1.7 : Carbides: Boron Carbide and Silicon Carbide
    # """""""""""""""""""""  

    # Natural B4C
    B4C_natural = openmc.Material(name="B4C_natural", temperature= params['common_temperature'])
    B4C_natural.add_element("boron", 4)
    B4C_natural.add_element("carbon", 1)
    B4C_natural.set_density("g/cm3", 2.52)

    SiC = openmc.Material(name='SiC') # Creates a new material named 'SiC'.
    SiC.set_density('g/cm3', 3.18)
    SiC.add_element('Si', 0.5) #  Adds silicon (Si) to the material with a fraction of 0.5.
    SiC.add_element('C', 0.5) # Adds carbon (C) to the material with a fraction of 0.5.

    materials.extend([B4C_natural, SiC])
    materials_database.update({'B4C_natural': B4C_natural, 'SiC': SiC})

    # """""""""""""""""""""
    # Sec. 1.8 : Carbon Based Materials : Graphite (Buffer)  & pyrolytic carbon (PyC) 
    # """""""""""""""""""""
   
    # Graphite
    Graphite = openmc.Material(name='Graphite')
    Graphite.set_density('g/cm3', 1.7)
    Graphite.add_element('C', 1.0)
    #This adds thermal scattering data for graphite. The add_S(α,β) thermal scattering treatment for carbon in graphite form.
    Graphite.add_s_alpha_beta('c_Graphite')

    # Graphite of lower density (buffer)
    buffer = openmc.Material(name='Buffer') # This creates a new material named 'Buffer'.
    buffer.set_density('g/cm3', 0.95)
    buffer.add_element('C', 1.0) #This adds carbon (C) 
    #This adds thermal scattering data for graphite. The add_S(α,β) thermal scattering treatment for carbon in graphite form.
    buffer.add_s_alpha_beta('c_Graphite') 

    # pyrolytic carbon (PyC) 
    PyC = openmc.Material(name='PyC') # Creates a new material named 'PyC'.
    PyC.set_density('g/cm3', 1.9)
    PyC.add_element('C', 1.0) #  Adds carbon (C) to the material with a fraction of 1.0
    
    #This adds thermal scattering data for graphite. The add_S(α,β) thermal scattering treatment for carbon in graphite form.
    PyC.add_s_alpha_beta('c_Graphite') 
    materials.extend([Graphite, buffer, PyC ])
    materials_database.update({'Graphite' : Graphite, 'buffer' : buffer, 'PyC': PyC})

    # I am not sure if this is necessary but keeping it for now
    # materials.export_to_xml()

    
    # colors = {materials_database['Zr']: 'green',
    #         materials_database['SS304']: 'pink',
    #         materials_database['NaK']: 'blue',
    #         materials_database['TRIGA_fuel']: 'red',
    #         materials_database['ZrH']: 'orange',
    #         materials_database['Be']: 'moccasin',
    #         materials_database['BeO']: 'seagreen',
    #         materials_database['B4C_nat']: 'black'}
    
    return materials_database


