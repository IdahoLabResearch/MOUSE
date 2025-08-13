# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import numpy as np
def fuel_calculations(params):
    
    U_mass = (params['Mass U235'] + params['Mass U238']) / 1000 # The mass of Uranium only (in Kg)
    nat_u_consum = U_mass*(params['Enrichment'] -0.0025)/(0.0071-0.0025) # Kg
    tail_waste = nat_u_consum - U_mass # Kg

    # value functions
    f_val_fun = (1-2*params['Enrichment'])*np.log((1-params['Enrichment'])/params['Enrichment'])
    tail_waste_val_fun = 5.96
    nat_u_waste_val_fun = 4.87


    kg_SWU = (U_mass*f_val_fun+tail_waste*tail_waste_val_fun- nat_u_consum *nat_u_waste_val_fun)
    
    params['Natural Uranium Mass'] = nat_u_consum
    params['Fuel Tail Waste Mass'] = tail_waste # Kg
    params['SWU'] = kg_SWU
 