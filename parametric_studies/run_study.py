import numpy as np
import itertools
np.random.seed(23)
import pandas as pd
import joblib
from pathlib import Path
import os,sys
import string,random
import pathos
from shutil import copyfile
# Add relative path to sys.path to be able to call the external scripts as import
# ---
sys.path.append(os.path.dirname(os.path.abspath('.')))
#print(os.path.dirname(os.path.abspath('.')),sys.path)
#
from OpenMC_GCMR import OpenMC_GCMR
from OpenMC_HPMR import OpenMC_HPMR
from OpenMC_LTMR import OpenMC_LTMR

def rstring(length):
    characters = string.ascii_letters + string.digits
    random_string = ''.join(random.choice(characters) for i in range(length))
    return random_string

def sample_search_space(params,rundir):
    """
    Sample the search space with modified params
    ---

    :param params: (dict) dictionary of parameters to change
    :param rundir: (str) /path/to/run sub-directory in which to run all cases
    
    Returns
    0:
    """
    mode = sys.argv[2]
    if mode == 'GCMR':
        gcmr = OpenMC_GCMR(working_dir=Path(rundir))
        if params == None:
            dummy = gcmr.fitness(gcmr.params,openmc_run=True)
        else:
            dummy = gcmr.fitness(params,openmc_run=True)
    elif mode == 'HPMR':
        hpmr = OpenMC_HPMR(working_dir=Path(rundir))
        if params == None:
            dummy = hpmr.fitness(hpmr.params,openmc_run=True)
        else:
            dummy = hpmr.fitness(params,openmc_run=True)
    elif mode == 'LTMR':
        ltmr = OpenMC_LTMR(working_dir=Path(rundir))
        if params == None:
            dummy = ltmr.fitness(ltmr.params,openmc_run=True)
        else:
            dummy = ltmr.fitness(params,openmc_run=True)
    return 0

if __name__ == "__main__":

    study_type = sys.argv[1]
    assert study_type in ['reflector','booster','influence_yttrium_hydride','fuel'], '---error: study_type is not valid.'
    if study_type == 'reflector':
        if sys.argv[2] == 'GCMR':
            reflectors = ['Graphite']#['Graphite','BeO','Be']
            reflector_thickness = [27.393]#[2,12.5,27.393,30,40]#[2,5,7.5,10,12.5,15,20,23,27.393,30]
            params = list(itertools.product(reflectors,reflector_thickness))
            for count,_ in enumerate(params):
                sample_search_space({'Control Drum Reflector':_[0],'Reflector': _[0],'Reflector Thickness':_[1]},rundir=study_type+'_'+sys.argv[2])  
        elif sys.argv[2] == 'HPMR':
            reflectors = ['Graphite','BeO','Be']#['Graphite']
            reflector_thickness = [25]#[30,35,40,45,50,55,60]#[50]
            params = list(itertools.product(reflectors,reflector_thickness))
            for count,_ in enumerate(params):
                sample_search_space({'Control Drum Reflector':_[0],'Reflector': _[0],'Reflector Thickness':_[1]},rundir=study_type+'_'+sys.argv[2])      
        else:
            reflectors = ['Graphite','BeO','Be']
            reflector_thickness = [13,14,20,23,27.393,30,35,40,45,50]
            params = list(itertools.product(reflectors,reflector_thickness))
            for count,_ in enumerate(params):
                sample_search_space({'Control Drum Reflector':_[0],'Reflector': _[0],'Reflector Thickness':_[1]},rundir=study_type+'_'+sys.argv[2])  
    elif study_type == 'fuel':
        if sys.argv[2] == 'GCMR':
            fuel_types = ['UN','UO2']
            enrichments = np.linspace(0.05,0.10,5)
            power = [10,12.5,15]
            rings_per_assembly = [6]
            params = list(itertools.product(fuel_types,enrichments,power,rings_per_assembly))
            for _ in params:
                sample_search_space({'Fuel Pin Materials': [_[0], 'buffer_graphite', 'PyC', 'SiC', 'PyC'],'Fuel':_[0],'Enrichment': _[1],'Power MWt':_[2],'Assembly Rings':_[3]},rundir=study_type+'_'+sys.argv[2])
        elif sys.argv[2] == 'HPMR':
            enrichments = [0.10,0.125,0.15,0.17]#np.linspace(0.05,0.10,5)
            power = [5,7,9]
            rings_per_assembly = [6,8,10]#6
            reflector = 'Graphite'
            params = list(itertools.product(enrichments,power,rings_per_assembly))
            for _ in params:
                sample_search_space({'Enrichment': _[0],'Control Drum Reflector':reflector,'Reflector':reflector,'Power MWt':_[1],'Number of Rings per Assembly':_[2]},rundir=study_type+'_'+sys.argv[2])
        else:
            fuel_types = ['TRIGA_fuel']
            enrichments = [0.0875,0.10]
            power = [15,17.5,20]
            rings_per_assembly = [12]
            reflector_thickness = [14,20,30,40,50]
            u_met_wo = [0.30,0.35,0.40,0.45]
            params = list(itertools.product(fuel_types,enrichments,power,rings_per_assembly,reflector_thickness,u_met_wo))
            for _ in params:
                sample_search_space({'U_met_wo': _[5],'Fuel Pin Materials': ['Zr', None, _[0], None, 'SS304'],'Fuel':_[0],'Enrichment': _[1],'Power MWt':_[2],'Number of Rings per Assembly':_[3],'Reflector Thickness':_[4]},rundir=study_type+'_'+sys.argv[2])
    elif 'booster' in study_type:
        if sys.argv[2] == 'GCMR':
            booster_type = ['ZrH','YHx']
            booster_radius = [0.45,0.50,0.55,0.60,0.65]
            params = list(itertools.product(booster_type,booster_radius))
            for _ in params:
                sample_search_space({'Moderator Booster': _[0],'Moderator Booster Radius':_[1]},rundir=study_type+'_'+sys.argv[2])
        else:
            booster_type = ['ZrH','YHx']
            delta = 1.5875 - 1.5367
            booster_radius = [1.5367, 1.75 ,2.0]
            params = list(itertools.product(booster_type,booster_radius))
            for _ in params:
                sample_search_space({'Moderator': _[0],'Moderator Pin Materials': [_[0], 'SS304'],'Moderator Pin Inner Radius': _[1],'Moderator Pin Radii': [_[1], _[1]+delta]},rundir=study_type+'_'+sys.argv[2])
    elif study_type == 'influence_yttrium_hydride':
        if sys.argv[2] == 'GCMR':
            gcmr = OpenMC_GCMR(working_dir=Path(study_type+'_'+sys.argv[2]))
            path_to_cost = '/home/seurpr/Desktop/LDRD/LDRD_microreactor/MOUSE_dev/mouse/parametric_studies/var_cost_for_Yh/'
            copyfile("{}/Cost_Database_1.xlsx".format(path_to_cost),"{}/Cost_Database.xlsx".format(path_to_cost))
            gcmr.fitness(parameters={'Moderator Booster': 'YHx','cost database':'/home/seurpr/Desktop/LDRD/LDRD_microreactor/MOUSE_dev/mouse/parametric_studies/var_cost_for_Yh/'},openmc_run=True)
            for _ in range(1,10):
                folderindex_new = rstring(10) # update sample
                sdir = Path(gcmr.working_dir / f"sample_{folderindex_new}")
                sdir.mkdir(parents = True, exist_ok = True)
                gcmr.folderindex = folderindex_new
                copyfile("{}/Cost_Database_{}.xlsx".format(path_to_cost,_+1),"{}/Cost_Database.xlsx".format(path_to_cost))
                gcmr.fitness(parameters={'Number of Samples':500},openmc_run=False) # only perform cost sampling but do not run openmc
        else:
            ltmr = OpenMC_LTMR(working_dir=Path(study_type+'_'+sys.argv[2]))
            path_to_cost = '/home/seurpr/Desktop/LDRD/LDRD_microreactor/MOUSE_dev/mouse/parametric_studies/var_cost_for_Yh/'
            copyfile("{}/Cost_Database_1.xlsx".format(path_to_cost),"{}/Cost_Database.xlsx".format(path_to_cost))
            ltmr.fitness(parameters={'Moderator': 'YHx','Moderator Pin Materials': ['YHx', 'SS304'],'cost database':'/home/seurpr/Desktop/LDRD/LDRD_microreactor/MOUSE_dev/mouse/parametric_studies/var_cost_for_Yh/'},openmc_run=True)
            for _ in range(1,10):
                folderindex_new = rstring(10) # update sample
                sdir = Path(ltmr.working_dir / f"sample_{folderindex_new}")
                sdir.mkdir(parents = True, exist_ok = True)
                ltmr.folderindex = folderindex_new
                copyfile("{}/Cost_Database_{}.xlsx".format(path_to_cost,_+1),"{}/Cost_Database.xlsx".format(path_to_cost))
                ltmr.fitness(parameters={'Number of Samples':500},openmc_run=False) # only perform cost sampling but do not run openmc
