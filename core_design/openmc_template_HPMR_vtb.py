# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
# Importing libraries
import openmc
import random
import itertools
import numpy as np

from core_design.openmc_materials_database import collect_materials_data

def drum_parameters_to_eq(radius, coating_angle, drum_angle, p0):
    """convert control drum parameters to A, B, C & D intersecting plane parameters"""
    coating_angle = coating_angle*np.pi/180
    drum_angle = drum_angle*np.pi/180
    l = radius*np.cos(coating_angle / 2)
    point = np.array([[np.cos(drum_angle), -np.sin(drum_angle)],[np.sin(drum_angle), np.cos(drum_angle)]])@np.array([l, 0]) + p0
    norm = np.array([np.cos(drum_angle), np.sin(drum_angle)])
    A, B = norm
    C = 0
    D = A*point[0] + B*point[1]
    if p0[0] > 0 or (np.isclose(p0[0], 0) and p0[1] > 0): #takes care of positivity of plane normal
        return A,B,C,D
    else:
        return -A,-B,-C,-D



def _normalize_vtb_params(params):
    active_height = params.get("active_fuel_height", params.get("Active Height"))
    axial_reflector = params.get("axial_reflector_height", params.get("Axial Reflector Thickness"))
    core_radius = params.get("Core Radius")
    reflector_width = params.get("reflector_width")
    if reflector_width is None and core_radius is not None:
        reflector_width = 2.0 * core_radius

    vtb_params = {
        "compact_radius": params.get("compact_radius", params.get("Compact Fuel Radius", params.get("Fuel Pin Radii", [None])[0])),
        "moderator_radius": params.get("moderator_radius", params.get("Moderator Radius", params.get("Moderator Booster Raddi"))),
        "coating_angle": params.get("coating_angle", params.get("Coating Angle", 90.0)),
        "B10_at_frac_B": params.get("B10_at_frac_B", params.get("B10 at frac B", 0.95)),
        "flake_width": params.get("flake_width", params.get("Assembly FTF")),
        "pin_pitch": params.get("pin_pitch", params.get("Lattice Pitch")),
        "enrichment": params.get("enrichment", params.get("Enrichment")),
        "reflector_width": reflector_width,
        "active_fuel_height": active_height,
        "axial_reflector_height": axial_reflector,
    }

    missing = [key for key, value in vtb_params.items() if value is None]
    if missing:
        raise ValueError(f"HPMR VTB is missing required parameters: {missing}")

    return {key: float(value) for key, value in vtb_params.items()}


def build_openmc_model_HPMR_vtb(params):
    params.setdefault("Shutdown Margin Calc", False)
    params.setdefault("Isothermal Temperature Coefficients", False)
    params.setdefault("VTB Axial Divisions", 11)
    params.setdefault("VTB Tally Mesh Dimension", (200, 200, 20))
    params.setdefault("VTB Low Packing Fraction", False)
    params.setdefault("VTB Drum Rotation", 180.0)

    vtb_params = _normalize_vtb_params(params)
    particles = int(params.get("Particles", 1000))
    batches = int(params.get("Batches", 100))
    inactive = int(params.get("Inactive", 50))
    axial_divs = int(params["VTB Axial Divisions"])
    low_pf = bool(params["VTB Low Packing Fraction"])
    drum_rotation = 0.0 if params["Shutdown Margin Calc"] else float(params["VTB Drum Rotation"])
    tally_mesh_dimension = params.get("VTB Tally Mesh Dimension")

    _export_openmc_model_HPMR_vtb(
        vtb_params,
        params,
        particles,
        batches,
        inactive,
        axial_divs,
        low_pf,
        drum_rotation,
        tally_mesh_dimension,
    )


def _export_openmc_model_HPMR_vtb(parms, params, particles, batches, inactive, axial_divs, low_pf, drum_rotation, tally_mesh_dimension):
    materials_database = collect_materials_data(params)

    air = materials_database["vtb_air"]
    shell_mod = materials_database["vtb_shell_ss"]
    shell_hp = shell_mod
    shell_air_mod = materials_database["vtb_shell_air_mod"]
    shell_air_hp = materials_database["vtb_shell_air_hp"]
    shell_air_center = materials_database["vtb_shell_air_center"]
    coo_vap = materials_database["vtb_potassium_vapor"]
    coo_liq = materials_database["vtb_potassium_liquid"]
    wick = materials_database["vtb_wick"]
    hp_vp_liq_wick = materials_database["vtb_hp_vapor_liquid_wick"]
    moderator = materials_database["vtb_yh_moderator"]
    fuel = materials_database[params.get("VTB Fuel Material", params.get("Fuel", "UCO"))]
    buffer = materials_database["vtb_buffer"]
    PyC1 = materials_database["vtb_PyC"]
    PyC2 = PyC1
    SiC = materials_database["vtb_SiC"]
    matrix_pin = materials_database["vtb_matrix_graphite"]
    matrix = matrix_pin
    beryllium = materials_database["vtb_beryllium"]
    beryllium_drum = beryllium
    B4C_drum = materials_database["vtb_B4C_drum"]
    B4C_central = materials_database["vtb_B4C_central"]

    fuel.depletable = True
    base_materials = [
        air,
        shell_mod,
        shell_air_mod,
        shell_air_hp,
        shell_air_center,
        coo_vap,
        coo_liq,
        wick,
        hp_vp_liq_wick,
        moderator,
        fuel,
        buffer,
        PyC1,
        SiC,
        matrix_pin,
        beryllium,
        B4C_drum,
        B4C_central,
    ]
    all_mats = openmc.Materials(list(dict.fromkeys(base_materials)))
    openmc.Materials.cross_sections = params["cross_sections_xml_location"]


    #% --- Compact --------------------------------------------------------------------------------------
    #surf  1 hexyprism 0.0 0.0 13.376      0 200 %fueled core
    # surf_1 isn't used in the OpenMC model
    #surf_1 = openmc.model.HexagonalPrism(edge_length=13.376, orientation='x', origin=(0,0), boundary_type='transmission', corner_radius=0.0)
    #surf 91 hexyprism 0.0 0.0 130.0       0 200 %fueled core
    # OpenMC does not have a z-truncated hex prism
    # create an infinite hex prism, also note OpenMC uses edge length while Serpent uses half width
    surf_91 = openmc.model.HexagonalPrism(edge_length=parms["reflector_width"]/2 * 2 / np.sqrt(3), orientation='x', origin=(0,0), boundary_type='transmission', corner_radius=0.0)
    # Create another slightly larger hexprism and truncating planes
    # What is universe zero in the Serpent model will be put inside a cell made of these surfaces
    surf_91o = openmc.model.HexagonalPrism(edge_length=(surf_91.plane_max.y0 + 0.1) * 2 / np.sqrt(3), orientation='x', origin=(0,0), boundary_type='vacuum', corner_radius=0.0)
    surf_91l = openmc.ZPlane(z0=0, boundary_type='vacuum')
    surf_91u = openmc.ZPlane(z0=2*parms["axial_reflector_height"] + parms["active_fuel_height"], boundary_type='vacuum')

    #surf 95l pz  20
    surf_95l = openmc.ZPlane(z0=parms["axial_reflector_height"])
    #surf 95u pz  180
    surf_95u = openmc.ZPlane(z0=parms["axial_reflector_height"] + parms["active_fuel_height"])
    #surf 96l pz  18.0
    surf_96l = openmc.ZPlane(z0=parms["axial_reflector_height"] - 2.)
    #surf 96u pz  182.0
    surf_96u = openmc.ZPlane(z0=parms["axial_reflector_height"] + parms["active_fuel_height"] + 2.)

    #surf  2 cyl 0.0 0.0 1.00   0.0 200
    surf_2 = openmc.ZCylinder(0, 0, parms["compact_radius"])
    #surf  3 cyl 0.0 0.0 0.825  0.0 200
    surf_3 = openmc.ZCylinder(0, 0, parms["moderator_radius"])
    #surf 35 cyl 0.0 0.0 0.875  0.0 200
    surf_35 = openmc.ZCylinder(0, 0, 0.875)
    #surf 36 cyl 0.0 0.0 0.900  0.0 200
    surf_36 = openmc.ZCylinder(0, 0, 0.900)
    #surf 3g cyl 0.0 0.0 0.920  0.0 200
    surf_3g = openmc.ZCylinder(0, 0, parms["moderator_radius"] + 0.095)
    #surf 5g cyl 0.0 0.0 1.07   0.0 200
    surf_5g = openmc.ZCylinder(0, 0, 1.07)
    #surf  5 cyl 0.0 0.0 1.05   0.0 200
    surf_5 = openmc.ZCylinder(0, 0, 1.05)
    #surf 51 cyl 0.0 0.0 0.97   0.0 200
    surf_51 = openmc.ZCylinder(0, 0, 0.970)
    #surf 52 cyl 0.0 0.0 0.90   0.0 200
    surf_52 = openmc.ZCylinder(0, 0, 0.90)
    #surf 53 cyl 0.0 0.0 0.80   0.0 200
    surf_53 = openmc.ZCylinder(0, 0, 0.80)

    #% infinite cells defining material universes
    # OpenMC cells are made infinite by not assigning a region
    #cell 51 802 moderator  -inf
    cell_51 = openmc.Cell(fill=moderator)
    uni_802 = openmc.Universe(cells=[cell_51])
    #cell 52 33  matrix_pin -inf
    cell_52 = openmc.Cell(fill=matrix_pin)
    uni_33 = openmc.Universe(cells=[cell_52])
    #cell 53 803    matrix  -inf
    cell_53 = openmc.Cell(fill=matrix)
    uni_803 = openmc.Universe(cells=[cell_53])
    #cell 55 804       air  -inf
    cell_55 = openmc.Cell(fill=air)
    uni_804 = openmc.Universe(cells=[cell_55])
    #cell 56 805  beryllium -inf
    cell_56 = openmc.Cell(fill=beryllium)
    uni_805 = openmc.Universe(cells=[cell_56])
    #cell 57 806  shell_mod -inf
    cell_57 = openmc.Cell(fill=shell_mod)
    uni_806 = openmc.Universe(cells=[cell_57])
    #cell 64 813  shell_hp  -inf
    cell_64 = openmc.Cell(fill=shell_hp)
    uni_813 = openmc.Universe(cells=[cell_64])
    #cell 58 820  shell_air_center -inf
    cell_58 = openmc.Cell(fill=shell_air_center)
    uni_820 = openmc.Universe(cells=[cell_58])

    #% mixed wick and coo_vap and liqu
    #cell 65 815 hp_vp_liq_wick -inf
    cell_65 = openmc.Cell(fill=hp_vp_liq_wick)
    uni_815 = openmc.Universe(cells=[cell_65])

    #% mixed air and ss in moderator
    #cell 68 816 shell_air_mod -inf
    cell_68 = openmc.Cell(fill=shell_air_mod)
    uni_816 = openmc.Universe(cells=[cell_68])

    #% mixed air and ss in heatpipe
    #cell 72 817 shell_air_hp   -inf
    cell_72 = openmc.Cell(fill=shell_air_hp)
    uni_817 = openmc.Universe(cells=[cell_72])

    #%  control drums 
    #cell 61 810 beryllium_drum -inf
    cell_61 = openmc.Cell(fill=beryllium_drum)
    uni_810 = openmc.Universe(cells=[cell_61])
    #cell 62 811   B4C_drum -inf
    cell_62 = openmc.Cell(fill=B4C_drum)
    uni_811 = openmc.Universe(cells=[cell_62])
    #cell 63 812   B4C_central -inf
    cell_63 = openmc.Cell(fill=B4C_central)
    uni_812 = openmc.Universe(cells=[cell_63])

    #beginning of pin cells
    #% Yan simplified model removed caps for the HP and moderator pins extended into the Be reflectors 
    #% moderator
    #cell 13   2 fill 802    -3      95l -95u % active length
    cell_13 = openmc.Cell(fill=uni_802, region=(-surf_3 & +surf_95l & -surf_95u)) # moderator material
    #cell 14   2 fill 816     3 -3g  95l -95u
    cell_14 = openmc.Cell(fill=uni_816, region=(+surf_3 & -surf_3g & +surf_95l & -surf_95u)) #shell_air_mod:mixed air and ss in moderator, moderator clad
    #cell 15   2 fill 803     3g     95l -95u
    cell_15 = openmc.Cell(fill=uni_803, region=(+surf_3g & +surf_95l & -surf_95u)) #graphite matrix
    #cell 151  2 fill 805       -36  96l -95l % 0.05 cm shell
    cell_151 = openmc.Cell(fill=uni_805, region=(-surf_36 & +surf_96l & -surf_95l))
    #cell 152  2 fill 805       -36  95u -96u % 0.05 cm shell
    cell_152 = openmc.Cell(fill=uni_805, region=(-surf_36 & +surf_95u & -surf_96u))
    #cell 151a 2 fill 805   -3g  36  96l -95l % 0.05 cm shell
    cell_151a = openmc.Cell(fill=uni_805, region=(-surf_3g & +surf_36 & +surf_96l & -surf_95l))
    #cell 152a 2 fill 805   -3g  36  95u -96u % 0.05 cm shell
    cell_152a = openmc.Cell(fill=uni_805, region=(-surf_3g & +surf_36 & +surf_95u & -surf_96u))
    #cell 153  2 fill 805        3g  96l -95l % 0.05 cm shell
    cell_153 = openmc.Cell(fill=uni_805, region=(+surf_3g & +surf_96l & -surf_95l))
    #cell 154  2 fill 805        3g  95u -96u % 0.05 cm shell
    cell_154 = openmc.Cell(fill=uni_805, region=(+surf_3g & +surf_95u & -surf_96u))
    #cell 155  2 fill 805                 96u
    cell_155 = openmc.Cell(fill=uni_805, region=(+surf_96u))
    #cell 156  2 fill 805                -96l
    cell_156 = openmc.Cell(fill=uni_805, region=(-surf_96l))
    uni_2 = openmc.Universe(cells=[cell_13,cell_14,cell_15,cell_151,cell_152,cell_151a,cell_152a,cell_153,cell_154,cell_155,cell_156])

    #% heat pipe cell
    #cell 16i  1 fill 815      -51 95l
    cell_16i = openmc.Cell(fill=uni_815, region=(-surf_51 & +surf_95l))
    #cell 16o  1 fill 817   51 -5g 95l
    cell_16o = openmc.Cell(fill=uni_817, region=(+surf_51 & -surf_5g & +surf_95l))
    #cell 17   1 fill 803       5g 95l -95u
    cell_17 = openmc.Cell(fill=uni_803, region=(+surf_5g & +surf_95l & -surf_95u))
    #cell 175  1 fill 805       5g 95u
    cell_175 = openmc.Cell(fill=uni_805, region=(+surf_5g & +surf_95u))
    #cell 176  1 fill 805         -96l
    cell_176 = openmc.Cell(fill=uni_805, region=(-surf_96l))
    #cell 177  1 fill 805      -5  96l -95l
    cell_177 = openmc.Cell(fill=uni_805, region=(-surf_5 & +surf_96l & -surf_95l))
    #cell 178  1 fill 805       5g 96l -95l
    cell_178 = openmc.Cell(fill=uni_805, region=(+surf_5g & +surf_96l & -surf_95l))
    #cell 178a 1 fill 805   -5g 5  96l -95l
    cell_178a = openmc.Cell(fill=uni_805, region=(-surf_5g & +surf_5 & +surf_96l & -surf_95l))
    uni_1 = openmc.Universe(cells=[cell_16i,cell_16o,cell_17,cell_175,cell_176,cell_177,cell_178,cell_178a])

    cell_18 = openmc.Cell(fill=uni_820)
    uni_8 = openmc.Universe(cells=[cell_18])

    #% monolith filled cell
    #cell 19  9 fill 803     -inf
    cell_19 = openmc.Cell(fill=uni_803)
    uni_9 = openmc.Universe(cells=[cell_19])

    #% monolith filled cell
    #cell 20 10 fill 805           -inf
    cell_20 = openmc.Cell(fill=uni_805)
    uni_10 = openmc.Universe(cells=[cell_20])

    #% monolith filled cell
    #% cell 21 11 fill 804           -inf
    #cell 21 11 fill 820           -inf
    cell_21 = openmc.Cell(fill=uni_820)
    uni_11 = openmc.Universe(cells=[cell_21])

    #lat 80  2  0.0 0.0 1 1 26.752
    # 10  % Be
    # Single element lattice? Use an "infinite" cell for now
    cell_lat_80 = openmc.Cell(fill=uni_10)
    uni_80 = openmc.Universe(cells=[cell_lat_80])

    #lat 90  2  0.0 0.0 1 1 26.752
    #  8  % air
    # Another single element lattice, using an "infinite" cell
    cell_lat_90 = openmc.Cell(fill=uni_8)
    uni_90 = openmc.Universe(cells=[cell_lat_90])

    #lat 91  2  0.0 0.0 1 1 26.752
    #  11  % air o  B4C
    # Another single element lattice, using an "infinite" cell
    cell_lat_91 = openmc.Cell(fill=uni_11)
    uni_91 = openmc.Universe(cells=[cell_lat_91])

    #% fuel
    #cell 11   3 fill 801        -2  95l -95u % fuel compact
    # cell_11 is the fuel compact containing TRISO particles. Syntax for this doesn't compare 1:1 with Serpent.
    # The Serpent model defines the fuel compact using surf_2 (a RCC that spans the full z of the problem space)
    # and further truncating it with surf_95l and surf_95u. The OpenMC TRISO particle packing and lattice algorithms
    # need a single RCC without further truncations.
    packing_divs_height = (surf_95u.z0 - surf_95l.z0)/axial_divs
    packing_div_surf = openmc.model.RightCircularCylinder((0.0,0.0,-packing_divs_height/2), packing_divs_height, surf_2.r, axis='z')
    packing_div_reg = -packing_div_surf
    if low_pf == True:
        print("RUNNING LOW PF FOR DEBUG")
        triso_centers = openmc.model.pack_spheres(radius=4.2750e-02, region=packing_div_reg, pf=0.01) # Calculate TRISO locations using OpenMC functionality
    else:
        triso_centers = openmc.model.pack_spheres(radius=4.2750e-02, region=packing_div_reg, pf=0.40)
    uni_900_spheres = [openmc.Sphere(r=r) for r in [2.1250e-02,3.1250e-02,3.5250e-02,3.8750e-02,4.2750e-02]]

    #complications arise from depletion divisions
    adiv_fuelmat_list = [[0 for _ in range(axial_divs)] for _ in range(3)]
    adiv_univ_list = [[0 for _ in range(axial_divs)] for _ in range(3)]
    flake_univ_list = [0 for _ in range(3)]
    print("Making depletable_regions")
    for nflake in range(3): #30 flakes in core
        for adiv in range(axial_divs):
            print(nflake, "/", 2, adiv, "/", axial_divs - 1)
            adiv_fuelmat_list[nflake][adiv] = fuel.clone()
            uni_900_cells = [openmc.Cell(fill=adiv_fuelmat_list[nflake][adiv], region=-uni_900_spheres[0]),
                             openmc.Cell(fill=buffer, region=+uni_900_spheres[0] & -uni_900_spheres[1]),
                             openmc.Cell(fill=PyC1, region=+uni_900_spheres[1] & -uni_900_spheres[2]),
                             openmc.Cell(fill=SiC, region=+uni_900_spheres[2] & -uni_900_spheres[3]),
                             openmc.Cell(fill=PyC2, region=+uni_900_spheres[3] & -uni_900_spheres[4])]
            uni_900 = openmc.Universe(cells=uni_900_cells)
            trisos = [openmc.model.TRISO(uni_900_spheres[-1].r, uni_900, center) for center in triso_centers]
            lower_left, upper_right = packing_div_reg.bounding_box
            shape = (2, 2, int(upper_right[2] - lower_left[2]))
            pitch = (upper_right - lower_left)/shape
            triso_lattice = openmc.model.create_triso_lattice(trisos, lower_left, pitch, shape, matrix_pin)
            cell_triso_lattice = openmc.Cell(fill=triso_lattice)
            adiv_univ_list[nflake][adiv] = openmc.Universe(cells=[cell_triso_lattice])

        stacked_packing_division_lattice = openmc.RectLattice()
        stacked_packing_division_lattice.pitch = (1e3, 1e3, packing_divs_height)
        stacked_packing_division_lattice.lower_left = (-1e3/2, -1e3/2, surf_95l.z0 - packing_divs_height)#packing_div_reg.bounding_box[0]
        universes = [[[uni_80]]] + [[[adiv_univ_list[nflake][adiv]]] for adiv in range(axial_divs)] + [[[uni_80]]]
        stacked_packing_division_lattice.universes = universes
        stacked_packing_division_cell = openmc.Cell(fill = stacked_packing_division_lattice)
        uni_801 = openmc.Universe(cells=[stacked_packing_division_cell])
        cell_12 = openmc.Cell(fill=uni_803, region=(+surf_2 & +surf_95l & -surf_95u))
        cell_125 = openmc.Cell(fill=uni_805, region=(+surf_95u))
        cell_126 = openmc.Cell(fill=uni_805, region=(-surf_95l))
        cell_11 = openmc.Cell(fill=uni_801, region=-surf_2 & +surf_95l & -surf_95u)
        uni_3 = openmc.Universe(cells=[cell_11, cell_12, cell_125, cell_126])

        #% --- Assembly lattice
    
        #lat 20  2  0.0 0.0 19 19 2.3
        #% New Assembly Configuration - 2 time less moderator - more fuel
        # 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
        #  9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
        #   9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
        #    9 9 9 9 9 9 9 9 9 1 3 1 3 1 3 1 9 9 9
        #     9 9 9 9 9 9 9 9 3 3 2 3 2 3 2 3 9 9 9
        #      9 9 9 9 9 9 9 1 2 1 3 1 3 1 3 1 9 9 9 
        #       9 9 9 9 9 9 3 3 3 3 2 3 2 3 2 3 9 9 9
        #        9 9 9 9 9 1 2 1 2 1 3 1 3 1 3 1 9 9 9 
        #         9 9 9 9 3 3 3 3 3 3 2 3 2 3 2 3 9 9 9
        #          9 9 9 1 2 1 2 1 2 1 3 1 3 1 3 1 9 9 9
        #           9 9 9 3 3 3 3 3 3 2 3 2 3 2 3 9 9 9 9
        #            9 9 9 1 2 1 2 1 3 1 3 1 3 1 9 9 9 9 9 
        #             9 9 9 3 3 3 3 2 3 2 3 2 3 9 9 9 9 9 9
        #              9 9 9 1 2 1 3 1 3 1 3 1 9 9 9 9 9 9 9 
        #               9 9 9 3 3 2 3 2 3 2 3 9 9 9 9 9 9 9 9
        #                9 9 9 1 3 1 3 1 3 1 9 9 9 9 9 9 9 9 9
        #                 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 
        #                  9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
        #                   9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
        # OpenMC lattices are defined via rings.
        lat_20 = openmc.HexLattice()
        lat_20.pitch = [parms["pin_pitch"]]
        lat_20.outer = uni_9
        lat_20.center = (0.0, 0.0)
        lat_20.orientation = 'x'
        # ORIGINAL
        lat_20.universes = [[x for x in itertools.chain.from_iterable(itertools.repeat([uni_1, uni_3], 18))], # 18 * 2 = 36
                            [x for x in itertools.chain.from_iterable(itertools.repeat([uni_3, uni_2], 15))], # 15 * 2 = 30
                            [x for x in itertools.chain.from_iterable(itertools.repeat([uni_1, uni_3], 12))], # 12 * 2 = 24
                            [x for x in itertools.chain.from_iterable(itertools.repeat([uni_3, uni_2], 9))],  # 9 * 2 = 18
                            [x for x in itertools.chain.from_iterable(itertools.repeat([uni_1, uni_3], 6))],  # 6 * 2 = 12
                            [x for x in itertools.chain.from_iterable(itertools.repeat([uni_3, uni_2], 3))],  # 3 * 2 = 6
                            [uni_1]]
        cell_lat_20 = openmc.Cell(fill=lat_20)
        #trans 20 0.0 0.0 0.0 0.0 0.0 30.
        cell_lat_20.rotation= (0.0,0.0,-30.0)
        flake_univ_list[nflake] = openmc.Universe(cells=[cell_lat_20])
  


    lat_50 = openmc.HexLattice()
    lat_50.pitch = [parms["flake_width"]]
    lat_50.outer = uni_90
    lat_50.center = (0.0, 0.0)
    lat_50.orientation = 'x'
    #assembly-wide depletable regions sit in a ring
    lat_50.universes = [[x for x in itertools.chain.from_iterable(itertools.repeat([uni_90, uni_80], 12))],        # 12 * 2 = 24
                        [x for x in itertools.chain.from_iterable(itertools.repeat([uni_90, flake_univ_list[0], flake_univ_list[0]], 6))], # 6 * 3 = 18
                        [x for x in itertools.chain.from_iterable(itertools.repeat([flake_univ_list[1]], 12))],                # 12 * 1 = 12
                        [x for x in itertools.chain.from_iterable(itertools.repeat([flake_univ_list[2]], 6))],                 # 6 * 1 = 6
                        [uni_91]]
    cell_lat_50 = openmc.Cell(fill=lat_50)
    uni_50 = openmc.Universe(cells=[cell_lat_50])


    #% -----------DRUMS DEFINITTION---------------------------------------------------------------------------------------
    drum_outer_r = (parms["flake_width"] - 0.252)/2
    drum_inner_r = drum_outer_r - 1.

    #surf 731 cyl 80.2560   0.00 13.2500   0.00 200.00
    surf_731 = openmc.model.RightCircularCylinder((parms["flake_width"]*3,0.0,0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 831 cyl 80.2560   0.00 12.2500   0.00 200.00
    surf_831 = openmc.model.RightCircularCylinder((parms["flake_width"]*3,0.0,0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    surf_771 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 180 + drum_rotation, np.array([parms["flake_width"]*3, 0])))
    #cell 801  731 fill 810   -731   -771
    cell_801 = openmc.Cell(fill=uni_810, region=(-surf_731 & -surf_771))
    #cell 802  731 fill 811    831  -731  771
    cell_802 = openmc.Cell(fill=uni_811, region=(+surf_831 & -surf_731 & +surf_771))
    #cell 803  731 fill 810   -831   771
    cell_803 = openmc.Cell(fill=uni_810, region=(-surf_831 & +surf_771))
    uni_731 = openmc.Universe(cells=[cell_801, cell_802, cell_803])
    #surf 732 cyl 40.1280  69.50 13.2500   0.00 200.00
    surf_732 = openmc.model.RightCircularCylinder((parms["flake_width"]*3/2,1.5*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 832 cyl 40.1280  69.50 12.2500   0.00 200.00
    surf_832 = openmc.model.RightCircularCylinder((parms["flake_width"]*3/2,1.5*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    surf_772 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 240 + drum_rotation, np.array([parms["flake_width"]*3/2, 1.5*np.sqrt(3)*parms["flake_width"]])))
    #cell 806  732 fill 810   -732   -772
    cell_806 = openmc.Cell(fill=uni_810, region=(-surf_732 & -surf_772))
    #cell 807  732 fill 811    832  -732  772
    cell_807 = openmc.Cell(fill=uni_811, region=(+surf_832 & -surf_732 & +surf_772))
    #cell 808  732 fill 810   -832   772
    cell_808 = openmc.Cell(fill=uni_810, region=(-surf_832 & +surf_772))
    uni_732 = openmc.Universe(cells=[cell_806, cell_807, cell_808])
    #surf 733 cyl -40.1280  69.50 13.2500   0.00 200.00
    surf_733 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3/2,1.5*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 833 cyl -40.1280  69.50 12.2500   0.00 200.00
    surf_833 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3/2,1.5*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 773 plane 1.0000  -1.73 0.0 -177.84
    surf_773 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 300 + drum_rotation, np.array([-parms["flake_width"]*3/2, 1.5*np.sqrt(3)*parms["flake_width"]])))
    #cell 811  733 fill 810   -733   773
    cell_811 = openmc.Cell(fill=uni_810, region=(-surf_733 & +surf_773))
    #cell 812  733 fill 811    833  -733  -773
    cell_812 = openmc.Cell(fill=uni_811, region=(+surf_833 & -surf_733 & -surf_773))
    #cell 813  733 fill 810   -833   -773
    cell_813 = openmc.Cell(fill=uni_810, region=(-surf_833 & -surf_773))
    uni_733 = openmc.Universe(cells=[cell_811, cell_812, cell_813])
    #surf 734 cyl -80.2560   0.00 13.2500   0.00 200.00
    surf_734 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3,0.0,0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 834 cyl -80.2560   0.00 12.2500   0.00 200.00
    surf_834 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3,0.0,0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 774 plane 1.0000  -0.00 0.0 -88.92
    surf_774 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 0 + drum_rotation, np.array([-parms["flake_width"]*3, 0.0])))
    #cell 816  734 fill 810   -734   774
    cell_816 = openmc.Cell(fill=uni_810, region=(-surf_734 & +surf_774))
    #cell 817  734 fill 811    834  -734  -774
    cell_817 = openmc.Cell(fill=uni_811, region=(+surf_834 & -surf_734 & -surf_774))
    #cell 818  734 fill 810   -834   -774
    cell_818 = openmc.Cell(fill=uni_810, region=(-surf_834 & -surf_774))
    uni_734 = openmc.Universe(cells=[cell_816, cell_817, cell_818])
    #surf 735 cyl -40.1280 -69.50 13.2500   0.00 200.00
    surf_735 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3/2,-1.5*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 835 cyl -40.1280 -69.50 12.2500   0.00 200.00
    surf_835 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3/2,-1.5*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 775 plane 1.0000   1.73 0.0 -177.84
    surf_775 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 60 + drum_rotation, np.array([-parms["flake_width"]*3/2, -1.5*np.sqrt(3)*parms["flake_width"]])))
    #cell 821  735 fill 810   -735   775
    cell_821 = openmc.Cell(fill=uni_810, region=(-surf_735 & +surf_775))
    #cell 822  735 fill 811    835  -735  -775
    cell_822 = openmc.Cell(fill=uni_811, region=(+surf_835 & -surf_735 & -surf_775))
    #cell 823  735 fill 810   -835   -775
    cell_823 = openmc.Cell(fill=uni_810, region=(-surf_835 & -surf_775))
    uni_735 = openmc.Universe(cells=[cell_821, cell_822, cell_823])
    #surf 736 cyl 40.1280 -69.50 13.2500   0.00 200.00
    surf_736 = openmc.model.RightCircularCylinder((parms["flake_width"]*3/2,-1.5*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 836 cyl 40.1280 -69.50 12.2500   0.00 200.00
    surf_836 = openmc.model.RightCircularCylinder((parms["flake_width"]*3/2,-1.5*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 776 plane 1.0000  -1.73 0.0 177.84
    surf_776 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 120 + drum_rotation, np.array([parms["flake_width"]*3/2, -1.5*np.sqrt(3)*parms["flake_width"]])))
    #cell 826  736 fill 810   -736   -776
    cell_826 = openmc.Cell(fill=uni_810, region=(-surf_736 & -surf_776))
    #cell 827  736 fill 811    836  -736  776
    cell_827 = openmc.Cell(fill=uni_811, region=(+surf_836 & -surf_736 & +surf_776))
    #cell 828  736 fill 810   -836   776
    cell_828 = openmc.Cell(fill=uni_810, region=(-surf_836 & +surf_776))
    uni_736 = openmc.Universe(cells=[cell_826, cell_827, cell_828])
    #surf 737 cyl 80.2563  46.34 13.2500   0.00 200.00
    surf_737 = openmc.model.RightCircularCylinder((parms["flake_width"]*3,np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 837 cyl 80.2563  46.34 12.2500   0.00 200.00
    surf_837 = openmc.model.RightCircularCylinder((parms["flake_width"]*3,np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 777 plane 1.0000   0.58 0.0 117.01
    surf_777 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 210 + drum_rotation, np.array([parms["flake_width"]*3, np.sqrt(3)*parms["flake_width"]])))
    #cell 831  737 fill 810   -737   -777
    cell_831 = openmc.Cell(fill=uni_810, region=(-surf_737 & -surf_777))
    #cell 832  737 fill 811    837  -737  777
    cell_832 = openmc.Cell(fill=uni_811, region=(+surf_837 & -surf_737 & +surf_777))
    #cell 833  737 fill 810   -837   777
    cell_833 = openmc.Cell(fill=uni_810, region=(-surf_837 & +surf_777))
    uni_737 = openmc.Universe(cells=[cell_831, cell_832, cell_833])
    #surf 738 cyl 0.0000  92.67 13.2500   0.00 200.00
    surf_738 = openmc.model.RightCircularCylinder((0.0000,2*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 838 cyl 0.0000  92.67 12.2500   0.00 200.00
    surf_838 = openmc.model.RightCircularCylinder((0.0000,2*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 778 plane 0.0000   1.00 0.0 101.33
    surf_778 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 270 + drum_rotation, np.array([0.0, 2*np.sqrt(3)*parms["flake_width"]])))
    #cell 836  738 fill 810   -738   -778
    cell_836 = openmc.Cell(fill=uni_810, region=(-surf_738 & -surf_778))
    #cell 837  738 fill 811    838  -738  778
    cell_837 = openmc.Cell(fill=uni_811, region=(+surf_838 & -surf_738 & +surf_778))
    #cell 838  738 fill 810   -838   778
    cell_838 = openmc.Cell(fill=uni_810, region=(-surf_838 & +surf_778))
    uni_738 = openmc.Universe(cells=[cell_836, cell_837, cell_838])
    #surf 739 cyl -80.2563  46.34 13.2500   0.00 200.00
    surf_739 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3,np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 839 cyl -80.2563  46.34 12.2500   0.00 200.00
    surf_839 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3,np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 779 plane 1.0000  -0.58 0.0 -117.01
    surf_779 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 330 + drum_rotation, np.array([-parms["flake_width"]*3, np.sqrt(3)*parms["flake_width"]])))
    #cell 841  739 fill 810   -739   779
    cell_841 = openmc.Cell(fill=uni_810, region=(-surf_739 & +surf_779))
    #cell 842  739 fill 811    839  -739  -779
    cell_842 = openmc.Cell(fill=uni_811, region=(+surf_839 & -surf_739 & -surf_779))
    #cell 843  739 fill 810   -839   -779
    cell_843 = openmc.Cell(fill=uni_810, region=(-surf_839 & -surf_779))
    uni_739 = openmc.Universe(cells=[cell_841, cell_842, cell_843])
    #surf 740 cyl -80.2563 -46.34 13.2500   0.00 200.00
    surf_740 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3,-np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 840 cyl -80.2563 -46.34 12.2500   0.00 200.00
    surf_840 = openmc.model.RightCircularCylinder((-parms["flake_width"]*3,-np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 780 plane 1.0000   0.58 0.0 -117.01
    surf_780 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 30 + drum_rotation, np.array([-parms["flake_width"]*3, -np.sqrt(3)*parms["flake_width"]])))
    #cell 846  740 fill 810   -740   780
    cell_846 = openmc.Cell(fill=uni_810, region=(-surf_740 & +surf_780))
    #cell 847  740 fill 811    840  -740  -780
    cell_847 = openmc.Cell(fill=uni_811, region=(+surf_840 & -surf_740 & -surf_780))
    #cell 848  740 fill 810   -840   -780
    cell_848 = openmc.Cell(fill=uni_810, region=(-surf_840 & -surf_780))
    uni_740 = openmc.Universe(cells=[cell_846, cell_847, cell_848])
    #surf 741 cyl -0.0000 -92.67 13.2500   0.00 200.00
    surf_741 = openmc.model.RightCircularCylinder((0.0,-2*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 841 cyl -0.0000 -92.67 12.2500   0.00 200.00
    surf_841 = openmc.model.RightCircularCylinder((0.0,-2*np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 781 plane 0.0000   1.00 0.0 -101.33
    surf_781 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 90 + drum_rotation, np.array([0, -2*np.sqrt(3)*parms["flake_width"]])))
    #cell 851  741 fill 810   -741   781
    cell_851 = openmc.Cell(fill=uni_810, region=(-surf_741 & +surf_781))
    #cell 852  741 fill 811    841  -741  -781
    cell_852 = openmc.Cell(fill=uni_811, region=(+surf_841 & -surf_741 & -surf_781))
    #cell 853  741 fill 810   -841   -781
    cell_853 = openmc.Cell(fill=uni_810, region=(-surf_841 & -surf_781))
    uni_741 = openmc.Universe(cells=[cell_851, cell_852, cell_853])
    #surf 742 cyl 80.2563 -46.34 13.2500   0.00 200.00
    surf_742 = openmc.model.RightCircularCylinder((parms["flake_width"]*3,-np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_outer_r, axis='z')
    #surf 842 cyl 80.2563 -46.34 12.2500   0.00 200.00
    surf_842 = openmc.model.RightCircularCylinder((parms["flake_width"]*3,-np.sqrt(3)*parms["flake_width"],0.0), parms["active_fuel_height"] + 2*parms["axial_reflector_height"], drum_inner_r, axis='z')
    #surf 782 plane 1.0000  -0.58 0.0 117.01
    surf_782 = openmc.Plane(*drum_parameters_to_eq(12.25, parms["coating_angle"], 150 + drum_rotation, np.array([parms["flake_width"]*3, -np.sqrt(3)*parms["flake_width"]])))
    #cell 856  742 fill 810   -742   -782
    cell_856 = openmc.Cell(fill=uni_810, region=(-surf_742 & -surf_782))
    #cell 857  742 fill 811    842  -742  782
    cell_857 = openmc.Cell(fill=uni_811, region=(+surf_842 & -surf_742 & +surf_782))
    #cell 858  742 fill 810   -842   782
    cell_858 = openmc.Cell(fill=uni_810, region=(-surf_842 & +surf_782))
    uni_742 = openmc.Universe(cells=[cell_856, cell_857, cell_858])

    #% -----------FILL LATTICE AND DRUMS TO UNIVERSE 0---------------------------------------------------

    #cell 102 0 fill 50 -91 731 732 733 734 735 736 737 738 739 740 741 742
    cell_102 = openmc.Cell(fill=uni_50, region=(-surf_91 & +surf_731 & +surf_732 & +surf_733 &
                                                +surf_734 & +surf_735 & +surf_736 & +surf_737 &
                                                +surf_738 & +surf_739 & +surf_740 & +surf_741 &
                                                +surf_742))
    #cell 180 0 fill       731 -731
    cell_180 = openmc.Cell(fill=uni_731, region=(-surf_731))
    #cell 181 0 fill       732 -732
    cell_181 = openmc.Cell(fill=uni_732, region=(-surf_732))
    #cell 182 0 fill       733 -733
    cell_182 = openmc.Cell(fill=uni_733, region=(-surf_733))
    #cell 183 0 fill       734 -734
    cell_183 = openmc.Cell(fill=uni_734, region=(-surf_734))
    #cell 184 0 fill       735 -735
    cell_184 = openmc.Cell(fill=uni_735, region=(-surf_735))
    #cell 185 0 fill       736 -736
    cell_185 = openmc.Cell(fill=uni_736, region=(-surf_736))
    #cell 186 0 fill       737 -737
    cell_186 = openmc.Cell(fill=uni_737, region=(-surf_737))
    #cell 187 0 fill       738 -738
    cell_187 = openmc.Cell(fill=uni_738, region=(-surf_738))
    #cell 188 0 fill       739 -739
    cell_188 = openmc.Cell(fill=uni_739, region=(-surf_739))
    #cell 189 0 fill       740 -740
    cell_189 = openmc.Cell(fill=uni_740, region=(-surf_740))
    #cell 190 0 fill       741 -741
    cell_190 = openmc.Cell(fill=uni_741, region=(-surf_741))
    #cell 191 0 fill       742 -742
    cell_191 = openmc.Cell(fill=uni_742, region=(-surf_742))
    #cell 104 0 outside  91
    # The hexprisms in OpenMC are infinite in z unlike Serpent
    # where they have z limits. I kept cell_102 infinite in order
    # to keep the above section as close to a 1:1 Serpent conversion
    # as possible.
    # Define cell_104 as void outside surf_91
    cell_104 = openmc.Cell(fill=None, region=(+surf_91))
    # Assembly uni_0
    uni_0 = openmc.Universe(cells=[cell_102,cell_180,cell_181,cell_182,
                                   cell_183,cell_184,cell_185,cell_186,
                                   cell_187,cell_188,cell_189,cell_190,
                                   cell_191,cell_104])
    # Stick uni_0 inside a hex prism with truncating z planes
    cell_00 = openmc.Cell(fill=uni_0, region=(-surf_91o & -surf_91u & +surf_91l))
    uni_00 = openmc.Universe(cells=[cell_00])
    #% --------------------------------------------------------------------------------------------------

    #####################################
    ##. Export geometry and materials
    #####################################
    geom = openmc.Geometry(uni_00)

    #######################
    ##. Settings to Run
    ########################

    settings = openmc.Settings()
    # Particles and batches set to match reference Serpent input so runtime can be compared.
    # Number of particles needed to generate cross sections with adequate statistics is TBD.
    settings.particles = particles
    settings.batches = batches 
    settings.inactive = inactive
    settings.output      = {'tallies': False, # turns off text tally output, hdf5 still produced
                            'summary': False} # the summary.h5 file is taking hours to write due to complexity of model with TRISO
    settings.temperature = {'default': params['Common Temperature'],
                            'method':'interpolation',
                            'range':(300.0, 1200.0)}
    # Source definition for first batch
    # Default is a isotropic Watt fission spectrum at (0,0,0). Move it to a point at core midplane
    initial_ksource = openmc.IndependentSource()
    initial_ksource.space = openmc.stats.Point(xyz=(0.0,0.0,parms["axial_reflector_height"] + parms["active_fuel_height"]/2))
    settings.source = initial_ksource


    #######################
    ##. Tallies
    ########################
    tallies = openmc.Tallies()
    if tally_mesh_dimension:
        pin_power_tally = openmc.Tally(name="pin_power_kappa")
        pin_power_tally.scores = ["kappa-fission"]
        low_pinz = parms["axial_reflector_height"]
        upper_pinz = parms["axial_reflector_height"] + parms["active_fuel_height"]
        zmesh = openmc.RegularMesh()
        zmesh.dimension = tuple(int(value) for value in tally_mesh_dimension)
        low_x = 3*parms["flake_width"]
        low_y = 2*parms["flake_width"]/np.sqrt(3)/2*5.5
        zmesh.lower_left = (-low_x, -low_y, low_pinz)
        zmesh.upper_right = (low_x, low_y, upper_pinz)
        pin_axial_position_filter = openmc.MeshFilter(zmesh)
        pin_power_tally.filters = [pin_axial_position_filter]
        tallies.append(pin_power_tally)

    material_colors = {
        air: 'white',  # rgb(255, 255, 255)
        shell_mod: 'grey',  # rgb(133, 133, 133)
        shell_hp: 'darkgrey',  # rgb(169, 169, 169)
        shell_air_mod: 'silver',  # rgb(192, 192, 192)
        shell_air_hp: 'lightgrey',  # rgb(211, 211, 211)
        shell_air_center: 'gainsboro',  # rgb(220, 220, 220)
        coo_vap: 'purple',  # rgb(128, 0, 128)
        coo_liq: 'green',  # rgb(0, 128, 0)
        wick: 'yellow',  # rgb(255, 255, 0)
        hp_vp_liq_wick: 'khaki',  # rgb(240, 230, 140)
        moderator: 'blueviolet',  # rgb(138, 43, 226)
        fuel: 'red',  # rgb(255, 0, 0)
        buffer: 'cyan',  # rgb(0, 255, 255)
        PyC1: 'teal',  # rgb(0, 128, 128)
        PyC2: 'darkcyan',  # rgb(0, 139, 139)
        SiC: 'darkgreen',  # rgb(0, 100, 0)
        matrix_pin: 'lightgray',  # rgb(211, 211, 211)
        matrix: 'darkgray',  # rgb(169, 169, 169)
        beryllium: 'magenta',  # rgb(0, 255, 0)
        beryllium_drum: 'limegreen',  # rgb(50, 205, 50)
        B4C_drum: 'crimson',  # rgb(220, 20, 60)
        B4C_central: 'darkred',  # rgb(139, 0, 0)
    }
    allfuel_mats = list(itertools.chain.from_iterable(adiv_fuelmat_list))
    for m in allfuel_mats:
        material_colors[m] = (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
    for a in allfuel_mats:
        all_mats.append(a)

    fuel_vol_per_flake = 4/3*np.pi*(2.1250e-02)**3*len(trisos)*63 #fuel vol in axial depletion region in cm^3
    for adiv in range(axial_divs):
        adiv_fuelmat_list[0][adiv].volume = fuel_vol_per_flake*12
        adiv_fuelmat_list[1][adiv].volume = fuel_vol_per_flake*12
        adiv_fuelmat_list[2][adiv].volume = fuel_vol_per_flake*6
    fuel.depletable = False
    print("Total Fuel Volume", fuel_vol_per_flake*(12+12+6)*axial_divs)


    # View of entire reactor at z=21 to match Serpent plot
    plot1 = openmc.Plot()
    plot1.colors = material_colors
    plot1.basis = 'xy'
    plot1.origin = (0.0,0.0,100+packing_divs_height/5)
    #plot1.origin = (0.0,0.0,199.0)
    plot1.width = (300.0, 300.0)
    #plot1.pixels = (1000, 1000)   # low fidelity for shorter runtime
    plot1.pixels = (2000, 2000) # high fidelity to match Serpent
    plot1.color_by = 'material'

    # View of entire reactor at y=5 to match Serpent plot
    plot2 = openmc.Plot()
    plot2.basis = 'xz'
    plot2.origin = (0.0,5.0,100.0)
    plot2.width = (280.0, 280.0)
    #plot2.pixels = (1000, 1000)   # low fidelity for shorter runtime
    plot2.pixels = (2000, 2000) # high fidelity to match Serpent
    plot2.color_by = 'material'
    plot2.colors = material_colors

    # Views of fuel assembly
    plot3 = openmc.Plot()
    plot3.basis = 'xy'
    plot3.origin = (parms["flake_width"],0.0,100+packing_divs_height/5)
    plot3.width = (parms["flake_width"], parms["flake_width"])
    plot3.pixels = (500, 500)
    plot3.color_by = 'material'
    plot3.colors = material_colors

    # Views of fuel compacts
    plot4 = openmc.Plot()
    plot4.basis = 'xy'
    plot4.origin = (parms["flake_width"],0.0,100+packing_divs_height/5)
    plot4.width = (parms["pin_pitch"]*3, parms["pin_pitch"]*3)
    plot4.pixels = (750, 750)
    plot4.color_by = 'material'
    plot4.colors = material_colors

    plots = openmc.Plots([plot1, plot2, plot3, plot4])

    #volume calculation
    #low_x = 3*parms["flake_width"]
    #low_y = 2*parms["flake_width"]/np.sqrt(3)/2*5.5
    #ll = (-low_x, -low_y, low_pinz)
    #ur = (low_x, low_y, upper_pinz)
    #vol_calc = openmc.VolumeCalculation(allfuel_mats, int(1e8), ll, ur)
    #settings.volume_calculations = [vol_calc]

    model = openmc.model.Model(geometry = geom, materials = all_mats, settings = settings, plots = plots, tallies = tallies)
    model.export_to_xml()
    #calculate volumes
    #openmc.calculate_volumes(threads = 48)
    #print("VOLUME RUNTIME", time.time() - start)
    #print("ANA RESULTS")
    #for adiv in range(axial_divs):
    #    print(adiv_fuelmat_list[0][adiv].id, fuel_vol_per_flake*12)
    #    print(adiv_fuelmat_list[1][adiv].id, fuel_vol_per_flake*12)
    #    print(adiv_fuelmat_list[2][adiv].id, fuel_vol_per_flake*6)
