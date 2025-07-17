import glob
from IPython.display import Image
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import pandas as pd
from math import pi, sqrt, tan, sin, cos
import openmc
import openmc.model
import openmc.deplete

#=================================================================================================
#                                    Input File Parameters
#=================================================================================================
batches      = 100
inactive     = 20
particles    = 1000
theta        = pi/180.0
#=================================================================================================
#                                Depletion simulation parameters
#=================================================================================================
tf            = 24*60*60
tot_power     = (5.0e+06)/160
time_steps    = [  0.01*tf,   0.99*tf,   3.00*tf,   6.00*tf,  20.00*tf,  70.00*tf, 100.00*tf, 165.00*tf, 365.00*tf, 365.00*tf, 365.00*tf]
power         = [tot_power, tot_power, tot_power, tot_power, tot_power, tot_power, tot_power, tot_power, tot_power, tot_power, tot_power]
fuel_radius   = 1.0
fuel_xs_area  = pi*fuel_radius*fuel_radius
no_rods_ass   = 72
no_ass = 18
mid_length    = 1
vol_fuel_12   = fuel_xs_area * mid_length * no_ass * no_rods_ass
#=================================================================================================
#                                    materials.xml file
#=================================================================================================
# UO2 fuel definition 
uo2_12 = openmc.Material(2, name='uo2_12')
uo2_12.set_density('atom/b-cm', 8.08250295E-02)
uo2_12.temperature = 600
uo2_12.volume      = vol_fuel_12
uo2_12.add_nuclide(     'O16', 2.59371545E-03, 'ao')
uo2_12.add_nuclide(     'O17', 1.05004397E-06, 'ao')
uo2_12.add_nuclide(     'O18', 5.99797186E-06, 'ao')
uo2_12.add_nuclide(    'U235', 2.59885326E-04, 'ao')
uo2_12.add_nuclide(    'U238', 1.04049396E-03, 'ao')
uo2_12.add_nuclide(    'Si28', 2.76954169E-03, 'ao')
uo2_12.add_nuclide(    'Si29', 1.40694868E-04, 'ao')
uo2_12.add_nuclide(    'Si30', 9.28556098E-05, 'ao')
uo2_12.add_nuclide(     'C12', 7.31619752E-02, 'ao')
uo2_12.add_nuclide(     'C13', 7.58819416E-04, 'ao')
uo2_12.add_s_alpha_beta('c_Graphite')

monolith_12 = openmc.Material(8, 'monolith_12')
monolith_12.set_density('g/cm3', 1.63)
monolith_12.temperature = 600
monolith_12.add_nuclide('C12' , 0.9893, 'ao')
monolith_12.add_nuclide('C13' , 0.0107, 'ao')
monolith_12.add_s_alpha_beta('c_Graphite')

monolith_rr = openmc.Material(15, 'monolith_rr')
monolith_rr.set_density('g/cm3', 1.63)
monolith_rr.temperature = 600
monolith_rr.add_nuclide('C12' , 0.9893, 'ao')
monolith_rr.add_nuclide('C13' , 0.0107, 'ao')
monolith_rr.add_s_alpha_beta('c_Graphite')

monolith_cr = openmc.Material(16, 'monolith_cr')
monolith_cr.set_density('g/cm3', 1.63)
monolith_cr.temperature = 600
monolith_cr.add_nuclide('C12' , 0.9893, 'ao')
monolith_cr.add_nuclide('C13' , 0.0107, 'ao')
monolith_cr.add_s_alpha_beta('c_Graphite')

htpipe = openmc.Material(17, 'htpipe')
htpipe.set_density('atom/b-cm', 2.74917E-02)
htpipe.temperature = 600
htpipe.add_nuclide('Si28',  1.49701E-02, 'ao')
htpipe.add_nuclide('Si29',  7.60143E-04, 'ao')
htpipe.add_nuclide('Si30',  5.01090E-04, 'ao')
htpipe.add_nuclide('Cr50',  6.46763E-03, 'ao')
htpipe.add_nuclide('Cr52',  1.24724E-01, 'ao')
htpipe.add_nuclide('Cr53',  1.41423E-02, 'ao')
htpipe.add_nuclide('Cr54',  3.52029E-03, 'ao')
htpipe.add_nuclide('Mn55',  1.66133E-02, 'ao')
htpipe.add_nuclide('Fe54',  3.12186E-02, 'ao')
htpipe.add_nuclide('Fe56',  4.90061E-01, 'ao')
htpipe.add_nuclide('Fe57',  1.13180E-02, 'ao')
htpipe.add_nuclide('Fe58',  1.50617E-03, 'ao')
htpipe.add_nuclide('Ni58',  6.33738E-02, 'ao')
htpipe.add_nuclide('Ni60',  2.44119E-02, 'ao')
htpipe.add_nuclide('Ni61',  1.06115E-03, 'ao')
htpipe.add_nuclide('Ni62',  3.38338E-03, 'ao')
htpipe.add_nuclide('Ni64',  8.61654E-04, 'ao')
htpipe.add_nuclide('Mo92',  1.75699E-03, 'ao')
htpipe.add_nuclide('Mo94',  1.09514E-03, 'ao')
htpipe.add_nuclide('Mo95',  1.88484E-03, 'ao')
htpipe.add_nuclide('Mo96',  1.97478E-03, 'ao')
htpipe.add_nuclide('Mo97',  1.13066E-03, 'ao')
htpipe.add_nuclide('Mo98',  2.85681E-03, 'ao')
htpipe.add_nuclide('Mo100', 1.14011E-03, 'ao')
htpipe.add_nuclide('Na23',  1.79266E-01, 'ao')

absorber = openmc.Material(18, 'absorber')
absorber.set_density('g/cm3', 2.52 )
absorber.temperature = 600
absorber.add_nuclide(  'C12' , 0.19786, 'ao')
absorber.add_nuclide(  'C13' , 0.00214, 'ao')
absorber.add_nuclide(  'B10' , 0.15720, 'ao')
absorber.add_nuclide(  'B11' , 0.64280, 'ao')

gap = openmc.Material(19, 'gap')
gap.set_density('atom/b-cm', 1.0E-010)
gap.temperature = 600
gap.add_nuclide(  'He4' , 1.0, 'ao')

beryllium = openmc.Material(20, 'beryllium')
beryllium.set_density('g/cm3',  1.85)
beryllium.temperature = 600
beryllium.add_nuclide('Be9' , 1.0, 'ao')
beryllium.add_s_alpha_beta('c_Be')

materials_file = openmc.Materials([uo2_12, monolith_12, monolith_rr, monolith_cr, htpipe, 
                                   absorber, gap, beryllium])
materials_file.export_to_xml()
#=================================================================================================
#                   geometry.xml file
#=================================================================================================
fuel_1                  = openmc.ZCylinder(surface_id=1, x0=0.0, y0=0.0, r=1.00,  name='fuel_1')
fuel_2                  = openmc.ZCylinder(surface_id=2, x0=0.0, y0=0.0, r=1.05,  name='fuel_2')
hp_1                    = openmc.ZCylinder(surface_id=3, x0=0.0, y0=0.0, r=1.10,  name='hp_1')
hp_2                    = openmc.ZCylinder(surface_id=4, x0=0.0, y0=0.0, r=1.15,  name='hp_2')

rods_pitch              = 3.4
ass_rings               = 6 
l2                      = 32.00/sqrt(3.0)
l2                      = l2 - 0.4/sqrt(3.0)
c2                      = sqrt(3.)/3.
x2                      = 0.0
y2                      = 0.0
r_right                 = openmc.YPlane(surface_id=5,y0= y2 + sqrt(3.)/2.*l2,)
r_left                  = openmc.YPlane(surface_id=6,y0= y2 - sqrt(3.)/2.*l2,)
r_upper_right           = openmc.Plane(surface_id=7,  b=c2,  a=1., d= l2+y2*c2+x2, name='r_upper_right')  # y = -x/sqrt(3) + a
r_upper_left            = openmc.Plane(surface_id=8, b=-c2,  a=1., d= l2-y2*c2+x2, name='r_upper_left ')  # y = x/sqrt(3) + a
r_lower_right           = openmc.Plane(surface_id=9, b=-c2,  a=1., d=-l2-y2*c2+x2, name='r_lower_right')  # y = x/sqrt(3) - a
r_lower_left            = openmc.Plane(surface_id=10,  b=c2,  a=1., d=-l2+y2*c2+x2, name='r_lower_left ')  # y = -x/sqrt(3) - a
l2                      = l2 + 0.4/sqrt(3.0)
c2                      = sqrt(3.)/3.
x2                      = 0.0
y2                      = 0.0
g_right                 = openmc.YPlane(surface_id=11,y0= y2 + sqrt(3.)/2.*l2,)
g_left                  = openmc.YPlane(surface_id=12,y0= y2 - sqrt(3.)/2.*l2,)
g_upper_right           = openmc.Plane(surface_id=13,  b=c2,  a=1., d= l2+y2*c2+x2, name='g_upper_right')  # y = -x/sqrt(3) + a
g_upper_left            = openmc.Plane(surface_id=14, b=-c2,  a=1., d= l2-y2*c2+x2, name='g_upper_left ')  # y = x/sqrt(3) + a
g_lower_right           = openmc.Plane(surface_id=15, b=-c2,  a=1., d=-l2-y2*c2+x2, name='g_lower_right')  # y = x/sqrt(3) - a
g_lower_left            = openmc.Plane(surface_id=16,  b=c2,  a=1., d=-l2+y2*c2+x2, name='g_lower_left ')  # y = -x/sqrt(3) - a

lc                      = 75.0*2.0/sqrt(3.0) 
cr                      = sqrt(3.)/3.
x2                      = 0.0
y2                      = 0.0
c_right                 = openmc.XPlane(surface_id=31,x0= x2 + sqrt(3.)/2.*lc,)
c_left                  = openmc.XPlane(surface_id=32,x0= x2 - sqrt(3.)/2.*lc,)
c_upper_right           = openmc.Plane(surface_id=33,  a=cr,  b=1., d= lc+x2*cr+y2, name='c_upper_right') 
c_upper_left            = openmc.Plane(surface_id=34, a=-cr,  b=1., d= lc-x2*cr+y2, name='c_upper_left')
c_lower_right           = openmc.Plane(surface_id=35, a=-cr,  b=1., d=-lc-x2*cr+y2, name='c_lower_right')
c_lower_left            = openmc.Plane(surface_id=36,  a=cr,  b=1., d=-lc+x2*cr+y2, name='c_lower_left')

lc                      = 75.05*2.0/sqrt(3.0) 
cr                      = sqrt(3.)/3.
x2                      = 0.0
y2                      = 0.0
c_right_1               = openmc.XPlane(surface_id=301,x0= x2 + sqrt(3.)/2.*lc,)
c_left_1                = openmc.XPlane(surface_id=302,x0= x2 - sqrt(3.)/2.*lc,)
c_upper_right_1         = openmc.Plane(surface_id=303,  a=cr,  b=1., d= lc+x2*cr+y2, name='c_upper_right_1') 
c_upper_left_1          = openmc.Plane(surface_id=304, a=-cr,  b=1., d= lc-x2*cr+y2, name='c_upper_left_1')
c_lower_right_1         = openmc.Plane(surface_id=305, a=-cr,  b=1., d=-lc-x2*cr+y2, name='c_lower_right_1')
c_lower_left_1          = openmc.Plane(surface_id=306,  a=cr,  b=1., d=-lc+x2*cr+y2, name='c_lower_left_1')

cr_top                  = openmc.Plane(surface_id=41, a=-tan(60 * theta),b=1.0,  name='cr_top')
cr_bot                  = openmc.Plane(surface_id=42, a=-tan(-60 * theta), b=1.0,  name='cr_bot')
cr_in                   = openmc.ZCylinder(surface_id=43, x0=0.0, y0=0.0, r=14.00,  name='cr_in')
cr_out                  = openmc.ZCylinder(surface_id=44, x0=0.0, y0=0.0, r=15.00,  name='cr_out')
cr_gap                  = openmc.ZCylinder(surface_id=45, x0=0.0, y0=0.0, r=15.05,  name='cr_gap')

bottom_0                = openmc.ZPlane(surface_id=51,  z0=-0.5, name='bottom_0')
top_0                   = openmc.ZPlane(surface_id=52, z0=+0.5, name='top_0')

CR_000_180              = openmc.YPlane(surface_id=70, y0= 0.0, name='CR_000')
CR_030_210              = openmc.Plane(surface_id=71,   a=-tan(30 * theta),  b=1.0, d= 0, name='CR_030')
CR_060_240              = openmc.Plane(surface_id=72,   a=-tan(60 * theta),  b=1.0, d= 0, name='CR_060')
CR_090_270              = openmc.XPlane(surface_id=73, x0= 0.0, name='CR_090')
CR_120_300              = openmc.Plane(surface_id=74,  a=-tan(120 * theta),  b=1.0, d= 0, name='CR_120')
CR_150_330              = openmc.Plane(surface_id=75,  a=-tan(150 * theta),  b=1.0, d= 0, name='CR_150')

core_out_2              = openmc.ZCylinder(surface_id=82, x0=0.0, y0=0.0, r=112.0,  name='core_out_2')


core_out_2.boundary_type= 'vacuum'
top_0.boundary_type     = 'reflective'
bottom_0.boundary_type  = 'reflective'

fuel_12                 = openmc.Cell(cell_id=103, name='fuel_12')
fgrp_12                 = openmc.Cell(cell_id=104, name='fgrp_12')
fvod_10                 = openmc.Cell(cell_id=108, name='fvod_10')

hpco_12                 = openmc.Cell(cell_id=113, name='hpco_12')
hpgr_12                 = openmc.Cell(cell_id=114, name='hpgr_12')
hpvd_10                 = openmc.Cell(cell_id=118, name='hpvd_10')

cr_drum                 = openmc.Cell(cell_id=21, name='cr_drum')
cr_refl                 = openmc.Cell(cell_id=22, name='cr_refl')
cr_gpp                  = openmc.Cell(cell_id=23, name='cr_gap')
cr_ass                  = openmc.Cell(cell_id=24, name='cr_ass')

cr_000                  = openmc.Cell(cell_id=40, name='cr_000')
cr_030                  = openmc.Cell(cell_id=41, name='cr_030')
cr_060                  = openmc.Cell(cell_id=42, name='cr_060')
cr_090                  = openmc.Cell(cell_id=43, name='cr_090')
cr_120                  = openmc.Cell(cell_id=44, name='cr_120')
cr_150                  = openmc.Cell(cell_id=45, name='cr_150')
cr_180                  = openmc.Cell(cell_id=46, name='cr_180')
cr_210                  = openmc.Cell(cell_id=47, name='cr_210')
cr_240                  = openmc.Cell(cell_id=48, name='cr_240')
cr_270                  = openmc.Cell(cell_id=49, name='cr_270')
cr_300                  = openmc.Cell(cell_id=50, name='cr_300')
cr_330                  = openmc.Cell(cell_id=51, name='cr_330')

assembly_reg_1          = openmc.Cell(cell_id=332, name='assembly_reg_1')
assembly_gap_12         = openmc.Cell(cell_id=334, name='assembly_gap_12')
grp_cc_cnt_2            = openmc.Cell(cell_id=342, name='grp_cc_cnt')
core_reg                = openmc.Cell(cell_id=345, name='core_reg')
core_reg_out            = openmc.Cell(cell_id=346, name='core_reg_out')
out_univ_2              = openmc.Cell(cell_id=348, name='out_univ_2')

cr_01                   = openmc.Cell(cell_id=61, name='cr_01')
cr_02                   = openmc.Cell(cell_id=62, name='cr_02')
cr_03                   = openmc.Cell(cell_id=63, name='cr_03')
cr_04                   = openmc.Cell(cell_id=64, name='cr_04')
cr_05                   = openmc.Cell(cell_id=65, name='cr_05')
cr_06                   = openmc.Cell(cell_id=66, name='cr_06')
cr_07                   = openmc.Cell(cell_id=67, name='cr_07')
cr_08                   = openmc.Cell(cell_id=68, name='cr_08')
cr_09                   = openmc.Cell(cell_id=69, name='cr_09')
cr_10                   = openmc.Cell(cell_id=70, name='cr_10')
cr_11                   = openmc.Cell(cell_id=71, name='cr_11')
cr_12                   = openmc.Cell(cell_id=72, name='cr_12')

fuel_12.region          = -fuel_1 & +bottom_0 & -top_0
fgrp_12.region          = +fuel_2 & +bottom_0 & -top_0
fvod_10.region          = +fuel_1 & -fuel_2 & +bottom_0 & -top_0

hpco_12.region          = -hp_1 & +bottom_0 & -top_0
hpgr_12.region          = +hp_2 & +bottom_0 & -top_0
hpvd_10.region          = +hp_1 & -hp_2 & +bottom_0 & -top_0

assembly_reg_1.region   = +r_left & -r_right & -r_upper_right & -r_upper_left & +r_lower_right & +r_lower_left & +bottom_0 & -top_0
assembly_gap_12.region  = (-r_left | +r_right | +r_upper_right | +r_upper_left | -r_lower_right | -r_lower_left) & +g_left & -g_right & -g_upper_right & -g_upper_left & +g_lower_right & +g_lower_left & +bottom_0 & -top_0

grp_cc_cnt_2.region     = +g_left & -g_right & -g_upper_right & -g_upper_left & +g_lower_right & +g_lower_left & +bottom_0 & -top_0
core_reg.region         = +c_left & -c_right & -c_upper_right & -c_upper_left & +c_lower_right & +c_lower_left & +bottom_0 & -top_0
core_reg_out.region     = (-c_left | +c_right | +c_upper_right | +c_upper_left | -c_lower_right | -c_lower_left) & +c_left_1 & -c_right_1 & -c_upper_right_1 & -c_upper_left_1 & +c_lower_right_1 & +c_lower_left_1 & +bottom_0 & -top_0

cr_drum.region          =  +cr_in & -cr_out & +cr_bot & -cr_top & +bottom_0 & -top_0
cr_refl.region          =  (-cr_in | -cr_bot | +cr_top) & -cr_out & +bottom_0 & -top_0
cr_gpp.region           =  +cr_out & -cr_gap & +bottom_0 & -top_0
cr_ass.region           =  +cr_gap & +bottom_0 & -top_0

out_univ_2.region       =  +bottom_0 & -top_0

fuel_12.fill            = uo2_12
fgrp_12.fill            = monolith_12
fvod_10.fill            = gap

hpco_12.fill            = htpipe
hpgr_12.fill            = monolith_12
hpvd_10.fill            = gap

cr_drum.fill            = absorber
cr_refl.fill            = beryllium
cr_gpp.fill             = gap
cr_ass.fill             = beryllium

grp_cc_cnt_2.fill       = monolith_cr
assembly_gap_12.fill    = monolith_12

out_univ_2.fill         = monolith_rr

f1                      = openmc.Universe()
h1                      = openmc.Universe()
uo                      = openmc.Universe()
c1                      = openmc.Universe()
c3                      = openmc.Universe()
c4                      = openmc.Universe()
root                    = openmc.Universe(universe_id=0, name='root universe')


f1.add_cells([fuel_12, fgrp_12, fvod_10,])
h1.add_cells([hpco_12, hpgr_12, hpvd_10,])
uo.add_cells([out_univ_2,])

lattice_hex_1             = openmc.HexLattice(lattice_id=55)
lattice_hex_1.center      = (0., 0.,)
lattice_hex_1.pitch       = (3.4,)
lattice_hex_1.n_rings     = (6)
lattice_hex_1.orientation = 'x'
lattice_hex_1.universes   = [[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1], 
                             [h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1],
                             [f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1]+[f1], 
                             [h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1]+[h1]+[f1], 
                             [f1]+[f1]+[f1]+[f1]+[f1]+[f1],
                             [h1]]
lattice_hex_1.outer       = uo
assembly_reg_1.fill       = lattice_hex_1


c1.add_cells([assembly_reg_1,assembly_gap_12,])
c3.add_cells([grp_cc_cnt_2,])
c4.add_cells([cr_drum,cr_refl,cr_gpp,cr_ass,])

cr_000.fill             = c4
cr_030.fill             = c4
cr_060.fill             = c4
cr_090.fill             = c4
cr_120.fill             = c4
cr_150.fill             = c4
cr_180.fill             = c4
cr_210.fill             = c4
cr_240.fill             = c4
cr_270.fill             = c4
cr_300.fill             = c4
cr_330.fill             = c4

rotation_angle          = 0
cr_000.rotation         = [0,  0,    0 + rotation_angle]
cr_330.rotation         = [0,  0,    0 + rotation_angle]
cr_030.rotation         = [0,  0,   60 + rotation_angle]
cr_060.rotation         = [0,  0,   60 + rotation_angle]
cr_090.rotation         = [0,  0,  120 + rotation_angle]
cr_120.rotation         = [0,  0,  120 + rotation_angle]
cr_150.rotation         = [0,  0,  180 + rotation_angle]
cr_180.rotation         = [0,  0,  180 + rotation_angle]
cr_210.rotation         = [0,  0, -120 + rotation_angle]
cr_240.rotation         = [0,  0, -120 + rotation_angle]
cr_270.rotation         = [0,  0,  -60 + rotation_angle]
cr_300.rotation         = [0,  0,  -60 + rotation_angle]

r_0                     = (78 + 112)/2.0
r_angle                 = (0.0+30.0)/2
cr_000.translation      = [r_0*cos((r_angle + 30*0)*theta),  r_0*sin((r_angle + 30*0)*theta),  0]
cr_030.translation      = [r_0*cos((r_angle + 30*1)*theta),  r_0*sin((r_angle + 30*1)*theta),  0]
cr_060.translation      = [r_0*cos((r_angle + 30*2)*theta),  r_0*sin((r_angle + 30*2)*theta),  0]
cr_090.translation      = [r_0*cos((r_angle + 30*3)*theta),  r_0*sin((r_angle + 30*3)*theta),  0]
cr_120.translation      = [r_0*cos((r_angle + 30*4)*theta),  r_0*sin((r_angle + 30*4)*theta),  0]
cr_150.translation      = [r_0*cos((r_angle + 30*5)*theta),  r_0*sin((r_angle + 30*5)*theta),  0]
cr_180.translation      = [r_0*cos((r_angle + 30*6)*theta),  r_0*sin((r_angle + 30*6)*theta),  0]
cr_210.translation      = [r_0*cos((r_angle + 30*7)*theta),  r_0*sin((r_angle + 30*7)*theta),  0]
cr_240.translation      = [r_0*cos((r_angle + 30*8)*theta),  r_0*sin((r_angle + 30*8)*theta),  0]
cr_270.translation      = [r_0*cos((r_angle + 30*9)*theta),  r_0*sin((r_angle + 30*9)*theta),  0]
cr_300.translation      = [r_0*cos((r_angle + 30*10)*theta),  r_0*sin((r_angle + 30*10)*theta),  0]
cr_330.translation      = [r_0*cos((r_angle + 30*11)*theta),  r_0*sin((r_angle + 30*11)*theta),  0]

a1                      = openmc.Universe()
a2                      = openmc.Universe()
a3                      = openmc.Universe()
a4                      = openmc.Universe()
a5                      = openmc.Universe()
a6                      = openmc.Universe()
a7                      = openmc.Universe()
a8                      = openmc.Universe()
a9                      = openmc.Universe()
a10                     = openmc.Universe()
a11                     = openmc.Universe()
a12                     = openmc.Universe()

a1.add_cells([cr_000,])
a2.add_cells([cr_030,])
a3.add_cells([cr_060,])
a4.add_cells([cr_090,])
a5.add_cells([cr_120,])
a6.add_cells([cr_150,])
a7.add_cells([cr_180,])
a8.add_cells([cr_210,])
a9.add_cells([cr_240,])
a10.add_cells([cr_270,])
a11.add_cells([cr_300,])
a12.add_cells([cr_330,])

core_hex                = openmc.HexLattice(lattice_id=65)
core_hex.center         = (0., 0.,)
core_hex.pitch          = (32.000,)
core_hex.n_rings        = (4)
core_hex.orientation    = 'y'
core_hex.universes      = [[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo]+[uo], 
                           [c1]+[c1]+[c1]+[c1]+[c1]+[c1]+[c1]+[c1]+[c1]+[c1]+[c1]+[c1], 
                           [c1]+[c1]+[c1]+[c1]+[c1]+[c1],
                           [c3]]
core_hex.outer          = uo
core_reg.fill           = core_hex
core_reg_out.fill       = gap

cr_01.region            = +CR_000_180 & -CR_030_210 & +c_right_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_02.region            = +CR_030_210 & -CR_060_240 & +c_upper_right_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_03.region            = +CR_060_240 & +CR_090_270 & +c_upper_right_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_04.region            = -CR_090_270 & +CR_120_300 & +c_upper_left_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_05.region            = -CR_120_300 & +CR_150_330 & +c_upper_left_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_06.region            = -CR_150_330 & +CR_000_180 & -c_left_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_07.region            = +CR_030_210 & -CR_000_180 & -c_left_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_08.region            = +CR_060_240 & -CR_030_210 & -c_lower_left_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_09.region            = -CR_090_270 & -CR_060_240 & -c_lower_left_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_10.region            = +CR_090_270 & -CR_120_300 & -c_lower_right_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_11.region            = +CR_120_300 & -CR_150_330 & -c_lower_right_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex
cr_12.region            = +CR_150_330 & -CR_000_180 & +c_right_1 & -core_out_2 & +bottom_0 & -top_0 #& ~core_hex

cr_01.fill              = a1
cr_02.fill              = a2
cr_03.fill              = a3
cr_04.fill              = a4
cr_05.fill              = a5
cr_06.fill              = a6
cr_07.fill              = a7
cr_08.fill              = a8
cr_09.fill              = a9
cr_10.fill              = a10
cr_11.fill              = a11
cr_12.fill              = a12

root.add_cells([core_reg,core_reg_out, cr_01, cr_02, cr_03, cr_04, cr_05, cr_06, cr_07, cr_08, cr_09, cr_10, cr_11, cr_12])
geometry                = openmc.Geometry(root)
geometry.export_to_xml()
#=================================================================================================
#                                   settings.xml file
#=================================================================================================
settings_file               = openmc.Settings()
settings_file.batches       = batches
settings_file.inactive      = inactive
settings_file.particles     = particles
settings_file.temperature   = {'method': 'interpolation'}
bounds                      = [-60.0, -60.0, -0.5, 60.0, 60.0, 0.5]
uniform_dist                = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source        = openmc.source.Source(space=uniform_dist)
settings_file.output        = {'tallies': False}
settings_file.sourcepoint   = {'write': False}
settings_file.export_to_xml()
#=================================================================================================
#                                     plots.xml File
#=================================================================================================
p1  = openmc.Plot(plot_id=1)
p1.width    = [250.0, 250.0]
p1.pixels   = [3600, 3600]
p1.basis    = 'xy'
p1.color_by = 'material'
plot_file = openmc.Plots([p1])
plot_file.export_to_xml()
openmc.plot_geometry()
#=================================================================================================
#                                      Initialize and run depletion calculation
#=================================================================================================
chain_file = '/home/abderi/cross_sections/endfb-viii.0-hdf5/chain_endfb71_msr.xml'
chain      = openmc.deplete.Chain.from_xml(chain_file)
model      = openmc.Model(geometry=geometry, settings=settings_file)
operator   = openmc.deplete.Operator(model, chain_file)
integrator = openmc.deplete.CECMIntegrator(operator, time_steps, power)
integrator.integrate()