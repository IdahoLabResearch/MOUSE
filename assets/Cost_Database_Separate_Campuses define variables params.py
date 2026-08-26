import numpy as np
from cost.fleet_mode import (
    servicing_facility_allocation,
    servicing_facility_occ_learning_multipliers,
)

## User Inputs

params['Production Rate'] = 100 #Must be an integer.

params['Average Distance From Serv to GenSite'] = 1000 #miles (statutory miles)
params['Max Reactors Per Servicing Facility'] = 1000
params['Servicing Facility OCC Learning Rate'] = 0.30
params['Servicing Facility Learning Cap'] = 5
params['Cask Learning Rate'] = 0.15
params['Cask Learning Cap'] = 100

# ***
## END of User Inputs
# ***


## CONSTANTS

params['m3_to_kg_He_RT_atmospheric'] = 0.1661 #kg/m^3
params['Shift To Headcount'] = 5


## GENERATING SITE

params['CoolantInventoryRPV_Mass'] = 40.5125 #kg
params['Cycle Length'] = 3 #[years] This is a MOUSE ouput, and we're ignoring that value in order to be conservative.
params['Fuel Mass In Core'] = 408.890 #[kgU] This is a MOUSE output.
params['Water Supply Frequency'] = 4 #[year^-1]
params['Maintenance Visit Frequency'] = 1 / (params['Cycle Length'] / 2)
params['GenSite Downtime'] = 1/12

#To Discuss: Not Used in Cost Database at this time
params['GenSite Operators Per Shift'] = 1
params['GenSite Security Staff Per Shift'] = 1
params['GenSite Shifts Per Day'] = 2
params['GenSite Staff Rotation Working Fraction'] = 0.45


## FLEET

#ToDo: 9.231 and 10 should change with cycle length
params['Generating Sites Count'] = int(np.floor(9.231*params['Production Rate']))
params['Fleet'] = 10 * params['Production Rate']
#/ThisSectionToDo


## MANUFACTURING CAMPUS

params['Manufacturing Campus Area'] = 56.1651 +	1.2066 * params['Production Rate'] ** 0.3919

params['MFG Emergency Generator Power'] = -10851.2764 + 10851.2773 * (params['Production Rate'] ** 0.0001)
params['MFG Warehouse Building Area'] = 5598.3987 + 70.3374 * (params['Production Rate'] ** 0.5539)
params['MFG Administration Building Area'] = 3947.3106 + 77.6167 * (params['Production Rate'] ** 0.4768)
params['MFG Warehouse Staff'] = np.ceil( 0.8576 + 0.2539 * (params['Production Rate'] ** 0.6533) )
params['MFG Security Staff Per Shift'] = np.ceil( 2.9981 + 0.3340 * (params['Production Rate'] ** 0.4768) )
params['MFG Maintenance Staff Headcount'] = np.ceil( 2.9981 + 0.3340 * (params['Production Rate'] ** 0.4768) )
params['MFG Local Transport Vehicle Count'] = np.rint( -4340.5105 + 4340.5109 * (params['Production Rate'] ** 0.0001) )
params['MFG Utility Vehicle Count'] = np.rint( -4339.5105 + 4340.5109 * (params['Production Rate'] ** 0.0001) )
# params['MFG Guard Station Count'] = 1.0000 + 1e-15 * (params['Production Rate'] ** 5.0000)
if params['Production Rate'] <= 500:
    params['MFG Guard Station Count'] = 1

elif params['Production Rate'] < 2000:
    params['MFG Guard Station Count'] = 2

else:
    params['MFG Guard Station Count'] = np.rint(params['Production Rate']/1000)


params['MFG Campus Power'] = 1.295614 +	0.00813495 * params['Production Rate'] **  0.620945 #MWe
params['MFG Switchyard Rating'] = np.ceil(params['MFG Campus Power'] * 2)

params['MFG Testing Line Annual Rate'] = int(np.floor(8000/92))
params['MFG Testing Line Count'] = np.ceil(params['Production Rate'] / params['MFG Testing Line Annual Rate'])
params['MFG Testing Engineering Headcount'] = 6 + ((params['MFG Testing Line Count'] - 1)*2)
params['MFG Testing Engineering Headcount'] = 6 + ((params['MFG Testing Line Count'] - 1) * 2)
params['MFG Testing Coolant Allotment'] = 0.15 * (1 * 24.417 + 9 * 11.114) * 8.2402
params['MFG Testing Line Count'] = np.ceil(params['Production Rate'] / params['MFG Testing Line Annual Rate'])


#ToDo: Road Length and Site Perimeter shouldn't be the same.
params['MFG Road Length'] =  -1539388.0000975644 + 1540607.3153829477 * (params['Production Rate'] ** 8.588964624690402e-05)
params['MFG Site Perimeter'] = -1539388.0000975644 + 1540607.3153829477 * (params['Production Rate'] ** 8.588964624690402e-05)
params['MFG Protected Perimeter'] = 1421.7487567019996 + 60.41360643368738 * (params['Production Rate'] ** 0.22340396387855505)

params['MFG Security Camera Count'] = 200 + -157.0909090909111 * (params['Production Rate'] ** -0.037788560889399164)
params['MFG Motion Detector Count'] = 93.59999999999997 + 2.8444444444444628 * (params['Production Rate'] ** 0.35218251811136164)


## SERVICING CAMPUS

params['SER Construction Duration'] = 120 #months; carried over from the former central-facility assumption.
params['Servicing Rate'] = params['Fleet'] / params['Cycle Length']  # reactors/year
(
    params['Servicing Facility Count'],
    params['Servicing Facility Reactor Counts'],
    params['Servicing Facility Design Capacity'],
    params['Servicing Rate Per Facility'],
) = servicing_facility_allocation(
    params['Fleet'],
    params['Servicing Rate'],
    params['Max Reactors Per Servicing Facility'],
)
params['Servicing Facility OCC Learning Multipliers'] = (
    servicing_facility_occ_learning_multipliers(
        params['Servicing Facility Count'],
        params['Servicing Facility OCC Learning Rate'],
        params['Servicing Facility Learning Cap'],
    )
)

params['SER Campus Area'] = (760 - 160.5) * (params['Servicing Rate Per Facility'] - 30) / (300 - 30) + 160.5
params['SER Campus Land Area'] = params['SER Campus Area']

scale_var_SER = np.rint(params['Servicing Rate Per Facility'] / 3) #Value of 3 should remain hardcoded. Couldn't let production rate be the variable for the Servicing Campus directly, for flexibility, but the relationships were built based on production rate with our basic assumption of a 3:1 ratio between servicing rate and production rate. Simple fix was to define this variable (scale_var_SER).

params['SER Switchyard Rating'] = 0 + 12 * scale_var_SER ** 0.698970
params['SER Switchyard Average Power'] = params['SER Switchyard Rating']/2

params['Servicing Hot Cell Annual Rate'] = int(np.floor(365* (11/12) / 3 ))
params['Servicing Hot Cell Count'] = np.ceil( params['Servicing Rate Per Facility'] / params['Servicing Hot Cell Annual Rate'] )
params['Radioactive Waste Processing Hot Cell Count'] = 1
params['He Gas Replenishment Per Hot Cell'] = ((3*3*5) * 2 * params['Servicing Hot Cell Annual Rate'] * 2 + 0.1 * (10*30*7)*12) * params['m3_to_kg_He_RT_atmospheric']
params['He Gas Replenishment'] = (params['Servicing Hot Cell Count'] * params['He Gas Replenishment Per Hot Cell'] + params['Radioactive Waste Processing Hot Cell Count'] * params['He Gas Replenishment Per Hot Cell'] + params['CoolantInventoryRPV_Mass'] * params['Servicing Rate Per Facility'])

params['SER Number of Operators Per Shift'] = np.ceil( 0.0 + 5.625 * (scale_var_SER ** 0.426) )
params['SER Engineering Headcount'] = np.ceil( 0.0 + 20.0 * (scale_var_SER ** 0.301) )
params['SER Maintenance Staff Per Shift'] = 40  # REVIEW NEEDED: provisional value carried over from the old central-facility example.
params['SER Security Staff Per Shift'] = 1  # REVIEW NEEDED: provisional value carried over from the old reactor-site staffing input.

params['Roundtrip Time'] = 2*params['Average Distance From Serv to GenSite'] / 40 / (8760/2) #[years] based on average speed of 40 miles/hour, 12 hours of driving per day
params['Roundtrip Time Reactor Transport'] = 2*params['Average Distance From Serv to GenSite'] / 15 / (8760/2) #[years] based on average speed of 15 miles/hour, 12 hours of driving per day
params['Dwell Time GenSite'] = 1/365 #[years]
params['Dwell Time Serv'] = 1/365 #[years]
params['Dwell Time Reactor Transport GenSite'] = 2/365 #[years]
params['Dwell Time Reactor Transport Serv'] = 2/365 #[years]

#ToDo: Road Length and Site Perimeter shouldn't be the same.
params['SER Road Length'] = -1539388.0000975644 + 1540607.3153829477 * (scale_var_SER ** 8.588964624690402e-05)
params['SER Site Perimeter'] = -1539388.0000975644 + 1540607.3153829477 * (scale_var_SER ** 8.588964624690402e-05)
params['SER Controlled Perimeter'] = 0
params['SER Protected Perimeter'] = 1421.7487567019996 + 60.41360643368738 * (scale_var_SER ** 0.22340396387855505)

params['Used Fuel Storage Lifetime Capacity'] = params['Servicing Facility Design Capacity'] * 50 / (params['Cycle Length'] + params['GenSite Downtime']) * params['Fuel Mass In Core']

params['SER Security Building Area'] = 8775 / (3.2808 ** 2)
params['SER Administration Building Area'] = 258000 / (3.2808 ** 2)

params['SER Helium Flowrate'] = params['He Gas Replenishment'] / (0.9 * 8766 * 3600)

params['SER Local Transport Vehicle Count'] = 80  # REVIEW NEEDED: provisional value carried over from the old central-facility example.
params['SER Utility Vehicle Count'] = 100  # REVIEW NEEDED: provisional value based on the old general transport vehicle count.
params['Reactor Transport Vehicle Count'] = np.ceil(params['Servicing Rate Per Facility'] * (params['Roundtrip Time Reactor Transport']+params['Dwell Time Reactor Transport GenSite']+params['Dwell Time Reactor Transport Serv']) + 1)
params['Helium Transport Truck Count'] = np.ceil(params['Servicing Facility Design Capacity'] * params['Annual Coolant Supply Frequency'] * ((params['Roundtrip Time']+params['Dwell Time GenSite']+params['Dwell Time Serv']) * 1.05))
params['Water Tanker Truck Count'] = np.ceil(params['Servicing Facility Design Capacity'] * params['Water Supply Frequency'] * ((params['Roundtrip Time']+params['Dwell Time GenSite']+params['Dwell Time Serv']) * 1.05))
params['Maintenance Truck Count'] = np.ceil(params['Servicing Facility Design Capacity'] * params['Maintenance Visit Frequency'] * ((params['Roundtrip Time']+params['Dwell Time GenSite']+params['Dwell Time Serv']) * 1.05))

params['SER Radiation Monitor Count'] = np.ceil(params['SER Site Perimeter'] / 1000 * 1.2)
params['SER Security Camera Count'] = np.ceil(params['SER Campus Area'] / params['Manufacturing Campus Area'] * params['MFG Security Camera Count'])
params['SER Motion Detector Count'] = np.ceil(params['SER Campus Area'] / params['Manufacturing Campus Area'] * params['MFG Motion Detector Count'])

params['Reactor Transport Cask Count'] = round(params['Servicing Facility Design Capacity'] * ((4 / 12) / params['Cycle Length']) * 1.5, -1)
params['Annual Used Fuel Cask Consumption'] = 1 * params['Servicing Rate Per Facility']
params['Annual Reactor Cask Replacement'] = np.ceil(0.05 * params['Reactor Transport Cask Count'])
params['Annual Radwaste Cask Consumption'] = 0.5 * params['Servicing Rate Per Facility']

params['Servicing Hot Cell Building Area'] = 0.0 + 411.428571 * (scale_var_SER ** 0.942)
params['Helium Purification and Storage Building Area'] = 0.0 + 72.0 * (scale_var_SER ** 0.6198)
params['SER Local Control Building Area'] = 40 * scale_var_SER
params['SER Remote Control Building Area'] = 40 * scale_var_SER
params['SER Rad Waste Management Building Area'] = 0.0 + 40.5 * (scale_var_SER ** 0.7447)
# The source correlation returns ft^2; the cost database expects m^2.
params['Radwaste Storage Warehouses Area'] = (0.0 + 15422.94 * (scale_var_SER ** 0.994)) / (3.2808 ** 2)
params['SER Emergency Generator Power'] = 0.0 + 3.266667 * (scale_var_SER ** 0.632)
params['Parts Service Center and Warehouse Building Area'] = 0.0 + 675.0 * (scale_var_SER ** 0.6021)
params['Service Air Water Building Count'] = np.rint( 0.0 + 1.333333 * (scale_var_SER ** 0.1761) )
params['SER Fire Station Count'] = np.rint( 0.0 + 1.333333 * (scale_var_SER ** 0.1761) )
params['SER Guard Station Count'] = np.rint( 0.0 + 1.0 * (scale_var_SER ** 0.301) )
params['Helicopter Count'] = np.rint( 0.0 + 0.5 * (scale_var_SER ** 0.301) )
