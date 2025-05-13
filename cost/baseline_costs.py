

# **************************************************************************************************************************
#                                                Sec. 0 :Inflation
# **************************************************************************************************************************

inflation_multipliers = [(	1985	,	2.351031042	),
(	1986	,	2.286828438	),
(	1987	,	2.217579774	),
(	1988	,	2.136971722	),
(	1989	,	2.046266296	),
(	1990	,	1.978944528	),
(	1991	,	1.915910254	),
(	1992	,	1.873063226	),
(	1993	,	1.821361194	),
(	1994	,	1.792614126	),
(	1995	,	1.752269334	),
(	1996	,	1.712641248	),
(	1997	,	1.681896823	),
(	1998	,	1.65443041	),
(	1999	,	1.633719786	),
(	2000	,	1.617597754	),
(	2001	,	1.604918477	),
(	2002	,	1.608313122	),
(	2003	,	1.601902774	),
(	2004	,	1.572518998	),
(	2005	,	1.542595467	),
(	2006	,	1.513543155	),
(	2007	,	1.477500092	),
(	2008	,	1.445814024	),
(	2009	,	1.447493593	),
(	2010	,	1.458139423	),
(	2011	,	1.432270014	),
(	2012	,	1.417006213	),
(	2013	,	1.416450085	),
(	2014	,	1.403505732	),
(	2015	,	1.391321967	),
(	2016	,	1.382535	),
(	2017	,	1.370060387	),
(	2018	,	1.348643837	),
(	2019	,	1.31987516	),
(	2020	,	1.285643232	),
(	2021	,	1.183842667	),
(	2022	,	1.047174995	),
(	2023	,	1	)]
inflation_dict = dict(inflation_multipliers)

# **************************************************************************************************************************
#                                                Sec. 1 : Baseline  Fixed Costs (dollars)
# **************************************************************************************************************************

### Sec. 1.1 Costs From MARVEL
plant_studies = ("Plant Studies", 5210451, 2023) 
# account 212 (island civil structure: mainly pit preparation)
pit_prep = ("Pit Preparation", 2573470, 2023) 



### Sec. 1.2 Costs From other sources
# from the WNA. (2023). Economics of Nuclear Power. . World Nuclear Association
#https://world-nuclear.org/information-library/economic-aspects/economics-of-nuclear-power#Licensingcosts
licensing_cost = ("Licensing Cost", 60e6, 2023) 


# **************************************************************************************************************************
#                                                Sec. 2 : Baseline  Unit Costs (dollars per unit)
# **************************************************************************************************************************


### Sec. 2.1 Unit Costs From MARVEL



### Sec. 2.2 Unit Costs From other sources

# from the EBD report
#https://inldigitallibrary.inl.gov/sites/sti/sti/Sort_46104.pdf
land_unit_cost = ("Land Unit Cost", 3800, 2021) # $/acre  in 2021
building_construction_cost = ("Building Construction Cost", 300, 2021) # $/ft^2


# **************************************************************************************************************************
#                                                Sec. 3 : Add all the costs to the cost dictionary
# **************************************************************************************************************************

# add all the costs to the cost dictionary
costs_list = [land_unit_cost, licensing_cost, plant_studies, pit_prep, building_construction_cost]

adjusted_cost_list = []
for item in costs_list:
    cost_factor = inflation_dict.get(item[2]) 
    adjusted_cost = cost_factor*item[1]
    adjusted_cost_tuple = (item[0], adjusted_cost)
    
    adjusted_cost_list.append(adjusted_cost_tuple)
cost_dictionary_data = dict(adjusted_cost_list)

