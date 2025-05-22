import numpy as np
import pandas as pd
# A function for the cost inflation needs to be added here
# interest cost
def TCI_estimate(debt_to_equity_ratio ,OCC, construction_duration, interest_rate):
    
    # Interest rate from this equation (from Levi)
    B =(1+ np.exp((np.log(1+ interest_rate)) * construction_duration/12))
    C  =((np.log(1+ interest_rate)*(construction_duration/12)/3.14)**2+1)
    Interest_expenses = debt_to_equity_ratio*OCC*((0.5*B/C)-1)
    return Interest_expenses



def energy_cost_levelized( plant_lifetime_years, cap_cost, ann_cost, discount_rate, power_MWe, capacity_factor  ):
    
    sum_cost = 0 # initialization 
    sum_elec = 0
    
    for i in range(  1 + plant_lifetime_years) :
        
        if i == 0:
        # assuming that the cap cost is split between the cons years
            cap_cost_per_year = cap_cost
            annual_cost = 0
            elec_gen = 0
            
        elif i >0:
            cap_cost_per_year  = 0
            annual_cost = ann_cost
            elec_gen = power_MWe *capacity_factor * 365 * 24       # MW hour. 
        sum_cost +=  (cap_cost_per_year + annual_cost)/ ((1+ discount_rate)**i) 
        sum_elec += elec_gen/ ((1 + discount_rate)**i) 
    LCOE =  sum_cost/ sum_elec
    return LCOE







def prettify(database1, caption):
    database1['Estimated Cost (2024 $)'] = database1['Estimated Cost (2024 $)'].apply(lambda x: '$ {:,.0f}'.format(x))
 
    
    database1 = database1.fillna('') # remove nan
    database1= database1.replace({'$ nan':'', 'nan hrs':''}) # remove nan

    database1[["Account"]] = database1[["Account"]].astype(str) 
    df_ = database1.loc[database1.Account.isin(['10', '20', '30', '40', '50', '60', '70', '80'])] # highlight high level accounts
    df_2 = database1.iloc[-8:] # final results  # df_2 = database1.iloc[-10:]
    df3 = pd.concat([df_, df_2])
    df4 = df_2 = database1.iloc[-3:]
    slice_ = pd.IndexSlice[df3.index, df3.columns]
    slice_2 = pd.IndexSlice[df4.index, df4.columns]

    slice_ = pd.IndexSlice[df3.index, df3.columns]
    slice_2 = pd.IndexSlice[df4.index, df4.columns]
    
    database_styled = (database1.style.set_properties(**{'font-weight': 'bold','color': 'black','background-color': 'pink' }, subset=slice_).\
                           set_properties(**{'font-weight': 'bold','color': 'black','background-color': 'pink' }, subset=slice_2).\
                           set_caption(caption).set_table_styles([{'selector': 'caption','props': [('color', 'red'),('font-size', '20px')]}]))
            
    return database_styled.hide()




# Estimate the cost factor to ajust for inflation
# Read the Excel file into a DataFrame


# Find the row that matches the observation_date
# row = df[df['observation_date'] == 'cost']

    # # Check if the row exists and return the value, otherwise return None
    # if not row.empty:
    #     return row['Industrial Machinery Manufacturing'].values[0]
    # else:
    #     print(f"This dollar year {dollar_year} is not found")





    
# # # # All these costs need to be adjusted for inflation at some point (to 2023$)
# def adjust_for_inflation(excel_file, cost_and_dollar_year):
    
#     # Read the Excel file into a DataFrame
#     df_cost_factors = pd.read_excel(excel_file)
    
#     # Check if the input is a tuple
#     if not isinstance(cost_and_dollar_year, tuple):
#         raise TypeError("The cost parameter must be a tuple")
#     parameter_name, raw_cost, dollar_year  = cost_and_dollar_year
    
    
#     # Find the row that matches the observation_date
#     row = df_cost_factors[df_cost_factors['observation_date'] == dollar_year]

#     cost_factor = row['Industrial Machinery Manufacturing'].values[0]
    
#     adjusted_cost  = raw_cost* cost_factor # 2023 dollars
#     return (parameter_name, adjusted_cost)
def custom_round(number):
    if number >= 1:
        rounded_number = round(number, 1)
    else:
        # Find the first non-zero digit after the decimal point
        str_number = f"{number:.10f}"  # Convert to string with high precision
        non_zero_index = next(i for i, ch in enumerate(str_number[2:], start=2) if ch != '0')
        rounded_number = round(number, non_zero_index - 1)
    return rounded_number