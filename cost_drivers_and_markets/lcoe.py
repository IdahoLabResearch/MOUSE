import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
# Replace 'file_path' with the path to your Excel file
file_path = 'LMTR_CGMR_Summary.xlsx'

# Function to check if a value is a double digit integer
# Function to check if a value is a double digit integer
def is_double_digit_excluding_multiples_of_10(val):
    return isinstance(val, (int, float)) and val == int(val) and 10 <= int(val) <= 99 and int(val) % 10 != 0


def energy_cost_levelized( plant_lifetime_years, capital_cost, ann_cost, discount_rate, power_MWe, capacity_factor  ):
    
    sum_cost = 0 # initialization 
    sum_elec = 0
    
    for i in range(  1 + plant_lifetime_years) :
        
        if i == 0:
            cap_cost_per_year = capital_cost
            annual_cost = 0 # no annual cost at year 0
            elec_gen = 0 # no electricity generated
            
        elif i >0:
            cap_cost_per_year  = 0
            annual_cost = ann_cost
            elec_gen = power_MWe *capacity_factor * 365.25 * 24       # MW hour. 
        sum_cost +=  (cap_cost_per_year + annual_cost)/ ((1+ discount_rate)**i) 
        sum_elec += elec_gen/ ((1 + discount_rate)**i) 
    LCOE =  sum_cost/ sum_elec
    return LCOE


plant_lifetime_years = 60
discount_rate = 0.06
power_MWe_LMTR = 20*0.31
capacity_factor_LMTR = 0.93


def calculate_LCOE_LMTR_FOAK(row):
    if row['Account ID'] < 70:
        return energy_cost_levelized(plant_lifetime_years, row['LTMR [$]'],\
            0, discount_rate, power_MWe_LMTR, capacity_factor_LMTR)
    else:
        return energy_cost_levelized(plant_lifetime_years, 0, row['LTMR [$]'],\
            discount_rate, power_MWe_LMTR, capacity_factor_LMTR)

def calculate_LCOE_LMTR_NOAK(row):
    if row['Account ID'] < 70:
        return energy_cost_levelized(plant_lifetime_years, row['NOAK LRMT [$]'],\
            0, discount_rate, power_MWe_LMTR, capacity_factor_LMTR)
    else:
        return energy_cost_levelized(plant_lifetime_years, 0, row['NOAK LRMT [$]'],\
            discount_rate, power_MWe_LMTR, capacity_factor_LMTR)



power_MWe_CGMR = 15*0.4
capacity_factor_CGMR = 0.93

def calculate_LCOE_CGMR_FOAK(row):
    if row['Account ID'] < 70:
        return energy_cost_levelized(plant_lifetime_years, row['GCMR [$]'],\
            0, discount_rate, power_MWe_CGMR, capacity_factor_CGMR)
    else:
        return energy_cost_levelized(plant_lifetime_years, 0, row['GCMR [$]'],\
            discount_rate, power_MWe_CGMR, capacity_factor_CGMR)

def calculate_LCOE_CGMR_NOAK(row):
    if row['Account ID'] < 70:
        return energy_cost_levelized(plant_lifetime_years, row['NOAK GCMR [$]'],\
            0, discount_rate, power_MWe_CGMR, capacity_factor_CGMR)
    else:
        return energy_cost_levelized(plant_lifetime_years, 0, row['NOAK GCMR [$]'],\
            discount_rate, power_MWe_CGMR, capacity_factor_CGMR)






# Read the Excel file
df = pd.read_excel(file_path, sheet_name='Sheet1')  

# # Apply the function to filter the dataframe
filtered_df = df[df['Account ID'].apply(is_double_digit_excluding_multiples_of_10)]

# Assuming your DataFrame is named df
df1 = filtered_df[(filtered_df['LTMR [$]'] != 0) & (filtered_df['LTMR [$]'].notna())]
df1['LCOE FOAK LMTR'] = df1.apply(calculate_LCOE_LMTR_FOAK, axis=1)
print("\nSum of LCOE FOAK LMTR:",(df1['LCOE FOAK LMTR']).sum())

df1['LCOE NOAK LMTR'] = df1.apply(calculate_LCOE_LMTR_NOAK, axis=1)
print("\nSum of LCOE NOAK LMTR:",(df1['LCOE NOAK LMTR']).sum())



df1['LCOE FOAK CGMR'] = df1.apply(calculate_LCOE_CGMR_FOAK, axis=1)
print("\nSum of LCOE FOAK CGMR:",(df1['LCOE FOAK CGMR']).sum())

df1['LCOE NOAK CGMR'] = df1.apply(calculate_LCOE_CGMR_NOAK, axis=1)
print("\nSum of LCOE NOAK CGMR:",(df1['LCOE NOAK CGMR']).sum())


df2 = df1[~((df1['LCOE FOAK LMTR'] < 1) & (df1['LCOE NOAK LMTR'] < 1))]




# LMTR
df2_sorted = df2.sort_values(by='LCOE FOAK LMTR', ascending=False)

# df2_sorted['Title'] = df2_sorted['Title'].apply(lambda x: '\n'.join(x.split()))

# Assume df2 is your DataFrame
bar_width = 0.35  # Width of the bars

# Set the positions of the bars on the X-axis
r1 = np.arange(len(df2_sorted['Title']))
r2 = [x + bar_width for x in r1]

# Create the plot
plt.figure(figsize=(22, 20))
plt.bar(r1, df2_sorted['LCOE FOAK LMTR'], color='orangered',\
    width=bar_width, edgecolor='black', label='FOAK LMTR')
plt.bar(r2, df2_sorted['LCOE NOAK LMTR'], color='royalblue',\
    width=bar_width, edgecolor='black', label='NOAK LMTR')

plt.xticks([r+bar_width/2  for r in range(len(df2_sorted['Title']))],df2_sorted['Title'] ,\
    rotation = 45, ha='right')  # Rotate the x-axis tick labels
plt.xlim(-0.2, len(df2_sorted['Title']) - 0.5 )  # Adjust the limits to fit the bars properly

plt.grid(axis='y', which='both', color='grey', linestyle='dashed', linewidth=0.5)
plt.minorticks_on()

plt.legend(loc='upper right' ,fontsize =38, frameon=True, edgecolor='black')


plt.ylabel('LCOE ($/MWh)', fontsize= 40)
plt.xticks(fontsize = 24); 
plt.yticks(fontsize =32)
plt.savefig('lcoe_lmtr.png')







# CGMR

df_CGMR = df2.sort_values(by='LCOE FOAK CGMR', ascending=False)

# Set the positions of the bars on the X-axis
r1 = np.arange(len(df_CGMR['Title']))
r2 = [x + bar_width for x in r1]

# Create the plot
plt.figure(figsize=(22, 20))
plt.bar(r1, df_CGMR['LCOE FOAK CGMR'], color='orangered',\
    width=bar_width, edgecolor='black', label='FOAK CGMR')
plt.bar(r2, df_CGMR['LCOE NOAK CGMR'], color='royalblue',\
    width=bar_width, edgecolor='black', label='NOAK CGMR')


plt.xticks([r+bar_width/2  for r in range(len(df_CGMR['Title']))],df_CGMR['Title'] ,\
    rotation = 45, ha='right')  # Rotate the x-axis tick labels
plt.xlim(-0.2, len(df_CGMR['Title']) - 0.5 )  # Adjust the limits to fit the bars properly

plt.grid(axis='y', which='both', color='grey', linestyle='dashed', linewidth=0.5)
plt.minorticks_on()

plt.legend(loc='upper right' ,fontsize =38, frameon=True, edgecolor='black')


plt.ylabel('LCOE ($/MWh)', fontsize= 40)
plt.xticks(fontsize = 24); 
plt.yticks(fontsize =32)
plt.savefig('lcoe_cgmr.png')