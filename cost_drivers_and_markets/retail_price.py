# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import pandas as pd
import matplotlib.pyplot as plt

LMTR_FOAK_LCOE = 242.7				
CGMR_FOAK_LCOE = 377.5

LMTR_NOAK_LCOE = 131.5
CGMR_NOAK_LCOE = 191.5

# Replace 'file_path' with the path to your Excel file
file_path = 'retail_elec_price_2023.xlsx'

# Read the Excel file
df = pd.read_excel(file_path, sheet_name='Sheet1', usecols=['State', 'Average retail price (cents/kWh)'])  

df['retail price ($/MWh)'] = df['Average retail price (cents/kWh)'] * 10

#  Create a bar plot

fig, axs = plt.subplots(1, 1);         # Declaring subplots
sq_size = 16;                          # Declaring size
fig.set_size_inches(2*sq_size, 1*sq_size)   # Setting Size
fig.subplots_adjust(wspace=0.3, hspace = 0.6) # Spacing out plots
plt.rc('text', usetex=True)

axs=plt.subplot(1, 1,1)

# plt.bar(bins[:-1], relative_frequencies, color = colors_[y], align='edge',\
#         edgecolor='black', width=np.diff(bins)) # , ,

bars = plt.bar(df['State'], df['retail price ($/MWh)'], color='lawngreen',\
    align='center', edgecolor='black', label = "Electricity Retail Price, 2023")

# Adjust the x-axis labels to be above the x-axis
ax = plt.gca()
ax.xaxis.set_ticks_position('bottom')  # Move the x-axis labels to the top
ax.xaxis.set_label_position('bottom')  # Move the x-axis label to the top
ax.tick_params(axis='x', direction='in', pad=-5)  # Adjust the position of the labels inside the bars

plt.setp(ax.xaxis.get_majorticklabels(),rotation=90 , ha='center', va='bottom')


# # Center the labels above the bars
# for tick in ax.get_xticklabels():
#     tick.set_verticalalignment('bottom')
#     tick.set_x(0.5)  # Adjust this value to fine-tune the vertical position

plt.axhline(y=LMTR_FOAK_LCOE, color='rosybrown', linestyle='solid', linewidth=5,\
    label = "LMTR FOAK LCOE")
plt.axhline(y=LMTR_NOAK_LCOE, color='rosybrown', linestyle='dashed', linewidth=5,\
    label = "LMTR NOAK LCOE")
plt.axhline(y=CGMR_FOAK_LCOE, color='hotpink', linestyle='solid', linewidth=5,\
    label = "CGMR FOAK LCOE")
plt.axhline(y=CGMR_NOAK_LCOE, color='hotpink', linestyle='dashed', linewidth=5,\
    label = "CGMR NOAK LCOE")

# colors_ = ['rosybrown','tomato','orangered','sandybrown','darkorange',
#            'lawngreen','limegreen','turquoise','royalblue','mediumblue','mediumpurple','blueviolet','mediumorchid',
#            'violet','hotpink', 'crimson']

plt.grid(axis='y', which='both', color='grey', linestyle='dashed', linewidth=0.5)
plt.minorticks_on()
    
# Add labels and title
plt.xlabel('State', fontsize=38)
plt.ylabel('Price or LCOE ($/MWh)', fontsize= 40)
# Rotate the x-axis labels to avoid overlap
plt.xticks(fontsize = 28); 
plt.yticks(fontsize =36)


# Set the x-axis limits based on the indices of the labels
ax.set_xlim(-0.5, len(df['State']) - 0.5)

plt.legend(loc='upper right', bbox_to_anchor=(1, 0.94) ,fontsize =38, frameon=True, edgecolor='black')
# Save the plot as a PNG file
plt.savefig('retail_elec_dist.png')