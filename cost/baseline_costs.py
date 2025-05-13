
import pandas as pd
# **************************************************************************************************************************
#                                                Sec. 0 :Inflation
# **************************************************************************************************************************


def calculate_inflation_multiplier(file_path, dollar_year):
    # Read the Excel file
    df = pd.read_excel(file_path, sheet_name = "Inflation Adjustment")
    # print((dollar_year))
    # print(type(dollar_year))
    # Check if the target year exists in the dataframe
    if dollar_year in df['Year'].values:
        # Get the multiplier for the target year
        multiplier = df.loc[df['Year'] == dollar_year, 'Multiplier'].values[0]
        return multiplier
    else:
        return f"Year {dollar_year} not found in the Excel file."

# example
# calculate_inflation_multiplier("Cost_Database.xlsx", 1998)


# # **************************************************************************************************************************
# #                                                Sec. 1 : Baseline Costs (dollars)
# # **************************************************************************************************************************

def read_cost_database(file_name):
    """
    Reads an Excel file with a specified sheet name into a Pandas DataFrame.
    
    Parameters:
    file_name (str): The name of the Excel file.
    sheet_name (str): The name of the sheet to read from the Excel file.
    
    Returns:
    pd.DataFrame: A DataFrame containing the data from the specified sheet.
    """
    # Read the Excel file into a Pandas DataFrame
    df = pd.read_excel(file_name, sheet_name= "Cost Database")
    df['inflation_multiplier'] = df['Dollar Year'].apply(lambda x: calculate_inflation_multiplier(file_name, x))
    # Add the inflation-adjusted columns
    df['Fixed Cost - Inflation adjusted ($)'] = df['Fixed Cost ($)'] * df['inflation_multiplier']
    df['Unit Cost - Inflation adjusted ($)'] = df['Unit Cost ($)'] * df['inflation_multiplier']
    
    return df
# example
# print(read_cost_database("Cost_Database.xlsx"))









# ### Sec. 1.1 Costs From MARVEL
# plant_studies = ("Plant Studies", 5210451, 2023) 
# # account 212 (island civil structure: mainly pit preparation)
# pit_prep = ("Pit Preparation", 2573470, 2023) 



# ### Sec. 1.2 Costs From other sources
# # from the WNA. (2023). Economics of Nuclear Power. . World Nuclear Association
# #https://world-nuclear.org/information-library/economic-aspects/economics-of-nuclear-power#Licensingcosts
# licensing_cost = ("Licensing Cost", 60e6, 2023) 


# # **************************************************************************************************************************
# #                                                Sec. 2 : Baseline  Unit Costs (dollars per unit)
# # **************************************************************************************************************************


# ### Sec. 2.1 Unit Costs From MARVEL



# ### Sec. 2.2 Unit Costs From other sources

# # from the EBD report
# #https://inldigitallibrary.inl.gov/sites/sti/sti/Sort_46104.pdf
# land_unit_cost = ("Land Unit Cost", 3800, 2021) # $/acre  in 2021
# building_construction_cost = ("Building Construction Cost", 300, 2021) # $/ft^2


# # **************************************************************************************************************************
# #                                                Sec. 3 : Add all the costs to the cost dictionary
# # **************************************************************************************************************************

# # add all the costs to the cost dictionary
# costs_list = [land_unit_cost, licensing_cost, plant_studies, pit_prep, building_construction_cost]

# adjusted_cost_list = []
# for item in costs_list:
#     cost_factor = inflation_dict.get(item[2]) 
#     adjusted_cost = cost_factor*item[1]
#     adjusted_cost_tuple = (item[0], adjusted_cost)
    
#     adjusted_cost_list.append(adjusted_cost_tuple)
# cost_dictionary_data = dict(adjusted_cost_list)

