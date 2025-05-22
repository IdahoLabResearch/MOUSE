
import pandas as pd
from cost.cost_utils import *
# **************************************************************************************************************************
#                                                Sec. 0 :Inflation
# **************************************************************************************************************************


def calculate_inflation_multiplier(file_path, base_dollar_year, cost_type, escalation_year):
    # Read the Excel file
    df = pd.read_excel(file_path, sheet_name = "Inflation Adjustment")

    # Check if the target year exists in the dataframe
    if base_dollar_year in df['Year'].values:
        if  escalation_year in df['Year'].values:
            # Get the multiplier for the target year
            multiplier = df.loc[df['Year'] == base_dollar_year, cost_type].values[0] / df.loc[df['Year'] == escalation_year, cost_type].values[0]
            
            return multiplier
        else:  
            return f"Year {escalation_year} not found in the Excel file."  
    else:
        return f"Year {base_dollar_year} not found in the Excel file."

# example
# print(calculate_inflation_multiplier("Cost_Database.xlsx", 2021,'General', 2022))


# # **************************************************************************************************************************
# #                                                Sec. 1 : Baseline Costs (dollars)
# # **************************************************************************************************************************

def escalate_cost_database(file_name, escalation_year):
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
    df['inflation_multiplier'] = df.apply(lambda row: calculate_inflation_multiplier(file_name, row['Dollar Year'], row['Type'], escalation_year), axis=1)

    # Add the inflation-adjusted columns
    df['Adjusted Fixed Cost ($)'] = df['Fixed Cost ($)'] * df['inflation_multiplier']
    df['Adjusted Unit Cost ($)'] = df['Unit Cost ($)'] * df['inflation_multiplier']
    
    return df
# example
# a = escalate_cost_database("Cost_Database.xlsx", 2024)
# print(a['Adjusted Unit Cost ($)'])



def scale_cost(initial_database, params):
    scaled_cost = initial_database[['Account', 'Level','Account Title']]
    
    escalation_year = params['escalation_year']
    
        # Define the specific account values to be considered
    low_level_accounts = [11, 12, 13, 14, 15,\
        211.1, 211.2, 211.3,\
            212.1, 212.2, 212.3,\
                213.11, 213.12, 213.13,\
                    213.21, 213.22, 213.23,\
                        214.111, 214.112, 214.113,\
                            214.121, 214.122, 214.123,\
                                214.711, 214.712, 214.713, 214.72,\
                                    215.11, 215.12, 215.13,\
                                        215.41, 215.42, 215.43 ]


    # Iterate through each row in the DataFrame
    for index, row in initial_database.iterrows():
        scaling_variable_value = params[row['Scaling Variable']] if pd.notna(row['Scaling Variable']) else 0
        
        # Calculate the 'Estimated Cost
        fixed_cost = row['Adjusted Fixed Cost ($)'] if pd.notna(row['Fixed Cost ($)']) else 0
        unit_cost = row['Adjusted Unit Cost ($)'] if pd.notna(row['Unit Cost ($)']) else 0
        scaling_variable_ref_value  = row['Scaling Variable Ref Value']
        exponent = row['Exponent']
        
        if row['Account'] in low_level_accounts:
            # Check if the Account value is in the specific accounts list
            if row['Scaling Variable Ref Value'] > 0:
                estimated_cost = fixed_cost + unit_cost * pow(scaling_variable_value,exponent) /(pow(scaling_variable_ref_value,exponent-1))

            else:
                # Calculate the 'Estimated Cost
                estimated_cost = fixed_cost + unit_cost * scaling_variable_value

        

            # Assign the calculated value to the corresponding row in the new DataFrame
            scaled_cost.at[index, f'Estimated Cost (${escalation_year } M$))'] = np.round(estimated_cost/ 1000000,2)
    return scaled_cost  


# # A method to update the high level costs when the low level costs change
# def update_high_level_costs(params, db):
    

#     cost_account = f'Estimated Cost ({params['escalation_year']} $)'  
        
#     #update total costs for accounts 10
#     # power is in kW
#     power = params['Power MWe']  * 1000

#     db.loc[db.Account == 211,  cost_account] =\
#         db.loc[db['Account'].isin(['211A', '211B' ,'211C']),  cost_account].sum() 
     
#     db.loc[db.Account == 212,  cost_account] =\
#         db.loc[db['Account'].isin(['212A', '212B' ,'212C']),  cost_account].sum() 
#     db.loc[db.Account == 214,  cost_account] =\
#         db.loc[db['Account'].isin([214.7]),  cost_account].sum()  
        
    # db.loc[db.Account == 21,  cost_account] =\
    #     db.loc[db['Account'].isin([211, 212, 214]),  cost_account].sum() 
     
    
    # db.loc[db.Account == 221.21,  cost_account] =\
    #     db.loc[db['Account'].isin(['221.21A','221.21B', '221.21C', '221.21D' ]),  cost_account].sum()     
    
    # db.loc[db.Account == 221.1,  cost_account] =\
    #     db.loc[db['Account'].isin([221.11, 221.12, 221.13]),  cost_account].sum() 
    # db.loc[db.Account == 221.2,  cost_account] =\
    #     db.loc[db['Account'].isin([221.21]),  cost_account].sum()  
    
    # db.loc[db.Account == 221.32,  cost_account] =\
    #     db.loc[db['Account'].isin(['221.32A','221.32B' ]),  cost_account].sum()
        
    # db.loc[db.Account == 221.3,  cost_account] =\
    #     db.loc[db['Account'].isin([221.31, 221.32]),  cost_account].sum() 
    
    # db.loc[db.Account == 221,  cost_account] =\
    #     db.loc[db['Account'].isin([221.1, 221.2, 221.3]),  cost_account].sum()  
    

        
    # db.loc[db.Account == 222,  cost_account] =\
    #     db.loc[db['Account'].isin([222.1, 222.2, 222.3]),  cost_account].sum() 
     
            
    # db.loc[db.Account == 223,  cost_account] =\
    #     db.loc[db['Account'].isin([223.2]),  cost_account].sum()   
    

    # db.loc[db.Account == 22,  cost_account] =\
    #     db.loc[db['Account'].isin([221, 222, 223, 227]),  cost_account].sum()
      
    
    
             
    # db.loc[db.Account == 232,  cost_account] =\
    #     db.loc[db['Account'].isin([232.1]),  cost_account].sum()

    # db.loc[db.Account == 23,  cost_account] =\
    #     db.loc[db['Account'].isin([232]),  cost_account].sum() 
    
    # db.loc[db.Account == 254,  cost_account] =\
    #     db.loc[db['Account'].isin(['254A', '254B']),  cost_account].sum()  
        
    # db.loc[db.Account == 25,  cost_account] =\
    #     db.loc[db['Account'].isin([251, 252, 253, 254]),  cost_account].sum()  

    # # account 10  
    # db.loc[db.Account == 10,  cost_account] =\
    #     db.loc[db['Account'].isin([11, 12, 13, 14, 15]),  cost_account].sum()
        
    # # account 20  
    # db.loc[db.Account == 20,  cost_account] =\
    #     db.loc[db['Account'].isin([21, 22, 23, 24, 25, 26]),  cost_account].sum()    
        
    
    # # account 30 
    # db.loc[db.Account == 30,  cost_account] =\
    #     db.loc[db['Account'].isin([31, 32, 33, 34, 35, 36]),  cost_account].sum()       
                      
        
    # # account 60
    # db.loc[db.Account == 60,  cost_account] =\
    #     db.loc[db['Account'].isin([62]),  cost_account].sum()  
    
    # # account 70
    
    # db.loc[db.Account == 71,  cost_account] =\
    #     db.loc[db['Account'].isin([711, 712, 713]),  cost_account].sum() 
        
    # db.loc[db.Account == 70,  cost_account] =\
    #     db.loc[db['Account'].isin([71, 75, 78]),  cost_account].sum()   
        
        
        
    # # account 80
    # db.loc[db.Account == 80,  cost_account] =\
    #     db.loc[db['Account'].isin([81, 82, 83]),  cost_account].sum()        

    # # # update final results
    # db.loc[db.Account == 'Overnight Capital Cost (OCC)',  cost_account] =\
    #     db.loc[db['Account'].isin([10, 20, 30]),  cost_account].sum()
        
    # db.loc[db.Account == 'OCC ($/kW)',  cost_account] = db.loc[db.Account == 'Overnight Capital Cost (OCC)',  cost_account].values[0]/power 
    
    # db.loc[db.Account == 'Total Capital Investment (TCI)',  cost_account] = db.loc[db.Account == 'Overnight Capital Cost (OCC)',\
    #     cost_account].values[0] +\
    #     db.loc[db.Account == 60,  cost_account].values[0]
     
    # db.loc[db.Account == 'TCI ($/kW)',  cost_account] = db.loc[db.Account == 'Total Capital Investment (TCI)',\
    #     cost_account].values[0]/power 
       
    # db.loc[db.Account == 'Annualized Cost ($/ year)',  cost_account] = db.loc[db.Account == 70,  cost_account].values[0] +\
    #     db.loc[db.Account == 80,  cost_account].values[0]
    
    
    	
    # return db
    
def update_high_level_costs(df):
    
    cost_column = next(col for col in df.columns if col.startswith("Estimated Cost"))

    for level in range(4, -1, -1):
        
        # Process each level from 4 to 0
        for idx, row in df[df['Level'] == level].iterrows():
            
            if pd.isna(row[cost_column]):
                # Sum costs of accounts below the current account
                below_costs = 0
                for sub_idx, sub_row in df[idx+1:].iterrows():
                    if sub_row['Level'] <= level:
                        break
                    if sub_row['Level'] == level + 1 and not pd.isna(sub_row[cost_column]):
                        below_costs += sub_row[cost_column]
                df.at[idx, cost_column] = below_costs
              

    return df
    
    

def bottom_up_cost_estimate(cost_database_filename, params):
    escalated_cost = escalate_cost_database(cost_database_filename, params['escalation_year'])
    scaled_cost = scale_cost(escalated_cost, params)
    updated_cost = update_high_level_costs(scaled_cost)
    return updated_cost 



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

