import pandas as pd

def remove_irrelevant_account(df, params):
    # Filter the dataframe
    indices_to_drop = []
    
    # Iterate over the rows and print the "Optional Variable" if it's not NaN
    for index, row in df.iterrows():
        if not pd.isna(row['Optional Variable']):
            if row['Optional Variable'] in params and params[row['Optional Variable']] == row['Optional Value']:
                    print("\n")
                    print(f"For the cost of the Account {row['Account']}: {row['Account Name']}, the {row['Optional Variable']} is selected to be {row['Optional Value']}")
            else:
                indices_to_drop.append(index)   
    
    # Drop the rows with collected indices (irrelevant accounts)
    df.drop(indices_to_drop, inplace=True)

    return df 



def find_children_accounts(df):
    # Find the column name that starts with "Estimated Cost"
    estimated_cost_column = [col for col in df.columns if col.startswith("FOAK Estimated Cost")][0]

    # Initialize a list for children accounts
    children_accounts = [None] * len(df)
    
    for target_level in range(4, -1, -1):
        source_level = target_level + 1
        # Iterate over the dataframe
        for i in range(len(df)):
            if df.iloc[i]['Level'] == target_level and pd.isna(df.iloc[i][estimated_cost_column]):
                children = []
                for j in range(i + 1, len(df)):
                    if df.iloc[j]['Level'] == source_level:
                        children.append(str(df.iloc[j]['Account']))  # Convert to string
                    elif df.iloc[j]['Level'] < source_level:
                        break
                # Convert the list to a comma-separated string
                children_str = ','.join(children) if children else None
                children_accounts[i] = children_str

    # Assign the list to the DataFrame
    df['Children Accounts'] = children_accounts
    return df    


def get_estimated_cost_column(df, option):
    if option == 'F':
        for col in df.columns:
            if col.startswith("FOAK Estimated Cost"):
                return col
    elif option == 'N'   :
        for col in df.columns:
            if col.startswith("NOAK Estimated Cost"):
                return col         
    return None