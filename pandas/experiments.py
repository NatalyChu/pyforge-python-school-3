import pandas as pd

# Load the example data
wells_df = pd.read_csv('pandas/wells.csv')
plates_df = pd.read_csv('pandas/plates.csv')
experiments_df = pd.read_csv('pandas/experiments.csv')

# Pivot the data for each table so that property names become columns
wells_pivot = wells_df.pivot(index=['well_id', 'well_row', 'well_column', 'plate_id'], 
                             columns='property_name', 
                             values='property_value').reset_index()

plates_pivot = plates_df.pivot(index='plate_id', 
                               columns='property_name', 
                               values='property_value').reset_index()

experiments_pivot = experiments_df.pivot(index='experiment_id', 
                                         columns='property_name', 
                                         values='property_value').reset_index()

# Merge wells with plates to get experiment_id
wells_pivot = wells_pivot.merge(plates_df[['plate_id', 'experiment_id']], on='plate_id', how='left')

# Merge the tables
merged = wells_pivot.merge(plates_pivot, on='plate_id', how='left', suffixes=('_well', '_plate'))
merged = merged.merge(experiments_pivot, on='experiment_id', how='left', suffixes=('', '_exp'))

# Function to prioritize well properties over plate and experiment
def get_value(row, col_base):
    well_val = row.get(f"{col_base}_well")
    plate_val = row.get(f"{col_base}_plate")
    exp_val = row.get(f"{col_base}_exp")
    return well_val if pd.notnull(well_val) else plate_val if pd.notnull(plate_val) else exp_val

# Get all property columns
property_columns = [col.split('_')[0] for col in merged.columns if '_well' in col]

# Apply the prioritization logic for each property
for col in property_columns:
    merged[col] = merged.apply(get_value, col_base=col, axis=1)

# Drop the suffix columns (well/plate/exp versions of properties)
cols_to_drop = [col for col in merged.columns if col.endswith('_well') or col.endswith('_plate') or col.endswith('_exp')]
merged.drop(columns=cols_to_drop, inplace=True)

# Save the result to an Excel file
merged.to_csv('pandas/result.csv', index=False)

# Print the result for inspection
print(merged)
