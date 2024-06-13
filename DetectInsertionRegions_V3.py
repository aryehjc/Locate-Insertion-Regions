import pandas as pd
from pathlib import Path
import re
import numpy as np
# Still detects in relatively same region after amended logic, so it's possibly needing a new reference for circular representation!
# Assign path
position_counts = {}
dir = '/home/aryeh/example_directory'
csv_files = [f for f in Path(dir).glob('*.csv')]

# Function to find all dot group positions in a sequence
def find_dot_groups(string):
    dot_groups = [m.start() for m in re.finditer(r'\.+', string)]
    return dot_groups

# Function to calculate the positions of dot groups
def calculate_dot_group_positions(dot_groups):
    positions = []
    start_pos = 1
    for group in dot_groups:
        positions.append(start_pos)
        start_pos = group + 1  # Start position of next dot group
    return positions

# Function to add start position to list of dot group positions
def add_start_to_list(row):
    try:
        return [val + row['Start_Position'] for val in row['Dot_Groups']]
    except TypeError:
        return []

# Function to filter positions within specified range
def filter_positions(row):
    return [val for val in row if 1 < val < 2697195] # Upper bound total genome length on NCBI.

# Process each CSV file
for csv in csv_files:
    print(csv)
    df = pd.read_csv(csv, sep=',', header=None)
    
    # Add column names
    df.columns = ["Header", "Sequence"]
    
    # Clean Header column
    df["Header"] = df["Header"].str.replace('^.*:|\($','',regex=True)
    df["Header"] = df["Header"].str.split('(').str[0]
    
    # Find start positions
    df['Start_Position'] = df.Header.str.split('-', expand=True)[0].astype(float)
    
    # Expand header to ranges
    df["Header"] = df.Header.str.split('\s*-\s*').apply(lambda x: range(int(x[0]), int(x[-1]) + 1))
    
    # Filter ranges within bounds
    lower_bound = 1
    upper_bound = 2697195
    df["Res"] = df['Header'].apply(lambda x: any(lower_bound <= val <= upper_bound for val in x))
    df = df[df.Res]
    del df['Res']
    
    # Find dot group positions
    df['Dot_Groups'] = df['Sequence'].apply(find_dot_groups)
    df['Dot_Group_Positions'] = df['Dot_Groups'].apply(calculate_dot_group_positions)
    df['All_Positions'] = df.apply(add_start_to_list, axis=1)
    df['All_Positions'] = df['All_Positions'].apply(filter_positions)
    df['Total_Insertions'] = df['All_Positions'].apply(lambda x: len(set(x)))
    
    # Sum indel regions
    df['Indel Regions'] = df['Dot_Groups'].apply(len)
    Total_Indels = df['Indel Regions'].sum()
    df["Total_Indels"] = Total_Indels
    
    # Prepare for saving
    Full_Coordinates = df[['All_Positions']].copy()
    Insertion_Total_Count = len(set(df['All_Positions'].explode().dropna()))
    
    # Print for debugging
    print(Full_Coordinates)
    print(df)
    print(Total_Indels)
    print(Insertion_Total_Count)
    print(f'{csv.name} saved.')
    
    # Update position counts
    unique_positions = set([item for sublist in df['All_Positions'] for item in sublist])
    for pos in unique_positions:
        position_counts[pos] = position_counts.get(pos, 0) + 1
    
    # Save individual CSV file result
    output_filename = f'{csv.stem}_processed.csv'
    df.to_csv(Path(dir, output_filename), index=False)

    # Append to summary file
    with open("TOTAL_Jun_13_12pm_Insertion_Regions_S_AUREUS.txt", "a") as f:
        print({csv.name}, Insertion_Total_Count, file=f)

# Create a DataFrame from the position_counts dictionary
unique_counts_df = pd.DataFrame(list(position_counts.items()), columns=['Position', 'Frequency'])
unique_counts_df = unique_counts_df.sort_values(by='Position')

# Save the frequency DataFrame
unique_counts_df.to_csv('TOTAL_Jun_13_12pm_INSERTION_FREQUENCIES_S_AUREUS.csv', sep=',', index=False)
print(unique_counts_df)
