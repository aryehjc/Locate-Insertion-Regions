import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Step 1: Read the CSV file into a DataFrame
unique_counts_df = pd.read_csv('example.csv')

# Step 2: Calculate bin width using Freedman-Diaconis rule
# Calculate interquartile range (IQR)
Q1 = unique_counts_df['Position'].quantile(0.25)
Q3 = unique_counts_df['Position'].quantile(0.75)
IQR = Q3 - Q1

# Calculate optimal bin width using the Freedman-Diaconis rule
bin_width = 2 * IQR / np.power(len(unique_counts_df['Position']), 1/3)

# Step 3: Calculate optimal number of bins
# Calculate data range
data_range = unique_counts_df['Position'].max() - unique_counts_df['Position'].min()

# Determine number of bins
num_bins = int(data_range / bin_width)

# Step 4: Plot histogram and KDE
plt.figure(figsize=(10, 6))

# Plot histogram with optimal bins and weights
plt.hist(unique_counts_df['Position'], bins=num_bins, weights=unique_counts_df['Frequency'], color='skyblue', edgecolor='blue', density=True)

# Overlay KDE plot
sns.kdeplot(unique_counts_df['Position'], color='red', linestyle='-', linewidth=2)

# Add labels and title
plt.xlabel('Position')
plt.ylabel('Frequency')
plt.title('Insertion Frequencies by Position (Fixed-width Bins)')
plt.grid(True)

# Show plot
plt.show()

# Print calculated bin width and number of bins
print(f'Calculated Bin Width: {bin_width}')
print(f'Optimal Number of Bins: {num_bins}')
