import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr  # Import pearsonr function for Pearson correlation coefficient

# Read the datasets
unique_counts_df1 = pd.read_csv('Plot1.csv')
unique_counts_df2 = pd.read_csv('Plot2.csv')

# Ensure both datasets have the same length
min_length = min(len(unique_counts_df1), len(unique_counts_df2))
unique_counts_df1_trimmed = unique_counts_df1.head(min_length)['Frequency']  # Extracting just the Frequency column
unique_counts_df2_trimmed = unique_counts_df2.head(min_length)['Frequency']  # Extracting just the Frequency column

# Calculate Pearson correlation coefficient
r, _ = pearsonr(unique_counts_df1_trimmed, unique_counts_df2_trimmed)

# Create scatter plot of y-values (Frequency)
plt.figure(figsize=(8, 6))

plt.scatter(unique_counts_df1_trimmed, unique_counts_df2_trimmed, color='green', alpha=0.6, label='M. Abs')
plt.scatter(unique_counts_df2_trimmed, unique_counts_df1_trimmed, color='blue', alpha=0.6, label='M. Mas')

plt.xlabel('Frequency in Plot1')
plt.ylabel('Frequency in Plot2')
plt.title('Frequency Comparison between Plot1 and Plot2')
plt.grid(True)

plt.text(0.05, 0.87, f"Pearson correlation coefficient (r): {r:.4f}", fontsize=10, family='FreeMono', weight='bold', color='red', transform=plt.gca().transAxes)

plt.legend()

plt.show()

print(f"Pearson correlation coefficient between unique_counts_df1 and unique_counts_df2: {r:.4f}")
