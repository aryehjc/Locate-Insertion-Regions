import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr  # Import pearsonr function for Pearson correlation coefficient

# First dataset
unique_counts_df1 = pd.read_csv('Plot1.csv')

plt.figure(figsize=(10, 6)) # Adjust the figure size if needed
plt.hist(unique_counts_df1['Position'], bins=510, weights=unique_counts_df1['Frequency'], color='skyblue', edgecolor='green', alpha=0.7, density=True)  # Alpha for transparency
sns.kdeplot(unique_counts_df1['Position'], color='red', linestyle='-', linewidth=2, bw_adjust=0.2)
plt.xlabel('Position')
plt.ylabel('Frequency')
plt.title('Insertion Frequencies by Position (Plot1 Vs Plot2)')
plt.grid(True)

# Second dataset
unique_counts_df2 = pd.read_csv('Plot2.csv')

plt.hist(unique_counts_df2['Position'], bins=500, weights=unique_counts_df2['Frequency'], color='skyblue', edgecolor='blue', alpha=0.7, density=True)  # Alpha for transparency
sns.kdeplot(unique_counts_df2['Position'], color='orange', linestyle='-', linewidth=2, bw_adjust=0.2)
plt.title('Insertion Frequencies by Position (Plot1 Vs Plot2)')  # Adjust the title for the second dataset
plt.grid(True)

plt.legend(['KDE Plot1', 'KDE Plot2'])  # Add legend to distinguish between the datasets


plt.text(0.05, 0.85, f"KL Divergence: {0.6219545736565859:.4f}", fontsize=10, family='FreeMono', weight='bold', color='red', transform=plt.gca().transAxes)



# First dataset
unique_counts_df1 = pd.read_csv('Plot1.csv')

# Second dataset
unique_counts_df2 = pd.read_csv('Plot2.csv')

# Ensure both datasets have the same length
min_length = min(len(unique_counts_df1), len(unique_counts_df2))
unique_counts_df1_trimmed = unique_counts_df1.head(min_length)
unique_counts_df2_trimmed = unique_counts_df2.head(min_length)

# Calculate Pearson correlation coefficient
r, _ = pearsonr(unique_counts_df1_trimmed['Position'], unique_counts_df2_trimmed['Position'])

print(f"Pearson correlation coefficient between unique_counts_df1 and unique_counts_df2: {r:.4f}")
plt.text(0.05, 0.87, f"Pearson correlation coefficient (r): {r:.4f}", fontsize=10, family='FreeMono', weight='bold', color='red', transform=plt.gca().transAxes)


plt.show()
