import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr  # Import pearsonr function for Pearson correlation coefficient
#Where plot 1 and 2 are insertion frequencies for both genomes in comparison respectively
# First dataset
unique_counts_df1 = pd.read_csv('Plot1.csv')

plt.figure(figsize=(10, 6)) # Adjust the figure size if needed
plt.hist(unique_counts_df1['Position'], bins=510, weights=unique_counts_df1['Frequency'], color='skyblue', edgecolor='green', alpha=0.7, density=True)  # Alpha for transparency
sns.kdeplot(unique_counts_df1['Position'], color='red', linestyle='-', linewidth=2, bw_adjust=0.2)
plt.xlabel('Position')
plt.ylabel('Frequency')
plt.title('Plot1 Vs Plot2')
plt.grid(True)

# Calculate Pearson correlation coefficient 'r' for first dataset
r1, _ = pearsonr(unique_counts_df1['Position'], unique_counts_df1['Frequency'])


# Second dataset
unique_counts_df2 = pd.read_csv('Plot2.csv')

plt.hist(unique_counts_df2['Position'], bins=500, weights=unique_counts_df2['Frequency'], color='skyblue', edgecolor='blue', alpha=0.7, density=True)  # Alpha for transparency
sns.kdeplot(unique_counts_df2['Position'], color='orange', linestyle='-', linewidth=2, bw_adjust=0.2)
plt.title('Plot1 Vs Plot2')  # Adjust the title for the second dataset
plt.grid(True)

# Calculate Pearson correlation coefficient 'r' for second dataset
r2, _ = pearsonr(unique_counts_df2['Position'], unique_counts_df2['Frequency'])


plt.legend(['KDE M. Abs', 'KDE M. Mas'])  # Add legend to distinguish between the datasets


plt.text(0.05, 0.85, f"KL Divergence: {0.6219545736565859:.4f}", fontsize=10, family='FreeMono', weight='bold', color='red', transform=plt.gca().transAxes)
plt.text(0.05, 0.79, f"Pearson correlation coefficient M.Abs (r): {r1:.4f}", fontsize=10, family='FreeMono', weight='bold', color='green', transform=plt.gca().transAxes)
plt.text(0.05, 0.82, f"Pearson correlation coefficient M.Mas (r): {r2:.4f}", fontsize=10, family='FreeMono', weight='bold', color='blue', transform=plt.gca().transAxes)

plt.show()
