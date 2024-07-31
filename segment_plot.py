import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the dataset
unique_counts_df = pd.read_csv('test.csv') #from doubleplot kde pearson 

# Determine segment boundaries
num_rows = len(unique_counts_df)
segment_size = num_rows // 4
boundaries = [unique_counts_df['Position'].sort_values().iloc[i * segment_size] for i in range(1, 4)]

# Plot histogram and KDE
plt.figure(figsize=(10, 6))  # Adjust the figure size if needed
plt.hist(unique_counts_df['Position'], bins=510, weights=unique_counts_df['Frequency'], color='skyblue', edgecolor='blue', density=True)  # Must be True for KDE
sns.kdeplot(unique_counts_df['Position'], color='red', linestyle='-', linewidth=2) #change edgecolor to green etc for different genome

# Add vertical lines to show segment boundaries
for boundary in boundaries:
    plt.axvline(x=boundary, color='black', linestyle='--', linewidth=1.5, label=f'Segment Boundary')

plt.xlabel('Position')
plt.ylabel('Frequency')
plt.title('M. Mas Insertion Frequencies by Position (Fixed-width Bins)')
plt.grid(True)
plt.legend()
plt.show()
