import pandas as pd
import matplotlib.pyplot as plt

# Load the data for M. Abs
unique_counts_df1 = pd.read_csv('Plot1.csv')

plt.figure(figsize=(10, 6))  # Adjust the figure size if needed
plt.hist(unique_counts_df1['Position'], bins=510, weights=unique_counts_df1['Frequency'], color='skyblue', edgecolor='green', alpha=0.7, density=True)

# Load the data for M. Mas
unique_counts_df2 = pd.read_csv('Plot2.csv')

plt.hist(unique_counts_df2['Position'], bins=500, weights=unique_counts_df2['Frequency'], color='skyblue', edgecolor='blue', alpha=0.7, density=True)

plt.xlabel('Position')
plt.ylabel('Frequency')
plt.title('Insertion Frequencies by Position')

# Create a legend based on edge colors
plt.legend(['M. Abs', 'M. Mas'], facecolor='white')  # Use labels based on edge colors

plt.grid(True)
plt.show()


#Change edgecolor for each plot to produce 2 different sets of graphs for overlay
