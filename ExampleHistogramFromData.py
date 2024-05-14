import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the CSV file
unique_counts_df = pd.read_csv('Insertion_Frequencies.csv')

# Plot the data as a bar plot with fixed-width bins
plt.figure(figsize=(10, 6))  # Adjust the figure size if needed
plt.hist(unique_counts_df['Position'], bins=30, weights=unique_counts_df['Frequency'], color='skyblue', edgecolor='black')
plt.xlabel('Position')
plt.ylabel('Frequency')
plt.title('Insertion Frequencies by Position (Fixed-width Bins)')
plt.grid(True)
plt.show()
