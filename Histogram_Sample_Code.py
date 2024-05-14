import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the CSV file
unique_counts_df = pd.read_csv('Insertion_Frequency.csv')

# Find the position and frequency of the highest y-value
max_frequency_index = unique_counts_df['Frequency'].idxmax()
max_frequency_position = unique_counts_df.loc[max_frequency_index, 'Position']
max_frequency = unique_counts_df.loc[max_frequency_index, 'Frequency']

# Plot the data
plt.figure(figsize=(10, 6))  # Adjust the figure size if needed
plt.plot(unique_counts_df['Position'], unique_counts_df['Frequency'])

# Mark the highest y-value on the plot
plt.scatter(max_frequency_position, max_frequency, color='red', label=f'Highest: ({max_frequency_position}, {max_frequency})')

plt.xlabel('Position')
plt.ylabel('Frequency')
plt.title('Insertion Frequencies by Position')
plt.grid(True)
plt.legend()
plt.show()
