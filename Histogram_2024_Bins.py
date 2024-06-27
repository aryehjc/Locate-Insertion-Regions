import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


unique_counts_df = pd.read_csv('example.csv')


plt.figure(figsize=(10, 6)) # Adjust the figure size if needed
plt.hist(unique_counts_df['Position'], bins=270, weights=unique_counts_df['Frequency'], color='skyblue', edgecolor='blue', density=True)  #Must be True for KDE , # can change bin number
sns.kdeplot(unique_counts_df['Position'], color='red', linestyle='-', linewidth=2)
plt.xlabel('Position')
plt.ylabel('Frequency')
plt.title('S. Aureus Insertion Frequencies by Position (Fixed-width Bins)')
plt.grid(True)
plt.show()
