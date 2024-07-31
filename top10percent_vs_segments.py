import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums

# Read the updated CSV file
df = pd.read_csv('test.csv')

# Step 1: Identify the top 10% of frequency values
top_10_percent_threshold = df['Frequency'].quantile(0.9)
top_10_percent_df = df[df['Frequency'] >= top_10_percent_threshold]

# Step 2: Calculate GC content for the top 10% frequency values
gc_content_top_10_percent = top_10_percent_df['GC_Content'].dropna()

# Step 3: Split the dataset into four equal parts
df_sorted = df.sort_values(by='Frequency').reset_index(drop=True)

# Number of rows for each segment
num_rows = len(df_sorted)
segment_size = num_rows // 4

# Create segments
segments = {
    'Segment 1': df_sorted.iloc[:segment_size]['GC_Content'].dropna(),
    'Segment 2': df_sorted.iloc[segment_size:2*segment_size]['GC_Content'].dropna(),
    'Segment 3': df_sorted.iloc[2*segment_size:3*segment_size]['GC_Content'].dropna(),
    'Segment 4': df_sorted.iloc[3*segment_size:]['GC_Content'].dropna(),
}

# Combine data for plotting
combined_data = [gc_content_top_10_percent] + list(segments.values())
labels = ['Top 10% Frequency'] + list(segments.keys())

# Perform Wilcoxon rank-sum tests and store statistics
test_results = []
for segment in segments.values():
    statistic, p_value = ranksums(gc_content_top_10_percent, segment)
    test_results.append((statistic, p_value))

# Plotting
plt.figure(figsize=(12, 6))
boxplot = plt.boxplot(combined_data, vert=False, patch_artist=True, labels=labels)
plt.title('GC Content Comparison: Top 10% Frequency vs Four Segments')
plt.xlabel('GC Content')
plt.grid(True)

# Annotate the plot with p-values and test statistics
for i, (statistic, p_value) in enumerate(test_results, start=1):
    plt.text(max(gc_content_top_10_percent), i + 0.5, f'Statistic = {statistic:.2f}\n p = {p_value:.3f}', 
             horizontalalignment='center', color='black', fontsize=10, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

plt.show()
