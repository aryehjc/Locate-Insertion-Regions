import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

def calculate_gc_content(sequence):
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    total_bases = len(sequence)
    if total_bases == 0:
        return 0.0
    return ((g_count + c_count) / total_bases) * 100

def process_fasta(file_path, chunk_size=10000):
    gc_contents = []
    
    with open(file_path, 'r') as file:
        seq_record = SeqIO.read(file, "fasta")
        sequence = str(seq_record.seq)
        
        sequence = ''.join(filter(lambda x: x in 'GATC', sequence))
        
        for i in range(0, len(sequence), chunk_size):
            chunk = sequence[i:i + chunk_size]
            gc_content = calculate_gc_content(chunk)
            gc_contents.append(gc_content)
    
    return gc_contents

def get_gc_content_around_positions(sequence, positions, window_size=10000):
    gc_contents = []
    seq_length = len(sequence)
    
    for pos in positions:
        start = max(0, pos - window_size // 2)
        end = min(seq_length, pos + window_size // 2)
        
        region = sequence[start:end]
        gc_content = calculate_gc_content(region)
        gc_contents.append(gc_content)
    
    return gc_contents

def get_top_10_percent(df):
    num_top_rows = max(1, int(len(df) * 0.1)) 
    sorted_df = df.sort_values(by='Frequency', ascending=False)
    top_10_percent_df = sorted_df.head(num_top_rows)
    return top_10_percent_df

def plot_overlayed_boxplots(gc_contents_all, gc_contents_top_10, output_pdf_path):
    sns.set(style="whitegrid")
    plt.figure(figsize=(12, 8))
    
    # Combine data into a DataFrame for plotting
    data = {
        'GC Content': gc_contents_all + gc_contents_top_10,
        'Type': ['All Positions'] * len(gc_contents_all) + ['Top 10% Positions'] * len(gc_contents_top_10)
    }
    
    df = pd.DataFrame(data)
    
    # Plotting boxplots for both datasets
    sns.boxplot(x='Type', y='GC Content', data=df, palette={'All Positions': 'skyblue', 'Top 10% Positions': 'salmon'})
    
    # Adding plot details
    plt.title('GC Content Distribution and Top 10% Frequency Positions')
    plt.ylabel('GC Content (%)')
    
    # Perform the Wilcoxon rank-sum test
    stat, p_value = mannwhitneyu(gc_contents_all, gc_contents_top_10, alternative='two-sided')
    
    # Add test results to the plot
    plt.text(0.5, min(min(gc_contents_all), min(gc_contents_top_10)) * 1.05, 
             f'Wilcoxon rank-sum test statistic: {stat:.2f}\nP-value: {p_value:.4f}', 
             fontsize=12, fontname='FreeMono', ha='center', va='bottom', color='red')

    plt.savefig(output_pdf_path, format='pdf')
    plt.close()


# File paths
csv_file_path = 'FASTA.csv'
fasta_file_path = 'FASTA.fasta'
output_pdf_path = 'gc_content_overlay_with_wilcoxon.pdf'

# Read the CSV file
df = pd.read_csv(csv_file_path)

# Get top 10% frequency values and their associated positions
top_10_percent_df = get_top_10_percent(df)
positions = top_10_percent_df['Position'].astype(int).tolist()

# Read the FASTA file and calculate GC content
sequence = SeqIO.read(fasta_file_path, "fasta").seq
sequence = ''.join(filter(lambda x: x in 'GATC', str(sequence)))

# Calculate GC content for each 10,000-character chunk
gc_contents_all_positions = process_fasta(fasta_file_path)

# Calculate GC content around the top 10% positions
gc_contents_top_10_percent = get_gc_content_around_positions(sequence, positions)

# Plot and save the overlayed boxplots with Wilcoxon test results
plot_overlayed_boxplots(gc_contents_all_positions, gc_contents_top_10_percent, output_pdf_path)
