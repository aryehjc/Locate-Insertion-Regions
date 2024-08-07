import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_gc_content(sequence):
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    total_bases = len(sequence)
    if total_bases == 0:
        return 0.0
    return ((g_count + c_count) / total_bases) * 100
# change to 10000 and see if there's a major difference (there is, fatter 10% graph. I saved it separately)
def get_gc_content_around_positions(sequence, positions, window_size=10000):
    gc_contents = []
    seq_length = len(sequence)
    
    for pos in positions:
        # Ensure the position is within the valid range
        if pos < 0 or pos >= seq_length:
            print(f"Position {pos} is out of bounds.")
            continue
        
        # Calculate window boundaries
        start = max(0, pos - window_size // 2)
        end = min(seq_length, pos + window_size // 2)
        
        # Extract the region around the position
        region = sequence[start:end]
        
        # Calculate GC content for the region
        gc_content = calculate_gc_content(region)
        
        # Print for debugging
        if gc_content == 0:
            print(f"Zero GC Content detected at position {pos}: Sequence region: {region}")
        
        gc_contents.append(gc_content)
    
    return gc_contents



def get_top_10_percent(df):
    # Determine the number of top rows to select (10% of the total rows)
    num_top_rows = max(1, int(len(df) * 0.1))  # At least one row if data is very small
    
    # Sort the dataframe by 'Frequency' in descending order
    sorted_df = df.sort_values(by='Frequency', ascending=False)
    
    # Select the top 10% rows
    top_10_percent_df = sorted_df.head(num_top_rows)
    
    return top_10_percent_df

def plot_gc_content_boxplot(gc_contents, output_pdf_path):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.boxplot(y=gc_contents)
    plt.title('GC Content Around Top 10% Frequency Positions')
    plt.ylabel('GC Content (%)')
    
    # Save the plot to a PDF file
    plt.savefig(output_pdf_path, format='pdf')
    plt.close()

# File paths
csv_file_path = 'FASTA.csv'
fasta_file_path = 'FASTA.fasta'
output_pdf_path = 'BOXPLOT_10_PERCENT.pdf'

# Read the CSV file
df = pd.read_csv(csv_file_path)

# Get top 10% frequency values and their associated positions
top_10_percent_df = get_top_10_percent(df)
positions = top_10_percent_df['Position'].astype(int).tolist()

# Read the FASTA file
with open(fasta_file_path, 'r') as file:
    seq_record = SeqIO.read(file, "fasta")
    sequence = str(seq_record.seq)
    
    # Remove non-base characters (if any)
    sequence = ''.join(filter(lambda x: x in 'GATC', sequence))

# Calculate GC content around the top 10% positions
gc_contents = get_gc_content_around_positions(sequence, positions)

# Plot the GC content as a boxplot and save to PDF
plot_gc_content_boxplot(gc_contents, output_pdf_path)

