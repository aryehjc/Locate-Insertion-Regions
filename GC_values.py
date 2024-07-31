import pandas as pd
from Bio import SeqIO

# Step 1: Read the FASTA file and calculate GC content with a sliding window
def calculate_gc_content(fasta_file, window_size=100):
    # Read the sequence from the FASTA file
    record = SeqIO.read(fasta_file, "fasta")
    sequence = str(record.seq)
    
    gc_content = []
    seq_length = len(sequence)
    
    for i in range(seq_length):
        # Define the window bounds
        start = max(0, i - window_size // 2)
        end = min(seq_length, i + window_size // 2 + 1)
        
        # Extract the window and calculate GC content
        window = sequence[start:end]
        gc_count = window.count('G') + window.count('C')
        gc_content.append(gc_count / len(window))
    
    return gc_content

# Step 2: Integrate GC content with insertion data
def update_csv_with_gc_content(csv_file, fasta_file, output_file):
    # Read the CSV file with insertion data
    df = pd.read_csv(csv_file, sep=',')  # Adjust separator if needed
    
    # Calculate GC content
    gc_content = calculate_gc_content(fasta_file, window_size=100)
    
    # Add GC content to the DataFrame
    df['GC_Content'] = df['Position'].apply(lambda x: gc_content[x-1] if x-1 < len(gc_content) else None)
    
    # Write the updated DataFrame to a new CSV file
    df.to_csv(output_file, index=False)

# Paths to your files
csv_file = 'Test.csv' #csv generated from insertion script
fasta_file = 'Reference.fasta' #corresponding genome
output_file = 'updated_test.csv' #to run for box plot generation

# Run the update
update_csv_with_gc_content(csv_file, fasta_file, output_file)
