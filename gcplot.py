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

def process_fasta(file_path, chunk_size=10000):
    gc_contents = []
    
    # Read the FASTA file
    with open(file_path, 'r') as file:
        seq_record = SeqIO.read(file, "fasta")
        sequence = str(seq_record.seq)
        
        # Remove non-base characters (like spaces)
        sequence = ''.join(filter(lambda x: x in 'GATC', sequence))
        
        # Calculate GC content for each chunk
        for i in range(0, len(sequence), chunk_size):
            chunk = sequence[i:i + chunk_size]
            gc_content = calculate_gc_content(chunk)
            gc_contents.append(gc_content)
    
    return gc_contents

def save_gc_content_plot(gc_contents, output_pdf_path):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.boxplot(y=gc_contents)
    plt.title('GC Content Distribution')
    plt.ylabel('GC Content (%)')
    
    # Save the plot to a PDF file
    plt.savefig(output_pdf_path, format='pdf')
    plt.close()

# Replace 'your_file.fasta' with the path to your FASTA file
fasta_file_path = 'GENOME.fasta'
output_pdf_path = 'gc_content_distribution.pdf'
#i did this with NZ too
gc_contents = process_fasta(fasta_file_path)
save_gc_content_plot(gc_contents, output_pdf_path)
