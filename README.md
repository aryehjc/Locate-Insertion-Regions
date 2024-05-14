# Locate-Insertion-Regions
This program was written by myself, Aryeh Chiam, to find possible phage insertions in genome assembly. Must align reference genome with genome assemblies using MUMMer dnadiff, from there convert FASTAs to csv and put all CSVs in a directory to run this program. Updated as of May 14 2024

First, run scripts in proteomics-scripts repo in order: Align_Script.py, delta2fasta.py, rubyscript.rb.
Requires BioPython and MUMMer

Getting complete genomes from NCBI first
```
mkdir genome_name
cd genome_name/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
grep 'Genome Name' assembly_summary_genbank.txt     | awk 'BEGIN{FS="\t"}{if($12=="Complete Genome"){print $20}}'     | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.fna.gz"}'     > urls.txt
ls
gedit urls.txt 
IFS=$'\n'; for NEXT in $(cat urls.txt); do wget "$NEXT"; done
```
```
./Align_Script.py
for i in *.delta; do ./delta2fasta.py $i alignments_fasta/`basename $i .fasta.delta`_delta.fa ; done
```
```
for i in *.fa; do ./rubyscript < $i > $i.csv; done
```
Use output csv directory in this program and set filepath in code.

Now using the CSVs, run this program. Annotations to be added shortly

