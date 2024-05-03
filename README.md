# Locate-Insertion-Regions
Find possible phage insertions in genome assembly. Must align reference genome with genome assemblies using MUMMer dnadiff, from there convert FASTAs to csv and put all CSVs in a directory to run this program. 


Getting complete genomes from NCBI first
mkdir genome_name
cd genome_name/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
grep 'Genome Name' assembly_summary_genbank.txt     | awk 'BEGIN{FS="\t"}{if($12=="Complete Genome"){print $20}}'     | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.fna.gz"}'     > urls.txt
ls
gedit urls.txt 
IFS=$'\n'; for NEXT in $(cat urls.txt); do wget "$NEXT"; done

./Align_Script.py, for i in *.delta; do ./delta2fasta.py $i alignments_fasta/`basename $i .fasta.delta`_delta.fa ; done 
for i in *.fa; do ./rubyscript < $i > $i.csv; done
Use output csvs.

# Running the below in-house scripts, before using Location-Insertion-Regions, which I wrote.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
===Align_Script===
#!/usr/bin/env python3
import sys
import subprocess as sp
import os

reference_genome=[
'GENOME_NAME.fasta',
        ]
    

all_assemblies=[]
f=open("assemblies.txt").readlines()
for i in f:
    all_assemblies.append(i.strip())

# if theres a slash after -p %s/ below its weird for some reason
       
for n,rg in enumerate(reference_genome):
    refid=rg.split("_")[0]
    for m,strains in enumerate(all_assemblies):
        if rg != strains:
            strains_id=strains.split('.')[0]
            sp.call('dnadiff FASTA/%s FASTA/%s -p %s_%s_%s'%(rg,strains,refid,refid,strains_id),shell=True)
====================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
===delta2fasta===
#!/usr/bin/env python3
import sys
import os
import subprocess
import re

def delta_to_fasta(delta_filename, output_filename=None):
    if output_filename is None:
        assert len(delta_filename.rsplit('.')) == 2, 'Give the delta file a file extension so that a basename can be identified and used for the output filename, otherwise supply a filename of your choosing with the `output_filename` parameter.'
        delta_file_basename = delta_filename.rsplit('.')[0]

        fasta_filename = '{}.fasta'.format(delta_file_basename)
    else:
        fasta_filename = output_filename

#    assert os.path.isfile(fasta_filename) == False, 'The filename {} is already used on the system. Please move or rename that file, or choose a different filename using the `output_filename` parameter.'.format(fasta_filename)

    coords_file = subprocess.run(['show-coords', '-c', '-l', '-r', '-T', delta_filename],
                                  stdout=subprocess.PIPE).stdout.decode('utf-8')

    coords_file_contents = coords_file.split('\n')[4:-1]

    seq_names = None
    alignments = []

    fasta_strings = []

    for line in coords_file_contents:
        pct_identity = line.split('\t')[6]

        if seq_names != tuple(line.split('\t')[-2:]): # this if clause should also work correctly in base case
#            assert len(alignments) == 0
            seq_names = tuple(line.split('\t')[-2:])
            first_seq_name, second_seq_name = seq_names

            aligns_file = subprocess.run(['show-aligns', delta_filename, first_seq_name, second_seq_name],
                                        stdout=subprocess.PIPE).stdout.decode('utf-8')
            alignments = alignments_from_aligns_file(aligns_file)

        alignment = alignments.pop(0)
        fasta_strings.append(sequences_lines_from_alignment(alignment, first_seq_name, second_seq_name, pct_identity))

    with open(fasta_filename, 'w') as fasta_file:
        fasta_file.write('\n\n'.join(fasta_strings))
    print('Finished, aligned sequences with coordinates and percent identities of matches should have been written in FASTA format to {}'.format(fasta_filename))

def alignments_from_aligns_file(aligns_file):
    file_contents = '\n'.join(aligns_file.split('\n')[3:-3]) # get rid of cruft in first 3 and last 3 lines of output
    alignments = file_contents.split('\n-- BEGIN alignment ')[1:]
    return alignments

def sequences_lines_from_alignment(alignment_string, first_seq_name='first_sequence', second_seq_name='second_sequence', pct_identity=None):
    header = alignment_string.split('\n\n')[0]
    lines = alignment_string.split('\n\n')[1:-1]
    footer = alignment_string.split('\n\n')[-1]

    coordinates_regex = r'\[ ([\+\-])1 ([0-9]+) \- ([0-9]+) \| ([\+\-])1 ([0-9]+) \- ([0-9]+) \]'
    first_strand_dir, first_seq_start_coord, first_seq_end_coord, second_strand_dir, second_seq_start_coord, second_seq_end_coord = re.findall(coordinates_regex, header)[0]


    #lines = lines.split('\n')
    #lines = ["\n".join(lines[i:i+4]) for i in range(0, len(lines), 4)]
    for n in range(len(lines)):
        lines[n] = lines[n].replace('^','').lstrip().rstrip()


#    of=open('tmp_file','w')
#    for n,line in enumerate(lines):
#        of.write("%d:%s"%(n,line.lstrip().rstrip()))
#    of.close()

    

    first_sequence = ""
    second_sequence = ""

    for n,line in enumerate(lines):
#        print(line.split('\n'))
        first_sequence_part, second_sequence_part = line.split('\n')
        remove_coord_numbers_regex = r'[0-9]+\s+(.*)'
        first_sequence_part = re.search(remove_coord_numbers_regex, first_sequence_part).group(1)
        second_sequence_part = re.search(remove_coord_numbers_regex, second_sequence_part).group(1)

        first_sequence += first_sequence_part
        second_sequence += second_sequence_part

    if pct_identity is not None:
        pct_identity_string = ' {}% identity'.format(pct_identity)
    else:
        pct_identity_string = ''

    first_fasta_string = '>{}:{}-{}({}){}\n{}'.format(first_seq_name, first_seq_start_coord, first_seq_end_coord,
                                                      first_strand_dir, pct_identity_string, first_sequence)
    second_fasta_string = '>{}:{}-{}({}){}\n{}'.format(second_seq_name, second_seq_start_coord, second_seq_end_coord,
                                                       second_strand_dir, pct_identity_string, second_sequence)

    fasta_string = '\n'.join([first_fasta_string, second_fasta_string])
    return fasta_string

print(delta_to_fasta(sys.argv[1],sys.argv[2]))
=========================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
===rubyscript===
#!/usr/bin/ruby

first_line = true

while line = STDIN.gets
  line.chomp!

  if line =~ /^>/
    puts unless first_line
    print line[1..-1]
    print ","  # <-- Change this to "\t" and it's a convert-fasta-to-tab
  else
    print line
  end

  first_line = false
end
puts
==================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Now using the CSVs, run this program. Annotations to be added shortly

