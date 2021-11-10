# FRET-FISH

## Probe design

The folder for probe design contains 3 scripts that are numbered from 1 to 3 as is the order they are supposed to be run. Before running the first script download the sequence of interest at UCSC, like for example, http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment=chrX:101640630,101742858 . The first script has written the gene name of "Magix" and filename extension is ".txt", however this can be changed according to users preference. Once the oligo list file was generated, look for homologies in the reference genome on terminal. In case your machine cannot handle the whole file then split the file into smaller ones with:

split -l 20000 Oligos_List.fa

Create a directory with the whole reference genome of the organism of interest:

makeblastdb -in RefGenome/GCF_000001635.26_GRCm38.p6_genomic.fna -out RefGenome/whole_mouse_genome.fa -dbtype nucl

Then, run BLAST to look for homologies.

blastn -query /Users/anamota/FRET_probe_design/Oligos_List.txt -db RefGenome/whole_mouse_genome.fa -out /Users/anamota/FRET_probe_design/homologies.txt -evalue 10 -word_size 11 -ungapped -perc_identity 80 -qcov_hsp_perc 80 -num_threads 4

Filter the BLAST file for Identities only:

awk '$1 ~ /^Query=/ {print $2} /^ Identities/ {print $3}' < /Users/anamota/FRET_probe_design/homologies.txt > /Users/anamota/FRET_probe_design/homologiesFilt.txt

Filter again for the specific homologies wanted:

sed -i 's#60/60^$##g' homologiesFilt.txt

Concatenate all the files into one:

cat *.txt > homologiesFiltAll.txt

The second script will calculate the penalty total for each oligo. User must change line 5 for the gene of interest.

The third script will generate a good set of oligos for FRET-FISH but first the user should confirm the variables in lines 16-20 that the spacing and number of oligos in groups as explained in the article. A question will pop-up that is if the user is interested to include mouse barcodes in the flaps. If answered 1, the barcode file will be read. If answered 0, the oligo list will be generated without the barcodes. The barcode file can be replaced with other barcodes if the user has another list.


## Optimization of design: Linker sequence

Open LinkerDesign.m and in line 6 select which replicate is going to be shown. The data of each replicate is in "datasets" folder. The script will open which files correspond to each replicate.

## Optimization of design: Group and Spacing

Open SpacingGroupDesigns.m and in line 5 select which replicate is going to be shown. The data of each replicate is in "datasets" folder. The script will open which files correspond to each replicate.
