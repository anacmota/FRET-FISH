# FRET-FISH

## Probe design

The folder for probe design contains 3 scripts that are numbered from 1 to 3 as is the order they are supposed to be run. Before running the first script download the sequence of interest at UCSC, like for example, http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment=chrX:101640630,101742858 . Save files like "gene.txt". The first script "generate_oligos_1.m" will break the 100kb region into all possible oligos, name each oligo, evaluate the GC content and other parameters and then prints a list. The script now has the gene name of "Magix" and filename extension is ".txt", however this can be changed according to users preference. Once the oligo list file was generated from MATLAB, one should look for homologies in the reference genome (GRCm38 mouse in the case of "Magix") with BLAST command on Unix. In case your machine cannot handle the whole file without running out RAM then split the file into smaller ones with:

split -l 20000 Oligos_List.fa

Create a directory with the whole reference genome of the organism of interest:

makeblastdb -in RefGenome/GCF_000001635.26_GRCm38.p6_genomic.fna -out RefGenome/whole_mouse_genome.fa -dbtype nucl

Then, run BLAST to look for homologies.

blastn -query /Users/anamota/Oligos_List.txt -db RefGenome/whole_mouse_genome.fa -out /Users/anamota/homologies.txt -evalue 10 -word_size 11 -ungapped -perc_identity 80 -qcov_hsp_perc 80 -num_threads 4

Filter the BLAST file for Identities only:

awk '$1 ~ /^Query=/ {print $2} /^ Identities/ {print $3}' < /Users/anamota/homologies.txt > /Users/anamota/homologiesFilt.txt

Filter again for the specific homologies wanted:

sed -i 's#60/60^$##g' homologiesFilt.txt

Concatenate all the files into one:

cat *.txt > homologiesFiltAll.txt

The second MATLAB script "identities_penalty_2.m" will calculate the penalty total for each oligo. User must change line 5 for the gene of interest.

The third script "select_best_window_print_3.m" will generate a good set of oligos for FRET-FISH but first the user should confirm the variables in lines 16-20 regarding the spacing and number of oligos in groups as explained in the article. A question will pop-up whether the user is interested to include mouse barcodes in the flaps. If answered 1, the barcode file will be read. If answered 0, the oligo list will be generated without the barcodes. The barcode file "mouse barcodes.txt" can be replaced with other barcodes for instance when working with other organisms like human cells. The barcodes are intended to be used as primers for PCR amplification to selectively isolate a probe from an oligopool, therefore each gene/probe should have an unique combination of barcodes. The forward primer can be changed in line 162 with any number from 3 to 25 since 1 and 2 are already used for reverse primer to distiguish the donor and acceptor oligos respectively.


## Optimization of design: Linker sequence

Open LinkerDesign.m and in line 6 select which replicate is going to be shown. The data of each replicate is in "datasets" folder. The script will open which files correspond to each replicate.

## Optimization of design: Group and Spacing

Open SpacingGroupDesigns.m and in line 5 select which replicate is going to be shown. The data of each replicate is in "datasets" folder. The script will open which files correspond to each replicate.

## Nuclear Staining Measurement

The nd2 files coming from Nikon microscope are placed in a folder (due to their large size, we only provide the images by request) and the Fiji script fromND2_toINTENSITIES.ijm.ijm will measure the average of the staining intensity within the nucleus segmentation. The average of the intensity values per nuclei are exported in a csv file and the R script Treatment_nuclear_staining.R will process the data and show the comparison between the control and the treated cells in boxplots.

## Automatic dot picking

Open script named fp_main_gui.m . This script will open a window with some pre-setted values that can be changed. The resulting picked dots are saved in a newly created "data" folder.

setname: folder name where the images are located.

autoPick: the pipeline will select the 2 FRET pairs automatically if written "1", other number will consider already pre-selected dots.

xyres: resolution of the xy plane.

zres: resolution of the z plane.

maxdist_nm: threshold distance in nm in order to pair donor and acceptor molecules.

calcdolder: directory of the folder with segmented nuclei.

imagefolder: directory of the folder with the raw images in tif format.

donorchannel: donor channel name (e.g. a488).

acceptorchannel: acceptor channel name (e.g. a594).

fretchannel: FRET channel name (e.g. a488fret).

dapichannel: DAPI channel name (e.g. dapi).

## FRET-FISH analysis tools

**The user only has to run the following MATLAB scripts in the "FRET-FISH analysis" folder. For the scripts Passage_DAPI_boxplot.m and G1_G2_M_FRET_boxplot.m, the user needs to dowload the processed files at [KI OneDrive](https://kise-my.sharepoint.com/:f:/g/personal/ana_faustino_mota_ki_se/EvfboM9oEnlInmGG8N-mJg8BRat-T6TiX7LJMlvYJOFcOg?e=8XP8Is) and update line 34-37 with the correct directory.**

The analysis related to the ATP depletion treatment is found in the script ATP_depletion.m . The data used for this analysis is in the folder "data".

The prolonged cell culture analysis is shown in Passage_DAPI_boxplot.m and in Passage_FRET_distribution.m . For Passage_DAPI_boxplot.m , the analysis requires the processed images files that due to their large size is only provided by request. The remaining data is provided in the folder "data".

The cell cycle analysis is presented in the scripts G1_G2_FRET_distribution.m and in G1_G2_M_FRET_boxplot.m . In the first script shows the FRET probability distribution and the second script shows the boxplots for G1 High Hoechst, G1 Low Hoechst, G2 and M cell cycle phases. The data used for these analysis is in the folder "data".

The lamina distance influence in the chromatin density analysis is presented in the script Lamina_distance_boxplot.m . The output figure shows boxplot distributions of FRET efficiency for each concentric layer of the nucleus in relation to the lamina. The data used for this analysis is in the folder "data".

The analysis of the correlations in between the orthogonal techniques: FRET-FISH, ATAC-seq and Hi-C are shown in the script FRET_ATAC_HiC_corr.m . The data used for this analysis is in the folder "data".
