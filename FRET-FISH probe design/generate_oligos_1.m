clear all
clc
% To download sequence of interest go to UCSC, for example:
% http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment=chrX:101640630,101742858)

gene = 'Magix'; % filename is Magix.txt
geneSeq = fileread([gene '.txt']); % reference sequence from mouse
geneSeq = regexprep(geneSeq,'[^\w'']',''); % remove white spaces
geneSeq = upper(geneSeq); % upper case all nucleotides

oligoSize = 60;
Oligo = [];
mer = 1;

mkdir(gene) % create a directory/folder for the corresponding gene
cd(gene)
fileID = fopen([gene '_oligos_list.fa'], 'w'); % name for the oligos list


while length(geneSeq) > mer + oligoSize
    % Calculate the GC content
    Oligo = geneSeq((1:oligoSize)+mer);
    num_G = nnz(Oligo == 'G');
    num_C = nnz(Oligo == 'C');
    gc_percentage = round( (num_G + num_C) / oligoSize * 100 );
    
    if gc_percentage < 80 && gc_percentage > 35
        % Remove the oligos with 7 consecutive nucleic acid that are the same
        A = 'AAAAAAA';
        G = 'GGGGGGG';
        C = 'CCCCCCC';
        T = 'TTTTTTT';
        RemoveCons = contains(string(Oligo), A) + contains(string(Oligo), G) + contains(string(Oligo), C) + contains(string(Oligo), T);
        
        if RemoveCons == 0
            fprintf(fileID, '> mer_%d \n %s \n', mer, Oligo); % save the filtered oligomers
        else
            Identities(mer) = 1000; % defined penalty for 7 consecutive nucleotides
        end
        
    else
        Identities(mer) = 1000; % defined penalty for too high or too low GC content
    end
    mer=1+mer;
end

fclose(fileID);
Identities(length(Identities)+1:mer-1) = 0;
save('Identities', 'Identities')

% Run on terminal for mouse
% if your machine cannot handle the whole file then split the file into smaller ones with:
% split -l 20000 Oligos_List.fa
% Create a directory with the whole reference genome of the organism of interest:
% makeblastdb -in RefGenome/GCF_000001635.26_GRCm38.p6_genomic.fna -out RefGenome/whole_mouse_genome.fa -dbtype nucl
% blastn -query /Users/anamota/FRET_probe_design/Oligos_List.txt -db RefGenome/whole_mouse_genome.fa -out /Users/anamota/FRET_probe_design/homologies.txt -evalue 10 -word_size 11 -ungapped -perc_identity 80 -qcov_hsp_perc 80 -num_threads 4
% awk '$1 ~ /^Query=/ {print $2} /^ Identities/ {print $3}' < /Users/anamota/Desktop/FRET_probe_design/test/homologies.txt > /Users/anamota/Desktop/FRET_probe_design/test/homologiesFilt.txt

% sed -i 's#60/60^$##g' homologiesFilt.txt
% cat *.txt > homologiesFiltAll.txt