# Clade3-NosZ
Clade3 NosZ discovery codes

This repository contains code required for reproduce the analyses inlcuded for identification of clade III NosZ

# Desulfitobacterium genomes are available from NCBI database and are compiled in folder Desulfitobacterium genome.

# Functional annotation using emapper, making the annotation and gff file available.

for filename in *.fasta; do base=${filename}; emapper.py -m diamond --itype genome --no_file_comments -i $filename --output_dir emapper -o  $filename --decorate_gff yes --cpu 56 --override & done

# Using the gff_extract to withdraw the positions of 4 genes right before and 5 genes right after the target nosZ gene (K00376)

python gff_extract.py -i input_directory -o output directory --prefix NosZ --before 4 --after 5 --target K00376

