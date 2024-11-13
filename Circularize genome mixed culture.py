#This script documents how to generate a genome comprising a single contig and polished with Illumina short-reads
#1 Correct pacbio sequence
canu -correct -p corrected_WS2 genomeSize=4m -pacbio-raw WS2.fastq.gz -correct useGrid=false

#2 assembly nanopore or PacBio long reads with FLYE
python /lustre/isaac/scratch/ghe3/software/Flye/bin/flye --nano-raw /lustre/isaac/scratch/ghe3/Sporomusa_hybrid_assembly/Sporomusa/fastq/SRX23778948_SRR28148124.fastq.gz --out-dir /lustre/isaac/scratch/ghe3/Sporomusa_hybrid_assembly/Sporomusa/fastq/FLYE --threads 56 --iterations 5 --meta

# Close genome with circlator
circlator all --threads 48 --data_type pacbio-corrected flye_circular_SAB.fasta  /lustre/isaac/scratch/ghe3/A_PacBio_sequencing_three/Raw_data/fastq_file/Canu_SAB/corrected_SAB.correctedReads.fasta.gz Circlator

# Polish iterate 10 times
python Polypolish_iteration.py flye_circular_SAB.fasta /lustre/isaac/scratch/ghe3/Circularize_SAB_genome/13Sapos_S5_L001_R1_001.fastq.gz /lustre/isaac/scratch/ghe3/Circularize_SAB_genome/13Sapos_S5_L001_R2_001.fastq.gz 10 56


