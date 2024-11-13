#This script documents how to generate a genome comprising a single contig and polished with Illumina short-reads
canu -correct -p corrected_WS2 genomeSize=4m -pacbio-raw WS2.fastq.gz -correct useGrid=false
#1 assembly nanopore or PacBio long reads with FLYE
python /lustre/isaac/scratch/ghe3/software/Flye/bin/flye --nano-raw /lustre/isaac/scratch/ghe3/Sporomusa_hybrid_assembly/Sporomusa/fastq/SRX23778948_SRR28148124.fastq.gz --out-dir /lustre/isaac/scratch/ghe3/Sporomusa_hybrid_assembly/Sporomusa/fastq/FLYE --threads 56 --iterations 5

#2 Map the Illumina short-reads to FLYE generated genome using bwa
mamba activate assembly

bwa index polish1.fasta

bwa mem -t 56 -a polish1.fasta /lustre/isaac/scratch/ghe3/Sporomusa_hybrid_assembly/Sporomusa/fastq/SRX23778371_SRR28147534_R1.fastq.gz > align1.sam

bwa mem -t 56 -a polish1.fasta /lustre/isaac/scratch/ghe3/Sporomusa_hybrid_assembly/Sporomusa/fastq/SRX23778371_SRR28147534_R2.fastq.gz > align2.sam

#3 filter align.sam files
polypolish filter --in1 align1.sam --in2 align2.sam --out1 flter.align1.sam --out2 flter.align2.sam

#4 polish with the filtered sam files
polypolish polish --careful polish1.fasta flter.align1.sam flter.align2.sam > polish2.fasta

# clean index file to save disk space
rm *.amb *.ann *.bwt *.pac *.sa *.sam
