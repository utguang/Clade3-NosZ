import subprocess
import argparse

def run_bowtie2(genome_fasta, reads_fastq, output_sam, num_threads=1):
    # Bowtie2 Indexing
    subprocess.run(["bowtie2-build", genome_fasta, "genome_index"])

    # Bowtie2 Alignment with -p parameter
    subprocess.run(["bowtie2", "-x", "genome_index", "-U", reads_fastq, "-S", output_sam, "--no-unal", "-p", str(num_threads)])

def run_repeatexplorer(sam_input, db_output, mag_output, genome_fasta):
    # rpe Build
    subprocess.run(["rpe", "build", "-d", db_output, "-r", sam_input, "-g", genome_fasta, "--mag"])

    # rpe Plot
    subprocess.run(["rpe", "plot", "-d", db_output, "-w", "3000"])

def main():
    parser = argparse.ArgumentParser(description="Metagenomics Analysis Pipeline")
    parser.add_argument("--genome-fasta", required=True, help="Path to the genome fasta file")
    parser.add_argument("--reads-fastq", required=True, help="Path to the reads fastq file")
    parser.add_argument("--num-threads", type=int, default=1, help="Number of threads for Bowtie2 (default: 1)")
    args = parser.parse_args()

    output_sam = "output.sam"
    db_output = "culture.db"
    mag_output = "Culture_output.sam"

    run_bowtie2(args.genome_fasta, args.reads_fastq, output_sam, args.num_threads)
    run_repeatexplorer(output_sam, db_output, mag_output, args.genome_fasta)

if __name__ == "__main__":
    main()
