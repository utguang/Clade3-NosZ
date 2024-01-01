import subprocess
import argparse
import os
import glob

def get_file_prefix(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

def run_bowtie2(genome_fasta, reads_fastq, output_dir, num_threads=1, keep_intermediate=False):
    # Extract prefix from the genome file
    genome_prefix = get_file_prefix(genome_fasta)

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Convert paths to absolute paths
    genome_fasta = os.path.abspath(genome_fasta)
    reads_fastq = os.path.abspath(reads_fastq)

    # Bowtie2 Indexing
    genome_index = os.path.join(output_dir, f"{genome_prefix}_index")
    subprocess.run(["bowtie2-build", genome_fasta, genome_index], check=True)

    # Bowtie2 Alignment with -p parameter
    output_sam = os.path.join(output_dir, f"{genome_prefix}_output.sam")
    subprocess.run(["bowtie2", "-x", genome_index, "-U", reads_fastq, "-S", output_sam, "--no-unal", "-p", str(num_threads)], check=True)

    if not keep_intermediate:
        # Remove intermediate files starting with "genome_index"
        intermediate_files = glob.glob(os.path.join(output_dir, f"{genome_prefix}_index*"))
        for file in intermediate_files:
            os.remove(file)

def run_repeatexplorer(sam_input, db_output, mag_output, genome_fasta, output_dir, keep_intermediate=False):
    # Extract prefix from the genome file
    genome_prefix = get_file_prefix(genome_fasta)

    # Convert paths to absolute paths
    sam_input_path = os.path.join(output_dir, sam_input)
    genome_fasta = os.path.abspath(genome_fasta)

    # rpe Build
    db_output = os.path.join(output_dir, f"{genome_prefix}_{db_output}")
    mag_output = os.path.join(output_dir, f"{genome_prefix}_{mag_output}")
    subprocess.run(["rpe", "build", "-d", db_output, "-r", sam_input_path, "-g", genome_fasta, "--mag"], check=True)

    # rpe Plot
    subprocess.run(["rpe", "plot", "-d", db_output, "-w", "3000"], check=True)

    if not keep_intermediate:
        # Remove intermediate files starting with "genome_index"
        intermediate_files = glob.glob(os.path.join(output_dir, f"{genome_prefix}_index*"))
        for file in intermediate_files:
            os.remove(file)

        # Remove intermediate SAM file
        os.remove(sam_input_path)

def main():
    parser = argparse.ArgumentParser(description="Metagenomics Analysis Pipeline")
    parser.add_argument("--genome-fasta", required=True, help="Path to the genome fasta file")
    parser.add_argument("--reads-fastq", required=True, help="Path to the reads fastq file")
    parser.add_argument("--num-threads", type=int, default=1, help="Number of threads for Bowtie2 (default: 1)")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--keep-intermediate", action="store_true", help="Keep intermediate files")
    args = parser.parse_args()

    run_bowtie2(args.genome_fasta, args.reads_fastq, args.output_dir, args.num_threads, args.keep_intermediate)

    db_output = "culture.db"
    mag_output = "Culture_output.sam"
    run_repeatexplorer(f"{get_file_prefix(args.genome_fasta)}_output.sam", db_output, mag_output, args.genome_fasta, args.output_dir, args.keep_intermediate)

if __name__ == "__main__":
    main()
