import os
import sys
import subprocess

# Usage: python Polypolish_iteration.py <input_genome.fasta> <read1.fastq.gz> <read2.fastq.gz> <iterations> <threads>

if len(sys.argv) != 6:
    print("Usage: python Polypolish_iteration.py <input_genome.fasta> <read1.fastq.gz> <read2.fastq.gz> <iterations> <threads>")
    sys.exit(1)

input_genome = sys.argv[1]  # Initial genome assembly (e.g., polish1.fasta)
read1 = sys.argv[2]         # Read 1 fastq file
read2 = sys.argv[3]         # Read 2 fastq file
iterations = int(sys.argv[4])  # Number of iterations for polishing
threads = sys.argv[5]       # Number of threads to use

current_genome = input_genome

for i in range(1, iterations + 1):
    print(f"Iteration {i}: Polishing {current_genome}...")

    # Step 0: Index the genome assembly before alignment
    print(f"Indexing the genome: {current_genome}")
    index_cmd = f"bwa index {current_genome}"
    subprocess.run(index_cmd.split())

    # Step 1: Align Read 1 separately using BWA
    align1_cmd = f"bwa mem -t {threads} -a {current_genome} {read1} > align1.sam"
    subprocess.run(align1_cmd, shell=True)

    # Step 2: Align Read 2 separately using BWA
    align2_cmd = f"bwa mem -t {threads} -a {current_genome} {read2} > align2.sam"
    subprocess.run(align2_cmd, shell=True)

    # Check if both alignment files were created
    if not os.path.exists("align1.sam") or not os.path.exists("align2.sam"):
        print("Error: align1.sam or align2.sam file not created. Check BWA alignment.")
        sys.exit(1)

    # Step 3: Filter the aligned SAM files using Polypolish
    filter_cmd = f"polypolish filter --in1 align1.sam --in2 align2.sam --out1 filtered.align1.sam --out2 filtered.align2.sam"
    subprocess.run(filter_cmd.split())

    # Check if filtered SAM files were created
    if not os.path.exists("filtered.align1.sam") or not os.path.exists("filtered.align2.sam"):
        print("Error: filtered.align1.sam or filtered.align2.sam not created. Check Polypolish filter step.")
        sys.exit(1)

    # Step 4: Polish the genome using Polypolish
    next_genome = f"polish{i + 1}.fasta"
    polish_cmd = f"polypolish polish --careful {current_genome} filtered.align1.sam filtered.align2.sam > {next_genome}"
    os.system(polish_cmd)

    # Clean up intermediate files: BWA index files and SAM files
    print(f"Cleaning up intermediate files for iteration {i}...")
    os.remove("align1.sam")
    os.remove("align2.sam")
    os.remove("filtered.align1.sam")
    os.remove("filtered.align2.sam")

    # Remove BWA index files
    for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
        index_file = current_genome + ext
        if os.path.exists(index_file):
            os.remove(index_file)

    current_genome = next_genome

    print(f"Finished polishing: {current_genome}")

print(f"Polishing completed after {iterations} iterations.")
print(f"Final polished genome: {current_genome}")
