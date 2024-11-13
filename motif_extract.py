import re
import pandas as pd
import argparse

# Function to read FASTA file and return sequences as a dictionary
def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        seq_id = None
        seq = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    sequences[seq_id] = ''.join(seq)
                seq_id = line[1:]
                seq = []
            else:
                seq.append(line)
        if seq_id is not None:
            sequences[seq_id] = ''.join(seq)
    return sequences

# Function to write sequences to a FASTA file
def write_fasta(sequences, output_path):
    with open(output_path, 'w') as file:
        for seq_id, seq in sequences.items():
            file.write(f'>{seq_id}\n')
            file.write(f'{seq}\n')

# Function to find sequences matching all of the motifs
def find_matching_sequences(sequences, motif_patterns):
    matching_sequences = {}
    for seq_id, seq in sequences.items():
        if all(pattern.search(seq) for pattern in motif_patterns):
            matching_sequences[seq_id] = seq
    return matching_sequences

# Main function
def main(input_file, output_file, motifs):
    # Compile the regex patterns for the motifs
    motif_patterns = [re.compile(motif) for motif in motifs]

    # Read the sequences from the provided FASTA file
    sequences = read_fasta(input_file)

    # Find sequences matching all of the motifs
    matching_sequences = find_matching_sequences(sequences, motif_patterns)

    # Write matching sequences to the output FASTA file
    write_fasta(matching_sequences, output_file)

    # Display the extracted sequences in a dataframe
    df = pd.DataFrame(matching_sequences.items(), columns=['Sequence ID', 'Sequence'])
    print(df)

# Argument parser setup
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences containing all specific motifs from a FASTA file.")
    parser.add_argument("input_file", type=str, help="Input FASTA file path")
    parser.add_argument("output_file", type=str, help="Output FASTA file path")
    parser.add_argument("motifs", type=str, nargs='+', help="Motif sequence patterns (e.g., 'C..FC...H.EM' 'D.H.' 'PHG' 'GPLH' 'EPH')")
    
    args = parser.parse_args()
    main(args.input_file, args.output_file, args.motifs)
