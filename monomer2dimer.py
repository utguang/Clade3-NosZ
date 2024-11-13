import sys

def replicate_sequences(file_lines):
    formatted_sequences = []
    current_seq_id = ""
    current_sequence = []
    
    for line in file_lines:
        if line.startswith('>'):
            if current_seq_id and current_sequence:
                sequence = ''.join(current_sequence)
                formatted_sequences.append(current_seq_id)
                formatted_sequences.append(f"{sequence}:{sequence}")
            current_seq_id = line.strip()
            current_sequence = []
        else:
            current_sequence.append(line.strip())
    
    # Add the last sequence to the list
    if current_seq_id and current_sequence:
        sequence = ''.join(current_sequence)
        formatted_sequences.append(current_seq_id)
        formatted_sequences.append(f"{sequence}:{sequence}")
    
    return formatted_sequences

def main(input_file, output_file):
    # Read the content of the input file
    with open(input_file, 'r') as file:
        file_content = file.readlines()

    # Replicate and format sequences as required
    replicated_sequences = replicate_sequences(file_content)

    # Join the formatted sequences into a single string
    formatted_output = "\n".join(replicated_sequences)

    # Save the formatted output to a new file
    with open(output_file, 'w') as output_file:
        output_file.write(formatted_output)

    print(f"Formatted sequences have been saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        main(input_file, output_file)
