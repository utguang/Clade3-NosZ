import os
import sys
import argparse

def extract_lines(input_path, output_directory, output_prefix, before_forward, after_forward, before_reverse, after_reverse, target):
    if os.path.isdir(input_path):
        for filename in os.listdir(input_path):
            if filename.endswith(".gff"):
                file_path = os.path.join(input_path, filename)
                process_file(file_path, output_directory, output_prefix, before_forward, after_forward, before_reverse, after_reverse, target)
    elif os.path.isfile(input_path):
        process_file(input_path, output_directory, output_prefix, before_forward, after_forward, before_reverse, after_reverse, target)
    else:
        print("Invalid input. Please provide a valid directory or file path.")

def process_file(input_file_path, output_directory, output_prefix, before_forward, after_forward, before_reverse, after_reverse, target):
    with open(input_file_path, 'r') as input_file:
        lines = input_file.readlines()

    output_lines = []

    for i, line in enumerate(lines):
        if target in line:
            orientation = line.split("\t")[6]  # Assuming orientation is in the 7th column
            if orientation == "+":
                # Extract the specified number of genes before and after for forward orientation
                start_index = max(0, i - before_forward)
                output_lines.extend(lines[start_index:i])

                # Extract the target line
                output_lines.append(line)

                end_index = min(i + after_forward + 1, len(lines))
                output_lines.extend(lines[i+1:end_index])
            elif orientation == "-":
                # Extract the specified number of genes before and after for reverse orientation
                start_index = max(0, i - before_reverse)
                output_lines.extend(lines[start_index:i])

                # Extract the target line
                output_lines.append(line)

                end_index = min(i + after_reverse + 1, len(lines))
                output_lines.extend(lines[i+1:end_index])

    # Build the output file path
    output_file_name = f"{output_prefix}_{os.path.basename(input_file_path)}"
    output_file_path = os.path.join(output_directory, output_file_name)

    # Write the extracted lines to the output file
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(output_lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract lines from GFF files.')
    parser.add_argument('-i', '--input', help='Input directory or file containing GFF files', required=True)
    parser.add_argument('-o', '--output', help='Output directory for extracted files', required=True)
    parser.add_argument('--prefix', help='Prefix for the output file', required=True)
    parser.add_argument('--before_forward', type=int, default=4, help='Number of genes right before for forward orientation')
    parser.add_argument('--after_forward', type=int, default=5, help='Number of genes right after for forward orientation')
    parser.add_argument('--before_reverse', type=int, default=5, help='Number of genes right before for reverse orientation')
    parser.add_argument('--after_reverse', type=int, default=4, help='Number of genes right after for reverse orientation')
    parser.add_argument('--target', help='Target line to search for in each GFF file', required=True)

    args = parser.parse_args()

    input_path = args.input
    output_directory = args.output
    output_prefix = args.prefix
    before_forward = args.before_forward
    after_forward = args.after_forward
    before_reverse = args.before_reverse
    after_reverse = args.after_reverse
    target = args.target

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Process input (directory or file)
    extract_lines(input_path, output_directory, output_prefix, before_forward, after_forward, before_reverse, after_reverse, target)
