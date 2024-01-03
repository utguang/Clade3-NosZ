import os
import sys
import argparse

def extract_lines(input_file_path, output_directory, output_prefix, before, after, target):
    with open(input_file_path, 'r') as input_file:
        lines = input_file.readlines()

    output_lines = []

    for i, line in enumerate(lines):
        if target in line:
            # Extract the specified number of lines before
            start_index = max(0, i - before)
            output_lines.extend(lines[start_index:i])

            # Extract the target line
            output_lines.append(line)

            # Extract the specified number of lines after
            end_index = min(i + after + 1, len(lines))
            output_lines.extend(lines[i+1:end_index])

    # Build the output file path
    output_file_name = f"{output_prefix}_{os.path.basename(input_file_path)}"
    output_file_path = os.path.join(output_directory, output_file_name)

    # Write the extracted lines to the output file
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(output_lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract lines from GFF files.')
    parser.add_argument('-i', '--input', help='Input directory containing GFF files', required=True)
    parser.add_argument('-o', '--output', help='Output directory for extracted files', required=True)
    parser.add_argument('--prefix', help='Prefix for the output file', required=True)
    parser.add_argument('--before', type=int, default=4, help='Number of genes right before the target line')
    parser.add_argument('--after', type=int, default=5, help='Number of genes right after the target line')
    parser.add_argument('--target', help='Target line to search for in each GFF file', required=True)

    args = parser.parse_args()

    input_directory = args.input
    output_directory = args.output
    output_prefix = args.prefix
    before = args.before
    after = args.after
    target = args.target

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Process each file in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".gff"):
            input_file_path = os.path.join(input_directory, filename)
            extract_lines(input_file_path, output_directory, output_prefix, before, after, target)
