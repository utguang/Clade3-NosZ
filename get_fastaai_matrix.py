import os
import re
import argparse

def main(result_folder):
    # Create a dictionary to store the genome names and AAI values
    genome_aai_dict = {}

    # List all files in the directory and sort them
    files = sorted(os.listdir(result_folder))

    # Iterate through each file
    for file_name in files:
        # Check if the file is a text file
        if file_name.endswith(".txt"):
            file_path = os.path.join(result_folder, file_name)
            genome_name = os.path.splitext(file_name)[0]

            # Open the file for reading
            with open(file_path, "r") as file:
                # Read the lines in the file
                lines = file.readlines()

                # Extract the AAI value from the last column of each line
                aai_values = [re.sub(r'[><%]', '', line.strip().split("\t")[-1]) for line in lines[1:]]

                # Add the genome name and AAI values to the dictionary
                genome_aai_dict[genome_name] = aai_values

            print("Processed:", file_name)

    # Create a new file to store the joined output
    output_file_path = os.path.join(result_folder, "joined_output.phylip")
    with open(output_file_path, "w") as output_file:
        # Write the total number of files as the first line
        output_file.write(str(len(files)) + "\n")

        # Iterate through the sorted genome names
        for genome_name in sorted(genome_aai_dict.keys()):
            # Get the corresponding AAI values from the dictionary
            aai_values = genome_aai_dict[genome_name]

            # Create a new line for each genome with the AAI values joined together
            joined_line = genome_name + "\t" + "\t".join(aai_values) + "\n"
            output_file.write(joined_line)

    print("Joined output file created:", output_file_path)

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Join AAI results into a single file")
    parser.add_argument("result_folder", help="Path to the folder containing AAI result files")
    args = parser.parse_args()

    # Call the main function with the specified result folder
    main(args.result_folder)
