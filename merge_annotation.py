import pandas as pd
import os

# Assuming your annotation files are in the same directory
annotation_files = [file for file in os.listdir('.') if file.endswith('.annotations')]

# Create a Pandas Excel writer
excel_writer = pd.ExcelWriter('output.xlsx', engine='xlsxwriter')

# Iterate through annotation files and write each as a separate sheet
for file in annotation_files:
    sheet_name = os.path.splitext(file)[0][:31]  # Truncate sheet name to 31 characters
    df = pd.read_csv(file, delimiter='\t')  # Change delimiter if needed
    df.to_excel(excel_writer, sheet_name=sheet_name, index=False)

# Save the Excel file
excel_writer.save()

print("Excel file with separate sheets created.")
