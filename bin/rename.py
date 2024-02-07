import os
import sys
import gzip

def rename_fasta_headers(file_path, output_path):
    # Get the file name without extension
    file_name = os.path.splitext(os.path.basename(file_path))[0].replace('.fasta','')
    print(file_name)

    # Read the contents of the file
    with gzip.open(file_path, 'rt') as file:
        lines = file.readlines()

    # Rename the headers
    new_lines = []
    for line in lines:
        if line.startswith('>'):
            new_lines.append(f'>{file_name}_{line[1:]}')
        else:
            new_lines.append(line)

    # Write the modified contents back to the file
    with gzip.open(output_path, 'wt') as file:
        file.writelines(new_lines)
    return

# Usage example
file_path = sys.argv[1]
output_path = sys.argv[2]
files=os.listdir(file_path)
for file in files:
    fasta_file_path = os.path.join(file_path, file)
    rename_fasta_headers(fasta_file_path, f'{output_path}/{file}')
