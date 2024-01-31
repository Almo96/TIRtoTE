import os

def remove_output_file(output_file):
    os.remove(output_file)


def reverse_complement(input_fasta, output_path):
    # Convert the provided output_path to an absolute path
    output_path = os.path.abspath(output_path)

    # Extract filename from the path and remove extension
    file_name = os.path.basename(input_fasta)
    file_name_no_extension = os.path.splitext(file_name)[0]

    # Define the output FASTA file path
    output_fasta = os.path.join(output_path, file_name_no_extension + "_revcomp.fa")

    # Dictionary for complementing bases
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    reverse_complement_sequences = []

    # Read the input FASTA file and generate reverse complement sequences
    with open(input_fasta, 'r') as infile:
        sequence_id = ''
        sequence_lines = []

        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id:
                    sequence = ''.join(sequence_lines)
                    reverse_complement_sequence = ''.join(complement[base] for base in reversed(sequence))
                    reverse_complement_sequences.append((sequence_id, reverse_complement_sequence))
                    sequence_lines = []
                sequence_id = line
            else:
                sequence_lines.append(line)

        if sequence_id:
            sequence = ''.join(sequence_lines)
            reverse_complement_sequence = ''.join(complement[base] for base in reversed(sequence))
            reverse_complement_sequences.append((sequence_id, reverse_complement_sequence))

    # Write the reverse complement sequences to the output FASTA file
    with open(output_fasta, 'w') as outfile:
        for seq_id, seq in reverse_complement_sequences:
            outfile.write(f"{seq_id}_recomp\n")
            outfile.write(f"{seq}\n")

    return output_fasta