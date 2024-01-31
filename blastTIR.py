import subprocess
import os
import utils

def blast(TIR, TIR_2, genome, output, path):
    output_blast1 = os.path.join(path, output + "_blast1" +".txt")
    output_blast2 = os.path.join(path, output + "_blast2" +".txt")
    output_blast_sum = os.path.join(path, output + "_blastsum" +".txt")

    # Define the BLAST command
    blast_command1 = [
        "blastn",
        "-query", TIR,
        "-subject", genome,
        "-outfmt", "6 qseqid sseqid sstart send evalue sstrand",
        "-strand", "plus",
        "-out", output_blast1,
    ]

    # Run the BLAST command
    subprocess.run(blast_command1)

    # Run with the second TIR
    blast_command2 = [
        "blastn",
        "-query", TIR_2,
        "-subject", genome,
        "-outfmt", "6 qseqid sseqid sstart send evalue sstrand",
        "-strand", "plus",
        "-out", output_blast2,
    ]

    # Run the BLAST command
    subprocess.run(blast_command2)


# Concatenate the output files
    with open(output_blast1, 'r') as file1, open(output_blast2, 'r') as file2, open(output_blast_sum, 'w') as file_sum:
        for line in file1:
            file_sum.write(line)
        for line in file2:
            file_sum.write(line)


    print(f"BLAST results saved to {output_blast_sum}")
    if os.path.getsize(output_blast1) != 0 and os.path.getsize(output_blast2) != 0:
        TE_presence = True
    else:
        TE_presence = False

    utils.remove_output_file(output_blast1)
    utils.remove_output_file(output_blast2)

    return output_blast_sum, TE_presence