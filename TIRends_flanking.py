# From each insertion 8 regions:
# R1 = Genomic DNA flanking TIR-1
# R2 = TIR-1-left 
# R3 = TIR-1-right
# R4 = TE coding region flanking TIR-1
# R5 = TE coding region flanking TIR-2
# R6 = TIR-2-left 
# R7 = TIR-2-right
# R8 = Genomic DNA flanking TIR-2

import os
import pandas as pd
from Bio import SeqIO

def TIRe_fDNA(output, path, output_blast, minlenTE, maxlenTE):

    # Part 1: create a .bed file with the TIR positions of each TE
    df = pd.read_csv(output_blast, sep='\t', header=None)

    # Sort the DataFrame by Chromosome and TIR starting position
    df_sorted = df.sort_values(by=[1, 2], ascending=[True, True])
    df_TIR = pd.DataFrame(columns=["Chr", "Start", "End", "Score", "Strand"])

    # Iterate through df_sorted and check for the right distance between two TIRs
    for i in range(len(df_sorted) - 1):
        if df_sorted.iloc[i + 1, 1] == df_sorted.iloc[i, 1]:        # Check if they are on the same Chromosome
            if df_sorted.iloc[i + 1, 3] - df_sorted.iloc[i, 2] < maxlenTE:
                if df_sorted.iloc[i + 1, 3] - df_sorted.iloc[i, 2] > minlenTE:        # Check if they are within the specified range
                    new_row1 = df_sorted.iloc[i, [1, 2, 3, 4, 5]].tolist()
                    new_row2 = df_sorted.iloc[i + 1, [1, 2, 3, 4, 5]].tolist()
                    df_TIR.loc[len(df_TIR)] = new_row1  # Add TIR left
                    df_TIR.loc[len(df_TIR)] = new_row2  # Add TIR right
    
    df_TIR["Strand"] = df_TIR["Strand"].replace("plus", "+")
    df_TIR["Strand"] = df_TIR["Strand"].replace("minus", "-")

    #Part 2: use the start and end position of each TIR pair to extract the sequence of the 8 regions of interest
    x = 10 # Genomic region
    y = 10 # TIR region
    c = 10 # Coding region

    df_regions = pd.DataFrame(columns=["Chr", "Start", "End", "Region", "Score", "Strand"])

    for index, row in df_TIR.iterrows():

        if index % 2 == 0:
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - x - 1, row['Start'] - 1, 'R1', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
        
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - 1, row['Start'] + y - 1, 'R2', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
            
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - y, row['End'], 'R3', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
            
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'], row['End'] + c, 'R4', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])

        else:
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - c - 1, row['Start'] - 1, 'R5', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
        
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - 1, row['Start'] + y - 1, 'R6', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
            
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - y, row['End'], 'R7', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
            
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'], row['End'] + x, 'R8', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])


    output_df_regions = os.path.join(path, output + "_regions" + ".bed")
    df_regions.to_csv(output_df_regions, sep='\t', header=False, index=False)
    print(f"Boundary regions of the TIR saved to {output_df_regions}")

    return output_df_regions


def TIRe_fDNA_out(regions4_fasta, output, path):

    sequences = list(SeqIO.parse(regions4_fasta, "fasta"))
    tr = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N', 'a':'t', 't': 'a', 'c': 'g', 'g': 'c', 'n':'n'}


    if len(sequences) % 8 == 0:
        str_regions = []
        string_seq = ""
        counter_seq = 1
        neg_counter = 1
        for sequence in sequences:
            if sequence.id.endswith("(+)"):
                if counter_seq == 1:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "||"
                elif counter_seq == 2:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "..."
                elif counter_seq == 3:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "|"
                elif counter_seq == 4:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "......"
                elif counter_seq == 5:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "|"
                elif counter_seq == 6:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "..."
                elif counter_seq == 7:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "||"
                elif counter_seq == 8:
                    string_seq = string_seq + (str(sequence.seq).upper())
                    str_regions.append(string_seq)
                    counter_seq = 0
                    string_seq = ""
                counter_seq = counter_seq + 1
            else:
                if neg_counter == 1:
                    result = [tr[base] for base in reversed(sequence.seq)]  #bedtools provided the reverse complement but on the genomic regions
                    R1 = ''.join(result)                                    # it is not needed, so by performing it again I get the original sequence
                    R1 = R1.upper()
                elif neg_counter == 2:
                    R7 = str(sequence.seq).upper()
                elif neg_counter == 3:
                    R6 = str(sequence.seq).upper()
                elif neg_counter == 4:
                    R5 = str(sequence.seq).upper()
                elif neg_counter == 5:
                    R4 = str(sequence.seq).upper()
                elif neg_counter == 6:
                    R3 = str(sequence.seq).upper()
                elif neg_counter == 7:
                    R2 = str(sequence.seq).upper()
                elif neg_counter == 8:
                    result2 = [tr[base] for base in reversed(sequence.seq)]
                    R8 = ''.join(result2)
                    R8 = R8.upper()                    
                    R = R1 + "||" + R2 + "..." + R3 + "|" + R4 + "......" + R5 + "|" + R6 + "..." + R7 + "||" + R8
                    str_regions.append(R)
                    neg_counter = 0
                    string_seq = ""
                neg_counter = neg_counter + 1

        output_str_regions = os.path.join(path, output + "_regions" + ".txt")

        with open(output_str_regions, "w") as file:
            for item in str_regions:
                file.write(f"{item}\n")

        print(f"The genomic flanking, the TIR and the coding regions were saved to {output_str_regions}")

    else:
        print("The number of sequences in the regions.bed file should be a multiple of 4")

    return output_str_regions

# This function works, but if there are deletion inside the TIR the nucleotides are shifted
# A better approach is to take the Start and End of each TIR instead of the start and end of the TE
# def TIRe_fDNAv1(output_bed, TIR, output):
#     x = 10 # Genomic region
#     y = 10 # TIR region
#     c = 10 # Coding region

#     sequence = list(SeqIO.parse(TIR, "fasta"))
#     first_sequence = sequence[0]
#     TIR_length = len(first_sequence)

#     column_names = ["Chr", "Start", "End", "Label", "Score", "Strand"]
#     df = pd.read_csv(output_bed, sep='\t', header=None, names=column_names)
#     df_regions = pd.DataFrame(columns=["Chr", "Start", "End", "Region", "Score", "Strand"])
    
#     for _, row in df.iterrows():

#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - x - 1, row['Start'] - 1, 'R1']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'], row['Start'] + y, 'R2']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] + TIR_length - y, row['Start'] + TIR_length, 'R3']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] + TIR_length + 1, row['Start'] + TIR_length + c + 1, 'R4']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - TIR_length - c - 1, row['End'] - TIR_length - 1, 'R5']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - TIR_length , row['End'] - TIR_length + y, 'R6']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - y, row['End'], 'R7']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] + 1, row['End'] + x + 1, 'R8']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
    
#     path = os.path.dirname(TIR)
#     output_df_regions = os.path.join(path, output + "_regions" + ".bed")
#     df_regions.to_csv(output_df_regions, sep='\t', header=False, index=False)
#     print(f"Boundary regions of the TIR saved to {output_df_regions}")

#     return output_df_regions