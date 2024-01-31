# Author Almorò Scarpa
# Program to find DNA transposons in a genome using the TIRs

import argparse
import os
import InputOutput
import blastTIR
import TEfromTIR
import TSD
import TIRends_flanking
import utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script with mandatory and optional fields.")
    
    # Add the mandatory positional arguments
    parser.add_argument("TIR", help="TIR file (.fasta or .fa)")
    parser.add_argument("genome", help="Genome file (.fasta or .fa)")
    parser.add_argument("output_path", help="Output path (.fasta or .fa)")
    parser.add_argument("output", help="Output file (.fasta or .fa)")
    parser.add_argument("minlenTE", help="Minimum lenght of the TE (integer)")
    parser.add_argument("maxlenTE", help="Maximum lenght of the TE (integer)")

    # Add the optional argument with a default value of 20
    parser.add_argument("--tsd", help="Optional TSD field (integer)", type=int)
    args = parser.parse_args()

    TIR = args.TIR
    genome = args.genome
    output = args.output
    output_path = args.output_path
    minlenTE = int(args.minlenTE)
    maxlenTE = int(args.maxlenTE)
    tsd = args.tsd

    InputOutput.controlIO(TIR, genome, output, output_path, minlenTE, maxlenTE)

    TIR_2 = utils.reverse_complement(TIR, output_path)

    output_blast, TE_presence = blastTIR.blast(TIR, TIR_2, genome, output, output_path)
    
    if TE_presence == True:

        output_bed, bed, blast2 = TEfromTIR.TE_position(output, output_blast, minlenTE, maxlenTE, output_path)

        output_fasta = TEfromTIR.getfasta(genome, output_bed, output, output_path, "")

        regions4_bed = TIRends_flanking.TIRe_fDNA(output, output_path, blast2, minlenTE, maxlenTE)

        regions4_fasta = TEfromTIR.getfasta(genome, regions4_bed, output, output_path, "_regions")

        output_regions = TIRends_flanking.TIRe_fDNA_out(regions4_fasta, output, output_path)


        if tsd is not None:

            TSD.checkTSD(output_regions, output, output_path)

    else:
        utils.remove_output_file(output_blast)

        utils.remove_output_file(TIR_2)
        
        print("No TEs have been found")
