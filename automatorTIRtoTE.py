import os
import subprocess
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script with mandatory and optional fields.")
    
    # Add the mandatory positional arguments
    parser.add_argument("TIR", help="TIR file (.fasta or .fa)")
    parser.add_argument("input_path", help="Folder containing the input files (.fasta or .fa)")
    parser.add_argument("output_path", help="Output path (.fasta or .fa)")
    parser.add_argument("minlenTE", help="Minimum lenght of the TE (integer)")
    parser.add_argument("maxlenTE", help="Maximum lenght of the TE (integer)")

    # Add the optional argument with a default value of 20
    args = parser.parse_args()

    TIR = args.TIR
    input_path = args.input_path
    output_path = args.output_path
    minlenTE = int(args.minlenTE)
    maxlenTE = int(args.maxlenTE)


input_files = os.listdir(input_path)

# Iterate through the input files
for input_file in input_files:
    # Check if the file has a .fa extension
    if input_file.endswith('.fa') or input_file.endswith('.fna') or input_file.endswith('.fasta'):
        # Remove the .fa extension to get the output file name
        output_name = os.path.splitext(input_file)[0]
        print("Output name:", output_name)

        # Define the full paths for input files and the script
        input_file_path = os.path.join(input_path, input_file)
        script_path = os.path.join(os.path.dirname(__file__), "TIRtoTE.py")

        # Run the TIRtoTE.py script with the specified parameters
        cmd = f"python {script_path} {TIR} {input_file_path} {output_path} {output_name} {minlenTE} {maxlenTE}"
        
        # Execute the command
        subprocess.call(cmd, shell=True)