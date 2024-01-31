import os

def controlIO(TIR, genome, output, output_path, minlenTE, maxlenTE):
    # Check if the extensions are .fasta or .fa
    valid_extensions = ['.fasta', '.fa', '.fna']
    for arg_name, arg_value in [("TIR", TIR), ("genome", genome)]:
        if not arg_value.endswith(tuple(valid_extensions)):
            raise ValueError(f"{arg_name} must have either a .fasta, .fna or .fa extension.")


    print("TIR:", os.path.basename(TIR))
    print("Genome:", os.path.basename(genome))
    print("Output path:", output_path)
    print("Output:", output)
    print("Minimum length:", minlenTE)
    print("Maximum length:", maxlenTE)

