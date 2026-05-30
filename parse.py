import os


def process_sequences(input_filename, output_directory="sequences/binders"):
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    try:
        with open(input_filename, "r") as file:
            for line in file:
                line = line.strip()

                # Skip empty lines or lines that don't start with '>'
                if not line or not line.startswith(">"):
                    continue

                try:
                    # Split the line at the separator
                    header_part, sequence_part = line.split("| all aa:")

                    # Extract and clean the sequence name (remove '>' and whitespace)
                    # Example: ">design_1_edits_5  " becomes "design_1_edits_5"
                    seq_name = header_part.strip()[1:]

                    # Extract and clean the sequence
                    sequence = sequence_part.strip()

                    # Define the output file path
                    output_path = os.path.join(output_directory, f"{seq_name}.fasta")

                    # Write to the individual FASTA file
                    with open(output_path, "w") as fasta_file:
                        fasta_file.write(f">{seq_name}\n{sequence}\n")

                    print(f"Successfully generated: {output_path}")

                except ValueError:
                    print(f"Skipping malformed line: {line}")

    except FileNotFoundError:
        print(f"Error: The file '{input_filename}' was not found.")


import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, default="seq..txt")
    args = parser.parse_args()
    process_sequences(args.path)
