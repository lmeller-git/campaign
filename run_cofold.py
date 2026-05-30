import os
import argparse
import subprocess
import yaml

def run_boltz2_cofold(target_fasta_path, binder_fasta_path, output_dir, samples, device):
    """
    Runs a Boltz2 cofold prediction using the recommended YAML input format.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Read sequence from FASTA files
    with open(target_fasta_path, 'r') as f:
        target_sequence = "".join([line for line in f.read().splitlines() if not line.startswith('>')])
    with open(binder_fasta_path, 'r') as f:
        binder_sequence = "".join([line for line in f.read().splitlines() if not line.startswith('>')])

    # Create the YAML configuration for Boltz
    boltz_config = {
        'version': 1,
        'sequences': [
            {'protein': {
                'id': 'A',
                'sequence': target_sequence,
                'msa': 'empty' # Can use server or provide MSA if needed
            }},
            {'protein': {
                'id': 'B',
                'sequence': binder_sequence,
                'msa': 'empty' # Force single-sequence mode for binder screening
            }}
        ]
    }
    
    yaml_path = os.path.join(output_dir, 'input.yaml')
    with open(yaml_path, 'w') as f:
        yaml.dump(boltz_config, f)

    print(f"Generated Boltz2 config at {yaml_path}")

    # Construct the boltz command
    command = [
        "boltz", "predict", yaml_path,
        "--out_dir", output_dir,
        "--diffusion_samples", str(samples),
        "--accelerator", "gpu",
        "--devices", str(device),
        "--override" # Override existing results if any
    ]
    
    # The documentation suggests using this flag for MSA generation if not providing one.
    # For binder screening, starting with single-sequence ('empty' msa) is often sufficient.
    # If you have an MSA server, you could add '--use_msa_server' here.
    
    print(f"Running command: {' '.join(command)}")

    # Execute the command
    try:
        subprocess.run(command, check=True, text=True, capture_output=True)
        print("Boltz2 run completed successfully.")
    except subprocess.CalledProcessError as e:
        print("--- Boltz2 Error ---")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        print("--------------------")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a Boltz2 cofold prediction.")
    parser.add_argument("--target_fasta_path", required=True, help="Path to the target FASTA file.")
    parser.add_argument("--binder_fasta_path", required=True, help="Path to the binder FASTA file.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the prediction outputs.")
    parser.add_argument("--samples", type=int, default=3, help="Number of samples (diffusion_samples) to generate.")
    # The --device argument is now optional. On a cluster, the scheduler usually handles GPU allocation.
    parser.add_argument("--device", type=int, default=0, help="GPU device ID to use for the prediction.")
    
    args = parser.parse_args()

    run_boltz2_cofold(
        args.target_fasta_path,
        args.binder_fasta_path,
        args.output_dir,
        args.samples,
        args.device
    )
