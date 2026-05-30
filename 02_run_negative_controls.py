import os
import glob
import subprocess

# --- Configuration ---
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '.'))
SEQUENCES_DIR = os.path.join(BASE_DIR, "sequences/binders")
DECOYS_DIR = os.path.join(BASE_DIR, "decoys")
PREDICTIONS_DIR = os.path.join(BASE_DIR, "predictions")
RUN_COFOLD_SCRIPT = os.path.join(os.path.dirname(__file__), "run_cofold.py")
SAMPLES_PER_CANDIDATE = 1 # Usually fewer samples for decoys is fine
# ---------------------

def main():
    scfv_fastas = glob.glob(os.path.join(SEQUENCES_DIR, "*.fasta"))
    decoy_fastas = glob.glob(os.path.join(DECOYS_DIR, "*.fasta"))

    if not scfv_fastas or not decoy_fastas:
        print(f"FASTA files not found in {SEQUENCES_DIR} or {DECOYS_DIR}")
        return

    print(f"Found {len(scfv_fastas)} scFv candidates and {len(decoy_fastas)} decoys.")

    for scfv_fasta in scfv_fastas:
        candidate_name = os.path.splitext(os.path.basename(scfv_fasta))[0]

        for decoy_fasta in decoy_fastas:
            decoy_name = os.path.splitext(os.path.basename(decoy_fasta))[0]
            output_dir = os.path.join(PREDICTIONS_DIR, f"{candidate_name}_vs_{decoy_name}")

            print(f"\n--- Running negative control for {candidate_name} vs {decoy_name} ---")

            command = [
                "python3", RUN_COFOLD_SCRIPT,
                "--target_fasta_path", decoy_fasta,
                "--binder_fasta_path", scfv_fasta,
                "--output_dir", output_dir,
                "--samples", str(SAMPLES_PER_CANDIDATE)
            ]

            subprocess.run(command, check=True)

if __name__ == "__main__":
    main()
