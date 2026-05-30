import os
import glob
import subprocess

# --- Configuration ---
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '.'))
SEQUENCES_DIR = os.path.join(BASE_DIR, "sequences/binders")
PREDICTIONS_DIR = os.path.join(BASE_DIR, "predictions")
EGFR_FASTA = os.path.join(BASE_DIR, "sequences/targets/egfr.fasta")
RUN_COFOLD_SCRIPT = os.path.join(os.path.dirname(__file__), "run_cofold.py")
SAMPLES_PER_CANDIDATE = 3
# ---------------------

def main():
    os.makedirs(PREDICTIONS_DIR, exist_ok=True)

    scfv_fastas = glob.glob(os.path.join(SEQUENCES_DIR, "*.fasta"))
    if not scfv_fastas:
        print(f"No FASTA files found in {SEQUENCES_DIR}")
        return

    print(f"Found {len(scfv_fastas)} scFv candidates.")

    for scfv_fasta in scfv_fastas:
        candidate_name = os.path.splitext(os.path.basename(scfv_fasta))[0]
        output_dir = os.path.join(PREDICTIONS_DIR, candidate_name)

        print(f"\n--- Running prediction for {candidate_name} ---")

        command = [
            "python3", RUN_COFOLD_SCRIPT,
            "--target_fasta_path", EGFR_FASTA,
            "--binder_fasta_path", scfv_fasta,
            "--output_dir", output_dir,
            "--samples", str(SAMPLES_PER_CANDIDATE)
        ]

        subprocess.run(command, check=True)
        break

if __name__ == "__main__":
    main()
