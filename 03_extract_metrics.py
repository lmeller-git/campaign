import os
import glob
import json
import pandas as pd
from collections import defaultdict

# --- Configuration ---
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '.'))
PREDICTIONS_DIR = os.path.join(BASE_DIR, "predictions")
METRICS_DIR = os.path.join(BASE_DIR, "metrics")
OUTPUT_CSV = os.path.join(METRICS_DIR, "all_metrics.csv")
# ---------------------

def extract_best_metrics(prediction_dir):
    """Extracts metrics from all samples and returns the one with the best ipSAE_min."""
    # Note: Boltz2 often stores results under prediction_dir/boltz_results_*/...
    # Adjust the glob pattern if your confidence files are nested differently.
    metrics_files = glob.glob(os.path.join(prediction_dir, "predictions", "*", "confidence*.json"))

    if not metrics_files:
        return None

    best_metrics = None
    best_ipsae = -1

    for mf in metrics_files:
        with open(mf, 'r') as f:
            try:
                metrics = json.load(f)
                # Assuming protein_iptm is the relevant metric as per boltz output
                current_ipsae = metrics.get('protein_iptm', -1)

                if current_ipsae > best_ipsae:
                    best_ipsae = current_ipsae
                    best_metrics = {
                        "ipSAE_min": metrics.get("protein_iptm"),
                        "pDockQ": metrics.get("iptm"), # Mapping iptm to pDockQ as a proxy
                        "interface_contacts": None, # This metric is not available in the json
                        "interface_pLDDT": metrics.get("complex_iplddt"),
                    }
            except json.JSONDecodeError:
                print(f"Warning: Could not parse {mf}")
                continue

    if best_metrics:
        # Use the passed prediction_dir instead of the loop variable `mf`
        best_metrics['prediction_dir'] = os.path.basename(prediction_dir)

    return best_metrics

def main():
    os.makedirs(METRICS_DIR, exist_ok=True)

    prediction_dirs = glob.glob(os.path.join(PREDICTIONS_DIR, "*"))

    all_metrics = []
    for pred_dir in prediction_dirs:
        if os.path.isdir(pred_dir):
            metrics = extract_best_metrics(pred_dir)
            if metrics:
                all_metrics.append(metrics)

    if not all_metrics:
        print("No metrics found. Please check your file paths and JSON structure.")
        return

    df = pd.DataFrame(all_metrics)

    # Separate main predictions from decoy predictions
    df['candidate'] = df['prediction_dir'].apply(lambda x: x.split('_vs_')[0])
    df['is_decoy'] = df['prediction_dir'].str.contains('_vs_')

    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Metrics extracted and saved to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
