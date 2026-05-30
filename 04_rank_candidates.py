import pandas as pd
import os

# --- Configuration ---
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '.'))
METRICS_DIR = os.path.join(BASE_DIR, "metrics")
RANKED_DIR = os.path.join(BASE_DIR, "ranked")
INPUT_CSV = os.path.join(METRICS_DIR, "all_metrics.csv")
OUTPUT_CSV = os.path.join(RANKED_DIR, "ranked_candidates.csv")
# ---------------------

def normalize(series):
    return (series - series.min()) / (series.max() - series.min())

def main():
    os.makedirs(RANKED_DIR, exist_ok=True)

    try:
        df = pd.read_csv(INPUT_CSV)
    except FileNotFoundError:
        print(f"Error: Metrics file not found at {INPUT_CSV}")
        print("Please run 03_extract_metrics.py first.")
        return

    main_preds = df[~df['is_decoy']].copy()
    decoy_preds = df[df['is_decoy']].copy()

    # Calculate specificity score
    if not decoy_preds.empty:
        max_decoy_ipsae = decoy_preds.groupby('candidate')['ipSAE_min'].max().reset_index()
        max_decoy_ipsae.rename(columns={'ipSAE_min': 'max_decoy_ipSAE'}, inplace=True)
        main_preds = pd.merge(main_preds, max_decoy_ipsae, on='candidate', how='left')
        main_preds['max_decoy_ipSAE'].fillna(0, inplace=True)
        main_preds['specificity_score'] = main_preds['ipSAE_min'] - main_preds['max_decoy_ipSAE']
    else:
        main_preds['specificity_score'] = main_preds['ipSAE_min'] # No decoys, so specificity is just ipSAE

    # Normalize metrics for scoring
    main_preds['ipSAE_norm'] = normalize(main_preds['ipSAE_min'])
    main_preds['pDockQ_norm'] = normalize(main_preds['pDockQ'])
    # The original script assumed 'interface_contacts' which is not in the boltz output
    # We will proceed without it in the score for now.
    main_preds['contacts_norm'] = 0
    if 'interface_contacts' in main_preds.columns and not main_preds['interface_contacts'].isnull().all():
        main_preds['contacts_norm'] = normalize(main_preds['interface_contacts'])

    # Calculate weighted score
    main_preds['score'] = (
        0.6 * main_preds['ipSAE_norm'] +
        0.2 * main_preds['pDockQ_norm'] +
        0.2 * main_preds['contacts_norm']
    )

    # Rank
    ranked_df = main_preds.sort_values(by='score', ascending=False)

    ranked_df.to_csv(OUTPUT_CSV, index=False)
    print(f"Ranking complete. Results saved to {OUTPUT_CSV}")
    print("\nTop 5 candidates:")
    print(ranked_df[['candidate', 'score', 'ipSAE_min', 'specificity_score', 'pDockQ', 'interface_contacts']].head())

if __name__ == "__main__":
    main()
