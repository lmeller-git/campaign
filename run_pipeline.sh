#!/bin/bash
PYTHON_EXEC="uv run"
# This script automates the setup and execution of the Boltz screening pipeline.
# --- 4. Run the pipeline steps ---
echo "-------------------------------------"
echo "STEP 1: Running cofolding predictions against target..."
uv run 01_run_predictions.py
echo "STEP 1: Done."
echo "-------------------------------------"

echo "STEP 2: Running negative controls against decoys..."
uv run 02_run_negative_controls.py
echo "STEP 2: Done."
echo "-------------------------------------"

echo "STEP 3: Extracting metrics from all runs..."
uv run 03_extract_metrics.py
echo "STEP 3: Done."
echo "-------------------------------------"

echo "STEP 4: Ranking candidates..."
uv run 04_rank_candidates.py
echo "STEP 4: Done."
echo "-------------------------------------"

echo "🎉 Pipeline finished successfully!"
echo "Your final ranked candidates are in: campaign/ranked/ranked_candidates.csv"
