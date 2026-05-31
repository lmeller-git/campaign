import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def generate_validation_plots(
    csv_path: str, output_img_path: str = "design_stats_violin.png"
):
    # 1. Load the data
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: Could not find CSV file at {csv_path}")
        return

    # 2. Select the key structural/binding metrics
    metrics = ["ipSAE", "ipTM_af", "pDockQ2"]

    # Verify columns exist in the CSV
    missing = [m for m in metrics if m not in df.columns]
    if missing:
        print(f"Warning: Missing columns {missing} in CSV. Adjusting selection...")
        metrics = [m for m in metrics if m in df.columns]

    if not metrics:
        print("No valid metrics found to plot.")
        return

    # 3. Initialize the matplotlib figure (Side-by-side layout)
    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(
        1, len(metrics), figsize=(6 * len(metrics), 7), sharey=False
    )

    # Handle single-metric case gracefully if axes isn't an array
    if len(metrics) == 1:
        axes = [axes]

    # Custom color palette for the plots
    colors = ["#4C72B0", "#55A868", "#C44E52"]

    # 4. Populate each subplot with a customized violin plot
    for idx, metric in enumerate(metrics):
        ax = axes[idx]

        # Draw the violin plot
        sns.violinplot(
            data=df,
            y=metric,
            ax=ax,
            color=colors[idx % len(colors)],
            inner="box",  # Draws a mini boxplot inside to clearly show median/IQR
            linewidth=2,
            cut=0,  # Limit the violin to actual data range (no wild extrapolations)
        )

        # Overlay an explicit median point marker for extreme clarity
        median_val = df[metric].median()
        ax.axhline(median_val, color="black", linestyle="--", alpha=0.4, linewidth=1.5)

        # Titles and Formatting
        ax.set_title(f"Distribution of {metric}", fontsize=16, pad=15, weight="bold")
        ax.set_ylabel(metric, fontsize=14, weight="bold")
        ax.set_xlabel("")  # Clear out the x-axis label since it's a single distribution

        # Annotate median text onto the plot area
        ax.text(
            0.05,
            0.95,
            f"Median: {median_val:.3f}",
            transform=ax.transAxes,
            fontsize=12,
            verticalalignment="top",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", alpha=0.7, edgecolor="gray"
            ),
        )

    # 5. Final Polish & Save
    plt.tight_layout()

    # Crucial Fix: Push the subplots down slightly to create a dedicated roof for the title
    plt.subplots_adjust(top=0.85)

    # Keep y safely under 1.0 so it never clips the top edge of the window
    plt.suptitle(
        f"Structural Validation Profile (n={len(df)} Designs)",
        fontsize=20,
        weight="bold",
        y=0.95,
    )

    # bbox_inches="tight" ensures no trailing elements are clipped in the saved PNG
    plt.savefig(output_img_path, dpi=300, bbox_inches="tight")
    print(f"Success! Validation violin plots saved to: {output_img_path}")
    plt.show()


if __name__ == "__main__":
    # Change 'metric/summary.csv' to the actual path of your output file
    generate_validation_plots("metrics/summary.csv")
