#!/usr/bin/env python3
"""
validate_sequences.py
─────────────────────
Batch Boltz2 + ipSAE binding validation for designed scFv sequences vs EGFR.

Pipeline per sequence:
  1. Parse binder FASTA  →  binder sequence
  2. Build Boltz2 YAML   →  chain B (binder) + chain A (EGFR target)
  3. boltz predict        →  .cif + pae .npz
  4. ipsae.py             →  binding metrics .txt
  5. Parse metrics        →  metric/<name>.csv

Final output: metric/summary.csv  (all sequences ranked by ipSAE, descending)

Usage:
  python validate_sequences.py

Parallelise across N GPUs (split the 300 sequences into equal slices):
  CUDA_VISIBLE_DEVICES=0 python validate_sequences.py --start   0 --end 100 &
  CUDA_VISIBLE_DEVICES=1 python validate_sequences.py --start 100 --end 200 &
  CUDA_VISIBLE_DEVICES=2 python validate_sequences.py --start 200 --end 300 &
  wait
  python validate_sequences.py --merge-only   # collects individual CSVs → summary
"""

import argparse
import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml
import torch

torch.set_float32_matmul_precision('medium')

# ─────────────────────────────────────── logging ─────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ──────────────────────────────────────── defaults ───────────────────────────

PAE_CUTOFF = 15.0  # passed to ipsae.py; Adaptyv competition default
DIST_CUTOFF = 15.0

# Key metrics written to every CSV (others appended after)
FRONT_COLS = [
    "sequence_name",
    "ipSAE",
    "ipSAE_d0chn",
    "ipSAE_d0dom",
    "ipTM_af",
    "ipTM_d0chn",
    "pDockQ",
    "pDockQ2",
    "LIS",
    "length",
]


# ═══════════════════════════════════════════ helpers ═════════════════════════


def parse_fasta(path: Path):
    """
    Return (header_id, sequence) from a FASTA file.
    header_id = first whitespace-delimited token after '>'.
    Only the first entry is used.
    """
    header, parts = "", []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not header:
                    header = line[1:].split()[0]
            else:
                parts.append(line.upper().replace("*", ""))  # strip stop codons
    if not parts:
        raise ValueError(f"No sequence in {path}")
    seq = "".join(parts)
    nonstandard = set(seq) - set("ACDEFGHIKLMNPQRSTVWY")
    if nonstandard:
        log.warning("Non-standard residues in %s: %s", path.name, nonstandard)
    return header, seq


def safe_name(stem: str) -> str:
    """
    Sanitise a filename stem to a valid Boltz job name.
    Boltz uses the YAML stem as the run name, so keep to [A-Za-z0-9_], ≤60 chars.
    """
    return re.sub(r"[^A-Za-z0-9_]", "_", stem)[:60]


# ════════════════════════════════════ Boltz helpers ══════════════════════════


def write_boltz_yaml(
    binder_seq: str,
    target_seq: str,
    out_path: Path,
    target_msa: Optional[Path] = None,
) -> None:
    """
    Emit a Boltz2 YAML for a two-chain complex.
      chain B  →  designed binder (no MSA)
      chain A  →  EGFR target     (optional pre-built MSA)

    Chain order matters for ipSAE parsing later (B before A in YAML so
    Boltz assigns the labels correctly; ipSAE reads from the CIF).
    """
    data = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "B",
                    "sequence": binder_seq,
                    "msa": "empty",
                }
            },
            {
                "protein": {
                    "id": "A",
                    "sequence": target_seq,
                    "msa": str(target_msa) if target_msa else "empty",
                }
            },
        ],
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        yaml.dump(data, fh, default_flow_style=False)


def run_boltz(yaml_path: Path, out_dir: Path) -> int:
    """
    Run `boltz predict` and return the exit code.
    Boltz derives the run name from the YAML filename stem, creating:
        {out_dir}/boltz_results_{stem}/predictions/{stem}/{stem}_model_0.cif
        {out_dir}/boltz_results_{stem}/predictions/{stem}/pae_{stem}_model_0.npz
    """
    cmd = [
        "boltz",
        "predict",
        str(yaml_path),
        "--out_dir",
        str(out_dir),
        "--write_full_pae",
    ]
    log.info("  boltz: %s", " ".join(cmd))
    return subprocess.run(cmd, text=True).returncode


def find_boltz_outputs(out_dir: Path, job: str):
    """
    Locate the PAE .npz and structure .cif from a completed Boltz run.
    Raises FileNotFoundError if either is missing.
    """
    pred = out_dir / f"boltz_results_{job}" / "predictions" / job
    cif = pred / f"{job}_model_0.cif"
    pae = pred / f"pae_{job}_model_0.npz"
    missing = [p for p in (cif, pae) if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Boltz output(s) not found: {missing}")
    return pae, cif


# ════════════════════════════════════ ipSAE helpers ══════════════════════════


def _ipsae_fmt(x: float) -> str:
    """Zero-pad single-digit cutoffs to match ipsae.py naming."""
    s = str(int(x))
    return ("0" + s) if x < 10 else s


def expected_ipsae_txt(cif: Path, pae_cutoff: float, dist_cutoff: float) -> Path:
    """Mirror the output path that ipsae.py generates from the CIF path."""
    return (
        cif.parent
        / f"{cif.stem}_{_ipsae_fmt(pae_cutoff)}_{_ipsae_fmt(dist_cutoff)}.txt"
    )


def run_ipsae(
    ipsae_script: Path,
    pae: Path,
    cif: Path,
    pae_cutoff: float,
    dist_cutoff: float,
) -> Path:
    """
    Call ipsae.py as a subprocess.
    Returns the path of the produced .txt file on success.
    """
    cmd = [
        sys.executable,
        str(ipsae_script),
        str(pae),
        str(cif),
        str(pae_cutoff),
        str(dist_cutoff),
    ]
    log.info(
        "  ipsae: %s %s %s %s %s",
        ipsae_script.name,
        pae.name,
        cif.name,
        pae_cutoff,
        dist_cutoff,
    )
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        # Print last 1000 chars of stderr for diagnosis
        log.error("ipsae.py stderr:\n%s", (proc.stderr or "")[-1000:])
        raise RuntimeError(f"ipsae.py exited {proc.returncode}")

    txt = expected_ipsae_txt(cif, pae_cutoff, dist_cutoff)
    if not txt.exists():
        raise FileNotFoundError(f"ipsae output missing: {txt}")
    return txt


def parse_ipsae_txt(txt: Path) -> dict:
    """
    Read ipsae's output CSV and return a flat dict of metrics.

    ipsae.py writes:
        <blank line>
        Chn1,Chn2,PAE,Dist,Type,ipSAE,...,Model
        A,B,15,15,asym,...          ← asym A→B
        B,A,15,15,asym,...          ← asym B→A
        A,B,15,15,max,...           ← symmetric max  ← we want this row
        <blank line>

    The "max" row is the symmetric best-of-both-directions estimate.
    Chain pair is always A (EGFR) – B (binder) in our setup.
    """
    # skip_blank_lines handles the leading blank line and intra-block separators
    df = pd.read_csv(txt, skip_blank_lines=True)
    df.columns = df.columns.str.strip()

    # Grab the first "max" row (there's only one chain pair in our 2-chain complex)
    max_rows = df[df["Type"].str.strip() == "max"]
    if max_rows.empty:
        raise ValueError(f"No 'max' row in {txt}")
    row = max_rows.iloc[0]

    metric_cols = [
        "ipSAE",
        "ipSAE_d0chn",
        "ipSAE_d0dom",
        "ipTM_af",
        "ipTM_d0chn",
        "pDockQ",
        "pDockQ2",
        "LIS",
        "n0res",
        "n0chn",
        "n0dom",
        "d0res",
        "d0chn",
        "d0dom",
        "nres1",
        "nres2",
        "dist1",
        "dist2",
    ]
    out = {
        "chain_pair": f"{str(row['Chn1']).strip()}-{str(row['Chn2']).strip()}",
    }
    for col in metric_cols:
        if col in row.index:
            out[col] = row[col]
    return out


# ═══════════════════════════════════ per-sequence pipeline ═══════════════════


def process_one(
    seq_path: Path,
    target_seq: str,
    work_dir: Path,
    metric_dir: Path,
    ipsae_script: Path,
    target_msa: Optional[Path],
    pae_cutoff: float,
    dist_cutoff: float,
    keep_work: bool,
    force: bool,
) -> Optional[dict]:
    """
    Full Boltz2 → ipSAE pipeline for one binder sequence.

    Returns a metrics dict on success, None on failure.
    The per-sequence CSV metric/<name>.csv is written before returning.
    """
    name = safe_name(seq_path.stem)
    metric_file = metric_dir / f"{name}.csv"

    # ── resume: skip already-completed sequences ──────────────────────────────
    if metric_file.exists() and not force:
        log.info("  ↩  cached – loading %s", metric_file.name)
        try:
            return pd.read_csv(metric_file).iloc[0].to_dict()
        except Exception:
            log.warning("  cached file unreadable, rerunning")

    # ── parse binder FASTA ───────────────────────────────────────────────────
    _, binder_seq = parse_fasta(seq_path)
    log.info("  binder: %d aa", len(binder_seq))

    # ── set up per-sequence scratch space ────────────────────────────────────
    seq_work = work_dir / name
    boltz_dir = seq_work / "boltz_out"
    yaml_path = seq_work / f"{name}.yaml"
    seq_work.mkdir(parents=True, exist_ok=True)
    boltz_dir.mkdir(exist_ok=True)

    try:
        # 1 ── write YAML ─────────────────────────────────────────────────────
        write_boltz_yaml(binder_seq, target_seq, yaml_path, target_msa)

        # 2 ── run Boltz2 ─────────────────────────────────────────────────────
        rc = run_boltz(yaml_path, boltz_dir)
        if rc != 0:
            raise RuntimeError(f"boltz predict returned {rc}")

        # 3 ── locate prediction files ────────────────────────────────────────
        pae, cif = find_boltz_outputs(boltz_dir, name)
        log.info("  cif: %s", cif.name)

        # 4 ── run ipSAE ──────────────────────────────────────────────────────
        txt = run_ipsae(ipsae_script, pae, cif, pae_cutoff, dist_cutoff)

        # 5 ── parse & save metrics ───────────────────────────────────────────
        metrics = parse_ipsae_txt(txt)
        metrics["sequence_name"] = name
        metrics["sequence"] = binder_seq
        metrics["length"] = len(binder_seq)

        metric_dir.mkdir(parents=True, exist_ok=True)
        row_df = pd.DataFrame([metrics])
        ordered = [c for c in FRONT_COLS if c in row_df.columns] + [
            c for c in row_df.columns if c not in FRONT_COLS
        ]
        row_df[ordered].to_csv(metric_file, index=False)

        ipsae_val = metrics.get("ipSAE", float("nan"))
        iptm_val = metrics.get("ipTM_af", float("nan"))
        pdockq_val = metrics.get("pDockQ", float("nan"))
        log.info(
            "  ✓  ipSAE=%.4f  ipTM=%.3f  pDockQ=%.4f", ipsae_val, iptm_val, pdockq_val
        )

        # 6 ── optional cleanup ───────────────────────────────────────────────
        if not keep_work:
            shutil.rmtree(seq_work, ignore_errors=True)

        return metrics

    except Exception as exc:
        log.error("  ✗  %s", exc)
        return None


# ═══════════════════════════════════════ summary helper ══════════════════════


def write_summary(metric_dir: Path) -> None:
    """
    Collect all per-sequence CSVs in metric_dir into a single ranked summary.
    Called automatically at the end of a run, and also by --merge-only.
    """
    csvs = sorted(metric_dir.glob("*.csv"))
    # exclude previous summary
    csvs = [c for c in csvs if c.stem not in ("summary",)]
    if not csvs:
        log.warning("No per-sequence CSVs found in %s", metric_dir)
        return

    frames = []
    for c in csvs:
        try:
            frames.append(pd.read_csv(c))
        except Exception as e:
            log.warning("Could not read %s: %s", c.name, e)

    if not frames:
        log.error("No valid CSVs to summarise")
        return

    summary = pd.concat(frames, ignore_index=True)
    if "ipSAE" in summary.columns:
        summary = summary.sort_values("ipSAE", ascending=False, na_position="last")

    ordered = [c for c in FRONT_COLS if c in summary.columns] + [
        c for c in summary.columns if c not in FRONT_COLS
    ]
    out_path = metric_dir / "summary.csv"
    summary[ordered].to_csv(out_path, index=False)
    log.info("Summary written: %s  (%d sequences)", out_path, len(summary))

    if "ipSAE" in summary.columns:
        log.info("Top-5 by ipSAE:")
        for _, r in summary.head(5).iterrows():
            log.info(
                "  %-42s  ipSAE=%.4f  ipTM=%.3f  pDockQ=%.4f",
                r.get("sequence_name", "?")[:42],
                r.get("ipSAE", float("nan")),
                r.get("ipTM_af", float("nan")),
                r.get("pDockQ", float("nan")),
            )


# ═══════════════════════════════════════════════ main ════════════════════════


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Batch Boltz2 + ipSAE binding validation against EGFR.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "--sequences-dir",
        type=Path,
        default=Path("sequences/binders"),
        metavar="DIR",
        help="Directory of binder FASTA files",
    )
    ap.add_argument(
        "--target",
        type=Path,
        default=Path("sequences/targets/egfr.fasta"),
        metavar="FASTA",
        help="EGFR target FASTA",
    )
    ap.add_argument(
        "--metric-dir",
        type=Path,
        default=Path("metrics"),
        metavar="DIR",
        help="Output directory for per-sequence metrics",
    )
    ap.add_argument(
        "--work-dir",
        type=Path,
        default=Path("work"),
        metavar="DIR",
        help="Scratch space for Boltz runs",
    )
    ap.add_argument(
        "--ipsae-script",
        type=Path,
        default=Path("ipsae.py"),
        metavar="PY",
        help="Path to ipsae.py",
    )
    ap.add_argument(
        "--target-msa",
        type=Path,
        default=Path("./egfr.a3m"),
        metavar="A3M",
        help="Optional EGFR MSA (.a3m) for Boltz2",
    )
    ap.add_argument("--pae-cutoff", type=float, default=PAE_CUTOFF)
    ap.add_argument("--dist-cutoff", type=float, default=DIST_CUTOFF)
    ap.add_argument(
        "--keep-work",
        action="store_true",
        help="Keep intermediate Boltz work dirs (useful for debugging)",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Rerun sequences that already have a metric file",
    )
    ap.add_argument(
        "--merge-only",
        action="store_true",
        help="Skip all predictions; just rebuild metric/summary.csv "
        "from existing per-sequence CSVs",
    )
    ap.add_argument(
        "--start",
        type=int,
        default=0,
        metavar="N",
        help="Index of first sequence to process (0-based)",
    )
    ap.add_argument(
        "--end",
        type=int,
        default=None,
        metavar="N",
        help="Index of last sequence (exclusive); default = all",
    )
    return ap


def main() -> int:
    args = build_parser().parse_args()

    # ── merge-only mode ───────────────────────────────────────────────────────
    if args.merge_only:
        write_summary(args.metric_dir)
        return 0

    # ── validate paths ────────────────────────────────────────────────────────
    errors = []
    for p, label in [
        (args.sequences_dir, "--sequences-dir"),
        (args.target, "--target"),
        (args.ipsae_script, "--ipsae-script"),
    ]:
        if not p.exists():
            errors.append(f"{label} not found: {p}")
    if args.target_msa and not args.target_msa.exists():
        errors.append(f"--target-msa not found: {args.target_msa}")
    if errors:
        for e in errors:
            log.error(e)
        return 1

    # ── load EGFR target ──────────────────────────────────────────────────────
    target_id, target_seq = parse_fasta(args.target)
    log.info("Target: %s  (%d aa)", target_id, len(target_seq))
    if args.target_msa:
        log.info("Target MSA: %s", args.target_msa)

    # ── collect binder FASTAs ─────────────────────────────────────────────────
    seq_paths = sorted(
        p
        for p in args.sequences_dir.iterdir()
        if p.suffix.lower() in {".fasta", ".fa", ".faa"} and p.is_file()
    )
    if not seq_paths:
        log.error("No FASTA files found in %s", args.sequences_dir)
        return 1

    end = args.end if args.end is not None else len(seq_paths)
    batch = seq_paths[args.start : end]
    log.info(
        "Binders: %d total, processing indices %d–%d (%d sequences)",
        len(seq_paths),
        args.start,
        end - 1,
        len(batch),
    )

    args.metric_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)

    # ── main loop ─────────────────────────────────────────────────────────────
    results, failed = [], []

    for i, sp in enumerate(batch, 1):
        log.info("─" * 66)
        log.info("[%d/%d]  %s", i, len(batch), sp.name)

        m = process_one(
            seq_path=sp,
            target_seq=target_seq,
            work_dir=args.work_dir,
            metric_dir=args.metric_dir,
            ipsae_script=args.ipsae_script,
            target_msa=args.target_msa,
            pae_cutoff=args.pae_cutoff,
            dist_cutoff=args.dist_cutoff,
            keep_work=args.keep_work,
            force=args.force,
        )

        if m is not None:
            results.append(m)
        else:
            failed.append(sp.name)

    # ── write summary ─────────────────────────────────────────────────────────
    log.info("═" * 66)
    log.info("Complete.  Success: %d  |  Failed: %d", len(results), len(failed))

    write_summary(args.metric_dir)

    if failed:
        fail_log = args.metric_dir / "failed.txt"
        fail_log.write_text("\n".join(failed) + "\n")
        log.warning("%d sequence(s) failed — see %s", len(failed), fail_log)

    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
