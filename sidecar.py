"""JSON-RPC sidecar for Brimer-PLAST Electron app.

Reads one JSON request per line from stdin, writes one JSON response per line
to stdout. The Electron main process spawns this as a child process and
communicates via stdio.

Request format (one JSON line):
    {"id": 1, "command": "run_pipeline", "params": {...}}

Response formats (one JSON line each):
    {"id": 1, "status": "progress", "message": "...", "pct": 50}
    {"id": 1, "status": "ok", "result": {...}}
    {"id": 1, "status": "error", "message": "..."}
"""

from __future__ import annotations

import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from dataclasses import asdict

import brimer_plast  # ensure package is importable (no-op)
from brimer_plast.pipeline import run_pipeline
from brimer_plast.pdf_report import build_pdf_report


def send(obj: dict) -> None:
    """Write one JSON line to stdout and flush."""
    json.dump(obj, sys.stdout, default=str)
    sys.stdout.write("\n")
    sys.stdout.flush()


def handle_run_pipeline(rid: int, params: dict) -> None:
    """Execute the full pipeline for one target and return results."""

    # Parse params
    genome = Path(params["genome"])
    annotations = Path(params["annotations"])
    target_key = params["target_key"]
    target_type = params.get("target_type", "gene")
    primer_args = params.get("primer_args", {})
    max_amplicon = params.get("max_amplicon", 2000)
    pdf_output_dir = params.get("pdf_output_dir")

    # Determine target_gene / target_transcript for PDF report
    target_gene = target_key if target_type == "gene" else None
    target_transcript = target_key if target_type == "transcript" else None

    # ── Run pipeline ──────────────────────────────────────────
    send({"id": rid, "status": "progress", "message": "Starting pipeline...", "pct": 5})

    try:
        result = run_pipeline(
            genome=genome,
            annotations=annotations,
            target_key=target_key,
            target_type=target_type,
            primer_args=primer_args,
            max_amplicon=max_amplicon,
        )
    except ValueError as e:
        send({"id": rid, "status": "error", "message": str(e)})
        return
    except (RuntimeError, FileNotFoundError) as e:
        send({"id": rid, "status": "error", "message": f"tnBLAST error: {e}"})
        return

    # ── Generate PDF ───────────────────────────────────────────
    pdf_path: str | None = None
    if pdf_output_dir and result.filtered_pairs:
        send({"id": rid, "status": "progress", "message": "Generating PDF report...", "pct": 90})
        try:
            out_dir = Path(pdf_output_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            timestamp = datetime.now().strftime("%Y-%m-%d")
            pdf_name = f"brimer-plast_{target_key}_{timestamp}.pdf"
            pdf_full = out_dir / pdf_name

            build_pdf_report(
                output_path=str(pdf_full),
                chains=result.chains,
                locus=result.locus,
                filtered_pairs=result.filtered_pairs,
                target_gene=target_gene,
                target_transcript=target_transcript,
                genome_path=str(genome),
                annotations_path=str(annotations),
            )
            pdf_path = str(pdf_full)
        except Exception as e:
            # PDF failure is non-fatal — still return primer results
            send({"id": rid, "status": "progress", "message": f"PDF generation failed: {e}", "pct": 95})

    # ── Return results ─────────────────────────────────────────
    pairs_data = []
    for pair in result.filtered_pairs:
        pairs_data.append({
            "pair_name": pair.pair_name,
            "chain_id": pair.chain_id,
            "forward_seq": pair.forward_seq,
            "reverse_seq": pair.reverse_seq,
            "forward_tm": pair.forward_tm,
            "reverse_tm": pair.reverse_tm,
            "forward_gc": pair.forward_gc,
            "reverse_gc": pair.reverse_gc,
            "product_size": pair.product_size,
            "pair_number": pair.pair_number,
        })

    send({
        "id": rid,
        "status": "ok",
        "result": {
            "filtered_pairs": pairs_data,
            "warnings": result.warnings,
            "pdf_path": pdf_path,
        },
    })


def main() -> None:
    """Read JSON requests from stdin and dispatch them."""
    # Suppress Brimer-PLAST's own logging so it doesn't mix with our JSON
    logging.getLogger("brimer_plast").setLevel(logging.WARNING)

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
        except json.JSONDecodeError as e:
            send({"id": None, "status": "error", "message": f"Invalid JSON: {e}"})
            continue

        rid = req.get("id")
        command = req.get("command", "")
        params = req.get("params", {})

        if command == "run_pipeline":
            handle_run_pipeline(rid, params)
        else:
            send({"id": rid, "status": "error", "message": f"Unknown command: {command}"})


if __name__ == "__main__":
    main()
