"""JSON-RPC sidecar for Brimer-PLAST Electron app.

Reads one JSON request per line from stdin, writes one JSON response per line
to stdout. The Electron main process spawns this as a child process and
communicates via stdio.

Request format (one JSON line):
    {"id": 1, "command": "run_pipeline", "params": {...}}

Response formats (one JSON line each):
    {"id": 1, "status": "progress", "message": "...", "pct": 50}
    {"id": 1, "status": "ok", "result": {...}}
    {"id": 1, "status": "error", "message": "...", "debug_zip": "/path/to/debug.zip"}
"""

from __future__ import annotations

import json
import logging
import os
import platform
import sys
import traceback
import zipfile
from datetime import datetime
from pathlib import Path

import brimer_plast  # ensure package is importable (no-op)
from brimer_plast.pipeline import run_pipeline
from brimer_plast.pdf_report import build_pdf_report


def send(obj: dict) -> None:
    """Write one JSON line to stdout and flush."""
    json.dump(obj, sys.stdout, default=str)
    sys.stdout.write("\n")
    sys.stdout.flush()


# ── Debug logging helpers ────────────────────────────────────────────────────


def _setup_debug_logging(debug_dir: str | None, rid: int, params: dict) -> Path | None:
    """Create a JSONL debug log file and write the pipeline_start event.

    Returns the log :class:`Path`, or ``None`` if *debug_dir* is falsy.
    """
    if not debug_dir:
        return None
    p = Path(debug_dir)
    p.mkdir(parents=True, exist_ok=True)
    log_path = p / f"pipeline_{rid}.jsonl"
    _write_debug_log(log_path, {
        "event": "pipeline_start",
        "timestamp": datetime.now().isoformat(),
        "target_key": params.get("target_key"),
        "target_type": params.get("target_type"),
        "primer_args": params.get("primer_args"),
        "tnblast_timeout": params.get("tnblast_timeout"),
        "genome_path": params.get("genome"),
        "gtf_path": params.get("annotations"),
    })
    return log_path


def _write_debug_log(log_path: Path | None, event: dict) -> None:
    """Append a JSON event line to the debug log, if a log path is set."""
    if log_path is None:
        return
    try:
        with open(log_path, "a") as f:
            f.write(json.dumps(event, default=str) + "\n")
    except OSError:
        pass  # best-effort logging


def _create_debug_zip(
    debug_dir: str | None,
    rid: int,
    log_path: Path | None,
    params: dict,
    extra_data: dict | None = None,
) -> str | None:
    """Bundle debug artifacts into a ZIP file.

    Returns the absolute path to the created ZIP, or ``None`` if the
    ZIP could not be created.
    """
    if not debug_dir:
        return None
    debug_path = Path(debug_dir)
    try:
        zip_name = f"debug_{rid}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"
        zip_path = debug_path / zip_name

        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            # JSONL log
            if log_path and log_path.exists():
                zf.write(log_path, "pipeline_log.jsonl")

            # Pipeline parameters
            zf.writestr("pipeline_params.json", json.dumps(params, default=str, indent=2))

            # System information
            sys_lines = [
                f"OS: {platform.system()} {platform.release()}",
                f"Python: {sys.version}",
                f"Platform: {platform.platform()}",
            ]
            zf.writestr("system_info.txt", "\n".join(sys_lines) + "\n")

            # Debug artifacts (assay file, tnBLAST outputs)
            for fname in ("assays.txt", "tntblast_genome.txt", "tntblast_transcriptome.txt"):
                fpath = debug_path / fname
                if fpath.exists():
                    zf.write(fpath, fname)

            # Extra error details
            if extra_data:
                zf.writestr("error_details.json", json.dumps(extra_data, default=str, indent=2))

        return str(zip_path)
    except (OSError, zipfile.BadZipFile) as exc:
        # ZIP creation failure is non-fatal — return None
        sys.stderr.write(f"Failed to create debug ZIP: {exc}\n")
        return None


# ── Helper: progress callback with debug logging ────────────────────────────


def _make_progress_callback(rid: int, log_path: Path | None):
    """Build a progress callback that sends IPC messages *and* writes to the debug log."""
    def cb(pct: int, msg: str) -> None:
        send({"id": rid, "status": "progress", "message": msg, "pct": pct})
        _write_debug_log(log_path, {"event": "progress", "pct": pct, "message": msg})
    return cb


def _make_debug_log_callback(log_path: Path | None):
    """Build a callback that writes structured events to the debug log."""
    def cb(event: dict) -> None:
        _write_debug_log(log_path, event)
    return cb


# ── Pipeline handler ─────────────────────────────────────────────────────────


def handle_run_pipeline(rid: int, params: dict) -> None:
    """Execute the full pipeline for one target and return results.

    All exceptions are caught, a debug ZIP is created, and the error
    response includes ``debug_zip`` pointing to it.
    """
    # ── Debug logging setup ──────────────────────────────────────────────
    debug_dir = params.get("debug_dir")
    log_path = _setup_debug_logging(debug_dir, rid, params)

    # Parse params
    genome = Path(params["genome"])
    annotations = Path(params["annotations"])
    target_key = params["target_key"]
    target_type = params.get("target_type", "gene")
    primer_args = params.get("primer_args", {})
    max_amplicon = params.get("max_amplicon", 2000)
    tnblast_timeout = params.get("tnblast_timeout", 1800)
    pdf_output_dir = params.get("pdf_output_dir")

    # Determine target_gene / target_transcript for PDF report
    target_gene = target_key if target_type == "gene" else None
    target_transcript = target_key if target_type == "transcript" else None

    # ── Run pipeline ────────────────────────────────────────────────────
    send({"id": rid, "status": "progress", "message": "Starting pipeline...", "pct": 5})
    _write_debug_log(log_path, {"event": "progress", "pct": 5, "message": "Starting pipeline..."})

    try:
        result = run_pipeline(
            genome=genome,
            annotations=annotations,
            target_key=target_key,
            target_type=target_type,
            primer_args=primer_args,
            max_amplicon=max_amplicon,
            tnblast_timeout=tnblast_timeout,
            debug_dir=debug_dir,
            debug_log_callback=_make_debug_log_callback(log_path),
            progress_callback=_make_progress_callback(rid, log_path),
        )
    except Exception as e:
        tb = traceback.format_exc()
        _write_debug_log(log_path, {"event": "pipeline_end", "status": "error", "traceback": tb})
        zip_path = _create_debug_zip(debug_dir, rid, log_path, params, {"traceback": tb})
        send({
            "id": rid,
            "status": "error",
            "message": str(e),
            "debug_zip": zip_path,
        })
        return

    # ── Generate PDF ─────────────────────────────────────────────────────
    pdf_path: str | None = None
    if pdf_output_dir and result.filtered_pairs:
        send({"id": rid, "status": "progress", "message": "Generating PDF report...", "pct": 90})
        _write_debug_log(log_path, {"event": "progress", "pct": 90, "message": "Generating PDF report..."})
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
            msg = f"PDF generation failed: {e}"
            send({"id": rid, "status": "progress", "message": msg, "pct": 95})
            _write_debug_log(log_path, {"event": "progress", "pct": 95, "message": msg})

    # ── Return results ───────────────────────────────────────────────────
    _write_debug_log(log_path, {"event": "pipeline_end", "status": "ok"})

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


# ── Main loop ────────────────────────────────────────────────────────────────


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