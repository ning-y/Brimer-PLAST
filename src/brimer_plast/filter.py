"""tnBLAST integration for Brimer-PLAST.

Wraps tnBLAST subprocess calls and filters primer pairs based on
predicted off-target amplification.
"""

import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from brimer_plast.models import PrimerPair


def write_assay_file(
    primer_pairs: list[PrimerPair],
    path: str | Path,
) -> None:
    """Write primer pairs as a tnBLAST tab-delimited assay file.

    Format:  name\\tforward_seq\\treverse_seq
    """
    with open(path, "w") as f:
        for i, pair in enumerate(primer_pairs, start=1):
            f.write(f"pair_{i}\t{pair.forward_seq}\t{pair.reverse_seq}\n")


def run_tnblast(
    assay_path: str | Path,
    database_path: str | Path,
    *,
    max_amplicon: int = 2000,
    min_tm: float = 55.0,
    max_tm: float = 65.0,
    salt_conc: float = 0.05,
    output_path: str | Path | None = None,
) -> dict[str, int]:
    """Run tnBLAST and return a dict mapping assay name to amplicon count.

    Args:
        assay_path: Path to the tnBLAST assay file.
        database_path: Path to the FASTA database (genome or transcriptome).
        max_amplicon: Maximum amplicon length for tnBLAST search.
        min_tm: Minimum primer Tm.
        max_tm: Maximum primer Tm.
        salt_conc: Salt concentration (M).
        output_path: Optional explicit path for tnBLAST output. If provided,
            the file is NOT cleaned up — the caller owns it. If None, a temp
            file is used and cleaned up automatically.

    Returns:
        { "pair_1": 3, "pair_2": 1, ... }
        The count is the number of predicted amplicons (1 = specific,
        >1 = off-target amplification, 0 = no amplification).

    Raises:
        RuntimeError: tnBLAST exited with a non-zero status.
        FileNotFoundError: tnBLAST is not installed or not on PATH.
    """
    owns_tmp = output_path is None
    if owns_tmp:
        tmp_dir = tempfile.mkdtemp()
        out_path = os.path.join(tmp_dir, "tntblast_output.txt")
    else:
        tmp_dir = None
        out_path = str(output_path)

    try:
        cmd = [
            "tntblast",
            "-i",
            str(assay_path),
            "-d",
            str(database_path),
            "-o",
            out_path,
            "-l",
            str(max_amplicon),
            "-e",
            str(min_tm),
            "-x",
            str(max_tm),
            "-s",
            str(salt_conc),
        ]
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                "tntblast: command not found. Ensure tnBLAST is installed "
                "and available on PATH."
            ) from None
        if result.returncode != 0:
            raise RuntimeError(
                f"tnBLAST failed (exit {result.returncode}):\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )

        return _parse_tnblast_output(out_path)
    finally:
        if owns_tmp and tmp_dir is not None:
            import shutil

            shutil.rmtree(tmp_dir, ignore_errors=True)


def filter_specific_pairs(
    primer_pairs: list[PrimerPair],
    genome_counts: dict[str, int],
    transcriptome_targets: dict[str, list[str]],
    *,
    target_gene: str,
    junction_mode: bool = True,
) -> list[PrimerPair]:
    """Return only primer pairs that are specific in both genome and transcriptome.

    In junction mode (default, at least one primer spans an exon-exon junction),
    a pair passes if:
      - genome count == 0   (junction-spanning primers don't match genome)
      - all transcriptome hits map to the target gene

    In non-junction mode (``--disable-junction-overlap``), a pair passes if:
      - genome count == 1   (single genomic target)
      - all transcriptome hits map to the target gene

    Args:
        primer_pairs: Candidate pairs from primer3.
        genome_counts: Pair -> amplicon count from tnBLAST against genome.
        transcriptome_targets: Pair -> list of gene names hit in transcriptome.
        target_gene: The gene the primers are designed for.
        junction_mode: Whether junction-spanning is enforced.
    """
    specific: list[PrimerPair] = []
    for i, pair in enumerate(primer_pairs, start=1):
        name = f"pair_{i}"
        gc = genome_counts.get(name, 0)
        targenes = transcriptome_targets.get(name, [])

        # All transcriptome hits must belong to the target gene
        tc_ok = len(targenes) > 0 and all(g == target_gene for g in targenes)

        if junction_mode:
            if gc == 0 and tc_ok:
                specific.append(pair)
        else:
            if gc == 1 and tc_ok:
                specific.append(pair)
    return specific


# ── helpers ──────────────────────────────────────────────────────────────────


def _parse_tnblast_output(path: str) -> dict[str, int]:
    """Parse tnBLAST output into assay-name -> amplicon-count mapping.

    Each predicted amplicon starts with a ``name = <assay_name>`` line.
    We count how many times each name appears.
    """
    counts: dict[str, int] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("name = "):
                name = line[len("name = ") :].strip()
                counts[name] = counts.get(name, 0) + 1
    return counts


@dataclass
class AmpliconHit:
    """A single tnBLAST amplicon prediction with genomic coordinates."""

    seqid: str
    amplicon_start: int  # 1-based
    amplicon_end: int  # 1-based


def _parse_tnblast_amplicons(path: str) -> dict[str, list[AmpliconHit]]:
    """Parse tnBLAST output into assay-name -> list of AmpliconHit.

    Extracts both the target sequence ID (from ``>seqid`` lines) and the
    amplicon range (from ``amplicon range = <start> .. <end>`` lines).

    tnBLAST outputs ``amplicon range`` *before* ``>seqid`` within each hit
    section, so this parser defers emitting each hit until the ``>`` line
    has been read — otherwise every hit would get the *previous* hit's seqid.
    """
    hits: dict[str, list[AmpliconHit]] = {}
    current_name: str | None = None
    pending_start: int | None = None
    pending_end: int | None = None

    def _emit() -> None:
        if (
            current_name is not None
            and pending_start is not None
            and pending_end is not None
            and current_seqid is not None
        ):
            hit = AmpliconHit(
                seqid=current_seqid,
                amplicon_start=pending_start,
                amplicon_end=pending_end,
            )
            hits.setdefault(current_name, []).append(hit)

    current_seqid: str | None = None
    with open(path) as f:
        for line in f:
            if line.startswith("name = "):
                # Flush any pending hit before switching to a new name
                _emit()
                pending_start = None
                pending_end = None
                current_name = line[len("name = ") :].strip()
            elif current_name is not None and line.startswith(">"):
                current_seqid = line[1:].strip()
                # The last ``amplicon range`` belongs to this seqid; emit now
                _emit()
                pending_start = None
                pending_end = None
            elif current_name is not None and line.startswith("amplicon range ="):
                # "amplicon range = 285 .. 494"
                m = re.search(r"= (\d+)\s*\.\.\s*(\d+)", line)
                if m:
                    pending_start = int(m.group(1))
                    pending_end = int(m.group(2))
    return hits
