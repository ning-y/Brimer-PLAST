"""Primer pair filtering for Brimer-PLAST.

Decides which candidate primer pairs pass the specificity filter after
tnBLAST has been run against both genome and transcriptome databases.
"""

from __future__ import annotations

from brimer_plast.models import PrimerPair


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