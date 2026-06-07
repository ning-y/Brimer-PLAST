"""Primer pair filtering for Brimer-PLAST.

Decides which candidate primer pairs pass the specificity filter after
tnBLAST has been run against both genome and transcriptome databases.
"""

from __future__ import annotations

from brimer_plast.models import PrimerPair, GeneLocus
from brimer_plast.tnblast import AmpliconHit


def filter_specific_pairs(
    primer_pairs: list[PrimerPair],
    genome_amplicons: dict[str, list[AmpliconHit]],
    transcriptome_targets: dict[str, list[str]],
    *,
    target_gene: str,
    target_locus: GeneLocus | None = None,
    junction_mode: bool = True,
    locus_padding: int = 1000,
) -> list[PrimerPair]:
    """Return only primer pairs that are specific in both genome and transcriptome.

    In junction mode (default, at least one primer spans an exon-exon junction),
    a pair passes if:
      - genome count == 0   (junction-spanning primers don't match genome; avoid gDNA and processed pseudogenes)
      - all transcriptome hits map to the target gene

    In non-junction mode (intron-spanning), a pair passes if:
      - all genomic hits (if any) are located within the target gene locus
      - all transcriptome hits map to the target gene

    Args:
        primer_pairs: Candidate pairs from primer3.
        genome_amplicons: Pair -> list of AmpliconHit from tnBLAST against genome.
        transcriptome_targets: Pair -> list of gene names hit in transcriptome.
        target_gene: The gene the primers are designed for.
        target_locus: Coordinates of the target gene for Mode B location check.
        junction_mode: Whether junction-spanning is enforced.
        locus_padding: Padding (bp) around target_locus for genomic hits.

    Each pair is looked up by its :attr:`PrimerPair.pair_name`.
    """
    specific: list[PrimerPair] = []
    for pair in primer_pairs:
        name = pair.pair_name
        if not name:
            continue
        hits = genome_amplicons.get(name, [])
        targenes = transcriptome_targets.get(name, [])

        # All transcriptome hits must belong to the target gene
        tc_ok = len(targenes) > 0 and all(g == target_gene for g in targenes)

        if junction_mode:
            if len(hits) == 0 and tc_ok:
                specific.append(pair)
        else:
            if not tc_ok:
                continue

            all_on_target = True
            for hit in hits:
                if target_locus:
                    on_locus = (
                        hit.seqid == target_locus.seqid
                        and hit.amplicon_start >= target_locus.min_start - locus_padding
                        and hit.amplicon_end <= target_locus.max_end + locus_padding
                    )
                    if not on_locus:
                        all_on_target = False
                        break
                else:
                    # If no locus info, we fallback to 1 hit max as a heuristic
                    if len(hits) > 1:
                        all_on_target = False
                    break

            if all_on_target:
                specific.append(pair)

    return specific
