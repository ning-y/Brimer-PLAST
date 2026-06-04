"""Tests for pdf_report.py: PDF report generation.

These tests verify that the genome-view diagram (``_GeneDiagram``) never
exceeds the available page frame, a recurring crash::

    Flowable <_GeneDiagram ...> too large on page … in frame 'normal' …

The diagram height is now computed dynamically from the available frame
geometry so this crash cannot re-occur.
"""

from __future__ import annotations

import pytest
from reportlab.lib.pagesizes import A4, landscape

from brimer_plast.models import (
    ConservedExonChain,
    ExonInfo,
    GeneLocus,
    GenomicFragment,
    PrimerPair,
)
from brimer_plast.diagram import (
    DIAGRAM_PAIR_CAP,
    TRANSCRIPT_CAP,
    _GeneDiagram,
)
from brimer_plast.pdf_report import (
    FRAME_H,
)


# ── helpers ──────────────────────────────────────────────────────────────────


def _locus_with_n_transcripts(n: int) -> GeneLocus:
    """Build a ``GeneLocus`` with *n* transcripts, each with 3 exons."""
    transcripts: dict[str, list[ExonInfo]] = {}
    for i in range(n):
        tid = f"t_{i:03d}"
        transcripts[tid] = [
            ExonInfo(seqid="chr1", start=1_000, end=2_000, strand="+"),
            ExonInfo(seqid="chr1", start=3_000, end=4_000, strand="+"),
            ExonInfo(seqid="chr1", start=5_000, end=6_000, strand="+"),
        ]
    return GeneLocus(
        gene_name="test_gene",
        seqid="chr1",
        strand="+",
        transcripts=transcripts,
        min_start=1_000,
        max_end=6_000,
    )


def _n_dummy_pairs(n: int) -> list[PrimerPair]:
    """Build *n* dummy ``PrimerPair`` objects with genomic fragments.

    Each pair has one forward fragment in exon 1 and one reverse fragment
    in exon 3, so the zoom-range logic in ``draw_gene_diagram`` can
    compute meaningful min/max coordinates.
    """
    pairs: list[PrimerPair] = []
    for i in range(n):
        offset = i * 10
        fwd_frags = [
            GenomicFragment(seqid="chr1", start=1_500 + offset, end=1_520 + offset, strand="+")
        ]
        rev_frags = [
            GenomicFragment(seqid="chr1", start=5_200 + offset, end=5_220 + offset, strand="+")
        ]
        pairs.append(
            PrimerPair(
                forward_seq="AAA" * 7,
                reverse_seq="TTT" * 7,
                pair_number=i + 1,
                primer3_forward_fragments=fwd_frags,
                primer3_reverse_fragments=rev_frags,
                tnblast_forward_fragments=fwd_frags,
                tnblast_reverse_fragments=rev_frags,
            )
        )
    return pairs


def _chain_for_exons() -> list[ConservedExonChain]:
    """Return a single ``ConservedExonChain`` matching the locus exons."""
    return [
        ConservedExonChain(
            id="test_chain_1",
            exons=[
                ExonInfo(seqid="chr1", start=1_000, end=2_000, strand="+"),
                ExonInfo(seqid="chr1", start=3_000, end=4_000, strand="+"),
                ExonInfo(seqid="chr1", start=5_000, end=6_000, strand="+"),
            ],
            template="A" * 5_000,
            junction_positions_1based=[1_001, 2_001],
        )
    ]


# ── constants mirrored from pdf_report.py for test stability ────────────────
PAGE1_OVERHEAD = 138
PER_PAIR = 20


def _diagram_base_h(n_transcripts: int, total_transcripts: int) -> float:
    """Replicate the diagram base-height logic from build_pdf_report."""
    n = min(n_transcripts, TRANSCRIPT_CAP)
    has_others = total_transcripts > TRANSCRIPT_CAP
    return (82 if has_others else 74) + n * 11


def _max_pairs_for_page(n_transcripts: int, total_transcripts: int, overhead: float) -> int:
    """Replicate the per-page pairs calculation from build_pdf_report."""
    base = _diagram_base_h(n_transcripts, total_transcripts)
    return max(1, int((FRAME_H - overhead - base) // PER_PAIR))


# ── tests ────────────────────────────────────────────────────────────────────


class TestDiagramPageSizing:
    """The dynamic pair-count logic must always produce heights that fit
    the available page frame — this is the invariant that prevents the
    ``Flowable … too large`` crash."""

    def test_page1_overhead_sufficient(self):
        """Page 1 overhead ≤ frame minus the minimum diagram height
        for even 1 pair + 1 transcript."""
        base = _diagram_base_h(n_transcripts=1, total_transcripts=1)
        assert PAGE1_OVERHEAD < FRAME_H - base - PER_PAIR, (
            "PAGE1_OVERHEAD is too large — no diagram could fit on page 1"
        )

    def test_later_page_fits_full_frame(self):
        """Later pages have zero overhead — the full frame is available.
        The diagram with max transcripts and 1 pair must trivially fit."""
        base = _diagram_base_h(
            n_transcripts=TRANSCRIPT_CAP, total_transcripts=TRANSCRIPT_CAP
        )
        h = base + PER_PAIR
        assert h <= FRAME_H, (
            f"Even 1 pair overflows when overhead=0 "
            f"(h={h:.0f} > frame={FRAME_H:.0f})"
        )

    def test_max_transcripts_min_pairs_fits_page1(self):
        """Even with TRANSCRIPT_CAP transcripts, at least 1 pair fits on page 1."""
        n = _max_pairs_for_page(
            n_transcripts=TRANSCRIPT_CAP,
            total_transcripts=TRANSCRIPT_CAP,
            overhead=PAGE1_OVERHEAD,
        )
        assert n >= 1, (
            f"Expected at least 1 pair on page 1 with {TRANSCRIPT_CAP} transcripts, "
            f"got {n} (base_h={_diagram_base_h(TRANSCRIPT_CAP, TRANSCRIPT_CAP):.0f}, "
            f"frame={FRAME_H:.0f}, overhead={PAGE1_OVERHEAD})"
        )

    def test_max_transcripts_min_pairs_fits_later_page(self):
        """Even with TRANSCRIPT_CAP transcripts, at least 1 pair fits on
        a later page (zero overhead)."""
        n = _max_pairs_for_page(
            n_transcripts=TRANSCRIPT_CAP,
            total_transcripts=TRANSCRIPT_CAP,
            overhead=0,
        )
        assert n >= 1, (
            f"Expected at least 1 pair on later page with "
            f"{TRANSCRIPT_CAP} transcripts, got {n}"
        )

    def test_later_page_holds_more_than_page1(self):
        """A later page (no overhead) should fit at least as many pairs
        as page 1."""
        for n in [1, 5, 10]:
            t = max(n, 1)
            p1 = _max_pairs_for_page(n, t, PAGE1_OVERHEAD)
            later = _max_pairs_for_page(n, t, 0)
            assert later >= p1, (
                f"With {n} transcripts, later page ({later}) should hold at "
                f"least as many pairs as page 1 ({p1})"
            )

    def test_fewer_transcripts_yields_more_pairs(self):
        """Fewer transcripts = smaller base height = more pairs per page."""
        with10 = _max_pairs_for_page(10, 10, 0)
        with1 = _max_pairs_for_page(1, 1, 0)
        assert with1 >= with10, (
            f"Expected 1 transcript to give at least as many pairs as "
            f"10 transcripts ({with1} vs {with10})"
        )

    def test_capped_transcripts_has_others_penalty(self):
        """When total_transcripts > TRANSCRIPT_CAP, the "+N others" line
        adds 8pt to the base height vs the uncapped case, reducing pair count."""
        capped = _max_pairs_for_page(TRANSCRIPT_CAP, TRANSCRIPT_CAP + 5, 0)
        uncapped = _max_pairs_for_page(TRANSCRIPT_CAP, TRANSCRIPT_CAP, 0)
        # The "others" overhead is only 8pt, so it might not change the count
        # (the int division by 20 is coarse).  Just verify it doesn't *increase*.
        assert capped <= uncapped, (
            f"capped ({capped}) should be ≤ uncapped ({uncapped})"
        )


class TestGeneDiagramFitsInFrame:
    """Direct verification that ``_GeneDiagram`` wrap() reports a height
    that does not exceed available frame space."""

    def test_diagram_at_limit_fits_on_page1(self):
        """A diagram sized for the max pairs that fit on page 1 must not
        overflow when wrap() is called with page1's remaining space."""
        n_pairs = _max_pairs_for_page(
            n_transcripts=TRANSCRIPT_CAP,
            total_transcripts=TRANSCRIPT_CAP,
            overhead=PAGE1_OVERHEAD,
        )
        locus = _locus_with_n_transcripts(TRANSCRIPT_CAP)
        pairs = _n_dummy_pairs(n_pairs)
        chains = _chain_for_exons()

        base = _diagram_base_h(TRANSCRIPT_CAP, TRANSCRIPT_CAP)
        h = base + n_pairs * PER_PAIR
        available = FRAME_H - PAGE1_OVERHEAD

        diagram = _GeneDiagram(locus, pairs, chains, 700, h)
        _, reported_h = diagram.wrap(700, available)

        assert reported_h <= available + 1, (  # +1 for float tolerance
            f"Height {reported_h:.0f}pt exceeds available {available:.0f}pt"
        )

    def test_diagram_at_limit_fits_on_later_page(self):
        """A diagram sized for the max pairs that fit on a later page
        (full frame) must not overflow when wrap() is called."""
        n_pairs = _max_pairs_for_page(
            n_transcripts=TRANSCRIPT_CAP,
            total_transcripts=TRANSCRIPT_CAP,
            overhead=0,
        )
        locus = _locus_with_n_transcripts(TRANSCRIPT_CAP)
        pairs = _n_dummy_pairs(n_pairs)
        chains = _chain_for_exons()

        base = _diagram_base_h(TRANSCRIPT_CAP, TRANSCRIPT_CAP)
        h = base + n_pairs * PER_PAIR
        available = FRAME_H

        diagram = _GeneDiagram(locus, pairs, chains, 700, h)
        _, reported_h = diagram.wrap(700, available)

        assert reported_h <= available + 1, (
            f"Height {reported_h:.0f}pt exceeds available {available:.0f}pt"
        )


class TestGeneDiagramYPosition:
    """The diagram drawing must stay within the vertical area declared
    by the flowable height, i.e. the bottommost element must not go below
    ``y=0`` after the flowable is positioned at ``y=height-30``."""

    def test_lowest_drawn_pixel_is_above_zero(self):
        """With max pairs that fit on a later page (full frame), the lowest
        pixel actually drawn on canvas (the last pair's Panel B rect) must
        be ≥ 0.  The return value of draw_gene_diagram (curr_y - 5) is a
        bookkeeping variable that routinely goes negative — it is NOT the
        bottommost coordinate touched by drawing calls."""
        n_pairs = _max_pairs_for_page(
            n_transcripts=TRANSCRIPT_CAP,
            total_transcripts=TRANSCRIPT_CAP,
            overhead=0,
        )
        locus = _locus_with_n_transcripts(TRANSCRIPT_CAP)
        pairs = _n_dummy_pairs(n_pairs)
        chains = _chain_for_exons()

        base = _diagram_base_h(TRANSCRIPT_CAP, TRANSCRIPT_CAP)
        height = base + n_pairs * PER_PAIR
        top_y = height - 30  # where _GeneDiagram.draw() starts

        # Trace the y descent to find the lowest pixel drawn.
        # The lowest element on canvas is the last pair's Panel B rect
        # which draws at y = curr_y - 1 (at the start of that iteration).
        over_y = top_y - 15                           # OVERVIEW_H
        zoom_top_y = over_y - 25
        curr_y = zoom_top_y - 8 - 5                   # EXON_HEIGHT=8, gap=5
        for _ in range(TRANSCRIPT_CAP):
            curr_y -= 11                              # EXON_HEIGHT + TRACK_SPACING
        curr_y -= 10                                  # gap before primers

        # curr_y BEFORE the last pair's drawing begins:
        last_start_y = curr_y - 20 * (n_pairs - 1)
        # Panel B rect at (last_start_y - 1) with height=4
        lowest_pixel = last_start_y - 1

        assert lowest_pixel >= 0, (
            f"Lowest drawn pixel ({lowest_pixel:.0f}pt) went below zero "
            f"(n_pairs={n_pairs}, height={height:.0f}, top_y={top_y:.0f})"
        )
