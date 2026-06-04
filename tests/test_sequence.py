"""Tests for sequence.py: exon ordering, coordinate mapping."""

import pytest

from brimer_plast.models import ExonInfo, GenomicFragment
from brimer_plast.sequence import (
    exons_in_template_order,
    genomic_range_to_fragments,
    reverse_complement,
    template_to_genomic,
)


# ── data helpers ─────────────────────────────────────────────────────────────

_POS_EXONS = [
    ExonInfo(seqid="chrI", start=101, end=350, strand="+"),
    ExonInfo(seqid="chrI", start=451, end=600, strand="+"),
]

_NEG_EXONS = [
    ExonInfo(seqid="chrI", start=500, end=600, strand="-"),
    ExonInfo(seqid="chrI", start=700, end=800, strand="-"),
]

_THREE_EXON_POS = [
    ExonInfo(seqid="chrI", start=100, end=199, strand="+"),
    ExonInfo(seqid="chrI", start=300, end=399, strand="+"),
    ExonInfo(seqid="chrI", start=500, end=599, strand="+"),
]


# ── exons_in_template_order ──────────────────────────────────────────────────


class TestExonsInTemplateOrder:
    def test_positive_strand_ascending(self):
        """Positive strand: exons sorted by start ascending."""
        result = exons_in_template_order(_POS_EXONS)
        assert result[0].start == 101
        assert result[1].start == 451

    def test_negative_strand_descending(self):
        """Negative strand: exons sorted by start descending (3' first)."""
        result = exons_in_template_order(_NEG_EXONS)
        assert result[0].start == 700  # downstream exon first
        assert result[1].start == 500

    def test_empty_returns_empty_list(self):
        """Empty input returns empty list."""
        assert exons_in_template_order([]) == []

    def test_single_exon_returns_that_exon(self):
        """Single exon returns itself unmodified."""
        ex = ExonInfo("chrI", 100, 200, "+")
        result = exons_in_template_order([ex])
        assert len(result) == 1
        assert result[0] == ex

    def test_three_exons_positive_strand(self):
        """Three exons on positive strand stay in genomic order."""
        result = exons_in_template_order(_THREE_EXON_POS)
        starts = [e.start for e in result]
        assert starts == [100, 300, 500]

    def test_three_exons_negative_strand(self):
        """Three exons on negative strand are reversed."""
        exons = [
            ExonInfo("chrI", 500, 599, "-"),
            ExonInfo("chrI", 300, 399, "-"),
            ExonInfo("chrI", 100, 199, "-"),
        ]
        result = exons_in_template_order(exons)
        starts = [e.start for e in result]
        assert starts == [500, 300, 100]

    def test_input_not_mutated(self):
        """The function should not mutate the input list."""
        original = [ExonInfo("chrI", 500, 600, "-"), ExonInfo("chrI", 100, 200, "-")]
        copy = list(original)
        exons_in_template_order(original)
        assert original == copy  # should not be sorted in-place


# ── template_to_genomic ──────────────────────────────────────────────────────


class TestTemplateToGenomic:
    """Maps 0-based template positions to genomic coordinates."""

    def test_single_exon_positive_strand(self):
        """Primer entirely within first exon (+, within bounds)."""
        # Exon1: 101-350 = 250bp. Template pos 0 = chrI:101.
        result = template_to_genomic(0, 20, _POS_EXONS)
        assert result == [GenomicFragment("chrI", 101, 120, "+")]

    def test_primer_in_second_exon(self):
        """Primer entirely within second exon (+, within bounds)."""
        # Exon1 = 250bp, template pos 250 = start of exon2 = chrI:451
        result = template_to_genomic(250, 30, _POS_EXONS)
        assert result == [GenomicFragment("chrI", 451, 480, "+")]

    def test_spanning_junction_two_fragments(self):
        """Primer spanning exon1→exon2 junction returns two fragments."""
        # Junction between exon1 (250bp) and exon2 is at template pos 250.
        # Primer at template pos 240, length 20:
        #   - 10bp in exon1 (pos 240-249)
        #   - 10bp in exon2 (pos 250-259)
        result = template_to_genomic(240, 20, _POS_EXONS)
        assert len(result) == 2
        # Exon1 fragment: g_start = 101 + 240 = 341, g_end = 341+10-1 = 350
        assert GenomicFragment("chrI", 341, 350, "+") in result
        # Exon2 fragment: g_start = 451, g_end = 451+10-1 = 460
        assert GenomicFragment("chrI", 451, 460, "+") in result

    def test_positive_last_exon_boundary(self):
        """Primer exactly at the end of the first exon (no spanning)."""
        # Template pos 249 = last base of exon1 = chrI:350
        result = template_to_genomic(249, 1, _POS_EXONS)
        assert result == [GenomicFragment("chrI", 350, 350, "+")]

    def test_single_exon_negative_strand(self):
        """Primer entirely within first exon (-, within bounds)."""
        # Template order: [Exon(700,800,-), Exon(500,600,-)]
        # Template pos 0 = 3' end of downstream exon = chrI:800
        result = template_to_genomic(0, 10, _NEG_EXONS)
        assert result == [GenomicFragment("chrI", 791, 800, "-")]

    def test_negative_strand_second_exon(self):
        """Primer entirely within second exon on negative strand."""
        # Exon1 (700-800) = 101bp of template
        # Template pos 101 = 3' end of exon2 (500-600) = chrI:600
        result = template_to_genomic(101, 15, _NEG_EXONS)
        assert result == [GenomicFragment("chrI", 586, 600, "-")]

    def test_spanning_junction_negative_strand(self):
        """Primer spanning exon2→exon1 junction on negative strand."""
        # Junction between exponential exons at template pos 101.
        # Primer at template pos 95, length 12:
        #   - 6bp in exon1 (pos 95-100)
        #   - 6bp in exon2 (pos 101-106)
        result = template_to_genomic(95, 12, _NEG_EXONS)
        assert len(result) == 2
        # Exon1 fragment: g_end = 800-95 = 705, g_start = 705-6+1 = 700
        assert GenomicFragment("chrI", 700, 705, "-") in result
        # Exon2 fragment: g_end = 600, g_start = 600-6+1 = 595
        assert GenomicFragment("chrI", 595, 600, "-") in result

    def test_spanning_three_exons(self):
        """Primer long enough to span two junctions (three exons)."""
        # Exons: [100-199] 100bp, [300-399] 100bp, [500-599] 100bp
        # Template: 0-99=exon1, 100-199=exon2, 200-299=exon3
        # Primer at template pos 90, length 30:
        #   10bp in exon1 (pos 90-99)
        #   20bp in exon2 (pos 100-119) — but remaining is only 20 after 10 in first
        # Wait, let me recalculate:
        #   offset_in_exon = 90, len_in_this = min(30, 100-90) = 10
        #   remaining = 20
        #   Next exon (300-399): len = min(20, 100) = 20
        #   remaining = 0
        # So only 2 exons.
        #
        # To span 3 exons: template pos 95, length 20
        #   offset = 95, len = min(20, 100-95) = 5, remaining = 15
        #   exon2: len = min(15, 100) = 15, remaining = 0
        # Nope. Wait:
        #   First exon: len = min(20, 5) = 5, remaining = 15
        #   Second exon: len = min(15, 100) = 15, remaining = 0
        # Only 2 exons.
        #
        # template pos 90, length 30
        #   First: len = min(30, 10) = 10, remaining = 20
        #   Second: len = min(20, 100) = 20, remaining = 0
        # Only 2 exons again.
        #
        # template pos 95, length 30
        #   First: len = min(30, 5) = 5, remaining = 25
        #   Second: len = min(25, 100) = 25, remaining = 0
        # Still 2.
        #
        # template pos 95, length 40
        #   First: len = min(40, 5) = 5, remaining = 35
        #   Second: len = min(35, 100) = 35, remaining = 0
        # Still 2. I need more length.
        #
        # The issue is the remaining primer after exon1 is only offset_in_exon = 5 bases available,
        # and exon2 has 100 bases. So I need remaining > 100 to span into exon3.
        # remaining_after_exon1 = 40 - 5 = 35, which is < 100. So 40 doesn't work.
        # I need remaining_after_exon1 > 100... that means total length - 5 > 100 → length > 105.
        # template pos 95, length 106:
        #   First: len = min(106, 5) = 5, remaining = 101
        #   Second: len = min(101, 100) = 100, remaining = 1
        #   Third: len = min(1, 100) = 1, remaining = 0
        # That gives 3 fragments!
        result = template_to_genomic(95, 106, _THREE_EXON_POS)
        assert len(result) == 3
        assert result[0] == GenomicFragment("chrI", 195, 199, "+")
        assert result[1] == GenomicFragment("chrI", 300, 399, "+")
        assert result[2] == GenomicFragment("chrI", 500, 500, "+")

    def test_empty_exon_list_returns_empty(self):
        """Empty exons list returns empty fragments list."""
        result = template_to_genomic(0, 20, [])
        assert result == []


# ── genomic_range_to_fragments ───────────────────────────────────────────────


class TestGenomicRangeToFragments:
    """Splits a genomic coordinate range across exons."""

    def test_range_within_single_exon(self):
        """Range entirely within one exon returns a single fragment."""
        result = genomic_range_to_fragments(200, 250, _POS_EXONS)
        assert result == [GenomicFragment("chrI", 200, 250, "+")]

    def test_range_spanning_two_exons(self):
        """Range crossing exon1→exon2 boundary returns two fragments."""
        # Exons: 101-350, 451-600. Range 300-500 spans both.
        result = genomic_range_to_fragments(300, 500, _POS_EXONS)
        assert len(result) == 2
        assert result[0] == GenomicFragment("chrI", 300, 350, "+")
        assert result[1] == GenomicFragment("chrI", 451, 500, "+")

    def test_range_no_overlap_returns_empty(self):
        """Range outside all exons returns empty list."""
        result = genomic_range_to_fragments(1, 50, _POS_EXONS)
        assert result == []

    def test_range_covers_entire_gene(self):
        """Range covering all exons returns fragments for each."""
        result = genomic_range_to_fragments(100, 600, _POS_EXONS)
        assert len(result) == 2
        assert result[0] == GenomicFragment("chrI", 101, 350, "+")
        assert result[1] == GenomicFragment("chrI", 451, 600, "+")

    def test_range_start_before_exact_boundary(self):
        """Range that starts before the first exon clamps to exon start."""
        result = genomic_range_to_fragments(50, 200, _POS_EXONS)
        assert result == [GenomicFragment("chrI", 101, 200, "+")]

    def test_range_end_after_exact_boundary(self):
        """Range that ends after the last exon clamps to exon end."""
        result = genomic_range_to_fragments(500, 700, _POS_EXONS)
        assert result == [GenomicFragment("chrI", 500, 600, "+")]

    def test_three_exon_range(self):
        """Range spanning three exons returns three fragments."""
        result = genomic_range_to_fragments(100, 599, _THREE_EXON_POS)
        assert len(result) == 3
        assert result[0] == GenomicFragment("chrI", 100, 199, "+")
        assert result[1] == GenomicFragment("chrI", 300, 399, "+")
        assert result[2] == GenomicFragment("chrI", 500, 599, "+")

    def test_empty_exons_returns_empty(self):
        """Empty exons list returns empty fragments list."""
        result = genomic_range_to_fragments(100, 200, [])
        assert result == []


# ── reverse_complement ───────────────────────────────────────────────────────


class TestReverseComplement:
    def test_reverses_and_complements(self):
        """Standard case."""
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        """Palindromic sequence is its own reverse complement."""
        assert reverse_complement("CGATCG") == "CGATCG"

    def test_empty_string(self):
        """Empty input returns empty."""
        assert reverse_complement("") == ""

    def test_handles_lowercase(self):
        """Lowercase input works correctly."""
        assert reverse_complement("atcg") == "cgat"
