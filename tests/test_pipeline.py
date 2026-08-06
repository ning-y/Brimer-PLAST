"""Tests for pipeline.py: core pipeline and dump_debug_info."""

import logging

from brimer_plast.models import (
    ConservedExonChain,
    ExonInfo,
    GeneLocus,
    GenomicFragment,
    PrimerPair,
)
from brimer_plast.pipeline import PipelineResult, dump_debug_info, make_pair_name


class TestPipelineResult:
    """PipelineResult dataclass correctness."""

    def test_empty_result_has_no_results(self):
        """Empty pipeline result should report has_results=False."""
        result = PipelineResult()
        assert not result.has_results

    def test_with_filtered_pairs_has_results(self):
        """Pipeline result with filtered pairs should report has_results=True."""
        result = PipelineResult(
            filtered_pairs=[PrimerPair(forward_seq="A", reverse_seq="T")]
        )
        assert result.has_results


class TestDumpDebugInfo:
    """Smoke tests for dump_debug_info — should not crash with real data."""

    def test_dump_empty_result(self, caplog):
        """dump_debug_info with empty result should not raise."""
        caplog.set_level(logging.DEBUG, logger="brimer_plast.pipeline")
        result = PipelineResult()
        dump_debug_info(result)
        assert "DEBUG DUMP" in caplog.text

    def test_dump_with_locus(self, caplog):
        """dump_debug_info with a locus should log transcript info."""
        caplog.set_level(logging.DEBUG, logger="brimer_plast.pipeline")
        locus = GeneLocus(
            gene_name="test_gene",
            seqid="chrI",
            strand="+",
            transcripts={
                "t1": [ExonInfo("chrI", 100, 199, "+"), ExonInfo("chrI", 300, 399, "+")],
            },
            min_start=100,
            max_end=399,
        )
        result = PipelineResult(locus=locus)
        dump_debug_info(result)
        assert "test_gene" in caplog.text
        assert "chrI" in caplog.text

    def test_dump_with_chains(self, caplog):
        """dump_debug_info with chains should log chain info."""
        caplog.set_level(logging.DEBUG, logger="brimer_plast.pipeline")
        chain = ConservedExonChain(
            id="test_chain_1",
            exons=[ExonInfo("chrI", 100, 199, "+"), ExonInfo("chrI", 300, 399, "+")],
            template="A" * 200,
            junction_positions_1based=[101],
            required_junction_positions_1based=[101],
        )
        result = PipelineResult(chains=[chain])
        dump_debug_info(result)
        assert "test_chain_1" in caplog.text
        assert "200 bp" in caplog.text

    def test_dump_with_filtered_pairs(self, caplog):
        """dump_debug_info with filtered pairs should log primer info."""
        caplog.set_level(logging.DEBUG, logger="brimer_plast.pipeline")
        pair = PrimerPair(
            forward_seq="ATCGATCGAT",
            reverse_seq="CGATCGATCG",
            forward_tm=60.0,
            reverse_tm=60.5,
            forward_gc=50.0,
            reverse_gc=50.0,
            product_size=150,
            pair_penalty=1.2,
            pair_number=1,
            primer3_forward_fragments=[
                GenomicFragment("chrI", 100, 109, "+"),
            ],
            primer3_reverse_fragments=[
                GenomicFragment("chrI", 300, 309, "+"),
            ],
        )
        result = PipelineResult(filtered_pairs=[pair])
        dump_debug_info(result)
        assert "Pair 1" in caplog.text
        assert "ATCGATCGAT" in caplog.text

    def test_dump_with_custom_logger(self, caplog):
        """A custom logger should be used when passed."""
        caplog.set_level(logging.DEBUG)
        logger = logging.getLogger("custom_dump_test")
        result = PipelineResult()
        dump_debug_info(result, log=logger)
        assert "DEBUG DUMP" in caplog.text

class TestMakePairName:
    """Pair naming: amplicon coordinates must match primer3 product size."""

    def test_name_amplicon_spans_matching_primer3_product_size(self):
        """PRIMER_RIGHT is the 5' base of the reverse primer, so the last
        amplicon base is reverse_start + 1 (1-based), not reverse_start + len."""
        # product_size == reverse_start - forward_start + 1 == 130
        name = make_pair_name("95388", 0, 645, 774)
        assert name == "95388:646-775"

    def test_name_honours_transcript_offset(self):
        name = make_pair_name("95388", 100, 645, 774)
        assert name == "95388:746-875"

    def test_name_matches_product_size_arithmetic(self):
        # amplicon length (end - start + 1) must equal reverse_start - forward_start + 1
        for fwd_start, rev_start in [(359, 545), (106, 274), (80, 179)]:
            name = make_pair_name("T", 0, fwd_start, rev_start)
            start, end = map(int, name.split(":")[1].split("-"))
            assert end - start + 1 == rev_start - fwd_start + 1
