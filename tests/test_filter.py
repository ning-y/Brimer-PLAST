"""Tests for filter.py: tnBLAST integration."""

import pytest

from brimer_plast.filter import (
    _parse_tnblast_output,
    filter_specific_pairs,
    filter_specific_pairs_single,
    run_tnblast,
    write_assay_file,
)
from brimer_plast.models import PrimerPair


class TestWriteAssayFile:
    def test_writes_tab_delimited(self, tmp_path):
        """Assay file should have three tab-separated columns."""
        pairs = [
            PrimerPair(forward_seq="ATCG", reverse_seq="CGAT"),
            PrimerPair(forward_seq="GCTA", reverse_seq="TAGC"),
        ]
        path = tmp_path / "assays.txt"
        write_assay_file(pairs, path)

        lines = path.read_text().strip().splitlines()
        assert len(lines) == 2
        for i, line in enumerate(lines):
            parts = line.split("\t")
            assert len(parts) == 3
            assert parts[0] == f"pair_{i + 1}"
            assert parts[1] == pairs[i].forward_seq
            assert parts[2] == pairs[i].reverse_seq

    def test_empty_pairs_writes_empty_file(self, tmp_path):
        """An empty list should produce an empty file."""
        path = tmp_path / "empty.txt"
        write_assay_file([], path)
        assert path.read_text() == ""


class TestParseTnblastOutput:
    def test_parses_single_hit(self, tmp_path):
        """A single amplicon should count as 1."""
        text = """name = pair_1
forward primer = ...
"""
        path = tmp_path / "tnt_single.txt"
        path.write_text(text)
        result = _parse_tnblast_output(str(path))
        assert result == {"pair_1": 1}

    def test_parses_multiple_hits(self, tmp_path):
        """Multiple amplicons for same assay should count >1."""
        text = """name = pair_1
forward primer = ...

name = pair_1
forward primer = ...

name = pair_2
forward primer = ...
"""
        path = tmp_path / "tnt_multi.txt"
        path.write_text(text)
        result = _parse_tnblast_output(str(path))
        assert result == {"pair_1": 2, "pair_2": 1}

    def test_empty_output(self, tmp_path):
        """Empty output should return empty dict."""
        path = tmp_path / "tnt_empty.txt"
        path.write_text("")
        result = _parse_tnblast_output(str(path))
        assert result == {}


class TestRunTnblast:
    def test_integration_on_mini_genome(self, mini_genome_fasta):
        """run_tnblast should find hits for our synthetic fixture."""
        import tempfile

        pairs = [
            PrimerPair(
                forward_seq="GCTAGCATCGCTACGTACGT",
                reverse_seq="TGCTAGCTACGTACGATCGC",
            ),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            assay_path = f"{tmp}/assays.txt"
            write_assay_file(pairs, assay_path)
            result = run_tnblast(
                assay_path,
                mini_genome_fasta,
                max_amplicon=1000,
                min_tm=55.0,
                max_tm=65.0,
            )
        # These primers are from the fixture template, so they should hit
        assert len(result) == 1
        assert result["pair_1"] >= 1

    def test_junk_primers_no_hits(self, mini_genome_fasta):
        """Completely unrelated primers should produce no hits."""
        import tempfile

        pairs = [
            PrimerPair(
                forward_seq="AAAAAAAAAAAAAAAAAAAA",
                reverse_seq="TTTTTTTTTTTTTTTTTTTT",
            ),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            assay_path = f"{tmp}/assays.txt"
            write_assay_file(pairs, assay_path)
            result = run_tnblast(
                assay_path,
                mini_genome_fasta,
                max_amplicon=200,
                min_tm=20.0,
                max_tm=40.0,
            )
        assert result == {}

    def test_tnblast_not_on_path(self, mini_genome_fasta, monkeypatch, tmp_path):
        """Missing tnBLAST should raise FileNotFoundError with a helpful message."""
        import os
        import tempfile

        pairs = [
            PrimerPair(forward_seq="ATCG", reverse_seq="CGAT"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            assay_path = os.path.join(tmp, "assays.txt")
            write_assay_file(pairs, assay_path)
            # Set PATH to an empty directory so tntblast is not found
            empty_dir = tmp_path / "empty_bin"
            empty_dir.mkdir()
            monkeypatch.setenv("PATH", str(empty_dir))
            with pytest.raises(FileNotFoundError, match="tntblast.*not found"):
                run_tnblast(assay_path, mini_genome_fasta)


class TestFilterSpecificPairs:
    """Tests for the dual-count filter (genome + transcriptome)."""

    def test_junction_mode_gc0_tc1_passes(self):
        """Junction mode: gc=0, all transcriptome hits on target should pass."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1")]
        result = filter_specific_pairs(
            pairs, {"pair_1": 0}, {"pair_1": ["target_gene"]},
            target_gene="target_gene",
        )
        assert len(result) == 1
        assert result[0] == pairs[0]

    def test_junction_mode_gc1_tc1_fails(self):
        """Junction mode: gc=1 should fail (expected gc=0)."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1")]
        result = filter_specific_pairs(
            pairs, {"pair_1": 1}, {"pair_1": ["target_gene"]},
            target_gene="target_gene",
        )
        assert len(result) == 0

    def test_junction_mode_gc0_tc0_fails(self):
        """Junction mode: no hits should fail."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1")]
        result = filter_specific_pairs(
            pairs, {"pair_1": 0}, {},
            target_gene="target_gene", junction_mode=True,
        )
        assert len(result) == 0

    def test_junction_mode_off_target_gene_fails(self):
        """Transcriptome hits on a different gene should fail."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1")]
        result = filter_specific_pairs(
            pairs, {"pair_1": 0}, {"pair_1": ["other_gene"]},
            target_gene="target_gene", junction_mode=True,
        )
        assert len(result) == 0

    def test_junction_mode_isoforms_pass(self):
        """Multiple transcriptome hits on the same target gene should pass."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1")]
        result = filter_specific_pairs(
            pairs, {"pair_1": 0}, {"pair_1": ["IL1B", "IL1B"]},
            target_gene="IL1B", junction_mode=True,
        )
        assert len(result) == 1

    def test_non_junction_mode_gc1_tc1_passes(self):
        """Non-junction mode: gc=1, tc=1 should pass."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1")]
        result = filter_specific_pairs(
            pairs,
            {"pair_1": 1},
            {"pair_1": ["target_gene"]},
            target_gene="target_gene",
            junction_mode=False,
        )
        assert len(result) == 1
        assert result[0] == pairs[0]

    def test_non_junction_mode_gc0_tc1_fails(self):
        """Non-junction mode: gc=0 should fail (expected gc=1)."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1")]
        result = filter_specific_pairs(
            pairs,
            {"pair_1": 0},
            {"pair_1": ["target_gene"]},
            target_gene="target_gene",
            junction_mode=False,
        )
        assert len(result) == 0

    def test_junction_mode_mixed_pairs(self):
        """Multiple pairs with mixed results."""
        pairs = [
            PrimerPair(forward_seq="F1", reverse_seq="R1"),  # passes
            PrimerPair(forward_seq="F2", reverse_seq="R2"),  # fails (gc=1)
            PrimerPair(forward_seq="F3", reverse_seq="R3"),  # fails (no tc hits)
            PrimerPair(forward_seq="F4", reverse_seq="R4"),  # passes
        ]
        result = filter_specific_pairs(
            pairs,
            {"pair_1": 0, "pair_2": 1, "pair_3": 0, "pair_4": 0},
            {"pair_1": ["target_gene"], "pair_4": ["target_gene", "target_gene"]},
            target_gene="target_gene",
            junction_mode=True,
        )
        assert len(result) == 2
        assert result[0] == pairs[0]
        assert result[1] == pairs[3]

    def test_empty_inputs(self):
        """Empty primer list and empty counts produce empty result."""
        assert filter_specific_pairs([], {}, {}, target_gene="x") == []


class TestFilterSpecificPairsSingle:
    """Legacy single-database tests for filter_specific_pairs_single."""

    def test_filters_out_off_target(self):
        """Pairs with amplicon count != 1 should be dropped."""
        pairs = [
            PrimerPair(forward_seq="F1", reverse_seq="R1"),
            PrimerPair(forward_seq="F2", reverse_seq="R2"),
            PrimerPair(forward_seq="F3", reverse_seq="R3"),
        ]
        counts = {"pair_1": 1, "pair_2": 3, "pair_3": 0}
        result = filter_specific_pairs_single(pairs, counts)
        assert len(result) == 1
        assert result[0] == pairs[0]

    def test_all_specific(self):
        """All pairs with exactly 1 hit should be kept."""
        pairs = [
            PrimerPair(forward_seq="F1", reverse_seq="R1"),
            PrimerPair(forward_seq="F2", reverse_seq="R2"),
        ]
        counts = {"pair_1": 1, "pair_2": 1}
        result = filter_specific_pairs_single(pairs, counts)
        assert result == pairs

    def test_none_specific(self):
        """No pairs with count 1 should produce empty list."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1")]
        counts = {"pair_1": 0}
        result = filter_specific_pairs_single(pairs, counts)
        assert result == []

    def test_empty_inputs(self):
        """Empty primer list and empty counts produce empty result."""
        assert filter_specific_pairs_single([], {}) == []
