"""Tests for filter.py: tnBLAST integration."""

import pytest

from brimer_plast.filter import filter_specific_pairs
from brimer_plast.tnblast import (
    AmpliconHit,
    _parse_tnblast_amplicons,
    _parse_tnblast_output,
    run_tnblast,
    write_assay_file,
)
from brimer_plast.models import PrimerPair, GeneLocus


def _hit_list(count: int, seqid: str = "chrI", start: int = 100, end: int = 200):
    """Helper to create a list of mock AmpliconHits."""
    return [AmpliconHit(seqid=seqid, amplicon_start=start, amplicon_end=end) for _ in range(count)]


class TestWriteAssayFile:
    def test_writes_tab_delimited(self, tmp_path):
        """Assay file should have three tab-separated columns."""
        pairs = [
            PrimerPair(forward_seq="ATCG", reverse_seq="CGAT", pair_name="test1"),
            PrimerPair(forward_seq="GCTA", reverse_seq="TAGC", pair_name="test2"),
        ]
        path = tmp_path / "assays.txt"
        write_assay_file(pairs, path)

        lines = path.read_text().strip().splitlines()
        assert len(lines) == 2
        assert lines[0] == f"test1\t{pairs[0].forward_seq}\t{pairs[0].reverse_seq}"
        assert lines[1] == f"test2\t{pairs[1].forward_seq}\t{pairs[1].reverse_seq}"

    def test_empty_pairs_writes_empty_file(self, tmp_path):
        """An empty list should produce an empty file."""
        path = tmp_path / "empty.txt"
        write_assay_file([], path)
        assert path.read_text() == ""


class TestParseTnblastOutput:
    def test_parses_single_hit(self, tmp_path):
        """A single amplicon should count as 1."""
        text = """name = 9746.1:45-199
forward primer = ...
"""
        path = tmp_path / "tnt_single.txt"
        path.write_text(text)
        result = _parse_tnblast_output(str(path))
        assert result == {"9746.1:45-199": 1}

    def test_parses_multiple_hits(self, tmp_path):
        """Multiple amplicons for same assay should count >1."""
        text = """name = 9746.1:45-199
forward primer = ...

name = 9746.1:45-199
forward primer = ...

name = 4567.1:331-470
forward primer = ...
"""
        path = tmp_path / "tnt_multi.txt"
        path.write_text(text)
        result = _parse_tnblast_output(str(path))
        assert result == {"9746.1:45-199": 2, "4567.1:331-470": 1}

    def test_empty_output(self, tmp_path):
        """Empty output should return empty dict."""
        path = tmp_path / "tnt_empty.txt"
        path.write_text("")
        result = _parse_tnblast_output(str(path))
        assert result == {}


class TestParseTnblastAmplicons:
    """Tests for _parse_tnblast_amplicons."""

    def test_parses_single_amplicon_with_coords(self, tmp_path):
        """A single amplicon with seqid and coordinate range."""
        text = """name = 9746.1:45-199
amplicon range = 100 .. 300
forward primer = ...

>chrI
"""
        path = tmp_path / "tnt.txt"
        path.write_text(text)
        result = _parse_tnblast_amplicons(str(path))
        assert "9746.1:45-199" in result
        assert len(result["9746.1:45-199"]) == 1
        hit = result["9746.1:45-199"][0]
        assert hit.seqid == "chrI"
        assert hit.amplicon_start == 100
        assert hit.amplicon_end == 300

    def test_multiple_amplicons_same_assay(self, tmp_path):
        """Multiple amplicons for the same assay in different genomic regions."""
        text = """name = 9746.1:45-199
amplicon range = 100 .. 200

>chrI
name = 9746.1:45-199
amplicon range = 1000 .. 1200

>chrII
name = 4567.1:331-470
amplicon range = 500 .. 600

>chrI
"""
        path = tmp_path / "tnt_multi.txt"
        path.write_text(text)
        result = _parse_tnblast_amplicons(str(path))
        assert len(result["9746.1:45-199"]) == 2
        assert len(result["4567.1:331-470"]) == 1

    def test_empty_output(self, tmp_path):
        """Empty file yields empty dict."""
        path = tmp_path / "empty.txt"
        path.write_text("")
        assert _parse_tnblast_amplicons(str(path)) == {}

    def test_amplicon_range_without_seqid(self, tmp_path):
        """Output with amplicon range but no seqid should not emit."""
        text = """name = 9746.1:45-199
amplicon range = 100 .. 200
forward primer = ...
"""
        path = tmp_path / "tnt_no_seqid.txt"
        path.write_text(text)
        result = _parse_tnblast_amplicons(str(path))
        assert result == {}

    def test_multiple_assays_with_results(self, tmp_path):
        """Multiple assays each with their own amplicons."""
        text = """name = 9746.1:45-199
amplicon range = 100 .. 200

>chrI
name = 9746.1:45-199
amplicon range = 300 .. 400

>chrI
name = 4567.1:331-470
amplicon range = 500 .. 600

>chrV
"""
        path = tmp_path / "tnt_multi_assay.txt"
        path.write_text(text)
        result = _parse_tnblast_amplicons(str(path))
        assert len(result["9746.1:45-199"]) == 2
        assert len(result["4567.1:331-470"]) == 1
        assert result["4567.1:331-470"][0].seqid == "chrV"


class TestRunTnblast:
    def test_integration_on_mini_genome(self, mini_genome_fasta):
        """run_tnblast should find hits for our synthetic fixture."""
        import tempfile

        pairs = [
            PrimerPair(
                forward_seq="GCTAGCATCGCTACGTACGT",
                reverse_seq="TGCTAGCTACGTACGATCGC",
                pair_name="test_1",
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
        assert result["test_1"] >= 1

    def test_junk_primers_no_hits(self, mini_genome_fasta):
        """Completely unrelated primers should produce no hits."""
        import tempfile

        pairs = [
            PrimerPair(
                forward_seq="AAAAAAAAAAAAAAAAAAAA",
                reverse_seq="TTTTTTTTTTTTTTTTTTTT",
                pair_name="junk_1",
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
            PrimerPair(forward_seq="ATCG", reverse_seq="CGAT", pair_name="test_1"),
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
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1")]
        result = filter_specific_pairs(
            pairs, {"p1": []}, {"p1": ["target_gene"]},
            target_gene="target_gene",
        )
        assert len(result) == 1
        assert result[0] == pairs[0]

    def test_junction_mode_gc1_tc1_fails(self):
        """Junction mode: gc=1 should fail (expected gc=0)."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1")]
        result = filter_specific_pairs(
            pairs, {"p1": _hit_list(1)}, {"p1": ["target_gene"]},
            target_gene="target_gene",
        )
        assert len(result) == 0

    def test_junction_mode_gc0_tc0_fails(self):
        """Junction mode: no hits should fail."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1")]
        result = filter_specific_pairs(
            pairs, {"p1": []}, {},
            target_gene="target_gene", junction_mode=True,
        )
        assert len(result) == 0

    def test_junction_mode_off_target_gene_fails(self):
        """Transcriptome hits on a different gene should fail."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1")]
        result = filter_specific_pairs(
            pairs, {"p1": []}, {"p1": ["other_gene"]},
            target_gene="target_gene", junction_mode=True,
        )
        assert len(result) == 0

    def test_junction_mode_isoforms_pass(self):
        """Multiple transcriptome hits on the same target gene should pass."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1")]
        result = filter_specific_pairs(
            pairs, {"p1": []}, {"p1": ["IL1B", "IL1B"]},
            target_gene="IL1B", junction_mode=True,
        )
        assert len(result) == 1

    def test_non_junction_mode_gc1_tc1_passes(self):
        """Non-junction mode: gc=1, tc=1 should pass."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1")]
        result = filter_specific_pairs(
            pairs,
            {"p1": _hit_list(1)},
            {"p1": ["target_gene"]},
            target_gene="target_gene",
            junction_mode=False,
        )
        assert len(result) == 1
        assert result[0] == pairs[0]

    def test_non_junction_mode_gc0_tc1_passes(self):
        """Non-junction mode: gc=0 should PASS (it means the intron is very large)."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1")]
        result = filter_specific_pairs(
            pairs,
            {"p1": []},
            {"p1": ["target_gene"]},
            target_gene="target_gene",
            junction_mode=False,
        )
        assert len(result) == 1
        assert result[0] == pairs[0]

    def test_non_junction_mode_off_locus_fails(self):
        """Non-junction mode: a hit on a different chromosome should fail."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1")]
        locus = GeneLocus("target_gene", "chrI", "+", {}, 1000, 2000)
        
        # Hit on chrII
        hits = [AmpliconHit(seqid="chrII", amplicon_start=1500, amplicon_end=1600)]
        result = filter_specific_pairs(
            pairs,
            {"p1": hits},
            {"p1": ["target_gene"]},
            target_gene="target_gene",
            target_locus=locus,
            junction_mode=False,
        )
        assert len(result) == 0

    def test_non_junction_mode_on_locus_passes(self):
        """Non-junction mode: a hit on the target chromosome and range should pass."""
        pairs = [PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1")]
        locus = GeneLocus("target_gene", "chrI", "+", {}, 1000, 2000)
        
        # Hit inside locus (with 1000bp padding)
        hits = [AmpliconHit(seqid="chrI", amplicon_start=1500, amplicon_end=1600)]
        result = filter_specific_pairs(
            pairs,
            {"p1": hits},
            {"p1": ["target_gene"]},
            target_gene="target_gene",
            target_locus=locus,
            junction_mode=False,
        )
        assert len(result) == 1

    def test_junction_mode_mixed_pairs(self):
        """Multiple pairs with mixed results."""
        pairs = [
            PrimerPair(forward_seq="F1", reverse_seq="R1", pair_name="p1"),   # passes
            PrimerPair(forward_seq="F2", reverse_seq="R2", pair_name="p2"),   # fails (gc=1)
            PrimerPair(forward_seq="F3", reverse_seq="R3", pair_name="p3"),   # fails (no tc hits)
            PrimerPair(forward_seq="F4", reverse_seq="R4", pair_name="p4"),   # passes
        ]
        result = filter_specific_pairs(
            pairs,
            {
                "p1": [],
                "p2": _hit_list(1),
                "p3": [],
                "p4": [],
            },
            {"p1": ["target_gene"], "p4": ["target_gene", "target_gene"]},
            target_gene="target_gene",
            junction_mode=True,
        )
        assert len(result) == 2
        assert result[0] == pairs[0]
        assert result[1] == pairs[3]

    def test_empty_inputs(self):
        """Empty primer list and empty counts produce empty result."""
        assert filter_specific_pairs([], {}, {}, target_gene="x") == []
