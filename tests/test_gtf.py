"""Tests for gtf.py: GTF parsing and attribute extraction."""

import pytest

from brimer_plast.gtf import (
    _parse_gtf_attributes,
    _find_gene_for_transcript,
    build_transcript_to_gene_map,
    parse_gtf,
    parse_gtf_all_transcripts,
    parse_gtf_grouped_by_transcript,
)


class TestParseGtfAttributes:
    """Tests for the core GTF attribute parser _parse_gtf_attributes."""

    def test_parses_gene_name_and_transcript_id(self):
        """Standard attributes with gene_name and transcript_id."""
        attrs = _parse_gtf_attributes(
            'gene_id "ENSG000001"; transcript_id "NM_001256799.1"; gene_name "GAPDH";'
        )
        assert attrs["gene_id"] == "ENSG000001"
        assert attrs["transcript_id"] == "NM_001256799.1"
        assert attrs["gene_name"] == "GAPDH"

    def test_parses_gene_only(self):
        """Attributes with just gene_id."""
        attrs = _parse_gtf_attributes('gene_id "WBGene00022277";')
        assert attrs["gene_id"] == "WBGene00022277"

    def test_handles_extra_whitespace(self):
        """Attributes with extra spaces around semicolons."""
        attrs = _parse_gtf_attributes(
            '  gene_id "ENSG001"  ;  transcript_id "NM_001" ;  '
        )
        assert attrs["gene_id"] == "ENSG001"
        assert attrs["transcript_id"] == "NM_001"

    def test_handles_attributes_with_underscores(self):
        """Attributes with underscore in the key name."""
        attrs = _parse_gtf_attributes(
            'gene_name "ACTB"; transcript_id "NM_001101.5";'
        )
        assert attrs["gene_name"] == "ACTB"
        assert attrs["transcript_id"] == "NM_001101.5"

    def test_handles_values_with_special_chars(self):
        """Values may contain dots, dashes, etc."""
        attrs = _parse_gtf_attributes(
            'gene_id "ENSG000001.1"; transcript_id "ENST000004.2";'
        )
        assert attrs["gene_id"] == "ENSG000001.1"
        assert attrs["transcript_id"] == "ENST000004.2"

    def test_returns_empty_dict_for_empty_string(self):
        """Empty attribute string should return empty dict."""
        assert _parse_gtf_attributes("") == {}

    def test_handles_multiple_gene_attributes(self):
        """GTF may have both gene_id and gene_name."""
        attrs = _parse_gtf_attributes(
            'gene_id "ENSG000001"; gene_name "GAPDH"; gene "GAPDH_sym";'
        )
        assert attrs["gene_id"] == "ENSG000001"
        assert attrs["gene_name"] == "GAPDH"
        assert attrs["gene"] == "GAPDH_sym"

    def test_only_matches_quoted_values(self):
        """Values without quotes should not be matched."""
        attrs = _parse_gtf_attributes('gene_id ENSG000001; transcript_id "NM_001";')
        assert "gene_id" not in attrs  # unquoted value not matched
        assert attrs["transcript_id"] == "NM_001"


class TestFindGeneForTranscript:
    """Tests for _find_gene_for_transcript."""

    def test_returns_gene_name_for_known_transcript(self, mini_genome_gtf):
        """test_transcript belongs to test_gene."""
        gene = _find_gene_for_transcript(mini_genome_gtf, "test_transcript")
        assert gene == "test_gene"

    def test_returns_none_for_missing_transcript(self, mini_genome_gtf):
        """Non-existent transcript should return None."""
        gene = _find_gene_for_transcript(mini_genome_gtf, "nonexistent")
        assert gene is None

    def test_finds_gene_in_multi_transcript_gtf(self, multi_genome_gtf):
        """transcript_A belongs to multi_test_gene."""
        gene = _find_gene_for_transcript(multi_genome_gtf, "transcript_A")
        assert gene == "multi_test_gene"

    def test_finds_gene_in_single_exon_gtf(self, single_exon_gtf):
        """single_exon_transcript belongs to single_exon_gene."""
        gene = _find_gene_for_transcript(single_exon_gtf, "single_exon_transcript")
        assert gene == "single_exon_gene"


class TestBuildTranscriptToGeneMap:
    """Tests for build_transcript_to_gene_map."""

    def test_returns_correct_mapping(self, multi_genome_gtf):
        """All three transcripts should map to multi_test_gene."""
        mapping = build_transcript_to_gene_map(multi_genome_gtf)
        assert mapping["transcript_A"] == "multi_test_gene"
        assert mapping["transcript_B"] == "multi_test_gene"
        assert mapping["transcript_C"] == "multi_test_gene"

    def test_includes_single_exon_transcripts(self, single_exon_gtf):
        """Single-exon transcripts should also appear."""
        mapping = build_transcript_to_gene_map(single_exon_gtf)
        assert mapping["single_exon_transcript"] == "single_exon_gene"

    def test_returns_empty_dict_for_empty_gtf(self, tmp_path):
        """A GTF with no exon lines should return empty dict."""
        empty_gtf = tmp_path / "empty.gtf"
        empty_gtf.write_text("##gff-version 3\n")
        assert build_transcript_to_gene_map(str(empty_gtf)) == {}
