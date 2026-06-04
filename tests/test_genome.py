"""Tests for genome.py: GTF parsing and sequence extraction."""

import pytest

from brimer_plast.genome import (
    compute_conserved_exon_chains,
    extract_sequence,
    get_target_information,
    parse_gtf,
    parse_gtf_grouped_by_transcript,
)
from brimer_plast.models import ExonInfo


# ── helpers ─────────────────────────────────────────────────────────────────


def _load_genome(fasta_path: str) -> str:
    with open(fasta_path) as f:
        lines = f.readlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


# ── parse_gtf (flat, backward-compat) ────────────────────────────────────────


class TestParseGTF:
    def test_returns_exons_for_gene_name(self, mini_genome_gtf):
        """parse_gtf should return exon coordinates for matching gene_name."""
        exons = parse_gtf(mini_genome_gtf, target_gene="test_gene")
        assert len(exons) == 2
        for exon in exons:
            assert exon.seqid == "chrI"
            assert exon.strand == "+"

    def test_returns_exons_for_transcript(self, mini_genome_gtf):
        """parse_gtf should return exon coordinates for matching transcript_id."""
        exons = parse_gtf(mini_genome_gtf, target_transcript="test_transcript")
        assert len(exons) == 2

    def test_raises_if_no_target_provided(self, mini_genome_gtf):
        """Must provide target_gene or target_transcript."""
        with pytest.raises(ValueError, match="either.*target.*gene.*transcript"):
            parse_gtf(mini_genome_gtf)

    def test_raises_if_both_targets_provided(self, mini_genome_gtf):
        """Cannot provide both target_gene and target_transcript."""
        with pytest.raises(ValueError, match="both"):
            parse_gtf(mini_genome_gtf, target_gene="x", target_transcript="y")

    def test_raises_if_gene_not_found(self, mini_genome_gtf):
        """Raise if target_gene does not appear in the GTF."""
        with pytest.raises(ValueError, match="not found"):
            parse_gtf(mini_genome_gtf, target_gene="nonexistent_gene")

    def test_raises_if_transcript_not_found(self, mini_genome_gtf):
        """Raise if target_transcript does not appear in the GTF."""
        with pytest.raises(ValueError, match="not found"):
            parse_gtf(mini_genome_gtf, target_transcript="nonexistent_transcript")

    def test_exon_coordinates_are_ordered(self, mini_genome_gtf):
        """Exons should be returned in genomic order (by start coordinate)."""
        exons = parse_gtf(mini_genome_gtf, target_gene="test_gene")
        starts = [e.start for e in exons]
        assert starts == sorted(starts)


# ── parse_gtf_grouped_by_transcript ──────────────────────────────────────────


class TestParseGTFGroupedByTranscript:
    def test_returns_dict_keyed_by_transcript_id(self, mini_genome_gtf):
        """Should return a dict with transcript_id keys and exon list values."""
        result = parse_gtf_grouped_by_transcript(
            mini_genome_gtf, target_gene="test_gene"
        )
        assert "test_transcript" in result
        assert len(result) == 1  # single transcript in mini_genome

    def test_each_entry_has_sorted_exons(self, mini_genome_gtf):
        """Exons for each transcript should be sorted by genomic position."""
        result = parse_gtf_grouped_by_transcript(
            mini_genome_gtf, target_gene="test_gene"
        )
        for tid, exons in result.items():
            starts = [e.start for e in exons]
            assert starts == sorted(starts), f"{tid} exons not sorted"

    def test_works_by_transcript_id_directly(self, mini_genome_gtf):
        """Filtering by target_transcript should return a single-entry dict."""
        result = parse_gtf_grouped_by_transcript(
            mini_genome_gtf, target_transcript="test_transcript"
        )
        assert "test_transcript" in result
        assert len(result["test_transcript"]) == 2

    def test_raises_if_no_target(self, mini_genome_gtf):
        with pytest.raises(ValueError, match="either.*target.*gene.*transcript"):
            parse_gtf_grouped_by_transcript(mini_genome_gtf)

    def test_raises_if_both_targets(self, mini_genome_gtf):
        with pytest.raises(ValueError, match="both"):
            parse_gtf_grouped_by_transcript(
                mini_genome_gtf, target_gene="x", target_transcript="y"
            )

    def test_raises_if_not_found(self, mini_genome_gtf):
        with pytest.raises(ValueError, match="not found"):
            parse_gtf_grouped_by_transcript(
                mini_genome_gtf, target_gene="nonexistent"
            )

    def test_multi_transcript_grouping(self, multi_genome_gtf):
        """multi_test_gene should have 3 transcripts."""
        result = parse_gtf_grouped_by_transcript(
            multi_genome_gtf, target_gene="multi_test_gene"
        )
        assert set(result.keys()) == {"transcript_A", "transcript_B", "transcript_C"}
        # A: 5 exons, B: 4 exons (skips exon4), C: 5 exons
        assert len(result["transcript_A"]) == 5
        assert len(result["transcript_B"]) == 4
        assert len(result["transcript_C"]) == 5


# ── parse_gtf_all_transcripts ────────────────────────────────────────────────


class TestParseGtfAllTranscripts:
    def test_returns_all_transcripts(self, multi_genome_gtf):
        """parse_gtf_all_transcripts should return ALL transcripts in the GTF."""
        from brimer_plast.genome import parse_gtf_all_transcripts

        result = parse_gtf_all_transcripts(multi_genome_gtf)
        assert "transcript_A" in result
        assert "transcript_B" in result
        assert "transcript_C" in result

    def test_does_not_raise_without_target(self, mini_genome_gtf):
        """parse_gtf_all_transcripts works without any target filter."""
        from brimer_plast.genome import parse_gtf_all_transcripts

        result = parse_gtf_all_transcripts(mini_genome_gtf)
        assert "test_transcript" in result

    def test_single_exon_transcripts_included(self, single_exon_gtf):
        """Single-exon transcripts should be included."""
        from brimer_plast.genome import parse_gtf_all_transcripts

        result = parse_gtf_all_transcripts(single_exon_gtf)
        assert "single_exon_transcript" in result
        assert len(result["single_exon_transcript"]) == 1


# ── build_transcriptome_fasta ────────────────────────────────────────────────


class TestBuildTranscriptomeFasta:
    def test_writes_fasta_with_spliced_sequences(self, mini_genome_fasta, mini_genome_gtf, tmp_path):
        """build_transcriptome_fasta should write a FASTA with spliced transcripts."""
        from brimer_plast.genome import build_transcriptome_fasta

        out_path = tmp_path / "transcriptome.fa"
        build_transcriptome_fasta(mini_genome_fasta, mini_genome_gtf, out_path)

        text = out_path.read_text()
        assert ">test_transcript" in text
        # The mini genome has one transcript with 2 exons (101-350, 451-600)
        # Spliced length = 250 + 150 = 400 bp
        lines = text.strip().splitlines()
        assert lines[0] == ">test_transcript"
        seq = "".join(lines[1:])
        assert len(seq) == 400

    def test_multi_transcript_fasta(self, multi_genome_fasta, multi_genome_gtf, tmp_path):
        """Should write all transcripts from the GTF."""
        from brimer_plast.genome import build_transcriptome_fasta

        out_path = tmp_path / "multi_transcriptome.fa"
        build_transcriptome_fasta(multi_genome_fasta, multi_genome_gtf, out_path)

        text = out_path.read_text()
        for tid in ("transcript_A", "transcript_B", "transcript_C"):
            assert f">{tid}" in text, f"Missing transcript {tid}"

        # Each transcript should have correct spliced length
        # transcript_A: 5 exons (150+150+150+100+100 = 650 bp)
        # transcript_B: 4 exons (150+150+100+100 = 500 bp)
        # transcript_C: 5 exons (150+150+150+100+100 = 650 bp)
        entries = text.strip().split(">")
        for entry in entries:
            if not entry:
                continue
            tid, _, seq_raw = entry.partition("\n")
            seq = seq_raw.replace("\n", "")
            expected = {"transcript_A": 650, "transcript_B": 500, "transcript_C": 650}
            assert len(seq) == expected[tid], (
                f"{tid} expected {expected[tid]} bp, got {len(seq)}"
            )

    def test_sequence_is_line_wrapped(self, mini_genome_fasta, mini_genome_gtf, tmp_path):
        """FASTA lines should be no longer than 60 characters."""
        from brimer_plast.genome import build_transcriptome_fasta

        out_path = tmp_path / "wrapped.fa"
        build_transcriptome_fasta(mini_genome_fasta, mini_genome_gtf, out_path)

        for line in out_path.read_text().splitlines():
            if line.startswith(">"):
                continue
            assert len(line) <= 60, f"Line too long ({len(line)} chars): {line}"


# ── compute_conserved_exon_chains ────────────────────────────────────────────


class TestComputeConservedExonChains:
    def test_single_transcript_returns_one_chain(self, mini_genome_gtf):
        """A single-transcript gene should produce one chain with all exons."""
        transcripts = parse_gtf_grouped_by_transcript(
            mini_genome_gtf, target_gene="test_gene"
        )
        chains = compute_conserved_exon_chains(list(transcripts.values()))
        assert len(chains) == 1
        assert len(chains[0].exons) == 2

    def test_two_conserved_chains(self, multi_genome_gtf):
        """multi_test_gene: B skips exon3, so conserved chains are [1,2] and [4,5]."""
        transcripts = parse_gtf_grouped_by_transcript(
            multi_genome_gtf, target_gene="multi_test_gene"
        )
        chains = compute_conserved_exon_chains(list(transcripts.values()))
        assert len(chains) == 2, f"Expected 2 conserved chains, got {len(chains)}"

        # Chain 1: exons 1,2 (starts 101, 301)
        chain1 = next(c for c in chains if c.exons[0].start == 101)
        assert len(chain1.exons) == 2
        assert chain1.exons[1].start == 301

        # Chain 2: exons 4,5 (starts 701, 901)
        chain2 = next(c for c in chains if c.exons[0].start == 701)
        assert len(chain2.exons) == 2
        assert chain2.exons[1].start == 901

    def test_raises_if_no_conserved_chains(self):
        """If transcripts share no exon adjacency, raise ValueError."""
        # Two transcripts with non-overlapping exon structures
        transcripts = {
            "t1": [
                ExonInfo("chrI", 1, 100, "+"),
                ExonInfo("chrI", 200, 300, "+"),
            ],
            "t2": [
                ExonInfo("chrI", 400, 500, "+"),
                ExonInfo("chrI", 600, 700, "+"),
            ],
        }
        with pytest.raises(
            ValueError, match="(?i)no conserved.*junction"
        ):
            compute_conserved_exon_chains(list(transcripts.values()))

    def test_raises_if_all_single_exon_transcripts(self):
        """A single-exon gene has no junctions at all."""
        transcripts = {
            "t1": [ExonInfo("chrI", 1, 100, "+")],
            "t2": [ExonInfo("chrI", 1, 100, "+")],
        }
        with pytest.raises(ValueError, match="(?i)at least 2 exons"):
            compute_conserved_exon_chains(list(transcripts.values()))


# ── extract_sequence ─────────────────────────────────────────────────────────


class TestExtractSequence:
    def test_concatenates_exons(self, mini_genome_fasta, mini_genome_gtf):
        """extract_sequence should concatenate exon sequences in order."""
        exons = parse_gtf(mini_genome_gtf, target_gene="test_gene")
        seq = extract_sequence(mini_genome_fasta, exons)
        # exon1: 101..350 (250 bp), exon2: 451..600 (150 bp) -> 400 bp total
        assert len(seq) == 400

    def test_empty_exons_returns_empty_string(self):
        """extract_sequence with empty exon list should return empty string."""
        seq = extract_sequence("dummy.fa", [])
        assert seq == ""

    def test_reverse_strand_reverse_complements(self, mini_genome_fasta):
        """Exons on - strand should be reverse-complemented."""
        exons = [ExonInfo(seqid="chrI", start=1, end=10, strand="-")]
        genome = _load_genome(mini_genome_fasta)
        forward = genome[0:10]
        seq = extract_sequence(mini_genome_fasta, exons)
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        expected = "".join(complement.get(b, b) for b in forward[::-1])
        assert seq == expected, (
            f"Expected reverse complement of {forward!r} = {expected!r}, got {seq!r}"
        )

    def test_raises_if_seqid_not_in_fasta(self, mini_genome_fasta):
        """Seqid not in FASTA should raise ValueError."""
        exons = [ExonInfo(seqid="nonexistent_chr", start=1, end=10, strand="+")]
        with pytest.raises(ValueError, match="nonexistent_chr"):
            extract_sequence(mini_genome_fasta, exons)

    def test_raises_if_exons_span_different_chromosomes(self, mini_genome_fasta):
        """Exons on different seqids should raise ValueError."""
        exons = [
            ExonInfo(seqid="chrI", start=1, end=10, strand="+"),
            ExonInfo(seqid="chrII", start=20, end=30, strand="+"),
        ]
        with pytest.raises(ValueError, match="different chromosomes"):
            extract_sequence(mini_genome_fasta, exons)


# ── get_target_information ───────────────────────────────────────────────────


class TestGetTargetInformation:
    def test_with_transcript(self, mini_genome_gtf, mini_genome_fasta):
        """--target-transcript should produce one chain with the transcript exons."""
        result = get_target_information(
            mini_genome_fasta, mini_genome_gtf, target_transcript="test_transcript"
        )
        assert len(result) == 1
        assert "test_transcript" in result[0].id
        assert len(result[0].template) == 400

    def test_with_transcript_has_junction_positions(self, mini_genome_gtf, mini_genome_fasta):
        """--target-transcript should produce 1-based junction positions."""
        result = get_target_information(
            mini_genome_fasta, mini_genome_gtf, target_transcript="test_transcript"
        )
        chain = result[0]
        # exon1=250bp, exon2=150bp, junction at 250 (0-indexed) → 251 (1-based)
        assert chain.junction_positions_1based == [251]

    def test_raises_if_neither(self, mini_genome_gtf, mini_genome_fasta):
        with pytest.raises(ValueError, match="either.*target.*gene.*transcript"):
            get_target_information(mini_genome_fasta, mini_genome_gtf)

    def test_raises_if_both(self, mini_genome_gtf, mini_genome_fasta):
        with pytest.raises(ValueError, match="both"):
            get_target_information(
                mini_genome_fasta, mini_genome_gtf, target_gene="x", target_transcript="y"
            )

    def test_multi_transcript_produces_two_chains(self, multi_genome_gtf, multi_genome_fasta):
        """multi_test_gene produces 2 conserved chains with templates."""
        result = get_target_information(
            multi_genome_fasta, multi_genome_gtf, target_gene="multi_test_gene"
        )
        assert len(result) == 2
        # Chain 1: exons [1,2] = 300bp, Chain 2: exons [4,5] = 200bp
        templates_by_len = {len(c.template): c for c in result}
        assert 300 in templates_by_len, f"Chain 1 template should be 300bp, got lens: {list(templates_by_len.keys())}"
        assert 200 in templates_by_len

    def test_target_transcript_unique_junctions(self, multi_genome_gtf, multi_genome_fasta):
        """transcript_B has unique junction B→D; only that should be required."""
        result = get_target_information(
            multi_genome_fasta, multi_genome_gtf, target_transcript="transcript_B"
        )
        assert len(result) == 1
        chain = result[0]
        # transcript_B exons: [101-250, 301-450, 701-800, 901-1000]
        # All junctions: [151, 301, 401] (1-based)
        assert chain.junction_positions_1based == [151, 301, 401], (
            f"Expected all 3 junctions, got {chain.junction_positions_1based}"
        )
        # Only B→D (301-450→701-800) is unique → position 301
        assert chain.required_junction_positions_1based == [301], (
            f"Expected only unique junction [301], got "
            f"{chain.required_junction_positions_1based}"
        )

    def test_target_transcript_no_unique_junctions(self, multi_genome_gtf, multi_genome_fasta):
        """transcript_A shares all junctions with siblings → no unique junctions."""
        result = get_target_information(
            multi_genome_fasta, multi_genome_gtf, target_transcript="transcript_A"
        )
        assert len(result) == 1
        chain = result[0]
        # All junctions populated, but required_junction_positions_1based is empty
        assert len(chain.junction_positions_1based) > 0
        assert chain.required_junction_positions_1based == [], (
            f"Expected empty required junctions, got {chain.required_junction_positions_1based}"
        )

    def test_target_transcript_single_transcript_gene_all_required(
        self, mini_genome_gtf, mini_genome_fasta
    ):
        """Single-transcript gene (test_gene): all junctions are trivially required."""
        result = get_target_information(
            mini_genome_fasta, mini_genome_gtf, target_transcript="test_transcript"
        )
        chain = result[0]
        # junction_positions = required_junction_positions for single-transcript genes
        assert chain.junction_positions_1based == [251]
        assert chain.required_junction_positions_1based == [251]

    def test_single_exon_transcript_raises_error(self, mini_genome_gtf, mini_genome_fasta):
        """A single-exon target has no junctions; should raise."""
        transcripts = {
            "t1": [ExonInfo("chrI", 1, 100, "+")],
        }
        with pytest.raises(ValueError, match="(?i)at least 2 exons"):
            compute_conserved_exon_chains(list(transcripts.values()))

    def test_single_exon_transcript_returns_junctionless_chain(self, single_exon_gtf, mini_genome_fasta):
        """A single-exon transcript should return a junctionless chain."""
        result = get_target_information(
            mini_genome_fasta, single_exon_gtf,
            target_transcript="single_exon_transcript",
        )
        assert len(result) == 1
        chain = result[0]
        assert len(chain.template) == 100  # 701..800 = 100 bp
        assert chain.junction_positions_1based == []
        assert chain.required_junction_positions_1based == []

    def test_single_exon_gene_returns_junctionless_chains(self, single_exon_gtf, mini_genome_fasta):
        """An all-single-exon gene returns junctionless chains."""
        result = get_target_information(
            mini_genome_fasta, single_exon_gtf,
            target_gene="single_exon_gene",
        )
        assert len(result) == 1
        chain = result[0]
        assert len(chain.template) == 100
        assert chain.junction_positions_1based == []
        assert chain.required_junction_positions_1based == []

    def test_mixed_exon_gene_returns_both_kinds(self, mixed_exon_gtf, mini_genome_fasta):
        """A gene with mixed single/multi-exon transcripts should
        produce conserved chains from multi-exon transcripts plus single-exon chains."""
        result = get_target_information(
            mini_genome_fasta, mixed_exon_gtf,
            target_gene="mixed_gene",
        )
        # mixed_gene has 2 transcripts: transcript_A (2 exons), transcript_B (1 exon)
        # -> 1 conserved chain (101-250, 301-450) + 1 single-exon chain
        assert len(result) == 2
        chain_with_junction = [c for c in result if c.junction_positions_1based]
        chain_without_junction = [c for c in result if not c.junction_positions_1based]
        assert len(chain_with_junction) == 1
        assert len(chain_without_junction) == 1
        assert chain_without_junction[0].required_junction_positions_1based == []

    def test_chain_ids_use_readable_format(self, multi_genome_gtf, multi_genome_fasta):
        """Chain IDs should use chain_N format, not concatenated exon starts."""
        result = get_target_information(
            multi_genome_fasta, multi_genome_gtf, target_gene="multi_test_gene"
        )
        for chain in result:
            assert "chain_" in chain.id, f"Expected chain_N format, got {chain.id}"
        assert result[0].id == "multi_test_gene_chain_1"
        assert result[1].id == "multi_test_gene_chain_2"

    def test_transcript_with_no_unique_junctions(self, multi_genome_gtf, multi_genome_fasta):
        """A transcript that shares all junctions with siblings returns
        a chain with empty required_junction_positions."""
        result = get_target_information(
            multi_genome_fasta, multi_genome_gtf,
            target_transcript="transcript_A",
        )
        assert len(result) == 1
        chain = result[0]
        # transcript_A has all 4 junctions shared with siblings
        # No unique junctions → required_junction_positions_1based is empty
        assert chain.required_junction_positions_1based == [], (
            f"Expected empty required junctions, got {chain.required_junction_positions_1based}"
        )
        # But junction_positions_1based should still be populated (Primer3 soft penalty)
        assert len(chain.junction_positions_1based) > 0


# ── chain template extraction ────────────────────────────────────────────────


class TestConservedChainTemplates:
    def test_template_matches_exon_concatenation(self, multi_genome_gtf, multi_genome_fasta):
        """Each chain template should be the concatenation of its exons."""
        result = get_target_information(
            multi_genome_fasta, multi_genome_gtf, target_gene="multi_test_gene"
        )
        for chain in result:
            expected = extract_sequence(multi_genome_fasta, chain.exons)
            assert chain.template == expected, (
                f"Template for {len(chain.exons)}-exon chain doesn't match extract_sequence"
            )

    def test_junction_positions_are_correct(self, multi_genome_gtf, multi_genome_fasta):
        """Junction positions should align with exon boundaries in the template."""
        result = get_target_information(
            multi_genome_fasta, multi_genome_gtf, target_gene="multi_test_gene"
        )
        for chain in result:
            cum_len = 0
            for i, exon in enumerate(chain.exons[:-1]):
                exon_len = exon.end - exon.start + 1
                cum_len += exon_len
                expected_1based = cum_len + 1
                assert expected_1based in chain.junction_positions_1based, (
                    f"Junction between exon {i+1} and {i+2} "
                    f"(1-based pos {expected_1based}) not in positions list"
                )


# ── _compute_junction_positions ──────────────────────────────────────────────


class TestComputeJunctionPositions:
    """Tests for the private _compute_junction_positions helper."""

    def test_two_exons(self):
        """Two exons of equal length produce one junction in the middle."""
        from brimer_plast.genome import _compute_junction_positions

        exons = [ExonInfo("chrI", 100, 299, "+"), ExonInfo("chrI", 400, 599, "+")]
        # Template order: + strand, so 100-299 then 400-599.
        # Exon1_len = 200, junction at 201 (1-based)
        result = _compute_junction_positions(exons)
        assert result == [201]

    def test_three_exons(self):
        """Three exons produce two junctions."""
        from brimer_plast.genome import _compute_junction_positions

        exons = [
            ExonInfo("chrI", 100, 199, "+"),
            ExonInfo("chrI", 300, 399, "+"),
            ExonInfo("chrI", 500, 599, "+"),
        ]
        # Exon1=100bp, junction at 101. Exon2=100bp, junction at 201.
        result = _compute_junction_positions(exons)
        assert result == [101, 201]

    def test_unequal_exon_lengths(self):
        """Junction positions are correct for unequal exon lengths."""
        from brimer_plast.genome import _compute_junction_positions

        exons = [ExonInfo("chrI", 100, 349, "+"), ExonInfo("chrI", 450, 599, "+")]
        # Exon1_len = 250, junction at 251
        result = _compute_junction_positions(exons)
        assert result == [251]

    def test_single_exon_returns_empty(self):
        """Single exon has no junctions."""
        from brimer_plast.genome import _compute_junction_positions

        exons = [ExonInfo("chrI", 100, 200, "+")]
        assert _compute_junction_positions(exons) == []

    def test_empty_exons_returns_empty(self):
        """Empty list returns empty."""
        from brimer_plast.genome import _compute_junction_positions

        assert _compute_junction_positions([]) == []

    def test_negative_strand_uses_template_order(self):
        """Negative-strand exons use template order (reversed)."""
        from brimer_plast.genome import _compute_junction_positions

        # Template order for negative strand: 500-599 first, then 300-399
        exons = [ExonInfo("chrI", 300, 399, "-"), ExonInfo("chrI", 500, 599, "-")]
        # The function takes exons already in template order.
        result = _compute_junction_positions(exons)
        # Exon1_len = 100 (500-599), junction at 101
        assert result == [101]


# ── get_gene_locus ───────────────────────────────────────────────────────────


class TestGetGeneLocus:
    """Tests for get_gene_locus."""

    def test_returns_locus_for_gene_name(self, mini_genome_gtf):
        """Should return a GeneLocus for a known gene."""
        from brimer_plast.genome import get_gene_locus

        locus = get_gene_locus(mini_genome_gtf, target_gene="test_gene")
        assert locus is not None
        assert locus.gene_name == "test_gene"
        assert locus.seqid == "chrI"
        assert locus.strand == "+"
        assert locus.min_start == 101
        assert locus.max_end == 600

    def test_returns_locus_for_transcript(self, mini_genome_gtf):
        """Should resolve transcript to its gene locus."""
        from brimer_plast.genome import get_gene_locus

        locus = get_gene_locus(mini_genome_gtf, target_transcript="test_transcript")
        assert locus is not None
        assert locus.gene_name == "test_gene"
        assert "test_transcript" in locus.transcripts

    def test_raises_for_missing_gene(self, mini_genome_gtf):
        """Non-existent gene raises ValueError (consistent with parse_gtf)."""
        from brimer_plast.genome import get_gene_locus

        with pytest.raises(ValueError, match="not found"):
            get_gene_locus(mini_genome_gtf, target_gene="nonexistent")

    def test_multi_transcript_locus_cover_full_range(self, multi_genome_gtf):
        """Locus range should span the full gene."""
        from brimer_plast.genome import get_gene_locus

        locus = get_gene_locus(multi_genome_gtf, target_gene="multi_test_gene")
        assert locus is not None
        assert locus.min_start == 101
        assert locus.max_end == 1000

    def test_deduplicates_identical_transcripts(self, multi_genome_gtf):
        """Transcripts with identical exon structure are deduplicated."""
        from brimer_plast.genome import get_gene_locus

        locus = get_gene_locus(multi_genome_gtf, target_gene="multi_test_gene")
        # transcript_A and transcript_C have identical exon structure
        # Only one of them should appear in unique_transcripts.
        # Transcript B has a different structure, so total = 2 unique.
        assert len(locus.transcripts) == 2, (
            f"Expected 2 unique transcript structures, got {len(locus.transcripts)}"
        )
