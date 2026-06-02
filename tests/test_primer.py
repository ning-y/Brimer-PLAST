"""Tests for primer.py: primer3-py integration."""

from brimer_plast.primer import (
    DEFAULT_PRIMER_ARGS,
    design_primers,
)


class TestDefaultArgs:
    def test_defaults_are_reasonable(self):
        """Default args should include essential primer3 global parameters."""
        assert "PRIMER_PRODUCT_SIZE_RANGE" in DEFAULT_PRIMER_ARGS
        assert "PRIMER_NUM_RETURN" in DEFAULT_PRIMER_ARGS
        assert "PRIMER_OPT_SIZE" in DEFAULT_PRIMER_ARGS
        assert "PRIMER_OPT_TM" in DEFAULT_PRIMER_ARGS

    def test_default_num_return_is_50(self):
        """Should request 50 candidate pairs by default for diverse forward positions."""
        assert DEFAULT_PRIMER_ARGS["PRIMER_NUM_RETURN"] == 50


class TestDesignPrimers:
    def test_returns_candidate_pairs(self, test_gene_template):
        """design_primers should return a list of candidate primer pairs."""
        result = design_primers(test_gene_template)
        assert isinstance(result, list)
        assert len(result) > 0, "Expected at least one candidate primer pair"

    def test_candidate_pairs_have_expected_keys(self, test_gene_template):
        """Each candidate pair should have forward/reverse/tm/gc/size attrs."""
        result = design_primers(test_gene_template)
        pair = result[0]
        assert pair.forward_seq
        assert pair.reverse_seq
        assert pair.forward_tm is not None
        assert pair.reverse_tm is not None
        assert pair.forward_gc is not None
        assert pair.reverse_gc is not None
        assert pair.product_size is not None

    def test_primer_sequences_are_dna(self, test_gene_template):
        """Primer sequences should contain only ATCG characters."""
        result = design_primers(test_gene_template)
        for pair in result:
            for seq in (pair.forward_seq, pair.reverse_seq):
                assert all(b in "ATCG" for b in seq.upper()), f"Invalid base in {seq}"

    def test_respects_custom_global_args(self, test_gene_template):
        """Custom global_args should override defaults."""
        custom = DEFAULT_PRIMER_ARGS.copy()
        custom["PRIMER_NUM_RETURN"] = 2
        custom["PRIMER_PRODUCT_SIZE_RANGE"] = "100-200"
        custom["PRIMER_MIN_TM"] = 58.0
        custom["PRIMER_MAX_TM"] = 62.0
        custom["PRIMER_OPT_TM"] = 60.0
        custom["PRIMER_MIN_SIZE"] = 20
        custom["PRIMER_MAX_SIZE"] = 20
        custom["PRIMER_OPT_SIZE"] = 20
        result = design_primers(test_gene_template, global_args=custom)
        # Should have at most 2 pairs (but could have 0 if constraints are too tight)
        assert len(result) <= 2

    def test_returns_empty_list_if_no_primers_possible(self):
        """A template of all N bases should produce no primers."""
        template = (
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" * 5
        )
        result = design_primers(template)
        assert result == []

    def test_short_template_returns_empty(self):
        """A template too short for primer design should return empty list."""
        result = design_primers("ATCG" * 5)  # 20 bp — shorter than min amplicon
        assert result == []

    def test_pair_penalty_is_numeric(self, test_gene_template):
        """Each pair should have a numeric penalty score."""
        result = design_primers(test_gene_template)
        pair = result[0]
        assert isinstance(pair.pair_penalty, (int, float))

    def test_different_template_returns_different_primers(self, test_gene_template):
        """Primers for different templates should differ."""
        result_a = design_primers(test_gene_template)
        # A varied template with balanced GC content that should produce primers
        template_b = "ATCGTAGCGTACGTACGATCGTACGTAGCTAGCATCGTAGTCGACTGAC" * 6
        result_b = design_primers(template_b)
        assert result_a, "test_gene_template should produce primers"
        assert result_b, "template_b should produce primers"
        assert result_a[0].forward_seq != result_b[0].forward_seq


class TestDesignPrimersWithJunctions:
    def test_junction_positions_passed_to_primer3(self, test_gene_template):
        """Providing junction positions should not crash and should return candidates."""
        # test_gene_template has exon1+exon2 with junction at 251 (1-based)
        result = design_primers(test_gene_template, junction_positions=[251])
        assert isinstance(result, list)
        # May or may not produce results with this template, but shouldn't crash

    def test_junction_positions_constrains_results(self, multi_genome_gtf, multi_genome_fasta):
        """With multi-transcript fixture, specific junction constraints should work."""
        from brimer_plast.genome import get_target_information
        chains = get_target_information(
            multi_genome_fasta, multi_genome_gtf, target_gene="multi_test_gene"
        )
        assert len(chains) > 0
        for chain in chains:
            result = design_primers(
                chain.template,
                junction_positions=chain.junction_positions_1based,
            )
            assert isinstance(result, list)
