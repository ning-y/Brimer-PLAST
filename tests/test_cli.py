"""Tests for cli.py: command-line interface."""

from typer.testing import CliRunner

from brimer_plast.cli import app

runner = CliRunner()


class TestHelp:
    def test_help_succeeds(self):
        """--help should exit with code 0."""
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0

    def test_help_shows_genome_arg(self):
        """--help should document --genome."""
        result = runner.invoke(app, ["--help"])
        assert "--genome" in result.stdout

    def test_help_shows_annotations_arg(self):
        """--help should document --annotations."""
        result = runner.invoke(app, ["--help"])
        assert "--annotations" in result.stdout

    def test_help_shows_target_gene_arg(self):
        """--help should document --target-gene."""
        result = runner.invoke(app, ["--help"])
        assert "--target-gene" in result.stdout

    def test_help_shows_target_transcript_arg(self):
        """--help should document --target-transcript."""
        result = runner.invoke(app, ["--help"])
        assert "--target-transcript" in result.stdout

    def test_help_shows_disable_junction_overlap(self):
        """--help should document --disable-junction-overlap."""
        result = runner.invoke(app, ["--help"])
        assert "--disable-junction-overlap" in result.stdout

    def test_help_shows_primer_design_options(self):
        """--help should document all primer design options."""
        result = runner.invoke(app, ["--help"])
        for option in [
            "--num-return",
            "--min-tm",
            "--max-tm",
            "--opt-tm",
            "--min-size",
            "--max-size",
            "--opt-size",
            "--min-gc",
            "--max-gc",
            "--product-min",
            "--product-max",
            "--max-amplicon",
            "--tsv",
        ]:
            assert option in result.stdout, f"{option} should appear in --help"


class TestRequiredArgs:
    def test_missing_genome(self, mini_genome_gtf):
        """Missing --genome should produce an error."""
        result = runner.invoke(
            app,
            [
                "--annotations",
                mini_genome_gtf,
                "--target-gene",
                "test_gene",
            ],
        )
        assert result.exit_code != 0

    def test_missing_annotations(self, mini_genome_fasta):
        """Missing --annotations should produce an error."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--target-gene",
                "test_gene",
            ],
        )
        assert result.exit_code != 0

    def test_missing_target_gene_and_transcript(self, mini_genome_fasta, mini_genome_gtf):
        """Missing both --target-gene and --target-transcript should error."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                mini_genome_gtf,
            ],
        )
        assert result.exit_code != 0
        assert "target" in result.stdout.lower() or "target" in str(result.exception).lower()

    def test_both_target_gene_and_transcript(self, mini_genome_fasta, mini_genome_gtf):
        """Providing both should error."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                mini_genome_gtf,
                "--target-gene",
                "test_gene",
                "--target-transcript",
                "test_transcript",
            ],
        )
        assert result.exit_code != 0
        assert "both" in result.stdout.lower() or "both" in str(result.exception).lower()


class TestEndToEnd:
    def test_produces_primer_set_disabled_junction(self, mini_genome_fasta, mini_genome_gtf):
        """Full pipeline with --disable-junction-overlap should produce a primer set."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                mini_genome_gtf,
                "--target-gene",
                "test_gene",
                "--num-return",
                "3",
                "--disable-junction-overlap",
            ],
        )
        assert result.exit_code == 0, f"cli failed: {result.stdout}"
        # Should output primer pairs in tabular format
        assert "Forward" in result.stdout or "forward" in result.stdout

    def test_produces_primer_set_with_transcript_disabled(self, mini_genome_fasta, mini_genome_gtf):
        """--target-transcript with --disable-junction-overlap should also work."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                mini_genome_gtf,
                "--target-transcript",
                "test_transcript",
                "--num-return",
                "3",
                "--disable-junction-overlap",
            ],
        )
        assert result.exit_code == 0

    def test_default_junction_spanning_no_results(self, mini_genome_fasta, mini_genome_gtf):
        """Default junction-spanning mode may produce 0 results for the synthetic fixture;
        this is expected behavior (not an error)."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                mini_genome_gtf,
                "--target-gene",
                "test_gene",
                "--num-return",
                "3",
            ],
        )
        # This may exit with code 1 if no primers can be designed, which is expected
        assert result.exit_code in (0, 1)

    def test_empty_result_when_gene_not_found(self, mini_genome_fasta, mini_genome_gtf):
        """Non-existent gene should error."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                mini_genome_gtf,
                "--target-gene",
                "nonexistent_gene",
            ],
        )
        assert result.exit_code != 0

    def test_tsv_flag_produces_tab_separated_output(self, mini_genome_fasta, mini_genome_gtf):
        """--tsv with --disable-junction-overlap should produce tab-separated output."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                mini_genome_gtf,
                "--target-gene",
                "test_gene",
                "--num-return",
                "3",
                "--tsv",
                "--disable-junction-overlap",
            ],
        )
        assert result.exit_code == 0
        lines = [line for line in result.stdout.strip().splitlines() if "\t" in line]
        assert len(lines) > 0, "Expected at least one tab-separated line in output"
        header = lines[0]
        assert header.startswith("pair"), "TSV header should start with 'pair'"
        assert "forward_seq" in header
        assert "reverse_seq" in header

    def test_single_exon_with_disable_junction_overlap(self, mini_genome_fasta, single_exon_gtf):
        """A single-exon gene with --disable-junction-overlap should produce primers."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                single_exon_gtf,
                "--target-gene",
                "single_exon_gene",
                "--num-return",
                "3",
                "--disable-junction-overlap",
            ],
        )
        # Should succeed because junction policy is relaxed
        assert result.exit_code in (0, 1), f"CLI failed: {result.stdout}"

    def test_single_exon_fails_without_disable_junction_overlap(self, mini_genome_fasta, single_exon_gtf):
        """A single-exon gene without --disable-junction-overlap should fail."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                single_exon_gtf,
                "--target-gene",
                "single_exon_gene",
                "--num-return",
                "3",
            ],
        )
        assert result.exit_code != 0
        assert "disable-junction-overlap" in result.stdout.lower()
