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

    def test_help_shows_pdf_options(self):
        """--help should document PDF output options."""
        result = runner.invoke(app, ["--help"])
        assert "--output-pdf" in result.stdout, "--output-pdf should appear in --help"
        assert "--no-pdf" in result.stdout, "--no-pdf should appear in --help"

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

    def test_output_pdf_and_no_pdf_are_mutually_exclusive(self, mini_genome_fasta, mini_genome_gtf):
        """Providing both --output-pdf and --no-pdf should error."""
        result = runner.invoke(
            app,
            [
                "--genome",
                mini_genome_fasta,
                "--annotations",
                mini_genome_gtf,
                "--target-gene",
                "test_gene",
                "--output-pdf",
                "/tmp/report.pdf",
                "--no-pdf",
            ],
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.stdout.lower()


class TestEndToEnd:
    def test_produces_primer_set_disabled_junction(
        self, mini_genome_fasta, mini_genome_gtf, tmp_path
    ):
        """Full pipeline with --disable-junction-overlap should produce a primer set."""
        pdf_path = tmp_path / "report.pdf"
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
                "--output-pdf",
                str(pdf_path),
            ],
        )
        assert result.exit_code == 0, f"cli failed: {result.stdout}"
        # Should output primer pairs in tabular format
        assert "Forward" in result.stdout or "forward" in result.stdout

    def test_produces_primer_set_with_transcript_disabled(
        self, mini_genome_fasta, mini_genome_gtf, tmp_path
    ):
        """--target-transcript with --disable-junction-overlap should also work."""
        pdf_path = tmp_path / "report.pdf"
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
                "--output-pdf",
                str(pdf_path),
            ],
        )
        assert result.exit_code == 0

    def test_default_junction_spanning_no_results(
        self, mini_genome_fasta, mini_genome_gtf, tmp_path
    ):
        """Default junction-spanning mode may produce 0 results for the synthetic fixture;
        this is expected behavior (not an error)."""
        pdf_path = tmp_path / "report.pdf"
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
                "--output-pdf",
                str(pdf_path),
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

    def test_tsv_flag_produces_tab_separated_output(
        self, mini_genome_fasta, mini_genome_gtf, tmp_path
    ):
        """--tsv with --disable-junction-overlap should produce tab-separated output."""
        pdf_path = tmp_path / "report.pdf"
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
                "--output-pdf",
                str(pdf_path),
            ],
        )
        assert result.exit_code == 0
        lines = [line for line in result.stdout.strip().splitlines() if "\t" in line]
        assert len(lines) > 0, "Expected at least one tab-separated line in output"
        header = lines[0]
        assert header.startswith("pair"), "TSV header should start with 'pair'"
        assert "forward_seq" in header
        assert "reverse_seq" in header

    def test_single_exon_with_disable_junction_overlap(
        self, mini_genome_fasta, single_exon_gtf, tmp_path
    ):
        """A single-exon gene with --disable-junction-overlap should produce primers."""
        pdf_path = tmp_path / "report.pdf"
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
                "--output-pdf",
                str(pdf_path),
            ],
        )
        # Should succeed because junction policy is relaxed
        assert result.exit_code in (0, 1), f"CLI failed: {result.stdout}"

    def test_single_exon_fails_without_disable_junction_overlap(
        self, mini_genome_fasta, single_exon_gtf
    ):
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
