"""Integration tests with real C. elegans genome (WBcel235)."""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from brimer_plast.cli import app
from brimer_plast.genome import get_target_information
from brimer_plast.primer import design_primers

CE11_DIR = Path(__file__).parent / "fixtures" / "ce11"
GENOME_FASTA = CE11_DIR / "genome.fna"
GTF_FILE = CE11_DIR / "annotations.gtf"

runner = CliRunner()


def _has_ce11() -> bool:
    return GENOME_FASTA.exists() and GTF_FILE.exists()


ce11 = pytest.mark.skipif(
    not _has_ce11(),
    reason="C. elegans data not downloaded. Run: bash tests/fixtures/download-ce11.sh",
)


# ── Target gene selection ──────────────────────────────────────────────────
# homt-1 (gene ID: WBGene00022277) on chromosome I: protein-coding, well-annotated,
# 5 exons, ~6kb gene — a good integration test candidate.


@ce11
class TestFindGenes:
    """Validate we can find known C. elegans genes in the GTF."""

    def test_find_homt1_by_gene_name(self):
        """homt-1 should be findable by gene name."""
        chains = get_target_information(
            GENOME_FASTA,
            GTF_FILE,
            target_gene="homt-1",
        )
        assert len(chains) > 0
        assert len(chains[0].template) > 0, "homt-1 template should be non-empty"
        assert "homt-1" in chains[0].id

    def test_find_homt1_by_transcript(self):
        """homt-1's NM transcript should be findable."""
        chains = get_target_information(
            GENOME_FASTA,
            GTF_FILE,
            target_transcript="NM_058260.8",
        )
        assert len(chains) == 1
        assert "NM_058260.8" in chains[0].id
        assert len(chains[0].template) > 0

    def test_template_is_spliced(self):
        """The spliced template should not contain intronic sequence."""
        chains = get_target_information(
            GENOME_FASTA,
            GTF_FILE,
            target_gene="homt-1",
        )
        template = chains[0].template
        # homt-1 has 5 exons spanning 4116..10230 (6.1 kb), spliced ~892 bp
        assert len(template) < 3000, (
            f"homt-1 template ({len(template)} bp) should be spliced, "
            f"much shorter than the 6 kb gene span"
        )
        assert len(template) > 200, f"homt-1 template ({len(template)} bp) seems too short"


@ce11
class TestPrimerDesignOnRealGenes:
    """Design primers for real C. elegans genes."""

    def test_homt1_produces_candidates(self):
        """homt-1 should produce at least one candidate primer pair from the first chain."""
        chains = get_target_information(
            GENOME_FASTA,
            GTF_FILE,
            target_gene="homt-1",
        )
        assert len(chains) > 0
        pairs = design_primers(
            chains[0].template,
            sequence_id=chains[0].id,
            junction_positions=chains[0].junction_positions_1based,
        )
        assert len(pairs) > 0, (
            f"Expected at least 1 candidate pair for homt-1 chain 1 "
            f"(template: {len(chains[0].template)} bp)"
        )

    def test_sesn1_produces_candidates(self):
        """sesn-1 should produce at least one candidate pair (no junction filter)."""
        chains = get_target_information(
            GENOME_FASTA,
            GTF_FILE,
            target_gene="sesn-1",
        )
        assert len(chains) > 0
        pairs = design_primers(
            chains[0].template,
            sequence_id=chains[0].id,
            # No junction constraints — just testing primer3 works on this gene
        )
        assert len(pairs) > 0

    def test_primers_are_dna(self):
        """All candidate primers should have only ATCG bases."""
        chains = get_target_information(
            GENOME_FASTA,
            GTF_FILE,
            target_gene="homt-1",
        )
        assert len(chains) > 0
        pairs = design_primers(
            chains[0].template,
            junction_positions=chains[0].junction_positions_1based,
        )
        for pair in pairs:
            for seq in (pair.forward_seq, pair.reverse_seq):
                assert all(b in "ATCG" for b in seq.upper()), f"Invalid base in {seq}"


@ce11
class TestEndToEnd:
    """Full pipeline via CLI."""

    def test_cli_homt1_produces_output_dual_mode(self, tmp_path):
        """brimer-plast --target-gene homt-1 with default dual mode."""
        pdf_path = tmp_path / "report.pdf"
        result = runner.invoke(
            app,
            [
                "--genome",
                str(GENOME_FASTA),
                "--annotations",
                str(GTF_FILE),
                "--target-gene",
                "homt-1",
                "--num-return",
                "5",
                "--output-pdf",
                str(pdf_path),
            ],
        )
        assert (
            result.exit_code == 0
        ), f"CLI failed:\nstdout: {result.stdout}"
        # Output should have primer table with columns
        assert "Forward" in result.stdout

    def test_cli_with_specific_tm(self, tmp_path):
        """Custom temperature range with dual mode."""
        pdf_path = tmp_path / "report.pdf"
        result = runner.invoke(
            app,
            [
                "--genome",
                str(GENOME_FASTA),
                "--annotations",
                str(GTF_FILE),
                "--target-gene",
                "homt-1",
                "--num-return",
                "3",
                "--min-tm",
                "58",
                "--max-tm",
                "62",
                "--opt-tm",
                "60",
                "--output-pdf",
                str(pdf_path),
            ],
        )
        assert (
            result.exit_code == 0
        ), f"CLI failed:\nstdout: {result.stdout}"

    def test_cli_rejects_nonexistent_gene(self):
        """Non-existent gene should exit with error."""
        result = runner.invoke(
            app,
            [
                "--genome",
                str(GENOME_FASTA),
                "--annotations",
                str(GTF_FILE),
                "--target-gene",
                "nonexistent_gene",
            ],
        )
        assert result.exit_code != 0
