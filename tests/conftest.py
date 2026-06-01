"""Pytest fixtures for Brimer-PLAST tests."""

import pytest


@pytest.fixture
def mini_genome_fasta() -> str:
    """Path to the synthetic mini genome FASTA."""
    from pathlib import Path

    return str(Path(__file__).parent / "fixtures" / "mini_genome.fa")


@pytest.fixture
def mini_genome_gtf() -> str:
    """Path to the synthetic mini genome GTF."""
    from pathlib import Path

    return str(Path(__file__).parent / "fixtures" / "mini_genome.gtf")


@pytest.fixture
def mini_genome_sequence(mini_genome_fasta: str) -> str:
    """The full chromosome sequence from the mini genome."""
    with open(mini_genome_fasta) as f:
        lines = f.readlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


@pytest.fixture
def test_gene_template() -> str:
    """The spliced template sequence for test_gene (exon1 + exon2)."""
    from pathlib import Path

    genome_path = Path(__file__).parent / "fixtures" / "mini_genome.fa"
    with open(genome_path) as f:
        lines = f.readlines()
    genome = "".join(line.strip() for line in lines if not line.startswith(">"))
    # exon1: 101..350, exon2: 451..600 (1-indexed GTF)
    exon1 = genome[100:350]
    exon2 = genome[450:600]
    return exon1 + exon2


@pytest.fixture
def multi_genome_fasta() -> str:
    """Path to the multi-transcript synthetic genome FASTA."""
    from pathlib import Path

    return str(Path(__file__).parent / "fixtures" / "multi_transcript_genome.fa")


@pytest.fixture
def multi_genome_gtf() -> str:
    """Path to the multi-transcript synthetic genome GTF."""
    from pathlib import Path

    return str(Path(__file__).parent / "fixtures" / "multi_transcript.gtf")


@pytest.fixture
def single_exon_gtf() -> str:
    """Path to a single-exon gene GTF."""
    from pathlib import Path

    return str(Path(__file__).parent / "fixtures" / "single_exon.gtf")


@pytest.fixture
def mixed_exon_gtf() -> str:
    """Path to a mixed (single + multi exon) gene GTF."""
    from pathlib import Path

    return str(Path(__file__).parent / "fixtures" / "mixed_exon.gtf")
