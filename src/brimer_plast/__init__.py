"""Brimer-PLAST: PCR primer design using primer3-py + tnBLAST."""

from __future__ import annotations

from ._version import __version__

# Single source of truth for the application title/branding.
# Both the PDF report and the Electron renderer derive from this.
APP_TITLE = "Brimer-PLAST by Wang Linfa Lab"

from brimer_plast.models import ConservedExonChain, ExonInfo, PrimerPair

__all__ = [
    "ConservedExonChain",
    "ExonInfo",
    "PrimerPair",
]
