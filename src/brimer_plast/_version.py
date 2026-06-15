"""Package version, derived from setuptools-scm at build time.

The version string is PEP 440 compliant.  When the package is installed,
it comes from the .dist-info metadata written by setuptools-scm.  When
running from source without install (e.g. via `python -m`), falls back
to a static placeholder.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _metadata_version

try:
    __version__ = _metadata_version("brimer-plast")
except (PackageNotFoundError, Exception):
    __version__ = "0.1.0"
