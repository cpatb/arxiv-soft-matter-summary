"""Cited fluid records: loaded from YAML, validated, envelope-enforced."""

from .loader import RECORDS_DIR, load_all, load_record
from .provenance import DOI_RE, ResolvedDOI, is_doi_syntax, resolve_doi
from .schema import (
    Envelope,
    EnvelopeError,
    FluidRecord,
    Licence,
    Parameter,
    ProvenanceError,
    Source,
)

__all__ = [
    "RECORDS_DIR", "load_all", "load_record",
    "DOI_RE", "ResolvedDOI", "is_doi_syntax", "resolve_doi",
    "Envelope", "EnvelopeError", "FluidRecord", "Licence", "Parameter",
    "ProvenanceError", "Source",
]
