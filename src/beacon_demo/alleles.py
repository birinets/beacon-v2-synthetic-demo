"""Allele representation: the same variant can be written with or without a shared anchor base.

A VCF writes a deletion of A after G as ``GA>G`` at the G; a query may write it as ``A>-`` one base later. Both
describe the same change. ``minimal`` trims the bases the two alleles share, from the end and then from the start,
and moves the start past any trimmed leading bases, so both spellings reduce to the same (start, ref, alt). The idea
follows the reference rule of genotype-pipeline-checks (https://github.com/birinets/genotype-pipeline-checks).
Trimming does not left-align a variant in a repeat; two spellings that differ by a shift in a repeat do not match.
"""

from __future__ import annotations

EMPTY = ("", "-", ".")


def clean(bases):
    """Upper-case a query allele and read "-" or "." as the empty allele."""
    if bases is None:
        return None
    bases = bases.strip().upper()
    return "" if bases in EMPTY else bases


def minimal(start, ref, alt):
    """(start, ref, alt) with shared trailing, then leading, bases removed; start is 0-based interbase."""
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, start = ref[1:], alt[1:], start + 1
    return start, ref, alt


def variant_type(min_ref, min_alt):
    if len(min_ref) == 1 and len(min_alt) == 1:
        return "SNP"
    if not min_alt and min_ref:
        return "DEL"
    if not min_ref and min_alt:
        return "INS"
    if len(min_ref) == len(min_alt):
        return "MNP"
    return "INDEL"


def length(min_ref, min_alt):
    """Number of bases the variant changes: 1 for a substitution, the inserted or deleted length for an indel."""
    return max(len(min_ref), len(min_alt))
