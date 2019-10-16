r"""
Constants (:mod:`dasi.constants`)
=============================

.. currentmodule:: dasi.constants

This module provides DASi constants.
"""


class Constants:
    """DASi constants."""

    FRAGMENT = (
        "PRE-MADE DNA FRAGMENT"
    )  #: an alignment that is generate from an already existing PCR product or fragment
    PCR_PRODUCT = (
        "PCR_PRODUCT"
    )  #: an alignment that is to be generated from a PCR product
    PCR_PRODUCT_WITH_PRIMERS = (
        "PCR_PRODUCT_WITH_PRIMERS"
    )  #: PCR product that can be produces from existing primers
    PCR_PRODUCT_WITH_LEFT_PRIMER = (
        "PCR_PRODUCT_WITH_LEFT_PRIMER"
    )  #: PCR product with existing left primer
    PCR_PRODUCT_WITH_RIGHT_PRIMER = (
        "PCR_PRODUCT_WITH_RIGHT_PRIMER"
    )  #: PCR product with existing right primer
    SHARED_FRAGMENT = (
        "FRAGMENT_SHARED_WITH_OTHER_QUERIES"
    )  #: A fragment alignment that is shared with other queries for potential reuse
    PRIMER_EXTENSION_PRODUCT = (
        "PRIMER_EXTENSION_PRODUCT"
    )  #: a pcr product that uses no template, but extends two primers
    PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER = (
        "PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER"
    )  #: a pcr product that uses no template, but extends an existing fwd primer.
    PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER = (
        "PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER"
    )  #: a pcr product that uses no template, but extends an existing rev primer.

    TEMPLATE = "TEMPLATE"  #: A template alignment group. Not an actual molecule.
    GAP = "GAP"  #: region that represents a gap that must be synthesized
    OVERLAP = "OVERLAP"  #: region that represents overlapping molecules
    MISSING = "__MISSING"  #: missing region

    PRIMER = "PRIMER"  #: a primer binding alignment

    PRIMER_MIN_BIND = 14  #: minimum primer binding for searching for primer alignments
    MIN_OVERLAP = 15  #: minimum overlap for searching for overlapping alignments
    MAX_HOMOLOGY = 100  #: maximum overlap for searching for overlapping alignments
    INF = 10.0 ** 6  #: almost infinity
