"""
Alignments are representative regions between two sequences.

.. module:: dasi.alignments

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    alignment
    alignment_container
"""

from .alignment import (
    Alignment,
    AlignmentGroup,
    ComplexAlignmentGroup,
    AlignmentGroupBase,
)
from .alignment_container import AlignmentContainer, AlignmentContainerFactory
