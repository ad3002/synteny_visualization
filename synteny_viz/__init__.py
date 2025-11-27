"""Synteny visualization tool for comparative genomics."""

from .models import Gene, HighlightGene, Species, Strand, StyleConfig, SyntenyData
from .visualizer import SyntenyVisualizer, render_synteny

__all__ = [
    "Gene",
    "HighlightGene",
    "Species",
    "Strand",
    "StyleConfig",
    "SyntenyData",
    "SyntenyVisualizer",
    "render_synteny",
]
