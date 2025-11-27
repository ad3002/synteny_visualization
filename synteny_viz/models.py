"""Data models for synteny visualization."""

from enum import Enum
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field


class Strand(str, Enum):
    PLUS = "+"
    MINUS = "-"
    UNKNOWN = "."


class Gene(BaseModel):
    """Single gene on a chromosome."""

    gene_id: str = Field(..., description="Unique gene identifier")
    name: str = Field(..., description="Gene name for display")
    chrom: str = Field(..., description="Chromosome/scaffold name")
    start: int = Field(..., description="Start position (0-based, can be negative if normalized)")
    end: int = Field(..., description="End position")
    strand: Strand = Field(default=Strand.UNKNOWN, description="Strand orientation")
    orthogroup: str | None = Field(default=None, description="Orthogroup ID for linking")

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def center(self) -> float:
        return (self.start + self.end) / 2


class Species(BaseModel):
    """Species with its genes in a genomic region."""

    species_id: str = Field(..., description="Unique species identifier")
    display_name: str = Field(..., description="Name shown on plot")
    genes: list[Gene] = Field(default_factory=list, description="Genes in the region")

    @property
    def region_start(self) -> int:
        if not self.genes:
            return 0
        return min(g.start for g in self.genes)

    @property
    def region_end(self) -> int:
        if not self.genes:
            return 0
        return max(g.end for g in self.genes)

    @property
    def region_length(self) -> int:
        return self.region_end - self.region_start


class HighlightGene(BaseModel):
    """Gene to highlight across all species."""

    name: str = Field(..., description="Gene name to highlight")
    orthogroup: str | None = Field(default=None, description="Orthogroup ID to highlight")
    color: str = Field(default="#e74c3c", description="Highlight color (hex)")
    label: str | None = Field(default=None, description="Custom label in legend")


class StyleConfig(BaseModel):
    """Visual style configuration."""

    # Colors
    gene_color: str = Field(default="#808080", description="Default gene color")
    gene_border_color: str = Field(default="#404040", description="Gene border color")
    link_color: str = Field(default="#cccccc", description="Orthology link fill color")
    link_alpha: float = Field(default=0.4, ge=0, le=1, description="Link transparency")
    highlight_link_alpha: float = Field(default=0.6, ge=0, le=1)
    chromosome_color: str = Field(default="#000000", description="Chromosome line color")

    # Sizes
    gene_height: float = Field(default=0.3, description="Gene rectangle height (relative)")
    chromosome_width: float = Field(default=1.0, description="Chromosome line width")
    species_spacing: float = Field(default=1.0, description="Vertical spacing between species")

    # Figure
    figure_width: float = Field(default=16, description="Figure width in inches")
    figure_height_per_species: float = Field(default=0.8, description="Height per species row")
    dpi: int = Field(default=150, description="Output resolution")

    # Fonts
    species_font_size: int = Field(default=10)
    gene_label_font_size: int = Field(default=8)
    scale_font_size: int = Field(default=9)

    # Layout
    show_scale_bar: bool = Field(default=True)
    show_gene_labels: bool = Field(default=True)
    gene_label_position: Literal["bottom", "inside", "none"] = Field(default="bottom")


class SpeciesGroup(BaseModel):
    """Group of species for bracket annotation."""

    name: str = Field(..., description="Group name (e.g., 'Mammals', 'Birds')")
    species_ids: list[str] = Field(..., description="Species IDs in this group")
    color: str = Field(default="#333333", description="Bracket color")
    children: list["SpeciesGroup"] = Field(
        default_factory=list, description="Nested subgroups for hierarchical brackets"
    )


class SyntenyData(BaseModel):
    """Complete synteny visualization data."""

    title: str | None = Field(default=None, description="Plot title")
    species: list[Species] = Field(..., description="Species in display order (top to bottom)")
    highlights: list[HighlightGene] = Field(
        default_factory=list, description="Genes to highlight"
    )
    groups: list[SpeciesGroup] = Field(
        default_factory=list, description="Species groups for bracket annotations"
    )
    style: StyleConfig = Field(default_factory=StyleConfig)

    def get_orthogroups(self) -> dict[str, list[tuple[str, Gene]]]:
        """Get all orthogroups with their genes."""
        groups: dict[str, list[tuple[str, Gene]]] = {}
        for sp in self.species:
            for gene in sp.genes:
                if gene.orthogroup:
                    if gene.orthogroup not in groups:
                        groups[gene.orthogroup] = []
                    groups[gene.orthogroup].append((sp.species_id, gene))
        return groups

    def is_highlighted(self, gene: Gene) -> HighlightGene | None:
        """Check if gene should be highlighted."""
        for h in self.highlights:
            # Match by gene name (case-insensitive)
            if gene.name.upper() == h.name.upper():
                return h
            # Match by orthogroup
            if h.orthogroup and gene.orthogroup and gene.orthogroup == h.orthogroup:
                return h
        return None

    def count_shared_orthogroups(self, species_ids: list[str]) -> int:
        """Count orthogroups shared by all species in the list."""
        if not species_ids:
            return 0

        # Get orthogroups for each species
        og_by_species: list[set[str]] = []
        for sp_id in species_ids:
            sp = next((s for s in self.species if s.species_id == sp_id), None)
            if sp:
                og_set = {g.orthogroup for g in sp.genes if g.orthogroup}
                og_by_species.append(og_set)

        if not og_by_species:
            return 0

        # Intersection of all
        shared = og_by_species[0]
        for og_set in og_by_species[1:]:
            shared = shared & og_set

        return len(shared)

    @classmethod
    def from_json_file(cls, path: str | Path) -> "SyntenyData":
        """Load from JSON file."""
        import json

        with open(path) as f:
            data = json.load(f)
        return cls.model_validate(data)

    def to_json_file(self, path: str | Path, indent: int = 2) -> None:
        """Save to JSON file."""
        with open(path, "w") as f:
            f.write(self.model_dump_json(indent=indent))
