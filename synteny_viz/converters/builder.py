"""Fluent builder for constructing SyntenyData from various sources."""

from pathlib import Path

from ..models import Gene, HighlightGene, Species, StyleConfig, SyntenyData
from .bed import load_genes_from_bed
from .gff import load_genes_from_gff
from .orthology import load_orthogroups, merge_orthogroups_into_species


class SyntenyBuilder:
    """Fluent builder for SyntenyData.

    Example:
        data = (
            SyntenyBuilder()
            .add_species_from_bed("human", "Human", "human_genes.bed", chrom="chr4")
            .add_species_from_bed("mouse", "Mouse", "mouse_genes.bed", chrom="chr5")
            .load_orthogroups("orthogroups.tsv")
            .highlight("TIGD4", color="#e74c3c")
            .build()
        )
    """

    def __init__(self, title: str | None = None):
        self._title = title
        self._species: list[Species] = []
        self._highlights: list[HighlightGene] = []
        self._style = StyleConfig()
        self._orthogroups: dict[str, str] = {}

    def add_species(
        self,
        species_id: str,
        display_name: str,
        genes: list[Gene],
    ) -> "SyntenyBuilder":
        """Add a species with pre-loaded genes."""
        species = Species(
            species_id=species_id,
            display_name=display_name,
            genes=genes,
        )
        self._species.append(species)
        return self

    def add_species_from_bed(
        self,
        species_id: str,
        display_name: str,
        bed_path: str | Path,
        chrom: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> "SyntenyBuilder":
        """Add a species by loading genes from a BED file."""
        genes = load_genes_from_bed(
            bed_path,
            species_id=species_id,
            chrom_filter=chrom,
            start_filter=start,
            end_filter=end,
        )
        return self.add_species(species_id, display_name, genes)

    def add_species_from_gff(
        self,
        species_id: str,
        display_name: str,
        gff_path: str | Path,
        feature_type: str = "gene",
        chrom: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> "SyntenyBuilder":
        """Add a species by loading genes from a GFF/GTF file."""
        genes = load_genes_from_gff(
            gff_path,
            species_id=species_id,
            feature_type=feature_type,
            chrom_filter=chrom,
            start_filter=start,
            end_filter=end,
        )
        return self.add_species(species_id, display_name, genes)

    def load_orthogroups(
        self,
        path: str | Path,
        gene_col: int = 0,
        orthogroup_col: int = 1,
    ) -> "SyntenyBuilder":
        """Load orthogroup assignments from a file."""
        og = load_orthogroups(path, gene_col=gene_col, orthogroup_col=orthogroup_col)
        self._orthogroups.update(og)
        return self

    def set_orthogroups(self, orthogroups: dict[str, str]) -> "SyntenyBuilder":
        """Set orthogroup assignments directly."""
        self._orthogroups.update(orthogroups)
        return self

    def highlight(
        self,
        gene_name: str,
        color: str = "#e74c3c",
        label: str | None = None,
    ) -> "SyntenyBuilder":
        """Add a gene to highlight across all species."""
        self._highlights.append(
            HighlightGene(name=gene_name, color=color, label=label)
        )
        return self

    def style(self, **kwargs) -> "SyntenyBuilder":
        """Update style configuration."""
        self._style = StyleConfig(**{**self._style.model_dump(), **kwargs})
        return self

    def build(self) -> SyntenyData:
        """Build the final SyntenyData object."""
        # Apply orthogroups to all species
        if self._orthogroups:
            merge_orthogroups_into_species(self._species, self._orthogroups)

        return SyntenyData(
            title=self._title,
            species=self._species,
            highlights=self._highlights,
            style=self._style,
        )
