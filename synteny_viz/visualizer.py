"""Core synteny visualization engine using matplotlib."""

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

from .models import Gene, HighlightGene, Species, SyntenyData, StyleConfig


class SyntenyVisualizer:
    """Renders synteny plots from SyntenyData."""

    def __init__(self, data: SyntenyData):
        self.data = data
        self.style = data.style
        self._fig: plt.Figure | None = None
        self._ax: plt.Axes | None = None

        # Computed layout
        self._x_scale: float = 1.0
        self._x_offset: float = 0.0
        self._species_y: dict[str, float] = {}

    def render(self, output_path: str | Path | None = None) -> plt.Figure:
        """Render the synteny plot."""
        self._compute_layout()
        self._create_figure()
        self._draw_scale_bar()
        self._draw_orthology_links()
        self._draw_species()
        self._draw_legend()
        self._finalize()

        if output_path:
            self._fig.savefig(
                output_path,
                dpi=self.style.dpi,
                bbox_inches="tight",
                facecolor="white",
                edgecolor="none",
            )

        return self._fig

    def _compute_layout(self) -> None:
        """Compute coordinate transforms and positions."""
        # Find global coordinate range
        all_starts = []
        all_ends = []
        for sp in self.data.species:
            if sp.genes:
                all_starts.append(sp.region_start)
                all_ends.append(sp.region_end)

        if not all_starts:
            self._x_scale = 1.0
            self._x_offset = 0.0
            return

        global_start = min(all_starts)
        global_end = max(all_ends)
        global_range = global_end - global_start

        # Scale to fit in plot width (0 to 1)
        self._x_scale = 1.0 / global_range if global_range > 0 else 1.0
        self._x_offset = -global_start * self._x_scale

        # Compute Y positions for each species (top to bottom)
        for i, sp in enumerate(self.data.species):
            self._species_y[sp.species_id] = -i * self.style.species_spacing

    def _create_figure(self) -> None:
        """Create matplotlib figure and axes."""
        n_species = len(self.data.species)
        height = n_species * self.style.figure_height_per_species + 1.5

        self._fig, self._ax = plt.subplots(figsize=(self.style.figure_width, height))
        self._ax.set_xlim(-0.05, 1.1)
        y_min = -(n_species - 1) * self.style.species_spacing - 0.5
        self._ax.set_ylim(y_min, 1.0)
        self._ax.set_aspect("auto")
        self._ax.axis("off")

    def _draw_scale_bar(self) -> None:
        """Draw scale bar at top."""
        if not self.style.show_scale_bar:
            return

        # Determine nice scale bar length
        all_lengths = [sp.region_length for sp in self.data.species if sp.genes]
        if not all_lengths:
            return

        max_length = max(all_lengths)
        scale_bp = self._nice_scale_value(max_length / 4)
        scale_width = scale_bp * self._x_scale

        # Draw scale bar
        y_pos = 0.7
        self._ax.plot([0, scale_width], [y_pos, y_pos], "k-", linewidth=2)
        self._ax.plot([0, 0], [y_pos - 0.02, y_pos + 0.02], "k-", linewidth=2)
        self._ax.plot(
            [scale_width, scale_width], [y_pos - 0.02, y_pos + 0.02], "k-", linewidth=2
        )

        # Label
        label = self._format_bp(scale_bp)
        self._ax.text(
            scale_width / 2,
            y_pos + 0.05,
            label,
            ha="center",
            va="bottom",
            fontsize=self.style.scale_font_size,
        )

    def _draw_orthology_links(self) -> None:
        """Draw polygons connecting orthologous genes between adjacent species."""
        orthogroups = self.data.get_orthogroups()

        for og_id, genes_list in orthogroups.items():
            # Group genes by species
            genes_by_species: dict[str, Gene] = {}
            for sp_id, gene in genes_list:
                genes_by_species[sp_id] = gene

            # Draw links between adjacent species
            for i in range(len(self.data.species) - 1):
                sp1 = self.data.species[i]
                sp2 = self.data.species[i + 1]

                if sp1.species_id not in genes_by_species:
                    continue
                if sp2.species_id not in genes_by_species:
                    continue

                gene1 = genes_by_species[sp1.species_id]
                gene2 = genes_by_species[sp2.species_id]

                self._draw_gene_link(sp1.species_id, gene1, sp2.species_id, gene2)

    def _draw_gene_link(
        self, sp1_id: str, gene1: Gene, sp2_id: str, gene2: Gene
    ) -> None:
        """Draw a single orthology link between two genes."""
        y1 = self._species_y[sp1_id]
        y2 = self._species_y[sp2_id]
        half_h = self.style.gene_height / 2

        # Gene positions in plot coordinates
        x1_left = gene1.start * self._x_scale + self._x_offset
        x1_right = gene1.end * self._x_scale + self._x_offset
        x2_left = gene2.start * self._x_scale + self._x_offset
        x2_right = gene2.end * self._x_scale + self._x_offset

        # Polygon vertices (connect bottom of top gene to top of bottom gene)
        vertices = [
            (x1_left, y1 - half_h),
            (x1_right, y1 - half_h),
            (x2_right, y2 + half_h),
            (x2_left, y2 + half_h),
        ]

        # Check if highlighted
        highlight = self.data.is_highlighted(gene1)
        if highlight:
            color = highlight.color
            alpha = self.style.highlight_link_alpha
        else:
            color = self.style.link_color
            alpha = self.style.link_alpha

        polygon = Polygon(vertices, closed=True, facecolor=color, alpha=alpha, edgecolor="none")
        self._ax.add_patch(polygon)

    def _draw_species(self) -> None:
        """Draw all species with their chromosomes and genes."""
        for sp in self.data.species:
            y = self._species_y[sp.species_id]
            self._draw_chromosome_line(sp, y)
            self._draw_genes(sp, y)
            self._draw_species_label(sp, y)

    def _draw_chromosome_line(self, species: Species, y: float) -> None:
        """Draw the chromosome backbone line."""
        if not species.genes:
            return

        x_start = species.region_start * self._x_scale + self._x_offset
        x_end = species.region_end * self._x_scale + self._x_offset

        # Add some padding
        padding = 0.02
        self._ax.plot(
            [x_start - padding, x_end + padding],
            [y, y],
            color=self.style.chromosome_color,
            linewidth=self.style.chromosome_width,
            solid_capstyle="round",
        )

    def _draw_genes(self, species: Species, y: float) -> None:
        """Draw gene rectangles for a species."""
        half_h = self.style.gene_height / 2

        # Collect gene names for bottom labels
        gene_labels = []

        for gene in species.genes:
            x = gene.start * self._x_scale + self._x_offset
            width = gene.length * self._x_scale

            # Check highlight
            highlight = self.data.is_highlighted(gene)
            if highlight:
                facecolor = highlight.color
            else:
                facecolor = self.style.gene_color

            rect = FancyBboxPatch(
                (x, y - half_h),
                width,
                self.style.gene_height,
                boxstyle="round,pad=0,rounding_size=0.02",
                facecolor=facecolor,
                edgecolor=self.style.gene_border_color,
                linewidth=0.5,
            )
            self._ax.add_patch(rect)

            gene_labels.append((gene.center * self._x_scale + self._x_offset, gene.name))

        # Draw gene labels at bottom for last species
        if species == self.data.species[-1] and self.style.show_gene_labels:
            for x_center, name in gene_labels:
                self._ax.text(
                    x_center,
                    y - half_h - 0.15,
                    name,
                    ha="center",
                    va="top",
                    fontsize=self.style.gene_label_font_size,
                    rotation=45,
                    style="italic",
                )

    def _draw_species_label(self, species: Species, y: float) -> None:
        """Draw species name on the right."""
        self._ax.text(
            1.02,
            y,
            species.display_name,
            ha="left",
            va="center",
            fontsize=self.style.species_font_size,
        )

    def _draw_legend(self) -> None:
        """Draw legend for highlighted genes."""
        if not self.data.highlights:
            return

        legend_elements = []
        for h in self.data.highlights:
            label = h.label or f"Gene: {h.name}"
            patch = mpatches.Patch(color=h.color, label=label)
            legend_elements.append(patch)

        # Add "other" entry
        other_patch = mpatches.Patch(color=self.style.gene_color, label="other")
        legend_elements.append(other_patch)

        self._ax.legend(
            handles=legend_elements,
            loc="upper right",
            frameon=False,
            fontsize=self.style.gene_label_font_size,
        )

    def _finalize(self) -> None:
        """Final plot adjustments."""
        if self.data.title:
            self._ax.set_title(self.data.title, fontsize=12, pad=20)
        plt.tight_layout()

    @staticmethod
    def _nice_scale_value(value: float) -> int:
        """Round to a nice scale bar value."""
        if value <= 0:
            return 100000

        # Find order of magnitude
        import math

        magnitude = 10 ** math.floor(math.log10(value))
        normalized = value / magnitude

        if normalized <= 1:
            nice = 1
        elif normalized <= 2:
            nice = 2
        elif normalized <= 5:
            nice = 5
        else:
            nice = 10

        return int(nice * magnitude)

    @staticmethod
    def _format_bp(bp: int) -> str:
        """Format base pairs for display."""
        if bp >= 1_000_000:
            return f"{bp / 1_000_000:.1f} Mbp"
        elif bp >= 1_000:
            return f"{bp / 1_000:.0f} Kbp"
        else:
            return f"{bp} bp"


def render_synteny(
    data: SyntenyData | str | Path,
    output: str | Path | None = None,
) -> plt.Figure:
    """Convenience function to render synteny plot.

    Args:
        data: SyntenyData object or path to JSON file
        output: Optional output file path

    Returns:
        matplotlib Figure
    """
    if isinstance(data, (str, Path)):
        data = SyntenyData.from_json_file(data)

    viz = SyntenyVisualizer(data)
    return viz.render(output)
