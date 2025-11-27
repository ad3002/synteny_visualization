"""Synteny rearrangement analysis and visualization."""

from dataclasses import dataclass
from enum import Enum

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch
import numpy as np

from .models import Species, SyntenyData, Gene


class RearrangementType(str, Enum):
    CONSERVED = "conserved"  # Same order, same orientation
    INVERSION = "inversion"  # Same genes, reversed order
    TRANSLOCATION = "translocation"  # Gene moved to different position
    INSERTION = "insertion"  # Gene present only in one species
    DELETION = "deletion"  # Gene absent in one species


@dataclass
class SyntenyBlock:
    """A block of genes with conserved or rearranged synteny."""
    orthogroups: list[str]
    start_idx_sp1: int
    end_idx_sp1: int
    start_idx_sp2: int
    end_idx_sp2: int
    rearrangement_type: RearrangementType


@dataclass
class GeneRearrangement:
    """Rearrangement status for a single gene/orthogroup."""
    orthogroup: str
    gene_sp1: Gene | None
    gene_sp2: Gene | None
    position_sp1: int  # Position index in species 1 (-1 if absent)
    position_sp2: int  # Position index in species 2 (-1 if absent)
    rearrangement_type: RearrangementType


def get_ordered_orthogroups(species: Species) -> list[tuple[str, Gene]]:
    """Get orthogroups ordered by gene position."""
    genes_with_og = [(g.orthogroup, g) for g in species.genes if g.orthogroup]
    return sorted(genes_with_og, key=lambda x: x[1].center)


def analyze_pairwise_rearrangements(
    sp1: Species,
    sp2: Species
) -> list[GeneRearrangement]:
    """Analyze rearrangements between two species."""

    # Get ordered genes
    og_order1 = get_ordered_orthogroups(sp1)
    og_order2 = get_ordered_orthogroups(sp2)

    # Create position maps
    pos1 = {og: i for i, (og, _) in enumerate(og_order1)}
    pos2 = {og: i for i, (og, _) in enumerate(og_order2)}

    gene_map1 = {og: g for og, g in og_order1}
    gene_map2 = {og: g for og, g in og_order2}

    # All orthogroups
    all_ogs = set(pos1.keys()) | set(pos2.keys())

    rearrangements = []

    for og in all_ogs:
        p1 = pos1.get(og, -1)
        p2 = pos2.get(og, -1)
        g1 = gene_map1.get(og)
        g2 = gene_map2.get(og)

        if p1 == -1:
            rtype = RearrangementType.INSERTION
        elif p2 == -1:
            rtype = RearrangementType.DELETION
        else:
            # Check if position is conserved relative to neighbors
            rtype = classify_gene_rearrangement(og, pos1, pos2, og_order1, og_order2)

        rearrangements.append(GeneRearrangement(
            orthogroup=og,
            gene_sp1=g1,
            gene_sp2=g2,
            position_sp1=p1,
            position_sp2=p2,
            rearrangement_type=rtype,
        ))

    return rearrangements


def classify_gene_rearrangement(
    og: str,
    pos1: dict[str, int],
    pos2: dict[str, int],
    og_order1: list[tuple[str, Gene]],
    og_order2: list[tuple[str, Gene]],
) -> RearrangementType:
    """Classify rearrangement type for a gene based on neighbor conservation."""

    p1 = pos1[og]
    p2 = pos2[og]

    # Get neighbors in species 1
    left1 = og_order1[p1 - 1][0] if p1 > 0 else None
    right1 = og_order1[p1 + 1][0] if p1 < len(og_order1) - 1 else None

    # Check if neighbors are conserved in species 2
    left1_p2 = pos2.get(left1, -1) if left1 else -1
    right1_p2 = pos2.get(right1, -1) if right1 else -1

    # Check conservation patterns
    neighbors_conserved = 0
    inversion_signal = 0

    if left1_p2 != -1:
        if left1_p2 == p2 - 1:
            neighbors_conserved += 1
        elif left1_p2 == p2 + 1:
            inversion_signal += 1

    if right1_p2 != -1:
        if right1_p2 == p2 + 1:
            neighbors_conserved += 1
        elif right1_p2 == p2 - 1:
            inversion_signal += 1

    if neighbors_conserved >= 1:
        return RearrangementType.CONSERVED
    elif inversion_signal >= 1:
        return RearrangementType.INVERSION
    else:
        return RearrangementType.TRANSLOCATION


def detect_synteny_blocks(
    sp1: Species,
    sp2: Species,
    min_block_size: int = 2,
) -> list[SyntenyBlock]:
    """Detect blocks of conserved/rearranged synteny."""

    rearrangements = analyze_pairwise_rearrangements(sp1, sp2)

    # Sort by position in sp1
    shared = [r for r in rearrangements if r.position_sp1 >= 0 and r.position_sp2 >= 0]
    shared.sort(key=lambda r: r.position_sp1)

    if not shared:
        return []

    blocks = []
    current_block_ogs = [shared[0].orthogroup]
    current_type = shared[0].rearrangement_type
    block_start_p1 = shared[0].position_sp1
    block_start_p2 = shared[0].position_sp2

    for i in range(1, len(shared)):
        r = shared[i]
        prev = shared[i-1]

        # Check if this gene continues the current block
        same_type = r.rearrangement_type == current_type
        consecutive_p1 = r.position_sp1 == prev.position_sp1 + 1

        if same_type and consecutive_p1:
            current_block_ogs.append(r.orthogroup)
        else:
            # Save current block if large enough
            if len(current_block_ogs) >= min_block_size:
                blocks.append(SyntenyBlock(
                    orthogroups=current_block_ogs,
                    start_idx_sp1=block_start_p1,
                    end_idx_sp1=prev.position_sp1,
                    start_idx_sp2=block_start_p2,
                    end_idx_sp2=prev.position_sp2,
                    rearrangement_type=current_type,
                ))

            # Start new block
            current_block_ogs = [r.orthogroup]
            current_type = r.rearrangement_type
            block_start_p1 = r.position_sp1
            block_start_p2 = r.position_sp2

    # Don't forget last block
    if len(current_block_ogs) >= min_block_size:
        blocks.append(SyntenyBlock(
            orthogroups=current_block_ogs,
            start_idx_sp1=block_start_p1,
            end_idx_sp1=shared[-1].position_sp1,
            start_idx_sp2=block_start_p2,
            end_idx_sp2=shared[-1].position_sp2,
            rearrangement_type=current_type,
        ))

    return blocks


# Colors for rearrangement types
REARRANGEMENT_COLORS = {
    RearrangementType.CONSERVED: "#27ae60",      # Green
    RearrangementType.INVERSION: "#e74c3c",      # Red
    RearrangementType.TRANSLOCATION: "#f39c12",  # Orange
    RearrangementType.INSERTION: "#3498db",      # Blue
    RearrangementType.DELETION: "#9b59b6",       # Purple
}


def render_rearrangement_plot(
    data: SyntenyData,
    output_path: str | None = None,
    reference_species_idx: int = 0,
) -> plt.Figure:
    """Render a plot showing rearrangements relative to a reference species.

    Shows each species as a horizontal track with genes colored by their
    rearrangement status relative to the reference.
    """

    if len(data.species) < 2:
        raise ValueError("Need at least 2 species for rearrangement analysis")

    ref_species = data.species[reference_species_idx]

    # Analyze rearrangements for all species vs reference
    all_rearrangements: dict[str, list[GeneRearrangement]] = {}
    for i, sp in enumerate(data.species):
        if i != reference_species_idx:
            all_rearrangements[sp.species_id] = analyze_pairwise_rearrangements(ref_species, sp)

    # Create figure
    n_species = len(data.species)
    fig_height = n_species * 0.8 + 2
    fig, ax = plt.subplots(figsize=(16, fig_height))

    # Get reference gene order for x-axis
    ref_og_order = get_ordered_orthogroups(ref_species)
    ref_positions = {og: i for i, (og, _) in enumerate(ref_og_order)}
    n_genes = len(ref_og_order)

    # Draw each species track
    y_positions = {}
    for i, sp in enumerate(data.species):
        y = -i * 1.0
        y_positions[sp.species_id] = y

        # Draw track background
        ax.axhline(y=y, color='#cccccc', linewidth=0.5, zorder=0)

        # Species label
        ax.text(-1.5, y, sp.display_name, ha='right', va='center',
                fontsize=9, style='italic')

        if i == reference_species_idx:
            # Reference species - just show genes in gray
            for og, gene in ref_og_order:
                x = ref_positions[og]
                rect = FancyBboxPatch(
                    (x - 0.4, y - 0.2), 0.8, 0.4,
                    boxstyle="round,pad=0,rounding_size=0.1",
                    facecolor='#808080',
                    edgecolor='#404040',
                    linewidth=0.5,
                    zorder=2,
                )
                ax.add_patch(rect)
        else:
            # Other species - color by rearrangement type
            rearrangements = all_rearrangements[sp.species_id]
            rearr_map = {r.orthogroup: r for r in rearrangements}

            # Get this species' gene order
            sp_og_order = get_ordered_orthogroups(sp)

            for og, gene in sp_og_order:
                if og in ref_positions:
                    # Position on x-axis based on reference
                    x = ref_positions[og]

                    rearr = rearr_map.get(og)
                    if rearr:
                        color = REARRANGEMENT_COLORS.get(
                            rearr.rearrangement_type, '#808080'
                        )
                    else:
                        color = '#808080'

                    rect = FancyBboxPatch(
                        (x - 0.4, y - 0.2), 0.8, 0.4,
                        boxstyle="round,pad=0,rounding_size=0.1",
                        facecolor=color,
                        edgecolor='#404040',
                        linewidth=0.5,
                        zorder=2,
                    )
                    ax.add_patch(rect)

                    # Draw inversion arrow for inverted genes
                    if rearr and rearr.rearrangement_type == RearrangementType.INVERSION:
                        ax.annotate('', xy=(x + 0.3, y), xytext=(x - 0.3, y),
                                   arrowprops=dict(arrowstyle='->', color='white', lw=1.5),
                                   zorder=3)

    # Gene labels at bottom
    for og, gene in ref_og_order:
        x = ref_positions[og]
        ax.text(x, -n_species * 1.0 - 0.5, gene.name,
               ha='center', va='top', fontsize=7, rotation=45, style='italic')

    # Legend
    legend_elements = [
        mpatches.Patch(color=REARRANGEMENT_COLORS[RearrangementType.CONSERVED],
                      label='Conserved'),
        mpatches.Patch(color=REARRANGEMENT_COLORS[RearrangementType.INVERSION],
                      label='Inversion'),
        mpatches.Patch(color=REARRANGEMENT_COLORS[RearrangementType.TRANSLOCATION],
                      label='Translocation'),
        mpatches.Patch(color='#808080', label='Reference'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=False, fontsize=8)

    # Title
    ax.set_title(f'Rearrangements relative to {ref_species.display_name}',
                fontsize=11, pad=10)

    # Axis settings
    ax.set_xlim(-2, n_genes + 1)
    ax.set_ylim(-n_species * 1.0 - 1.5, 1)
    ax.axis('off')

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight',
                   facecolor='white', edgecolor='none')

    return fig


def render_breakpoint_plot(
    data: SyntenyData,
    output_path: str | None = None,
) -> plt.Figure:
    """Render a dot plot showing synteny breakpoints between adjacent species."""

    n_comparisons = len(data.species) - 1
    if n_comparisons < 1:
        raise ValueError("Need at least 2 species")

    fig, axes = plt.subplots(1, n_comparisons, figsize=(5 * n_comparisons, 5))
    if n_comparisons == 1:
        axes = [axes]

    for i, ax in enumerate(axes):
        sp1 = data.species[i]
        sp2 = data.species[i + 1]

        # Get gene orders
        og_order1 = get_ordered_orthogroups(sp1)
        og_order2 = get_ordered_orthogroups(sp2)

        pos1 = {og: idx for idx, (og, _) in enumerate(og_order1)}
        pos2 = {og: idx for idx, (og, _) in enumerate(og_order2)}

        # Find shared orthogroups
        shared = set(pos1.keys()) & set(pos2.keys())

        # Analyze rearrangements
        rearrangements = analyze_pairwise_rearrangements(sp1, sp2)
        rearr_map = {r.orthogroup: r for r in rearrangements}

        # Plot dots
        for og in shared:
            x = pos1[og]
            y = pos2[og]
            rearr = rearr_map.get(og)
            if rearr:
                color = REARRANGEMENT_COLORS.get(rearr.rearrangement_type, '#808080')
            else:
                color = '#808080'
            ax.scatter(x, y, c=color, s=50, edgecolors='white', linewidths=0.5, zorder=2)

        # Draw diagonal for reference
        max_pos = max(len(og_order1), len(og_order2))
        ax.plot([0, max_pos], [0, max_pos], 'k--', alpha=0.3, linewidth=1, zorder=1)

        ax.set_xlabel(sp1.display_name, fontsize=9, style='italic')
        ax.set_ylabel(sp2.display_name, fontsize=9, style='italic')
        ax.set_aspect('equal')
        ax.set_title(f'{sp1.display_name} vs {sp2.display_name}', fontsize=10)

    # Shared legend
    legend_elements = [
        mpatches.Patch(color=REARRANGEMENT_COLORS[RearrangementType.CONSERVED],
                      label='Conserved'),
        mpatches.Patch(color=REARRANGEMENT_COLORS[RearrangementType.INVERSION],
                      label='Inversion'),
        mpatches.Patch(color=REARRANGEMENT_COLORS[RearrangementType.TRANSLOCATION],
                      label='Translocation'),
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=3,
              frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.02))

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight',
                   facecolor='white', edgecolor='none')

    return fig


def summarize_rearrangements(data: SyntenyData) -> dict:
    """Generate summary statistics of rearrangements across all species pairs."""

    summary = {
        'pairwise': [],
        'total_conserved': 0,
        'total_inversions': 0,
        'total_translocations': 0,
    }

    for i in range(len(data.species) - 1):
        sp1 = data.species[i]
        sp2 = data.species[i + 1]

        rearrangements = analyze_pairwise_rearrangements(sp1, sp2)

        counts = {
            RearrangementType.CONSERVED: 0,
            RearrangementType.INVERSION: 0,
            RearrangementType.TRANSLOCATION: 0,
            RearrangementType.INSERTION: 0,
            RearrangementType.DELETION: 0,
        }

        for r in rearrangements:
            counts[r.rearrangement_type] += 1

        summary['pairwise'].append({
            'species1': sp1.display_name,
            'species2': sp2.display_name,
            'conserved': counts[RearrangementType.CONSERVED],
            'inversions': counts[RearrangementType.INVERSION],
            'translocations': counts[RearrangementType.TRANSLOCATION],
            'insertions': counts[RearrangementType.INSERTION],
            'deletions': counts[RearrangementType.DELETION],
        })

        summary['total_conserved'] += counts[RearrangementType.CONSERVED]
        summary['total_inversions'] += counts[RearrangementType.INVERSION]
        summary['total_translocations'] += counts[RearrangementType.TRANSLOCATION]

    return summary
