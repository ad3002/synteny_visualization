"""Converter for CSV synteny format (attributes + region files)."""

import csv
from pathlib import Path

from ..models import Gene, Species, SyntenyData, HighlightGene, StyleConfig


def load_orthogroups_from_attributes(path: str | Path) -> dict[str, str]:
    """Load orthogroups from synteny_attributes.csv.

    Format: gene_ID, group, strand, genome
    Returns: dict mapping gene_ID -> group
    """
    gene_to_group: dict[str, str] = {}

    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id = row['gene_ID'].strip().strip('"')
            group = row['group'].strip().strip('"')
            gene_to_group[gene_id] = group

    return gene_to_group


def load_species_from_region_csv(
    path: str | Path,
    species_id: str,
    display_name: str,
    orthogroups: dict[str, str] | None = None,
) -> Species:
    """Load species from synteny_*_region.csv.

    Format: Chr, start, end, gene_ID
    """
    genes = []

    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_name = row['gene_ID'].strip()
            chrom = str(row['Chr']).strip()
            start = int(row['start'])
            end = int(row['end'])

            # Get orthogroup
            orthogroup = None
            if orthogroups:
                orthogroup = orthogroups.get(gene_name)

            gene = Gene(
                gene_id=f"{species_id}:{gene_name}",
                name=gene_name,
                chrom=chrom,
                start=start,
                end=end,
                orthogroup=orthogroup,
            )
            genes.append(gene)

    return Species(
        species_id=species_id,
        display_name=display_name,
        genes=genes,
    )


def flip_species_coordinates(species: Species) -> Species:
    """Flip/mirror gene coordinates (multiply by -1).

    Used to minimize crossings when synteny is inverted.
    """
    if not species.genes:
        return species

    flipped_genes = []
    for gene in species.genes:
        flipped_genes.append(Gene(
            gene_id=gene.gene_id,
            name=gene.name,
            chrom=gene.chrom,
            start=-gene.end,  # Swap and negate
            end=-gene.start,
            strand=gene.strand,
            orthogroup=gene.orthogroup,
        ))

    return Species(
        species_id=species.species_id,
        display_name=species.display_name,
        genes=flipped_genes,
    )


def count_crossings(sp1: Species, sp2: Species) -> int:
    """Count orthology link crossings between two species.

    A crossing occurs when two links intersect.
    """
    # Build list of ortholog pairs with their positions
    pairs = []
    for g1 in sp1.genes:
        if not g1.orthogroup:
            continue
        for g2 in sp2.genes:
            if g2.orthogroup == g1.orthogroup:
                pairs.append((g1.center, g2.center))
                break

    # Count crossings
    crossings = 0
    for i, (x1a, x1b) in enumerate(pairs):
        for x2a, x2b in pairs[i + 1:]:
            # Crossing if one pair "crosses over" the other
            if (x1a < x2a and x1b > x2b) or (x1a > x2a and x1b < x2b):
                crossings += 1

    return crossings


def normalize_species_coordinates(
    species: Species,
    anchor_orthogroup: str | None = None,
    max_distance: int | None = None,
    compress_gaps: int | None = 500000,  # Compress gaps larger than this to fixed size
    compressed_gap_size: int = 100000,  # Size of compressed gap
) -> Species:
    """Normalize gene coordinates to align species.

    If anchor_orthogroup is provided, aligns by that gene's center.
    Otherwise, aligns by region center.
    If max_distance is set, filters out genes beyond that distance from anchor.
    If compress_gaps is set, large gaps are compressed to compressed_gap_size.
    """
    if not species.genes:
        return species

    # Sort genes by position first
    sorted_genes = sorted(species.genes, key=lambda g: g.start)

    # Compress large gaps if enabled
    if compress_gaps is not None and len(sorted_genes) > 1:
        compressed_genes = [sorted_genes[0]]
        offset = 0  # Cumulative compression offset

        for i in range(1, len(sorted_genes)):
            prev_gene = sorted_genes[i - 1]
            curr_gene = sorted_genes[i]
            gap = curr_gene.start - prev_gene.end

            if gap > compress_gaps:
                # Compress this gap
                compression = gap - compressed_gap_size
                offset += compression

            # Create gene with adjusted coordinates
            compressed_genes.append(Gene(
                gene_id=curr_gene.gene_id,
                name=curr_gene.name,
                chrom=curr_gene.chrom,
                start=curr_gene.start - offset,
                end=curr_gene.end - offset,
                strand=curr_gene.strand,
                orthogroup=curr_gene.orthogroup,
            ))

        sorted_genes = compressed_genes

    # Find anchor point
    anchor_center = None
    if anchor_orthogroup:
        for gene in sorted_genes:
            if gene.orthogroup == anchor_orthogroup:
                anchor_center = (gene.start + gene.end) / 2
                break

    # Fall back to region center if anchor not found
    if anchor_center is None:
        region_start = min(g.start for g in sorted_genes)
        region_end = max(g.end for g in sorted_genes)
        anchor_center = (region_start + region_end) / 2

    # Normalize all genes relative to anchor
    normalized_genes = []
    for gene in sorted_genes:
        normalized_start = int(gene.start - anchor_center)
        normalized_end = int(gene.end - anchor_center)

        # Filter by distance if specified
        if max_distance is not None:
            gene_center = (normalized_start + normalized_end) / 2
            if abs(gene_center) > max_distance:
                continue

        normalized_genes.append(Gene(
            gene_id=gene.gene_id,
            name=gene.name,
            chrom=gene.chrom,
            start=normalized_start,
            end=normalized_end,
            strand=gene.strand,
            orthogroup=gene.orthogroup,
        ))

    return Species(
        species_id=species.species_id,
        display_name=species.display_name,
        genes=normalized_genes,
    )


def get_gene_order(species: Species) -> dict[str, int]:
    """Get gene order by orthogroup (position index in sorted order)."""
    sorted_genes = sorted(
        [g for g in species.genes if g.orthogroup],
        key=lambda g: g.center
    )
    return {g.orthogroup: i for i, g in enumerate(sorted_genes)}


def calculate_synteny_similarity(sp1: Species, sp2: Species) -> float:
    """Calculate synteny similarity based on gene order correlation.

    Uses Spearman-like rank correlation: how similar is the order of shared genes.
    Returns value between -1 (inverted) and 1 (identical order).
    """
    order1 = get_gene_order(sp1)
    order2 = get_gene_order(sp2)

    # Find shared orthogroups
    shared = set(order1.keys()) & set(order2.keys())
    if len(shared) < 2:
        return 0.0

    # Get ranks for shared genes
    ranks1 = [order1[og] for og in shared]
    ranks2 = [order2[og] for og in shared]

    # Calculate Spearman correlation
    n = len(shared)
    mean1 = sum(ranks1) / n
    mean2 = sum(ranks2) / n

    numerator = sum((r1 - mean1) * (r2 - mean2) for r1, r2 in zip(ranks1, ranks2))
    denom1 = sum((r - mean1) ** 2 for r in ranks1) ** 0.5
    denom2 = sum((r - mean2) ** 2 for r in ranks2) ** 0.5

    if denom1 == 0 or denom2 == 0:
        return 0.0

    return numerator / (denom1 * denom2)


def build_distance_matrix(species_list: list[Species]) -> list[list[float]]:
    """Build pairwise distance matrix based on synteny dissimilarity."""
    n = len(species_list)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            # Similarity can be negative (inverted), so use abs for distance
            sim = calculate_synteny_similarity(species_list[i], species_list[j])
            # Convert to distance: 1 - |similarity| (0 = identical, 1 = unrelated)
            dist = 1.0 - abs(sim)
            matrix[i][j] = dist
            matrix[j][i] = dist

    return matrix


def hierarchical_cluster_species(
    species_list: list[Species],
) -> list[Species]:
    """Reorder species using hierarchical clustering by synteny similarity.

    Groups similar species together to minimize visual rearrangements.
    Uses simple agglomerative clustering with average linkage.
    """
    if len(species_list) <= 2:
        return species_list

    n = len(species_list)
    dist_matrix = build_distance_matrix(species_list)

    # Track which species belong to which cluster
    # Initially each species is its own cluster
    clusters: list[list[int]] = [[i] for i in range(n)]
    active = set(range(n))  # Active cluster indices

    # Agglomerative clustering
    while len(active) > 1:
        # Find closest pair of clusters
        min_dist = float('inf')
        merge_i, merge_j = -1, -1

        active_list = sorted(active)
        for i_idx, i in enumerate(active_list):
            for j in active_list[i_idx + 1:]:
                # Average linkage: mean distance between all pairs
                total_dist = 0.0
                count = 0
                for si in clusters[i]:
                    for sj in clusters[j]:
                        total_dist += dist_matrix[si][sj]
                        count += 1
                avg_dist = total_dist / count if count > 0 else float('inf')

                if avg_dist < min_dist:
                    min_dist = avg_dist
                    merge_i, merge_j = i, j

        # Merge clusters
        clusters[merge_i] = clusters[merge_i] + clusters[merge_j]
        active.remove(merge_j)

    # Get final order from the single remaining cluster
    final_cluster_idx = list(active)[0]
    ordered_indices = clusters[final_cluster_idx]

    return [species_list[i] for i in ordered_indices]


def optimize_species_order(
    species_list: list[Species],
    anchor_species: str | None = None,
) -> list[Species]:
    """Optimize species order to minimize total crossings.

    1. Hierarchically cluster species by synteny similarity
    2. Within clusters, order to minimize adjacent crossings
    3. Optionally anchor a specific species at top

    Args:
        species_list: List of species to reorder
        anchor_species: Species ID to keep at top (e.g., reference species)
    """
    if len(species_list) <= 1:
        return species_list

    # First, cluster hierarchically
    clustered = hierarchical_cluster_species(species_list)

    # If anchor species specified, rotate to put it first
    if anchor_species:
        anchor_idx = None
        for i, sp in enumerate(clustered):
            if sp.species_id == anchor_species:
                anchor_idx = i
                break
        if anchor_idx is not None:
            clustered = clustered[anchor_idx:] + clustered[:anchor_idx]

    # Now optimize local order: for each adjacent pair, check if swapping reduces crossings
    # Use greedy optimization with multiple passes
    improved = True
    max_passes = 10
    passes = 0

    while improved and passes < max_passes:
        improved = False
        passes += 1

        for i in range(1, len(clustered) - 1):
            # Try swapping species[i] with species[i+1]
            # Calculate total crossings before and after

            # Before: crossings(i-1, i) + crossings(i, i+1)
            cross_before = (
                count_crossings(clustered[i - 1], clustered[i]) +
                count_crossings(clustered[i], clustered[i + 1])
            )

            # After swap: crossings(i-1, i+1) + crossings(i+1, i)
            cross_after = (
                count_crossings(clustered[i - 1], clustered[i + 1]) +
                count_crossings(clustered[i + 1], clustered[i])
            )

            if cross_after < cross_before:
                clustered[i], clustered[i + 1] = clustered[i + 1], clustered[i]
                improved = True

    return clustered


def optimize_strand_orientation(species_list: list[Species]) -> list[Species]:
    """Optimize strand orientation to minimize crossings.

    For each species (except first), checks if flipping reduces crossings
    with the previous species. Propagates flips through the list.
    """
    if len(species_list) < 2:
        return species_list

    optimized = [species_list[0]]

    for i in range(1, len(species_list)):
        prev_sp = optimized[-1]
        curr_sp = species_list[i]

        # Count crossings with current orientation
        crossings_normal = count_crossings(prev_sp, curr_sp)

        # Count crossings with flipped orientation
        flipped_sp = flip_species_coordinates(curr_sp)
        crossings_flipped = count_crossings(prev_sp, flipped_sp)

        # Choose orientation with fewer crossings
        if crossings_flipped < crossings_normal:
            optimized.append(flipped_sp)
        else:
            optimized.append(curr_sp)

    return optimized


def get_conserved_orthogroups(species_list: list[Species]) -> set[str]:
    """Find orthogroups present in all species."""
    if not species_list:
        return set()

    # Get orthogroups for each species
    og_by_species = []
    for sp in species_list:
        og_set = {g.orthogroup for g in sp.genes if g.orthogroup}
        og_by_species.append(og_set)

    # Intersection of all
    conserved = og_by_species[0]
    for og_set in og_by_species[1:]:
        conserved = conserved & og_set

    return conserved


# Pleasant color palette for conserved genes (colorblind-friendly)
CONSERVED_PALETTE = [
    "#4e79a7",  # blue
    "#59a14f",  # green
    "#9c755f",  # brown
    "#76b7b2",  # teal
    "#edc949",  # yellow
    "#af7aa1",  # purple
    "#ff9da7",  # pink
    "#f28e2c",  # orange
    "#b07aa1",  # lavender
    "#86bcb6",  # cyan
]


def load_synteny_from_csv_folder(
    folder: str | Path,
    species_order: list[tuple[str, str, str]] | None = None,
    highlight_genes: list[str] | None = None,
    highlight_orthogroups: list[str] | None = None,
    normalize: bool = True,
    anchor_orthogroup: str | None = None,
    max_distance: int | None = None,
    optimize_strands: bool = True,
    optimize_order: bool = False,
    anchor_species: str | None = None,
    color_conserved: bool = False,
) -> SyntenyData:
    """Load complete synteny data from a folder with CSV files.

    Args:
        folder: Path to folder containing synteny_attributes.csv and synteny_*_region.csv
        species_order: List of (file_suffix, species_id, display_name) tuples.
                      If None, auto-detected from files.
        highlight_genes: Gene names to highlight (special color)
        highlight_orthogroups: Orthogroup IDs to highlight (special color)
        normalize: If True, normalize coordinates to align species
        anchor_orthogroup: Orthogroup to use as anchor for alignment
        max_distance: Maximum distance from anchor to include genes (bp)
        optimize_strands: If True, flip species to minimize link crossings
        optimize_order: If True, reorder species by synteny similarity (hierarchical clustering)
        anchor_species: Species ID to keep at top when optimizing order
        color_conserved: If True, color genes present in all species

    Returns:
        SyntenyData ready for visualization
    """
    folder = Path(folder)

    # Load orthogroups
    attributes_file = folder / "synteny_attributes.csv"
    orthogroups = {}
    if attributes_file.exists():
        orthogroups = load_orthogroups_from_attributes(attributes_file)

    # Default species order matching the reference image
    if species_order is None:
        species_order = [
            ("danio_rerio", "zebrafish", "Zebrafish"),
            ("xenopus_laevis_S", "xenopus_s", "Xenopus laevis S"),
            ("xenopus_laevis_L", "xenopus_l", "Xenopus laevis L"),
            ("xenopus_tropicales", "xenopus_t", "Xenopus tropicalis"),
            ("chiken", "chicken", "Chicken"),
            ("quail", "quail", "Quail"),
            ("zebra_finch", "zebra_finch", "Zebra finch"),
            ("canary", "canary", "Canary"),
            ("anolis", "anolis", "Anolis lizard"),
            ("mouse", "mouse", "Mouse"),
            ("human", "human", "Human"),
        ]

    # Load species
    species_list = []
    for file_suffix, species_id, display_name in species_order:
        region_file = folder / f"synteny_{file_suffix}_region.csv"
        if region_file.exists():
            species = load_species_from_region_csv(
                region_file, species_id, display_name, orthogroups
            )
            if normalize:
                species = normalize_species_coordinates(species, anchor_orthogroup, max_distance)
            species_list.append(species)

    # Optimize species order by synteny similarity (hierarchical clustering)
    if optimize_order and len(species_list) > 2:
        species_list = optimize_species_order(species_list, anchor_species)

    # Optimize strand orientations to minimize crossings
    if optimize_strands and len(species_list) > 1:
        species_list = optimize_strand_orientation(species_list)

    # Create highlights
    highlights = []

    # Color conserved genes if requested
    if color_conserved:
        conserved = get_conserved_orthogroups(species_list)
        # Get representative gene names for conserved orthogroups (from first species)
        og_to_name = {}
        for sp in species_list:
            for g in sp.genes:
                if g.orthogroup and g.orthogroup not in og_to_name:
                    og_to_name[g.orthogroup] = g.name

        # Sort orthogroups by position in first species for consistent ordering
        first_sp = species_list[0]
        og_positions = {}
        for g in first_sp.genes:
            if g.orthogroup:
                og_positions[g.orthogroup] = g.center

        sorted_conserved = sorted(conserved, key=lambda og: og_positions.get(og, 0))

        for i, og in enumerate(sorted_conserved):
            color = CONSERVED_PALETTE[i % len(CONSERVED_PALETTE)]
            gene_name = og_to_name.get(og, og)
            highlights.append(HighlightGene(
                name=gene_name,
                orthogroup=og,
                color=color,
                label=gene_name.upper() if gene_name else og
            ))

    # Add specific gene highlights (override conserved colors)
    if highlight_genes:
        for gene in highlight_genes:
            og = orthogroups.get(gene)
            # Remove existing highlight for this orthogroup
            highlights = [h for h in highlights if h.orthogroup != og]
            highlights.append(HighlightGene(
                name=gene,
                orthogroup=og,
                color="#e74c3c",
                label=f"Gene: {gene}"
            ))

    if highlight_orthogroups:
        for og in highlight_orthogroups:
            # Remove existing highlight for this orthogroup
            highlights = [h for h in highlights if h.orthogroup != og]
            highlights.append(HighlightGene(
                name=og,
                orthogroup=og,
                color="#e74c3c",
                label=f"Orthogroup: {og}"
            ))

    return SyntenyData(
        species=species_list,
        highlights=highlights,
        style=StyleConfig(
            figure_width=16,
            figure_height_per_species=0.7,
            gene_height=0.25,
        ),
    )
