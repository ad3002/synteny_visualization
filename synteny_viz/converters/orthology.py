"""Orthology data converter."""

from pathlib import Path

from ..models import Gene, Species


def load_orthogroups(
    path: str | Path,
    gene_col: int = 0,
    orthogroup_col: int = 1,
    species_col: int | None = None,
    delimiter: str = "\t",
    skip_header: bool = True,
) -> dict[str, str]:
    """Load orthogroup assignments from a tabular file.

    Expected format (tab-separated):
        gene_name    orthogroup_id    [species]

    Or OrthoFinder format:
        Orthogroup    Gene1, Gene2, Gene3...

    Args:
        path: Path to orthogroup file
        gene_col: Column index for gene name
        orthogroup_col: Column index for orthogroup ID
        species_col: Optional column for species (not used, for filtering)
        delimiter: Column delimiter
        skip_header: Skip first line

    Returns:
        Dict mapping gene_name -> orthogroup_id
    """
    gene_to_og: dict[str, str] = {}
    path = Path(path)

    with open(path) as f:
        for i, line in enumerate(f):
            if skip_header and i == 0:
                continue

            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split(delimiter)

            # OrthoFinder format: Orthogroup    Gene1, Gene2, Gene3
            if len(parts) == 2 and "," in parts[1]:
                og_id = parts[0]
                genes = [g.strip() for g in parts[1].split(",")]
                for gene in genes:
                    if gene:
                        gene_to_og[gene] = og_id
            # Simple format: gene    orthogroup
            elif len(parts) > max(gene_col, orthogroup_col):
                gene = parts[gene_col].strip()
                og_id = parts[orthogroup_col].strip()
                if gene and og_id:
                    gene_to_og[gene] = og_id

    return gene_to_og


def merge_orthogroups_into_species(
    species_list: list[Species],
    orthogroups: dict[str, str],
    match_by: str = "name",
) -> None:
    """Assign orthogroup IDs to genes in species.

    Modifies Species objects in place.

    Args:
        species_list: List of Species objects
        orthogroups: Dict mapping gene identifier -> orthogroup_id
        match_by: Which gene attribute to match ("name" or "gene_id")
    """
    for species in species_list:
        for gene in species.genes:
            key = gene.name if match_by == "name" else gene.gene_id
            if key in orthogroups:
                # Create new gene with orthogroup (genes are immutable after creation)
                gene.orthogroup = orthogroups[key]


def load_ensembl_orthologs(
    path: str | Path,
    source_species: str,
    target_species: str,
    delimiter: str = "\t",
) -> dict[str, str]:
    """Load orthology from Ensembl BioMart export.

    Expected columns:
        Gene stable ID, Gene name, [Target] Gene stable ID, [Target] Gene name

    Args:
        path: Path to Ensembl export file
        source_species: Source species ID (for orthogroup naming)
        target_species: Target species ID
        delimiter: Column delimiter

    Returns:
        Dict mapping gene_name -> orthogroup_id (using source gene as group name)
    """
    gene_to_og: dict[str, str] = {}
    path = Path(path)

    with open(path) as f:
        header = next(f)  # Skip header

        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split(delimiter)
            if len(parts) < 4:
                continue

            source_id = parts[0]
            source_name = parts[1]
            target_id = parts[2]
            target_name = parts[3]

            if not source_name or not target_name:
                continue

            # Use source gene name as orthogroup ID
            og_id = source_name

            gene_to_og[source_name] = og_id
            gene_to_og[target_name] = og_id

    return gene_to_og
