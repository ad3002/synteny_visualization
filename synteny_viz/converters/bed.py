"""BED format converter."""

from pathlib import Path

from ..models import Gene, Strand


def load_genes_from_bed(
    path: str | Path,
    species_id: str | None = None,
    name_column: int = 3,
    chrom_filter: str | None = None,
    start_filter: int | None = None,
    end_filter: int | None = None,
) -> list[Gene]:
    """Load genes from a BED file.

    Standard BED format:
        chrom, start, end, name, score, strand

    Extended BED (BED12) also supported - extra columns ignored.

    Args:
        path: Path to BED file
        species_id: Optional species ID to prefix gene IDs
        name_column: Column index for gene name (0-based), default 3
        chrom_filter: Only include genes from this chromosome
        start_filter: Only include genes after this position
        end_filter: Only include genes before this position

    Returns:
        List of Gene objects
    """
    genes = []
    path = Path(path)

    with open(path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith("#") or line.startswith("track"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])

            # Apply filters
            if chrom_filter and chrom != chrom_filter:
                continue
            if start_filter is not None and end < start_filter:
                continue
            if end_filter is not None and start > end_filter:
                continue

            # Get name
            if len(parts) > name_column:
                name = parts[name_column]
            else:
                name = f"gene_{line_num}"

            # Get strand
            strand = Strand.UNKNOWN
            if len(parts) > 5:
                strand_char = parts[5]
                if strand_char == "+":
                    strand = Strand.PLUS
                elif strand_char == "-":
                    strand = Strand.MINUS

            # Create gene ID
            gene_id = f"{species_id}:{name}" if species_id else name

            gene = Gene(
                gene_id=gene_id,
                name=name,
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
            )
            genes.append(gene)

    return genes


def genes_to_bed(genes: list[Gene], output_path: str | Path) -> None:
    """Export genes to BED format.

    Args:
        genes: List of Gene objects
        output_path: Output file path
    """
    with open(output_path, "w") as f:
        for gene in genes:
            strand = gene.strand.value if gene.strand != Strand.UNKNOWN else "."
            line = f"{gene.chrom}\t{gene.start}\t{gene.end}\t{gene.name}\t0\t{strand}\n"
            f.write(line)
