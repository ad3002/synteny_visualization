"""GFF/GTF format converter."""

import re
from pathlib import Path

from ..models import Gene, Strand


def parse_gff_attributes(attr_string: str, is_gtf: bool = False) -> dict[str, str]:
    """Parse GFF3 or GTF attribute string.

    GFF3: ID=gene001;Name=BRCA1;biotype=protein_coding
    GTF:  gene_id "BRCA1"; gene_name "BRCA1"; gene_biotype "protein_coding";
    """
    attrs = {}

    if is_gtf:
        # GTF format: key "value";
        pattern = r'(\w+)\s+"([^"]+)"'
        for match in re.finditer(pattern, attr_string):
            attrs[match.group(1)] = match.group(2)
    else:
        # GFF3 format: key=value;
        for pair in attr_string.split(";"):
            pair = pair.strip()
            if "=" in pair:
                key, value = pair.split("=", 1)
                attrs[key] = value

    return attrs


def load_genes_from_gff(
    path: str | Path,
    species_id: str | None = None,
    feature_type: str = "gene",
    name_attr: str | None = None,
    chrom_filter: str | None = None,
    start_filter: int | None = None,
    end_filter: int | None = None,
) -> list[Gene]:
    """Load genes from a GFF3 or GTF file.

    GFF format (tab-separated):
        seqid, source, type, start, end, score, strand, phase, attributes

    Args:
        path: Path to GFF/GTF file
        species_id: Optional species ID to prefix gene IDs
        feature_type: Feature type to extract (default: "gene")
        name_attr: Attribute to use as gene name. If None, tries:
                   GFF3: Name, gene_name, ID
                   GTF: gene_name, gene_id
        chrom_filter: Only include genes from this chromosome
        start_filter: Only include genes after this position
        end_filter: Only include genes before this position

    Returns:
        List of Gene objects
    """
    genes = []
    path = Path(path)

    # Detect format from extension
    is_gtf = path.suffix.lower() in (".gtf", ".gtf.gz")

    # Priority order for name attribute
    if name_attr:
        name_priority = [name_attr]
    elif is_gtf:
        name_priority = ["gene_name", "gene_id"]
    else:
        name_priority = ["Name", "gene_name", "ID"]

    with open(path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            seqid = parts[0]
            feat_type = parts[2]
            start = int(parts[3]) - 1  # GFF is 1-based, convert to 0-based
            end = int(parts[4])
            strand_char = parts[6]
            attributes = parts[8]

            # Filter by feature type
            if feat_type != feature_type:
                continue

            # Apply filters
            if chrom_filter and seqid != chrom_filter:
                continue
            if start_filter is not None and end < start_filter:
                continue
            if end_filter is not None and start > end_filter:
                continue

            # Parse attributes
            attrs = parse_gff_attributes(attributes, is_gtf)

            # Get name
            name = None
            for attr in name_priority:
                if attr in attrs:
                    name = attrs[attr]
                    break
            if not name:
                name = f"gene_{line_num}"

            # Get strand
            strand = Strand.UNKNOWN
            if strand_char == "+":
                strand = Strand.PLUS
            elif strand_char == "-":
                strand = Strand.MINUS

            # Create gene ID
            gene_id = attrs.get("ID", attrs.get("gene_id", f"gene_{line_num}"))
            if species_id:
                gene_id = f"{species_id}:{gene_id}"

            gene = Gene(
                gene_id=gene_id,
                name=name,
                chrom=seqid,
                start=start,
                end=end,
                strand=strand,
            )
            genes.append(gene)

    return genes
