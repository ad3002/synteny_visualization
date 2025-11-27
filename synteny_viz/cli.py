"""Command-line interface for synteny visualization."""

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.table import Table

from .models import SyntenyData
from .visualizer import render_synteny

app = typer.Typer(
    name="synteny-viz",
    help="Synteny visualization tool for comparative genomics",
    no_args_is_help=True,
)
console = Console()


@app.command()
def render(
    input_file: Annotated[
        Path, typer.Argument(help="Input JSON file with synteny data")
    ],
    output: Annotated[
        Optional[Path],
        typer.Option("-o", "--output", help="Output image file (PNG, PDF, SVG)"),
    ] = None,
    show: Annotated[
        bool, typer.Option("--show", help="Display plot interactively")
    ] = False,
    dpi: Annotated[
        int, typer.Option("--dpi", help="Output resolution")
    ] = 150,
    width: Annotated[
        float, typer.Option("--width", help="Figure width in inches")
    ] = 16,
) -> None:
    """Render synteny plot from JSON configuration."""
    import matplotlib.pyplot as plt

    if not input_file.exists():
        console.print(f"[red]Error:[/red] File not found: {input_file}")
        raise typer.Exit(1)

    try:
        data = SyntenyData.from_json_file(input_file)
    except Exception as e:
        console.print(f"[red]Error parsing JSON:[/red] {e}")
        raise typer.Exit(1)

    # Override style if specified
    if dpi != 150:
        data.style.dpi = dpi
    if width != 16:
        data.style.figure_width = width

    output_path = output or input_file.with_suffix(".png")

    console.print(f"Rendering {len(data.species)} species...")
    fig = render_synteny(data, output_path)
    console.print(f"[green]Saved to:[/green] {output_path}")

    if show:
        plt.show()


@app.command()
def validate(
    input_file: Annotated[
        Path, typer.Argument(help="Input JSON file to validate")
    ],
) -> None:
    """Validate a synteny JSON configuration file."""
    if not input_file.exists():
        console.print(f"[red]Error:[/red] File not found: {input_file}")
        raise typer.Exit(1)

    try:
        data = SyntenyData.from_json_file(input_file)
    except Exception as e:
        console.print(f"[red]Validation failed:[/red] {e}")
        raise typer.Exit(1)

    # Print summary
    table = Table(title="Synteny Data Summary")
    table.add_column("Property", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Title", data.title or "(none)")
    table.add_row("Species", str(len(data.species)))
    table.add_row("Highlights", str(len(data.highlights)))

    total_genes = sum(len(sp.genes) for sp in data.species)
    table.add_row("Total genes", str(total_genes))

    orthogroups = data.get_orthogroups()
    table.add_row("Orthogroups", str(len(orthogroups)))

    console.print(table)

    # Species details
    sp_table = Table(title="Species Details")
    sp_table.add_column("ID", style="cyan")
    sp_table.add_column("Display Name")
    sp_table.add_column("Genes", justify="right")
    sp_table.add_column("Region", style="dim")

    for sp in data.species:
        region = f"{sp.region_start:,}-{sp.region_end:,}" if sp.genes else "-"
        sp_table.add_row(sp.species_id, sp.display_name, str(len(sp.genes)), region)

    console.print(sp_table)
    console.print("[green]Validation passed![/green]")


@app.command()
def init(
    output: Annotated[
        Path, typer.Option("-o", "--output", help="Output JSON file")
    ] = Path("synteny.json"),
) -> None:
    """Create a template JSON configuration file."""
    from .models import Gene, HighlightGene, Species, StyleConfig

    # Create example data
    example = SyntenyData(
        title="Example Synteny Plot",
        species=[
            Species(
                species_id="species1",
                display_name="Species 1",
                genes=[
                    Gene(
                        gene_id="sp1:gene1",
                        name="GeneA",
                        chrom="chr1",
                        start=0,
                        end=50000,
                        orthogroup="OG001",
                    ),
                    Gene(
                        gene_id="sp1:gene2",
                        name="GeneB",
                        chrom="chr1",
                        start=100000,
                        end=150000,
                        orthogroup="OG002",
                    ),
                    Gene(
                        gene_id="sp1:target",
                        name="TargetGene",
                        chrom="chr1",
                        start=200000,
                        end=250000,
                        orthogroup="OG003",
                    ),
                ],
            ),
            Species(
                species_id="species2",
                display_name="Species 2",
                genes=[
                    Gene(
                        gene_id="sp2:gene1",
                        name="GeneA",
                        chrom="chrX",
                        start=10000,
                        end=60000,
                        orthogroup="OG001",
                    ),
                    Gene(
                        gene_id="sp2:gene2",
                        name="GeneB",
                        chrom="chrX",
                        start=120000,
                        end=170000,
                        orthogroup="OG002",
                    ),
                    Gene(
                        gene_id="sp2:target",
                        name="TargetGene",
                        chrom="chrX",
                        start=220000,
                        end=270000,
                        orthogroup="OG003",
                    ),
                ],
            ),
        ],
        highlights=[
            HighlightGene(name="TargetGene", color="#e74c3c", label="Gene: TargetGene"),
        ],
        style=StyleConfig(),
    )

    example.to_json_file(output)
    console.print(f"[green]Created template:[/green] {output}")
    console.print("Edit the file and run: synteny-viz render synteny.json")


@app.command("convert-bed")
def convert_bed(
    bed_files: Annotated[
        list[Path], typer.Argument(help="BED files to convert (format: species_id:path)")
    ],
    output: Annotated[
        Path, typer.Option("-o", "--output", help="Output JSON file")
    ] = Path("synteny.json"),
    orthogroups: Annotated[
        Optional[Path], typer.Option("--orthogroups", "-g", help="Orthogroups TSV file")
    ] = None,
    highlight: Annotated[
        Optional[list[str]], typer.Option("--highlight", "-H", help="Gene names to highlight")
    ] = None,
) -> None:
    """Convert BED files to synteny JSON.

    Example:
        synteny-viz convert-bed human:human.bed mouse:mouse.bed -g ortho.tsv -H TIGD4
    """
    from .converters import SyntenyBuilder

    builder = SyntenyBuilder()

    for bed_spec in bed_files:
        # Parse species_id:path format
        spec_str = str(bed_spec)
        if ":" in spec_str:
            species_id, path = spec_str.split(":", 1)
            display_name = species_id.replace("_", " ").title()
        else:
            path = spec_str
            species_id = Path(path).stem
            display_name = species_id.replace("_", " ").title()

        bed_path = Path(path)
        if not bed_path.exists():
            console.print(f"[red]Error:[/red] BED file not found: {bed_path}")
            raise typer.Exit(1)

        builder.add_species_from_bed(species_id, display_name, bed_path)
        console.print(f"Loaded: {species_id} from {bed_path}")

    if orthogroups:
        if not orthogroups.exists():
            console.print(f"[red]Error:[/red] Orthogroups file not found: {orthogroups}")
            raise typer.Exit(1)
        builder.load_orthogroups(orthogroups)
        console.print(f"Loaded orthogroups from: {orthogroups}")

    if highlight:
        for gene in highlight:
            builder.highlight(gene)

    data = builder.build()
    data.to_json_file(output)
    console.print(f"[green]Created:[/green] {output}")


if __name__ == "__main__":
    app()
