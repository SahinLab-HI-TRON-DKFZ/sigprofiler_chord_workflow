# scripts/run_sigprofiler.py
import os, sys
from pathlib import Path
from SigProfilerAssignment import Analyzer as Analyze

# Snakemake params
input_dir  = Path(snakemake.input[0])
output_dir = Path(snakemake.params.outdir)
genome     = str(snakemake.params.genome_build)
cosmic_ver = float(snakemake.params.cosmic_version)
exome      = bool(snakemake.params.exome)

output_dir.mkdir(parents=True, exist_ok=True)

if not input_dir.exists() or not any(input_dir.iterdir()):
    raise FileNotFoundError(f"SigProfilerAssignment input directory is missing or empty: {input_dir}")

Analyze.cosmic_fit(
    samples=str(input_dir),
    output=str(output_dir),
    input_type="vcf",            # expects folder with VCF files
    context_type="96",
    genome_build=genome,
    cosmic_version=cosmic_ver,
    exome=exome,
    collapse_to_SBS96=True,
    make_plots=True,
    verbose=True
)

# Marker to signal completion
Path(snakemake.output[0]).touch()
