# --- Snakemake: Unified SigProfilerAssignment + CHORD workflow ---
from pathlib import Path
import os, re, json

configfile: "config/config.yaml"

# ---- Config ----
REF_BUILD      = config.get("sigprofiler", {}).get("genome_build", "GRCh38")
COSMIC_VERSION = config.get("sigprofiler", {}).get("cosmic_version", 3.4)
EXOME_DATA     = config.get("sigprofiler", {}).get("exome", False)
SP_INPUT_DIR   = config.get("sigprofiler", {}).get("input_vcf_dir", "input_vcfs")
SP_OUTPUT_DIR  = config.get("sigprofiler", {}).get("output_dir", "sigprofiler_output")

CH_VCF_DIR     = config.get("chord", {}).get("vcf_dir", "vcf")
CH_SV_CALLER   = config.get("chord", {}).get("sv_caller", "manta")
CH_RECURSIVE   = config.get("chord", {}).get("recursive", False)
CH_REF_UCSC    = config.get("chord", {}).get("ref_ucsc", "hg38")
CH_OUTPUT_DIR  = config.get("chord", {}).get("output_dir", "chord_result")

# ---- Targets ----
rule all:
    input:
        # SigProfiler outputs: marker file + directory
        directory(SP_OUTPUT_DIR),
        f"{SP_OUTPUT_DIR}/.done",
        # CHORD outputs
        f"{CH_OUTPUT_DIR}/merged_contexts.txt",
        f"{CH_OUTPUT_DIR}/chord_pred.txt"

# ---- SigProfilerAssignment rule ----
rule sigprofiler_assignment:
    input:
        # use the directory as an input wildcard (presence checked in script)
        SP_INPUT_DIR
    output:
        touch(f"{SP_OUTPUT_DIR}/.done")
    params:
        outdir = SP_OUTPUT_DIR,
        genome_build = REF_BUILD,
        cosmic_version = COSMIC_VERSION,
        exome = EXOME_DATA
    conda:
        "envs/sigprofiler.yaml"
    threads: 2
    log:
        "logs/sigprofiler.log"
    script:
        "scripts/run_sigprofiler.py"

# ---- CHORD rule ----
rule chord_predict:
    input:
        vcf_dir = CH_VCF_DIR
    output:
        merged = f"{CH_OUTPUT_DIR}/merged_contexts.txt",
        pred    = f"{CH_OUTPUT_DIR}/chord_pred.txt"
    params:
        outdir    = CH_OUTPUT_DIR,
        sv_caller = CH_SV_CALLER,
        recursive = CH_RECURSIVE,
        ref_ucsc  = CH_REF_UCSC
    conda:
        "envs/chord.yaml"
    threads: 2
    log:
        "logs/chord.log"
    script:
        "scripts/run_chord.R"
