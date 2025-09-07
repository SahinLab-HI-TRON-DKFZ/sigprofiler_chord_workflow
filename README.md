#SigProfiler + CHORD Workflow (Snakemake)

Run **SigProfilerAssignment** (Python) and **CHORD** (R) in one go.

## Layout
```
.
├── Snakefile
├── config/
│   └── config.yaml
├── envs/
│   ├── sigprofiler.yaml
│   └── chord.yaml
└── scripts/
    ├── run_sigprofiler.py
    └── run_chord.R
```

## Configure
Edit `config/config.yaml`:
- `sigprofiler.input_vcf_dir`: VCF folder for SigProfilerAssignment.
- `sigprofiler.genome_build`: `GRCh38` or `GRCh37`; `cosmic_version`: e.g., `3.4`; `exome`: `true/false`.
- `chord.vcf_dir`: folder with `*snv_indel.vcf(.gz)` and `*sv.vcf(.gz)` per sample.
- `chord.ref_ucsc`: `hg38` (default) or `hg19`; `sv_caller`: `manta` by default.
- `chord.recursive`: set `true` to search subfolders.

## Run
```bash
snakemake -n -q                 # dry-run
snakemake --use-conda --cores 8 # run
```

## Outputs
- `sigprofiler_output/` (SigProfiler results; `.done` file marks completion)
- `chord_result/merged_contexts.txt`
- `chord_result/chord_pred.txt`

Author
Sakshi Singh
Email: sakshi.singh@dkfz-heidelberg.de
Developer and pipeline designer at SahinLab-HI-TRON-DKFZ
