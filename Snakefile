# =====================================================================================================================
# Mouse vestibular (saccule and utricle of BXD strains) project: main workflow.
# Each analysis in analyses/ has its own Snakefile and can be run alone; this file runs them in order,
# then rebuilds the manuscript outputs (manuscript/Snakefile: paper figures + tables, enrichment report, lab deck).
#   snakemake --profile profiles/slurm -n                  # everything (dry run)
#   snakemake --profile profiles/slurm 05_dge              # one analysis
#   snakemake -s analyses/05_dge/Snakefile --profile profiles/slurm   # the same, directly
# Shared paths and parameters: config/config.yaml.
# =====================================================================================================================
# ---- configuration: shared settings, then lab-specific overrides if present
configfile: "config/config.yaml"
from pathlib import Path
if Path("config/config_lab.yaml").exists():   # lab-specific settings, not in git (overrides config.yaml)
    configfile: "config/config_lab.yaml"

# ---- Snakemake used to run each analysis's own Snakefile
SMK = f"{config['conda_prefix']}/{config['envs']['pheast']}/bin/snakemake"      # Snakemake 9

# ---- steps: the analyses in analyses/, then the manuscript outputs
ANALYSES = ["01_rna_genotype_matching", "02_gwas_sumstats", "03_pantry", "04_twas", "05_dge", "06_enrichment",
            "07_single_cell_rna", "08_deconvolution", "09_interaction_qtl_mapping", "10_gwas"]
STEPS = ANALYSES + ["manuscript"]

# ---- run order: a step starts once the steps it depends on have finished (each writes logs/steps/<step>.done)
DEPENDS_ON = {"03_pantry": ["01_rna_genotype_matching", "02_gwas_sumstats"], "04_twas": ["02_gwas_sumstats", "03_pantry"],
              "05_dge": ["01_rna_genotype_matching", "03_pantry"], "06_enrichment": ["05_dge"],
              "08_deconvolution": ["03_pantry", "07_single_cell_rna"], "09_interaction_qtl_mapping": ["03_pantry", "08_deconvolution"],
              "manuscript": ["05_dge", "06_enrichment"]}

# ---- command for each step: run its Snakefile with the Slurm profile
CMD = {s: f"{SMK} -s {f} --profile {config['profile']}" for s, f in
       [(a, f"analyses/{a}/Snakefile") for a in ANALYSES] + [("manuscript", "manuscript/Snakefile")]}

# ---- these rules only launch the step workflows, so they run in this controller job, not as Slurm jobs
localrules: all


# ---- target: every step finished
rule all:
    input: [f"logs/steps/{s}.done" for s in STEPS]


# ---- one rule per step: wait for the steps it depends on
for _s in STEPS:
    rule:
        name: _s
        localrule: True
        input: [f"logs/steps/{d}.done" for d in DEPENDS_ON.get(_s, [])]
        output: touch(f"logs/steps/{_s}.done")
        params: cmd=CMD[_s]
        shell: "{params.cmd}"


### END
