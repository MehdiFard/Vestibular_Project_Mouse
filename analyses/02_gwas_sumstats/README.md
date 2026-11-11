# 02 GWAS summary statistics and LD reference

Mouse vestibular GWAS (six traits: aged VsEP females / males, young beam test, young VsEP P1 amplitude, P2 latency,
threshold; received, in `../../data/gwas_results`) cleaned to the genotyped variants and converted to LDSC
`.sumstats` for FUSION TWAS; LD reference = BXD genotypes, REF/ALT from GRCm38, one PLINK set per chromosome.

```
analyses/02_gwas_sumstats/
  Snakefile   ld_reference -> clean_gwas -> munge (LDSC) -> to_pheast
  scripts/    clean_gwas_to_make_sumstat.R; develop_notes.sh = interactive notes of the original run
  src/ldsc/   bulik/ldsc clone (munge_sumstats.py, Python 2.7 env `ldsc`)
  input/      gwas_results/*.gwas (copies as received, from ../../data/gwas_results); the BXD genotypes are
              read in place from ../../data/genotype/mouse_genotype_mm10.{tped,tfam}
  results/    gwas_results/ (cleaned .gwas.cln, sumstats_format/*.sumstats, gwas_metadata.tsv),
              genotype/ (pgen sets, LDREF/mouse_mm10_<chr>.{bed,bim,fam})
```

Both outputs are copied into `../03_pantry/pheast/input/{gwas_sumstat,LDREF}` (checked identical, October 2026),
where pheast and `../04_twas` read them. N is set to 100 for every trait (as in the original run).

Run: `snakemake -s analyses/02_gwas_sumstats/Snakefile --profile profiles/slurm` (from the project root).
