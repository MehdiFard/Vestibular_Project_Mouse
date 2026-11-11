# Saccule vs Utricle — PCA and differential expression

Expression comparison of the two mouse vestibular otolith organs (saccule, utricle)
at 8 and 22 weeks, built on the counts that the Pantry phenotyping pipeline already
produced.

Run everything with the analysis Snakefile (from the project root; DE steps go to Slurm):

```bash
snakemake -s analyses/05_dge/Snakefile --profile profiles/slurm
```

(`bash scripts/run_all.sh` is the original launcher and still works.) Before October 2026 this folder was
`Pantry/sac_vs_utr/`; its pathway analysis is now `../06_enrichment`, the manuscript, enrichment report and
lab-meeting deck are in `../../manuscript/`, and `rebuild_291/` (291-library sensitivity run) is still here.
Sample metadata (`src/`) was copied from the 2024 project; the genotype match is
`../01_rna_genotype_matching/results/ase_counter_result/rna_to_geno_summary.tsv`.

Note (October 2026): checks 06, 07 and 09 are older than the final `data/sample_table.tsv` (14 Sep), so
`snakemake -n` lists them for rerun; their inputs may have changed slightly since.

The deck lands at `report/saccule_vs_utricle_report.pptx`.

## Data

| | |
|---|---|
| Libraries | 260 (Saccule w8 73, Saccule w22 55, Utricle w8 74, Utricle w22 58) |
| Counts | kallisto `est_counts` from `../03_pantry/phenotyping/intermediate/<dataset>/expression/`, summed transcript→gene with **tximport** |
| Genes | 55,401 annotated → 27,488 protein-coding + lincRNA → 18,273 after the expression filter |
| Annotation | `Mus_musculus.GRCm38.102.chr.gtf` |
| Organ / age | taken from which Pantry dataset a library belongs to |
| Batch | `pheast/input/<dataset>/<dataset>.batches.tsv` (plates 1–4) |
| Sex | old project's `metadata_mouse_all_plates_with_matched_geno_ids.csv` (stops at FTL_945) |
| Strain | `vestibulart_project_mouse/result/ase_counter_result/rna_to_geno_summary.tsv` — RNA-to-genotype concordance |

Pantry's own `output/*/expression.bed.gz` is **quantile-normalised TPM**, and its
`unnorm/expression.bed` is raw TPM — neither is a valid DESeq2 input, which is why
counts are re-imported from the kallisto runs.

## Design

This is a BXD-style strain panel. Saccule and utricle come from the **same animal**. The animal is taken from the
submission workbook (`data/fastq_bulk/bulk_rnaseq_metadata_mouse_b1_to_b4.xlsx`: the two organs of an animal have
consecutive library IDs, same strain and age); where the workbook has no such pair, strain + age + sex is used.
Until October 2026 the pair was strain + age + sex only, which joined two different BXD128 w22 mice (FTL_1031,
FTL_1036) into one pair; fixed in `scripts/01_build_inputs.R` (same rule as `rebuild_291/`; contrast B: 104 -> 102
libraries, 860 -> 863 DEGs).

- 70 animal pairs (both organs) at w8, **51** at w22
- **32** strains (saccule) / **36** (utricle) sampled at **both** ages

| Organ | Age | Libraries | In a pair | Unpaired | Female | Male |
|---|---|---|---|---|---|---|
| Saccule | w8 | 73 | 70 | 3 | 73 | 0 |
| Saccule | w22 | 55 | 51 | 4 | 30 | 25 |
| Utricle | w8 | 74 | 70 | 4 | 74 | 0 |
| Utricle | w22 | 58 | 51 | 7 | 32 | 26 |
| **All** | | **260** | **242** | **18** | **209** | **51** |

A pair is one animal with both organs present. Within any one organ/age group each strain appears exactly once.

| Batch | w8 | w22 | Female | Male | Median reads (M) | Mapped % |
|---|---|---|---|---|---|---|
| 1 | 88 | 0 | 88 | 0 | 28.5 | 81.0 |
| 2 | 48 | 28 | 73 | 3 | 30.8 | 82.9 |
| 3 | 11 | 61 | 24 | 48 | 28.0 | 80.8 |
| 4 | 0 | 24 | 24 | 0 | 17.9 | 74.8 |

Organs are balanced inside every batch (44/44, 38/38, 35/37, 11/13). Batch 4 is the
later "oldmouse" run: it adds 10 w22 pairs but is the shallowest and worst-mapping
data in the set, which is why adding it raised the w22 pair count while slightly
lowering the w22 DEG count.

Consequences, all of which are visible in the PCA covariate figure:

- **Organ is clean.** It is balanced inside every batch and every strain.
- **Age is partly confounded with batch.** Batch 1 is w8-only, batch 4 is w22-only;
  the age effect is identified from batches 2 and 3, which span both.
- **Sex is confounded with age.** Every w8 library is female, so age and sex cannot
  be separated using all libraries. Sex itself is now known for **all 260** — recovered
  from expression, see below — which makes a clean female-only age contrast possible.

## Where the strain labels come from, and why they can be trusted

Strain is taken from the **RNA-to-genotype concordance table**, not the metadata
spreadsheet. That table compares each library's ASEReadCounter allelic fractions to
every strain in the genotype VCF and takes the best match. It is the same source the
Pantry pipeline itself used to attach genotypes (`11.*_Rscript06_make_fasta_list_make_geno_*.R`
line ~90 pulls genotypes by `Top_Simil_ID_ByMixAna`), so the strain label and the
genotype behind each library agree by construction.

Audit of that table (`scripts/01_build_inputs.R` records `strain_simil`, `strain_gap`
and `strain_source` for every library so any call can be re-checked):

| Check | Result |
|---|---|
| Pantry libraries covered | **260 / 260** |
| Minimum concordance | 0.914 (the pipeline's build filter was ≥ 0.90) |
| Libraries with low concordance **and** an ambiguous runner-up | **0** |
| Self-consistency: both organs of one animal, matched independently, get the same strain | **112 / 113 = 99.1%** |
| Agrees with the metadata sheet (after normalising `BXD45` vs `BXD45/RwwJ`) | 216 |
| Disagrees with the metadata sheet | 20 |
| Absent from the metadata sheet entirely | 24 |

**The 20 disagreements are label swaps, not match failures.** The clearest example:

| Library | Metadata says | Genotype says | Concordance |
|---|---|---|---|
| FTL_733 / FTL_734 (w8 pair) | BXD68 | **BXD202** | 0.965 |
| FTL_735 / FTL_736 (w8 pair) | BXD202 | **BXD68** | 0.997 |
| FTL_884 / FTL_916 (w22) | BXD68 | BXD68 ✓ | 0.997 |

Two w8 pairs have their spreadsheet labels exchanged, while the same strain's w22
libraries match correctly. Finding these is exactly what the concordance analysis was
built to do, and Pantry already used the corrected assignment.

Two caveats on the source table:
- Its `is_correct` column is **meaningless here** — it tests `RNA_ID == top_simil_ID`,
  i.e. a library ID against a strain name, so it is `False` for all 321 rows. It is not
  evidence of a failed match.
- The lowest-gap calls (0.001) are near-identical substrains — `C3HeB/FeJ` vs `C3H/HeJ`,
  `BALB/cJ` vs `BALB/cByJ`. Genetically near-equivalent, so harmless for pairing.

## Validated against the complete submission sheet

The sheet this project was built from
(`metadata_mouse_all_plates_with_matched_geno_ids.csv`) covers **FTL_588-FTL_945**, 285
libraries. Every batch-4 library is numbered FTL_1022-FTL_1057 and so falls outside it —
which is exactly why 24 libraries had no recorded strain. The sheet was written before
batch 4 existed and never updated.

The **complete four-batch workbook** turned up in the backup archive:

```
/data/rss/baker/mohammadi_p/data_backup/vestibular_project_data/mouse/
    rnaseq_fastq/bulk/bulk_rnaseq_metadata_mouse_b1_to_b4.xlsx
```

323 libraries, FTL_588-FTL_1057, one sheet per batch. `scripts/17_validate_against_backup_sheet.py`
checks every data-derived annotation against it (`backup_sheet_validation.tsv`).

| Check | Result |
|---|---|
| Libraries found in the sheet | **260 / 260** |
| **Age** agrees | **260 / 260** |
| **Sex** agrees where recorded | **232 / 232** |
| — of the 52 sexes inferred from *Xist*/Y genes, now checkable | 24, **all 24 confirmed** |
| **Organ** agrees | 258 / 260 |
| **Strain** agrees | 230 / 260 |
| Animal pairs whose two organs got different strain calls | **0 / 122** |

### The 24 batch-4 libraries: all confirmed

All 24 are in the sheet and **all 24 strain assignments made by genotype concordance match
the sheet exactly** — BXD125, BXD144, BXD32, BXD128, BXD63, BXD95, BXD160, BXD100, BXD38,
BXD156, BXD85, BXD177, BXD178, BXD157. Their sexes (all female) are confirmed too. The
recovery that took the w22 organ pairs from 43 to 52 is independently correct.

### The 30 strain disagreements are not 30 errors

| Type | n | Example |
|---|---|---|
| Spelling only | 8 | `BXD131a`→BXD131, `BXD147a`, `BXD183a`, `BXF213` (typo for BXD213), `C57BI6J` (capital I for L) |
| Near-tie, genuinely ambiguous | 4 | `C3H/H2J`: C3HeB/FeJ 0.9961 vs C3H/HeJ 0.9952, gap **0.0009**; BXD151/BXD152, gap 0.0035 |
| Genuine relabel | 18 | see below |

Two of the genuine relabels are **clean reciprocal swaps**, which is what a mix-up of two
tubes looks like:

- **BXD34 ↔ BXD60** — FTL_785/786 labelled BXD60 match BXD34 at 0.995; FTL_787/788 labelled
  BXD34 match BXD60 at 0.996
- **BXD68 ↔ BXD202** — FTL_733/734 labelled BXD68 match BXD202 at 0.965; FTL_735/736
  labelled BXD202 match BXD68 at 0.997

**None of this affects the differential expression results.** Strain enters the organ
contrasts only as a blocking factor identifying which two libraries came from one animal,
and no animal pair (0 of 122) received different strain calls for its two organs. It would
matter for QTL mapping — but Pantry attaches genotypes *by* the concordance match, so
labels and genotypes agree there by construction.

### The two organ disagreements: the data is right

`FTL_600` and `FTL_601` (KK/HlJ, w8) are the same animal with saccule and utricle swapped
relative to the sheet. Expression settles it: PC1 separates the organs at R² = 0.85, and
FTL_601 — which the sheet calls saccule — sits at **PC1 = -16.6**, deep in the utricle
cluster, while FTL_600 sits at +1.9 on the saccule side. **The analysis labels are correct
and the sheet row is wrong.**

Separately, seven libraries sit on the wrong side of the PC1 midpoint (FTL_596, 591, 622,
624, 722, 661, 782), all within PC1 -3.6 to +5.9. Their organ labels *agree* with the
sheet, so these are not labelling errors; they are transcriptionally intermediate, most
likely dissection cross-contamination between two adjacent organs.

## Were all samples used?

Yes, through the whole Pantry pipeline — 260 in, 260 out:

| Stage | saccule_w8 | saccule_w22 | utricle_w8 | utricle_w22 |
|---|---|---|---|---|
| phenotyping input | 73 | 55 | 74 | 58 |
| Pheast input | 73 | 55 | 74 | 58 |
| genotype `.fam` | 73 | 55 | 74 | 58 |
| QTL covariate matrix | 73 | 55 | 74 | 58 |

Nothing was lost between phenotyping and QTL mapping.

### From 323 sequenced libraries to 260 analysed

The backup workbook lists **323** libraries; the RNA-to-genotype concordance table has
**321** (FTL_752 and FTL_775 have FASTQs in the archive but never reached the matching
step); **260** are analysed here. The 61 dropped between 321 and 260 account exactly:

| Reason | n |
|---|---|
| Concordance below the 0.90 threshold — failed genotype QC | **30** |
| Same strain sampled twice within one organ-age group, deduplicated | **23** |
| Lost to an edge case in the deduplication rule (see below) | **8** |

The selection is done by `scripts/11.1`–`11.4_Rscript06_make_fasta_list_make_geno_*.R` in
the `vestibulart_project_mouse` project, which write
`data/genotype/<organ>_<age>_similarity0.9_unik.{tfam,tped}`. Those four files contain
73 / 55 / 74 / 58 libraries — exactly the 260 used here.

**The deduplication rule is more careful than it first looks.** When a strain appears twice
within an organ-age group it does *not* keep an arbitrary one: it keeps the library whose
*originally recorded* strain already matched the genotype call (`Original_ID == "YES"`) and
drops the one that became a duplicate only because the mix-up analysis reassigned it.
`dplyr::distinct(..., .keep_all = TRUE)` runs afterwards and only breaks remaining ties —
and it keeps the **first** row, not a random one as its comment claims.

⚠️ **The edge case: BXD187 and BXD195 were lost entirely.** If *every* library claiming a
strain got there by reassignment, they all score `"NO"` and the rule drops all of them,
leaving the strain unrepresented. That happened to two strains at w22, costing 8 libraries
(4 animals):

| Library | Organ | Recorded strain | Genotype call | Concordance |
|---|---|---|---|---|
| FTL_771 / FTL_772 | Sac / Utr | BXD195 | **BXD187** | 0.932 / 0.931 |
| FTL_783 / FTL_784 | Sac / Utr | BXD174a | **BXD187** | 0.934 / 0.931 |
| FTL_773 / FTL_774 | Sac / Utr | BXD187 | **BXD195** | 0.941 / 0.941 |
| FTL_931 / FTL_932 | Sac / Utr | BXD89 | **BXD195** | 0.942 / 0.943 |

Note the first and third rows: **BXD187 and BXD195 are swapped with each other** — a third
reciprocal swap, alongside BXD34↔BXD60 and BXD68↔BXD202. Their concordances (0.93–0.94) are
lower than the dataset median (0.991) but the runner-up strains sit near 0.70, so the calls
are unambiguous relative to the alternatives.

Recovering them would add **2 strains and up to 2 organ pairs** at w22 (one animal per
strain after deduplication). Worth doing if the w22 contrast is ever refitted; it does not
invalidate anything already computed.

## What the pairing buys, in numbers

Each strain contributed exactly **one saccule and one utricle library per age**, from
the same animal — so `strain` is the pairing key, the same role a subject ID plays in
a paired t-test. Across the w8 libraries, averaged per gene:

| Variance explained | mean R² |
|---|---|
| animal (genetic background), raw R² | 0.619 |
| animal, **adjusted** R² | **0.243** |
| animal, excess over a permuted 70-level grouping | 0.123 |
| organ (the effect being measured), raw R² | 0.104 |
| organ, **adjusted** R² | **0.097** |

⚠️ **Do not quote the raw 0.62 against 0.10 as "six times".** The animal factor has 69
parameters against organ's one, on 140 libraries, so it explains a lot by construction — a
*randomly permuted* 70-level grouping already scores R² = 0.496. After adjusting, the gap
is ~2.5×, and on an excess-over-permutation basis ~1.3×. Genetic background is genuinely the
larger nuisance, but not by the margin the raw numbers suggest.

The argument for pairing does not rest on that ratio anyway — it rests on the measured
outcome below. An unpaired model has to carry all
of that in the residual; `~ strain + organ` cancels it. That is worth ~1,400 extra
genes at FDR<0.05 (10,923 vs 9,537 at w8).

Cost of pairing: animals contributing only one organ are dropped. At w8 that is 7 of
147 libraries; at w22, 11 of 113.

The 24 "oldmouse" libraries that the metadata sheet lacks were recovered from the
concordance table (all ≥ 0.974, gap ≥ 0.113), which took the w22 organ contrast from
**43 pairs to 52** (51 after the October 2026 pairing fix). Note the trade-off: those new pairs are all batch 4, which was
sequenced to about **60% of the depth** of the others (median 17.9M reads at 74.8%
mapped, against 28–31M at 81–83%), so the fold changes now rest on more animals but
the significant-gene count fell slightly (926 → 860). Both versions are defensible;
the current one generalises better.

## Why the UNPAIRED run reports more genes than the paired one

`DE_summary.tsv` looks backwards at first glance: the unpaired sensitivity run reports
**773** genes at w8 against the paired run's 756, and **928** at w22 against 863. More
noise should mean fewer findings, not more. `scripts/15_paired_vs_unpaired_check.R`
refits both models on **exactly the same libraries** so the model and the sample set
can be told apart:

| Age | Model | Libraries | Resid. df | FDR < 0.05 | \|log2FC\| > 1 | Both | Median SE |
|---|---|---|---|---|---|---|---|
| w8 | `~ strain + organ` (paired) | 140 | 69 | **10,923** | 756 | 756 | 0.0359 |
| w8 | `~ batch + organ` (same libraries) | 140 | 136 | 9,448 | 796 | 796 | 0.0456 |
| w8 | `~ batch + organ` (all libraries) | 147 | 143 | 9,537 | 773 | 773 | 0.0447 |
| w22 | `~ strain + organ` (paired) | 102 | 50 | **10,934** | 863 | 863 | 0.0450 |
| w22 | `~ batch + organ` (same libraries) | 102 | 98 | 9,414 | 911 | 911 | 0.0567 |
| w22 | `~ batch + organ` (all libraries) | 113 | 109 | 9,550 | 928 | 928 | 0.0551 |

Three things to read out of it:

1. **The headline count is not a power measure.** `Both` equals `|log2FC| > 1` on every
   row — every gene clearing two-fold also clears FDR, in both models. For the organ
   contrasts the FDR test filters nothing; the headline number is counting effect sizes.
2. **Pairing's benefit is in the FDR column**, and it is large: on identical libraries
   strain adds **1,475** genes at w8 and **1,520** at w22, even though it spends 70 (w8)
   and 51 (w22) degrees of freedom doing it. Residual df falls from 136 to 69 and the
   model still wins, because the between-animal variance it removes is larger than the organ effect.
3. **Noise inflates an effect-size threshold.** Dropping strain makes each fold change
   ~20% less precise (median SE ratio 1.20 at w8, 1.17 at w22). A noisier estimate
   overshoots and undershoots equally often, but far more genes sit just *below*
   |log2FC| = 1 than just above, so random error promotes more genes than it demotes:

   | Same libraries | crosses into \|log2FC\|>1 | falls back out | net |
   |---|---|---|---|
   | w8 | 76 | 36 | **+40** |
   | w22 | 98 | 50 | **+48** |

   The extra unpaired genes are the threshold catching noise, not extra discovery.

### Why strain is not in the unpaired model

Because "unpaired" *means* "without strain" — strain **is** the pairing. The run exists
only to show what the same contrast looks like if you ignore that the two organs came
from one animal, so putting strain back would just reproduce analysis A/B.

There is also nothing to gain from a hybrid. In `~ strain + organ`, an animal that
contributed only one organ gets a strain coefficient that fits its single library
exactly: it adds one parameter and one observation, and contributes nothing to the
organ contrast. That is why A and B use 140 of 147 libraries at w8 and 102 of 113 at
w22 — the rest are not thrown away for tidiness, they carry no information once the
pairing is in the model.

## What the organ × age interaction actually tests

Analyses A–D ask "is this gene different?". Contrast E asks whether the **organ
difference itself** changes with age:

```
gap(w8)      = log2( Utricle / Saccule ) at 8 weeks
gap(w22)     = log2( Utricle / Saccule ) at 22 weeks
interaction  = gap(w22) - gap(w8)          <- what is tested
```

With design `~ batch + organ + age + organ:age`, this is exactly the `organUtricle.agew22`
coefficient, and the four group means are simple sums of coefficients — so the arithmetic
in the table below is the model's, not an approximation.

| Gene | gap at w8 | gap at w22 | Interaction | FDR | Shape |
|---|---|---|---|---|---|
| Gpnmb | +3.23 | +3.09 | −0.14 | 0.83 | no interaction — parallel lines |
| Kcne1 | +3.99 | +7.18 | **+3.19** | <0.001 | gap widens, utricle pulls ahead |
| Jchain | −0.99 | −3.96 | **−2.97** | <0.001 | gap widens, saccule pulls ahead |
| S100a8 | −4.63 | −2.33 | **+2.30** | 0.0037 | gap narrows |
| Slc14a1 | +1.86 | +0.55 | **−1.31** | 0.0044 | gap narrows |
| Cxcl10 | +1.31 | −0.42 | **−1.73** | 0.0023 | gap reverses sign |

Worked example, Kcne1: at w8 the utricle carries 2^3.99 ≈ 16× more than the saccule; by
w22 it carries 2^7.18 ≈ 145×. The gap grew by 3.19 log2 units. Both organs *lose* Kcne1
with age — the saccule simply loses far more of it. The interaction says nothing about
the direction of the age change, only about the gap between the organs.

Gpnmb is the useful control: a 9-fold organ difference that is the same at both ages, so
the interaction is ≈ 0. **A gene can be hugely organ-specific and score zero here.**

Only **122 of 18,273** genes have a significant interaction, against 756 and 863 for the
organ contrasts. The organ identity of these two tissues is largely age-stable, and the
interaction list is the least reproducible result in the project (see the limma
cross-check). `scripts/16_interaction_illustration.R`, figure `31`.

Five of those six genes are significant; **Gpnmb is deliberately not** — it is the
control showing that a large, constant organ difference scores zero on this test.

### The nine most significant interaction genes

Figure `32` draws the same panels for the top nine by adjusted p (`interaction_top9.tsv`):

| Gene | Mean count | gap w8 | gap w22 | Interaction | FDR |
|---|---|---|---|---|---|
| Cma1 | 11 | +2.56 | −3.10 | −5.66 | 2.2e-16 |
| Mcpt4 | 9 | +2.71 | −2.88 | −5.60 | 1.5e-15 |
| Trpm1 | 113 | +3.41 | +1.85 | −1.57 | 1.2e-13 |
| Gchfr | 87 | +2.96 | +1.27 | −1.69 | 1.2e-13 |
| Fam163b | 93 | −1.30 | −0.07 | +1.22 | 1.0e-11 |
| Slc5a7 | 15 | +0.87 | −1.73 | −2.60 | 4.1e-11 |
| Jchain | 127 | −0.99 | −3.96 | −2.97 | 1.8e-09 |
| Tpsb2 | 11 | +2.50 | −1.11 | −3.61 | 5.6e-09 |
| Slitrk2 | 34 | −0.50 | −1.72 | −1.22 | 1.2e-08 |

Two things stand out. **Cma1, Mcpt4 and Tpsb2 are all mast-cell proteases** (chymase,
mast cell protease 4, tryptase β2) and move together — utricle-enriched at w8, saccule-
enriched at w22 — which reads as a real difference in mast-cell content rather than three
independent findings. But all three have **mean counts of 9–11**: they top the ranking on
a large fold change carried by very few reads, so treat them as a lead, not a result.
`Trpm1`, `Gchfr`, `Fam163b` and `Jchain` are the best-supported entries on the list.

The `Mean count` column is in the figure's panel labels for exactly this reason.

## Analyses

| Key | Contrast | Model |
|---|---|---|
| A | Utricle vs Saccule at w8 | `~ strain + organ` (paired) |
| B | Utricle vs Saccule at w22 | `~ strain + organ` (paired) |
| C | w22 vs w8 in Saccule | `~ batch + age` |
| D | w22 vs w8 in Utricle | `~ batch + age` |
| E | Organ × age interaction | `~ batch + organ + age + organ:age` |
| S | Unpaired organ contrasts | `~ batch + organ` (sensitivity check) |

Positive `log2FoldChange` = higher in **Utricle** (A, B), higher at **w22** (C, D),
or the organ gap **widening** with age (E).

Significance: `padj < 0.05` **and** `|log2FC| > 1`. Fold changes are shrunken with
`ashr`; the unshrunken value is kept as `log2FC_raw`.

### Gene filtering

Two filters, in this order — both set in `00_config.R`:

1. **Biotype.** `KEEP_BIOTYPES <- c("protein_coding", "lincRNA")`. Everything else is
   dropped *before* testing so it never enters the FDR correction: TEC, antisense,
   pseudogenes, processed_transcript, miRNA, snoRNA and the rest.
2. **Expression.** ≥ 10 counts in ≥ 55 libraries (the size of the smallest group).

55,401 → 27,488 → 18,273. Each contrast re-applies both on its own subset (~17,000 tested).

⚠️ This GTF is **Ensembl 102**, where lncRNAs are still split across several biotypes.
`lincRNA` (1,186 genes tested) is the *intergenic* class only. If you want lncRNA in
the broad sense, add `antisense` (1,056), `processed_transcript` (425),
`bidirectional_promoter_lncRNA` (126) and `sense_intronic` (121) to `KEEP_BIOTYPES`
and re-run. Ensembl merged all of these into a single `lncRNA` biotype in later releases.

## Differences from the earlier `vestibular_mouse` scripts

| | Earlier | Here | Why |
|---|---|---|---|
| Counts | RSEM `expected_count` | kallisto + tximport | matches the Pantry quantification the rest of this project uses; tximport supplies the transcript-length offset |
| Gene filter | non-zero in ≥ 20 samples | protein-coding + lincRNA, then ≥ 10 counts in ≥ 55 samples | the old filter keeps genes that are near zero everywhere, which only inflates the multiple-testing burden |
| Organ model | `~ plate + age + tissue` | `~ strain + organ`, per age | the organs are paired within animal; ignoring that throws away the pairing and leaves genetic background in the residual |
| PCA input | `log2(rawCount + 1)` on TPM/counts | VST | behaves properly at low counts |
| Fold changes | raw | `ashr`-shrunken | stops noisy low-count genes dominating the ranking |
| Sex | in the sample table | recovered for all 260 from Xist/Y genes, then used in a female-only age contrast | 52 libraries had no recorded sex; sex is confounded with age across all libraries, so it is handled by restriction rather than as a covariate |
| Interaction | not tested | contrast E | answers "does the organ difference change with age" directly rather than eyeballing two lists |

The earlier volcano plots also thresholded on **raw** `pvalue` while labelling the
axis `padj`; here everything is on adjusted p.

## The age contrast, done within strain

37 strains were sampled at **both** ages (32 of them in saccule, 36 in utricle), so
their age effect can be measured within a genetic line rather than between animals. `scripts/10_pairing_and_paired_age.R`.

Batch turns out to be estimable alongside strain and age here — strains differ in
which batches their two libraries landed in (b1→b3, b2→b3, b1→b2, b2→b2 …), so the
model is `~ strain + batch + age`. This is the only contrast in the project that is
free of *both* genetic background and batch.

| Organ | Strains | Libraries | Main analysis (FDR+LFC) | Within-strain | Shared | Same direction |
|---|---|---|---|---|---|---|
| Saccule | 24 | 48 | 156 | **193** | 104 | 99.4% |
| Utricle | 25 | 50 | 141 | **200** | 102 | 100% |

Two things follow. First, 48–50 libraries find **more** significant genes than 128–132
do — that is how much noise the between-animal comparison carries. Second, **99–100% of
the main-analysis hits move in the same direction** in this cleaner fit, which is the
best evidence available here that the age results are real biology rather than batch.

## Why batch is in the design rather than corrected out first

`limma::removeBatchEffect` is used here **only for plotting** (figure 05). The DE
models keep `batch` as a term. `scripts/06_batch_correction_check.R` measures why,
on the saccule data (128 libraries, 17,328 genes):

| Check | batch in the design | removeBatchEffect, then test `~ group` |
|---|---|---|
| Real w22-vs-w8 contrast, genes at FDR<0.05 | 6,097 | 9,875 (**+62%**) |
| Null label **balanced** within batch — false positives | 0.2 | 0.3 |
| Null label **confounded** with batch — false positives | **0.0** | **50.5** (max 986) |

The usual "you lose degrees of freedom" argument is real but negligible at this n
(123 vs 126 residual df). The damage comes from **collinearity**: when batch overlaps
the factor of interest, the correct model widens that factor's standard error to
reflect how little independent information there is. A corrected matrix has no batch
term left to be collinear with, so the standard error comes out far too small.

That is exactly this study's situation — batch 1 is w8-only, batch 4 is w22-only — so
correcting first would have made the age results look far more certain than they are.
The limma man page says the same thing: *"This function is not intended to be used to
prepare data for linear modeling by `lmFit`."*

## Does dropping the other biotypes cost real separation?

No — `scripts/07_biotype_separation_check.R` measures it:

| | all biotypes (23,428) | protein-coding + lincRNA (18,273) |
|---|---|---|
| PC1 % variance | 22.2% | 28.7% |
| Organ separation on PC1 (standardised) | 4.50 | **4.71** |
| PC1 R² organ | 0.836 | **0.848** |
| PC2 R² age | 0.652 | 0.066 |
| PC2 R² **batch** | **0.781** | 0.067 |
| Best PC for age | PC2, R²=0.652 | **PC3, R²=0.701** |

Two things are going on:

1. **"% variance explained" is a ratio.** Adding 5,155 noisy genes inflates the
   denominator, so PC1's *share* falls (28.7% → 22.2%) even though the organ
   separation is unchanged. Judge separation with a scale-free measure, not PC1%.
2. **The extra biotypes carry batch, not biology.** Mean per-gene R²: the other
   biotypes track batch 0.191 vs organ 0.035 — they are low-count, multi-mapping
   features (pseudogenes, TEC) that are exquisitely sensitive to library prep and
   RNA integrity. 18.4% of the top-2000 variable genes are these, and in that top
   set they track batch (0.162) over organ (0.064), the reverse of the coding genes
   (organ 0.236, batch 0.081).

Because batch and age are confounded here, that batch axis *looks* like age
separation on PC2. After filtering, age is not lost — it moves to PC3 where it is
both stronger (0.701 vs 0.652) and finally beats batch (0.552) instead of losing to it.

## Cross-check: DESeq2 vs limma-voom

`scripts/08_limma_crosscheck.R` reruns every contrast with limma-voom, same
libraries, filters and designs:

| Contrast | DESeq2 | limma | Shared | % of DESeq2 recovered | r (log2FC) |
|---|---|---|---|---|---|
| Utricle vs Saccule, w8 | 754 | 846 | 729 | 96.7% | 0.981 |
| Utricle vs Saccule, w22 | 924 | 1,044 | 885 | 95.8% | 0.977 |
| w22 vs w8, Saccule | 167 | 358 | 150 | 89.8% | 0.870 |
| w22 vs w8, Utricle | 146 | 256 | 136 | 93.2% | 0.912 |
| Organ × age interaction | 90 | 51 | 33 | 36.7% | 0.798 |

The organ results are essentially engine-independent. The age results reproduce at
~90%, with limma calling roughly twice as many overall. **The interaction is the one
place the two disagree** (37% shared) — interaction terms carry the least
information, so the two shrinkage schemes pull apart. Treat that list as provisional.

## Two notes on the PCA

**"VST" = variance-stabilising transform** (`DESeq2::vst`), not counts — my earlier
shorthand "VST counts" was sloppy. It maps counts onto a log2-like scale on which a
gene's spread no longer depends on how highly it is expressed. Median per-gene SD
across expression octiles, each scaled to its own maximum (`scripts/09_vst_and_ntop_check.R`):

| Expression octile | Q1 (low) | Q3 | Q5 | Q8 (high) |
|---|---|---|---|---|
| normalised counts | 0.01 | 0.08 | 0.19 | **1.00** |
| log2(count + 1) | **1.00** | 0.58 | 0.30 | 0.27 |
| VST | 0.63 | 1.00 | 0.79 | 0.85 |

Raw counts let the handful of highest-expressed genes dominate any distance measure;
`log2(count+1)` overcorrects, making the *lowest*-expressed genes the noisiest (a gene
going 1 → 2 counts reads as a two-fold change). VST is the flat one. It converges to
log2 at high counts and shrinks toward a constant at low ones — e.g. a gene with mean
count 16 varies by 0.79 on the log2 scale but only 0.15 on VST, while a gene with mean
count 10,044 reads 13.16/13.46/13.60 on log2 and 13.18/13.48/13.61 on VST.

**The top-2000 applies to the PCA only.** Every DE test in script 03 uses all ~17,000
genes that survive the two filters — there is no variable-gene selection in the DE path.
And the choice of N barely matters:

| genes used | 250 | 500 | 1,000 | **2,000** | 5,000 | 10,000 | all 18,273 |
|---|---|---|---|---|---|---|---|
| PC1 % variance | 32.2 | 31.3 | 30.6 | 28.7 | 25.1 | 22.4 | 20.4 |
| PC1 R² organ | 0.822 | 0.838 | 0.847 | **0.848** | 0.838 | 0.823 | 0.813 |

PC1 tracks organ at R² = 0.81–0.85 across the whole range, so the conclusion does not
depend on the cutoff. (The equivalent number for age is unstable — 0.40 to 0.74 — because
age is a weaker signal smeared across several PCs, which is consistent with everything
else here.)

## Sex: recovered from expression

52 libraries had no sex in the metadata (24 absent from the sheet, 28 present with the
field blank), and no other file on disk carries it. It does not need a spreadsheet:
**Xist** is female-specific and the Y genes **Ddx3y / Uty / Eif2s3y / Kdm5d** are
male-specific. `scripts/01_build_inputs.R` makes the call; `scripts/13_infer_sex.R`
validates it and draws the figure.

| | median Xist CPM | median Y-gene CPM |
|---|---|---|
| known females (160) | 1653.8 | 4.7 |
| known males (48) | 3.4 | 275.4 |

The two groups are separated by a gap of 11.5 log2 units, and the rule is **100%
concordant** with all 208 libraries whose sex *was* recorded (160/160 F, 48/48 M).
Final counts: w8 = 147 F / 0 M, w22 = 62 F / 51 M.

### Did the missing sex harm the analyses?

**Organ contrasts: no.** Both organs come from one animal, so sex is constant within a
pair and is absorbed by the strain term exactly as batch is. Nothing to correct.

**Age contrasts: a little, and it is now fixable.** Only a handful of genes respond to
sex at all, and the age effect survives holding sex constant:

| Organ | Age DEGs (all libs) | Age DEGs (females only) | Shared | % recovered | r | Sex DEGs in w22 | Age DEGs that are also sex DEGs |
|---|---|---|---|---|---|---|---|
| Saccule | 168 | **380** | 151 | 89.9% | 0.876 | 31 | 11 |
| Utricle | 146 | **275** | 138 | 94.5% | 0.917 | 8 | 3 |

Only 11 of 168 and 3 of 146 age DEGs are also sex-responsive, so contamination is
~2–7%. The female-only contrast (`DE_age_females_*.tsv.gz`) recovers 90–95% of the
original DEGs *and finds more*, because dropping the males removes heterogeneity from
the w22 group. Treat it as the cleaner age result.

Of the 31 sex-responsive genes in w22 saccule, 6 are the expected X/Y markers
(Xist, Ddx3y, Eif2s3y, Kdm5d, Uty, Gm29650); the remaining 25 are autosomal with
modest effects.

### Sex on the PCA

`scripts/02_pca.R` also draws the PCA coloured by sex: `08_pca_by_sex.png` (PC1 vs PC2)
and `09_pca_pc1_pc3_by_sex.png` (PC1 vs PC3). Colour = sex, shape = organ. Each figure
has two panels — all 260 libraries, and the 22-week libraries alone, which are the only
ones where both sexes are present.

- **PC1 is clean.** Organ explains R² = 0.85 of it; sex explains **0.001** across all
  libraries and **0.000** within 22 weeks. The two sexes sit on top of each other inside
  each organ cluster, so the organ contrasts cannot be confounded by sex.
- **PC3 is where sex lives**, and it is crowded: age 0.70, batch 0.55, sex 0.45. Across
  all libraries these are not separable — every 8-week animal is female, so "sex" and
  "age" colour the same split.
- **Restricting to 22 weeks removes age but not batch.** There, sex explains R² = 0.36
  of PC3 and batch 0.31, because sex and batch are themselves confounded
  (`w22_sex_by_batch.tsv`):

| Batch (w22 only) | Female | Male |
|---|---|---|
| 2 | 25 | 3 |
| 3 | 13 | 48 |
| 4 | 24 | 0 |

This is the concrete reason sex is handled by **restriction** — a female-only age
contrast — rather than as a covariate: adding `sex` to `~ batch + age` would leave it
competing with batch for the same variance.

One correction this surfaced: the covariate R² table (`pca_covariate_r2.tsv`, figure
`06_pca_covariate_r2.png`) used the metadata `sex` column, which is blank for 52
libraries, so the sex row was computed on a biased subset of 208. It now uses the
recovered `sex_final`, complete for all 260. PC3 sex accordingly reads **0.45**, not the
0.64 the earlier table showed.

### A pairing bug this uncovered

These are inbred strains, so two different mice of the same strain are genetically
identical and the genotype match cannot tell them apart. One "pair" — BXD51 at w22 —
turned out to be a **female saccule and a male utricle**, i.e. two different animals.
`pair_id` is therefore `strain + age + sex`, not `strain + age`. That drops the w22
organ pairs from 53 to **52**, and the pairing is now correct by construction.

## Master metadata

`results/tables/master_metadata.tsv` (and a `.csv` copy) is the single place where
everything about the 260 libraries lives — 39 columns, one row per library.
`master_metadata_dictionary.tsv` describes every column. Built by
`scripts/12_build_master_metadata.R`.

| Group | Columns |
|---|---|
| Identity | `library_id`, `dataset`, `organ`, `age`, `group` |
| Genetics | `strain`, `strain_source`, `strain_simil`, `strain_gap`, `strain_runnerup`, `strain_runnerup_simil`, `strain_metadata_sheet`, `strain_label_from_rna_submission`, `geno_id_original` |
| Technical | `batch`, `sex`, `insert_size`, `read_length`, `paired_end`, `n_fastq_pairs`, `fastq_r1`, `fastq_r2` |
| Sequencing QC | `n_reads_processed`, `n_pseudoaligned`, `pct_pseudoaligned`, `pct_unique` |
| Expression QC | `total_gene_counts`, `genes_detected_min10`, `genes_detected_min10_pc_linc` |
| Design flags | `pair_id`, `has_organ_pair`, `sampled_at_both_ages`, `in_paired_organ_contrast`, `in_paired_age_contrast` |
| Coordinates | `PC1`–`PC5` |

Every strain call carries its own evidence, so disagreements stay visible rather than
being silently overwritten:

| library_id | organ | strain | strain_metadata_sheet | strain_simil | strain_source |
|---|---|---|---|---|---|
| FTL_733 | Saccule | BXD202/RwwJ | BXD68 | 0.965 | DISAGREES with metadata |
| FTL_734 | Utricle | BXD202/RwwJ | BXD68 | 0.966 | DISAGREES with metadata |
| FTL_735 | Saccule | BXD68/RwwJ | BXD202 | 0.997 | DISAGREES with metadata |
| FTL_736 | Utricle | BXD68/RwwJ | BXD202 | 0.997 | DISAGREES with metadata |

### QC summary from the master file

Mapping rate is healthy throughout — median 81.2%, minimum 62.8%. Batch 4 is the
weakest on both depth and mapping:

| Batch | n | Median reads (M) | Median % mapped |
|---|---|---|---|
| 1 | 88 | 28.5 | 81.0 |
| 2 | 76 | 30.8 | 82.9 |
| 3 | 72 | 28.0 | 80.8 |
| **4** | 24 | **17.9** | **74.8** |

That is the reason the w22 organ contrast gained 10 pairs but lost a few significant
genes — the added libraries are real animals, just shallower ones.

246 of 260 libraries are in the paired organ contrast; the 14 that are not are animals
where only one organ was sequenced.

## Layout

```
analyses/05_dge/
├── Snakefile
├── src/                               sample metadata (strain, sex, plate)
├── rebuild_291/                       291-library sensitivity re-analysis (own scripts/, report/)
├── scripts/
│   ├── 00_config.R                    paths, thresholds, colours
│   ├── 01_build_inputs.R              sample table, annotation, tximport matrix
│   ├── 02_pca.R                       PCA, scree, covariate R², batch and sex checks
│   ├── 03_differential_expression.R   the five DESeq2 analyses
│   ├── 04_figures.R                   volcano, concordance, top-gene figures
│   ├── 05_build_report.py             builds the PowerPoint
│   ├── 06_batch_correction_check.R    why batch stays in the design
│   ├── 07_biotype_separation_check.R  what the biotype filter costs
│   ├── 08_limma_crosscheck.R          limma-voom rerun of every contrast
│   ├── 09_vst_and_ntop_check.R        what VST does; PCA top-N sensitivity
│   ├── 10_pairing_and_paired_age.R    what pairing buys; within-strain age contrast
│   ├── 11_model_diagram.R             explains what each formula compares
│   ├── 12_build_master_metadata.R     the master metadata file
│   ├── 13_infer_sex.R                 validates the sex calls, draws the figure
│   ├── 14_sex_effect.R                sex effect and the female-only age contrast
│   ├── 15_paired_vs_unpaired_check.R  why the unpaired run reports more genes
│   ├── 16_interaction_illustration.R  what the organ x age interaction tests
│   ├── 17_validate_against_backup_sheet.py  checks every annotation against the
│   │                                   complete four-batch submission workbook
│   ├── run_de.slurm                   steps 3-5 as a SLURM job
│   ├── run_pca.slurm                  step 2 as a SLURM job
│   ├── run_paired_check.slurm         step 15 as a SLURM job
│   └── run_all.sh
├── data/                              count matrix, sample table, cached objects
├── results/tables/                    full + significant-only results (tsv)
├── results/figures/                   300 dpi PNGs
├── report/                            the deck
└── logs/
```

## Environment

`Env1` is the only conda env here with DESeq2, tximport and python-pptx together:

```
$(conda info --base)/envs/Env1      # envs/Env1.yml
```

DESeq2 fits one small GLM per gene, so the BLAS is pinned to a single thread and
parallelism is done over genes with `BiocParallel` — letting OpenBLAS spread each
tiny fit over all 128 cores is dramatically slower.

**Run the DE step on a compute node.** The login node's reaper kills it a few
minutes in, with no error in the log. `run_all.sh` submits `run_de.slurm`
(set your `--partition` and `--account` there) for exactly this reason; the whole job
takes about three minutes there.
