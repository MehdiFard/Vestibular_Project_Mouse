# Template for running FUSION TWAS using the weights produced by Pantry Pheast
# This example is based on GWAS data from https://zenodo.org/records/3629742

import pandas as pd

datasets = ['saccule_w8', 'saccule_w22','utricle_w8', 'utricle_w22']
# datasets = ['utricle_w22']
modalities = ['alt_polyA', 'alt_TSS', 'RNA_editing', 'expression', 'isoforms', 'splicing', 'stability']
# modalities = ['expression']
trait_df = pd.read_csv('../03_pantry/pheast/input/gwas_sumstat/gwas_metadata.tsv', sep='\t')
traits = trait_df['trait'].tolist()
gwas_n = {traits[i]: n for i, n in enumerate(trait_df['sample_size'].tolist())}
# chroms = [f'chr{i + 1}' for i in range(22)]
chroms = [f'{i + 1}' for i in range(20)]


wildcard_constraints:
    modality = "[^.]+",
    dataset = "[a-zA-Z0-9_-]+"

# rule all:
#     input:
#         expand('{dataset}/output/{modality}/{modality}.{trait}.tsv',
#         dataset=datasets, modality=modalities, trait=traits),
#         expand("{dataset}/output/{modality}/twas_hits.tsv", dataset=datasets, modality=modalities)
#         # expand('{dataset}/output/twas_hits.tsv', dataset=datasets)


# rule assoc_test_trait:
#     """Run all chromosomes for a trait together to submit fewer jobs"""
#     input:
#         sumstats = '../03_pantry/pheast/input/gwas_sumstat/{trait}.sumstats',
#         weights = lambda wildcards: f'../03_pantry/pheast/intermediate/{wildcards.dataset}/twas/{wildcards.modality}.pos',
#         ref_ld_chr = expand('../03_pantry/pheast/input/LDREF/mouse_mm10_{chrom}.{ext}', chrom=chroms, ext=['bed', 'bim', 'fam']),
#     output:
#         '{dataset}/intermediate/{modality}/{trait}/{modality}.{trait}.{chrom}.tsv'
#     params:
#         output_dir = lambda w: f"{w.dataset}/intermediate/{w.modality}/{w.trait}",
#         weights_dir = lambda w: f"../03_pantry/pheast/intermediate/{w.dataset}/twas/",
#         ref_ld_prefix = '../03_pantry/pheast/input/LDREF/mouse_mm10_',
#         chroms = ' '.join(chroms),
#         coloc_p = '5e-8', # Hits with p below this are tested for coloc, recommended to be TWAS threshold. A threshold of 5e-8 / 6 is recommended for 6 modalities, but this will run coloc for everything below 5e-8 in case the values are useful.
#         gwas_n = lambda w: gwas_n[w.trait],
#     group: 'twas'
#     threads: 16
#     resources:
#         mem_mb = 32000,
#     shell:
#         """
#         mkdir -p {params.output_dir}
#         parallel -j {threads} \
#             Rscript ../03_pantry/pheast/scripts/fusion_twas/FUSION.assoc_test.R \
#                 --sumstats {input.sumstats} \
#                 --weights {input.weights} \
#                 --weights_dir {params.weights_dir} \
#                 --ref_ld_chr {params.ref_ld_prefix} \
#                 --chr {{}} \
#                 --coloc_P {params.coloc_p} \
#                 --GWASN {params.gwas_n} \
#                 --out {params.output_dir}/{wildcards.modality}.{wildcards.trait}.{{}}.tsv \
#             ::: {params.chroms}
#         """


rule all:
    input:
        expand("{dataset}/intermediate/{modality}/{trait}/{modality}.{trait}.{chrom}.tsv",
            modality=modalities, trait=traits, chrom=chroms, dataset=datasets),
            expand("{dataset}/output/{modality}/twas_hits.tsv", dataset=datasets, modality=modalities),


rule assoc_test_chr:
    input:
        sumstats = "../03_pantry/pheast/input/gwas_sumstat/{trait}.sumstats",
        weights = "../03_pantry/pheast/intermediate/{dataset}/twas/{modality}.pos",
        ref_ld = lambda wildcards: expand("../03_pantry/pheast/input/LDREF/mouse_mm10_{chrom}.{ext}", chrom=[wildcards.chrom], ext=['bed', 'bim', 'fam'])
    output:
        "{dataset}/intermediate/{modality}/{trait}/{modality}.{trait}.{chrom}.tsv"
    params:
        weights_dir = "../03_pantry/pheast/intermediate/{dataset}/twas/",
        ref_ld_prefix = "../03_pantry/pheast/input/LDREF/mouse_mm10_",
        coloc_p = "5e-8",
        gwas_n = lambda w: gwas_n[w.trait]
    resources:
        mem_mb = 4000
    shell:
        """
        mkdir -p {wildcards.dataset}/intermediate/{wildcards.modality}/{wildcards.trait}
        Rscript ../03_pantry/pheast/scripts/fusion_twas/FUSION.assoc_test.R \
            --sumstats {input.sumstats} \
            --weights {input.weights} \
            --weights_dir {params.weights_dir} \
            --ref_ld_chr {params.ref_ld_prefix} \
            --chr {wildcards.chrom} \
            --coloc_P {params.coloc_p} \
            --GWASN {params.gwas_n} \
            --out {output}
        """


rule assemble_assoc_tests:
    # one group's chromosomes only (fixed October 2026: expanding over all datasets merged the four groups)
    input:
        expand('{{dataset}}/intermediate/{{modality}}/{{trait}}/{{modality}}.{{trait}}.{chrom}.tsv', chrom=chroms)
    output:
        out = '{dataset}/output/{modality}/{modality}.{trait}.tsv'
    params:
        output_dir = '{dataset}/output/{modality}',
    group: 'twas'
    shell:
        """
        mkdir -p {params.output_dir}
        head -n1 {input[0]} | cut -f3- > {output.out}
        for fname in {input}; do
            tail -n+2 $fname | cut -f3- >> {output.out}
        done
        """

rule assemble_hits_all_traits_modalities:
    # one group's tables only (same fix); column 18 = TWAS.P after cut -f3- (tab-separated; FUSION pads fields with
    # spaces, e.g. "     NA" for untested genes, so the value is trimmed and must be a number)
    input:
        expand('{{dataset}}/output/{modality}/{modality}.{trait}.tsv', modality=modalities, trait=traits)
    output:
        out = '{dataset}/output/{modality}/twas_hits.tsv'
    params:
        traits = ' '.join(traits),
        modalities = ' '.join(modalities),
        n_modalities = len(modalities),
    shell:
        """
        head -n1 {input[0]} \
            | sed 's/^/TRAIT\tMODALITY\t/' \
            > {output.out}
        for trait in {params.traits}; do
            for modality in {params.modalities}; do
                tail -n+2 {wildcards.dataset}/output/$modality/$modality.$trait.tsv \
                    | awk -F'\t' '{{p = $18; gsub(/ /, "", p)}} p ~ /^[-+0-9.eE]+$/ && p + 0 < 5e-8 / {params.n_modalities}' \
                    | sed "s/^/$trait\t$modality\t/" \
                    >> {output.out}
            done
        done
        """

### END