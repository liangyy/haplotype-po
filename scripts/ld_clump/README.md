# About

We perform LD clumping on GWAS.

Run on CRI to distribute the computation.

# Workflow

1. Using the SNP map table, extract the set of variants in GWAS that appear in our genotype cohort.
2. Have the LD panel prepared.
3. LD clump on the subset GWAS.

# Neale's lab UK Biobank GWAS

* Step 1: Subset GWAS

* Step 2: Prepare LD panel

To build reference LD panel, I reuse pipeline build [here](https://github.com/liangyy/ptrs-ukb/tree/master/pipeline/subset_bgen)
But I modified the snakemake so that it does not limit on variant.

```
screen -dmS subset_bgen bash run_subset_bgen.screen  
```

Clean up test.

```
export CHRNUM=4
screen -dmS cleanup_test bash run_cleanup.screen
```

Clean up run all chromosomes

```
screen -dmS cleanup_all bash submit_cleanup.sh
```

__Note__: For Neale's lab UK Biobank GWAS, I reuse pipeline built [here](https://github.com/liangyy/ptrs-ukb/tree/master/pipeline/ld_clump). 
But I copied and modified from the original to fix the missing `ref-first` so that we can get `chr:pos:ref:alt` which matches with Neale's lab variant ID.


* Step 3: LD clumping

For illustration, I run on one chromosome for test. 
Run on all Neale's lab GWASs and use chromosome 16.

```
screen -dmS ld_clump_chr bash -c "bash submit_ld_clump_one_chromosome.sh \
  <(ls /vol/bmd/yanyul/UKB/neale_lab_gwas/|sed 's#.tsv.bgz##g') \
  16" 
```

