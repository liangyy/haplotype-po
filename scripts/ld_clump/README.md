# About

We perform LD clumping on GWAS.

Run on CRI to distribute the computation.

# Workflow

1. Using the SNP map table, extract the set of variants in GWAS that appear in our genotype cohort.
2. Have the LD panel prepared.
3. LD clump on the subset GWAS.

# Neale's lab UK Biobank GWAS

For Neale's lab UK Biobank GWAS, I reuse pipeline built [here](https://github.com/liangyy/ptrs-ukb/tree/master/pipeline/ld_clump). 
But I copied and modified from the original to fix the missing `ref-first` so that we can get `chr:pos:ref:alt` which matches with Neale's lab variant ID.

