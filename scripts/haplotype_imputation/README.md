# Idea 1: using PRS for external GWAS

Given PRS results and observed phenotypes, the parent of origin of the haplotypes.
 
* Inputs:
    - PRS matrix (haplotype 1 and 2).
    - Observed phenotypes (paternal and maternal).
    - A map to link the PRS phenotypes and observed phenotypes.
*  Output:
    - Pr(Z) per-individual representing the probability of haplotype 1 is from father.

# Idea 2: fit PRS on-the-fly

* Inputs:
    - Genotype matrix.
    - A list of SNP to work with, they could be trait-specific.
    - Observed phenotypes.
    - A map to link the PRS phenotypes and observed phenotypes.
* Output:
    - Pr(Z) per-individual representing the probability of haplotype 1 is from father.
