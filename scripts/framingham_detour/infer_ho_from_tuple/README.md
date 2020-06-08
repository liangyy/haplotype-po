Here we infer haplotype origin from phased genotypes of family tuple.

* We first extract the phased genotypes of these individuals (using `bcftools`).
* And then we run a python script to calculate all possible haplotype origin configuration and select the "best" as the inferred one.

