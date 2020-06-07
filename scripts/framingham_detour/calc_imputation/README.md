Imputation will be performed chromosome by chromosome.

About the file format: 

* observed expression: indiv x gene (SampleID column records individual id)
* covariate: indiv x covariate (SampleID column records individual id)
* predicted expression: haplo_indiv x gene (FID, IID, columns records haplo-individual id) and the id with `_h1` is for haplotype 1 and the one with `_h2` is for haplotype 2
* pedigree: one row for one child-father-mother tuple (`shareid`, `fshare`, `mshare`)

