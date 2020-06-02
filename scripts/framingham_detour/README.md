Here we take a detour on analyzing Framingham data. The goal is to evaluate the performance of the imputation scheme.

For Framingham data, we will make use of the family data where we observe genotypes of the parents and offsprings. 
For phenotype, we will use the microarray based expression data.

For book keeping, here are some of the related data files.

* Imputed genotypes: `/gpfs/data/im-lab/nas40t2/Data/dbGaP/Transcriptome/Framingham/imputed_hrc1.1`
* Pedigree data: `/gpfs/data/im-lab/nas40t2/Data/Framingham/40031/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v23.p8.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v23.pht000183.v10.p8.Framingham_Pedigree.MULTI.txt.gz`
* Micorarray data (expression matrix): `/gpfs/data/im-lab/nas40t2/Data/dbGaP/Transcriptome/Framingham/apt-gene/rma-sketch.summary.txt`
* Linking microarray sample ID to the sample ID in Framingham: `/gpfs/data/im-lab/nas40t2/Data/Framingham/43832/PhenoGenotypeFiles/ChildStudyConsentSet_phs000363.Framingham.v12.p9.c2.HMB-IRB-NPU-MDS/ExpressionFiles/phe000002.v5.FHS_SABRe_project3.sample-info.MULTI/phe000002.v5_release_manifest.txt` 

