# Eukaryotic_5mC_Evolution
 Files and code associated to Sarre et al 2025: Repressive Cytosine Methylation is a marker of Viral Gene Transfer across distant eukaryotes. [Molecular Biology and Evolution](https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msaf176/8213644)



 The data and code used for the *Acanthamoeba* / *Naegleria* / *Cyanophora* methylome project.
 
 EM-seq and RNA-seq raw data can be found in GEO: [GSE249241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287846)
 
 This repository includes processed files and is divided into 3 folders:

## DNMT_phylogeny
 This folder contains the DNMT sequences, alignments and trees from EukProt, Giant Virus DB, Mirusvirus, Asgardarchaeota and other sequences.
 
## CytidineAnalogueTreatments
 This folder contains the DEseq2 differential expression results for the TElocal RNA-seq counts on the *N. gruberi* and *A. castellanii* cytidine analogue treatments. 
 
## LGTanalysis
 This folder contains the DIAMOND ten best NCBI non-redundant database hits for a set of eukaryotes (tax.tsv), their gene methylation levels (gene_meth.tsv), and their potential taxonomic origin category. 
 
