# Eukaryotic_5mC_Evolution
 Files and code associated to Sarre et al 2025



 The data and code used for the *Acanthamoeba* / *Naegleria* / *Cyanophora* methylome project.
 
 EM-seq and RNA-seq raw data can be found here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249241
 
 This repository includes processed files and is divided into 5 folders:

## DNMT_phylogeny
 This folder contains the DNMT sequences, alignments and trees from EukProt, Giant Virus DB, Mirusvirus, Asgardarchaeota and other sequences.
 
## CytidineAnalogueTreatments
 This folder contains the DEseq2 differential expression results for the TElocal RNA-seq counts on the *N. gruberi* and *A. castellanii* cytidine analogue treatments. 
 
## LGTanalysis
 This folder contains the DIAMOND ten best NCBI non-redundant database hits for a set of eukaryotes (tax.tsv), their gene methylation levels (gene_meth.tsv), and their potential taxonomic origin category. 
 