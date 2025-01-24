##### Load libraries and functions #####

source("loadFunctions.R")

##### Org info ######

CGMAP_PATH = "C_par/Cyanophora_EMseq.merged.formatted.CGmap.gz"
ORGANISM_FASTA_PATH = "C_par/Cyapar1_AssemblyScaffolds.fasta"
width_of_context = 5
FASTA_INDEX_PATH = "C_par/Cyapar1_AssemblyScaffolds_pUC19_lambda.fasta.fai"

CG_LOCI <- build_cg_loci_from_fasta_mod(ORGANISM_FASTA_PATH, width_of_context, positiveCntrlPath = "pUC19.fasta", negativeCntrlPath = "lambda.fasta")

SLENGTHS <- index2slengths(FASTA_INDEX_PATH)

methylationData <- loadMethylationData()

CGMAP_SUMMARY <- methylationData[[1]]
CG_METH_ARRAY <- methylationData[[2]]

CG_BS <- build_bsseq(SLENGTHS, CGMAP_SUMMARY, CG_METH_ARRAY)

org_CG_BS <- CG_BS[!seqnames(CG_BS) %in% c("M77789.2", "NC_001416.1")]
CG_BS_CG <- org_CG_BS[rowRanges(org_CG_BS)$dinucleotide == "CG",]
globalMethylation <- sum(assays(CG_BS_CG)$M, na.rm = TRUE) / sum(assays(CG_BS_CG)$Cov, na.rm = TRUE)
CG_BS_CHG <- org_CG_BS[substr(rowRanges(org_CG_BS)$context, 4, 4) != "g" & substr(rowRanges(org_CG_BS)$context, 5, 5) == "g",]
CG_BS_CHH <- org_CG_BS[substr(rowRanges(org_CG_BS)$context, 4, 4) != "g" & substr(rowRanges(org_CG_BS)$context, 5, 5) != "g",]

pUC19_BS <- CG_BS[seqnames(CG_BS) == "M77789.2" & rowRanges(CG_BS)$dinucleotide == "CG"]
lambda_BS <- CG_BS[seqnames(CG_BS) == "NC_001416.1"]

falsePositiveRate <- 100 * sum(assays(lambda_BS)$M, na.rm = TRUE) / sum(assays(lambda_BS)$Cov, na.rm = TRUE)
positiveControlRate <- 100 * sum(assays(pUC19_BS)$M, na.rm = TRUE) / sum(assays(pUC19_BS)$Cov, na.rm = TRUE)
CG_rate <- 100 * sum(assays(CG_BS_CG)$M, na.rm = TRUE) / sum(assays(CG_BS_CG)$Cov, na.rm = TRUE)
CHG_rate <- 100 * sum(assays(CG_BS_CHG)$M, na.rm = TRUE) / sum(assays(CG_BS_CHG)$Cov, na.rm = TRUE)
CHH_rate <- 100 * sum(assays(CG_BS_CHH)$M, na.rm = TRUE) / sum(assays(CG_BS_CHH)$Cov, na.rm = TRUE)


export_cg_map(CG_BS_CG, "C_par/C_par_CG.CGmap")
export_cg_map(CG_BS_CHG, "C_par/C_par_CHG.CGmap")
export_cg_map(CG_BS_CHH, "C_par/C_par_CHH.CGmap")

singleBaseMethDistribution <- as.data.frame(assays(CG_BS_CG[assays(CG_BS_CG)$Cov >= 10,])$Fraction)
singleBaseMethDistribution <- ggplot(singleBaseMethDistribution, mapping = aes(x = V1)) + 
  geom_histogram(aes(y = after_stat(density)), color = "black", fill ="grey", breaks = seq(0, 1, by = 1/24)) + 
  theme_classic() +
  xlab("mCG / CG") + ylab("Density")

ggsave("C_par/C_par_cytosine_meth_ratios.pdf", singleBaseMethDistribution, width = 3, height = 3, units = "in")

##### Make dinucleotide plots #####

globalDinucleotideLevels <- data.frame(Dinucleotide = c("CA", "CC", "CG", "CT"),
                                       mC = c(dinucleotideLevels(org_CG_BS, "CA"),
                                              dinucleotideLevels(org_CG_BS, "CC"),
                                              dinucleotideLevels(org_CG_BS, "CG"),
                                              dinucleotideLevels(org_CG_BS, "CT")))

dinucleotidePlot <- ggplot(globalDinucleotideLevels, mapping = aes(x = Dinucleotide, y = mC*100)) +
  geom_col(fill = "#000000") + theme_bw() + ylab("mC (%)") + geom_hline(yintercept = falsePositiveRate, linetype="dashed", color = "red")

ggsave("C_par/C_par_dinucleotide_mC.pdf", dinucleotidePlot, width = 3, height = 3, units = "in")


##### Run kmer analysis #####

CG_kmer_tree <- generate_kmer_tree("c", width_of_context)

CG_whole_genome_results <- extract_kmer_stats(kmer_tree=CG_kmer_tree, BS=org_CG_BS, width_of_context) # CG_whole_genome_results <- readRDS(file.path(SCRATCH_DIR, "CG_whole_genome_results.rds"))
CG_whole_genome_results$mCG <- CG_whole_genome_results$mCG*100
CG_whole_genome_results <- CG_whole_genome_results[order(CG_whole_genome_results$mCG, decreasing = TRUE),]

CG_whole_genome_results$plus1 <- substr(CG_whole_genome_results$context, 4, 4)
CG_whole_genome_results$plus2 <- substr(CG_whole_genome_results$context, 5, 5)
CG_whole_genome_results$minus1 <- substr(CG_whole_genome_results$context, 2, 2)
CG_whole_genome_results$minus2 <- substr(CG_whole_genome_results$context, 1, 1)

#Need to find if there is lambda reads in the CGmap
# falsePositiveRate

dinucleotide <- CG_whole_genome_results %>% group_by(plus1) %>% summarise(weighted_mCG = sum(mCG * n) / sum(n))

ggplot(dinucleotide, mapping = aes(x = plus1, y = weighted_mCG)) +
  geom_col() + 
  # geom_hline(yintercept = falsePositiveRate, color = "red", linetype = "dashed") + 
  theme_bw() + xlab("CX") + ylab("mC (%)")

CG_global_level <- sum(assays(CG_BS_CG)$M, na.rm = TRUE) / sum(assays(CG_BS_CG)$Cov, na.rm = TRUE)
CHG_global_level <- sum(assays(CG_BS_CHG)$M, na.rm = TRUE) / sum(assays(CG_BS_CHG)$Cov, na.rm = TRUE)
CHH_global_level <- sum(assays(CG_BS_CHH)$M, na.rm = TRUE) / sum(assays(CG_BS_CHH)$Cov, na.rm = TRUE)

contextLevels <- data.frame(context = c("CG", "CHG", "CHH"),
                            mC = c(CG_global_level, CHG_global_level, CHH_global_level))

contextPlot <- ggplot(contextLevels, mapping = aes(x = context, y = mC*100)) +
  geom_col() + xlab("Cytosine context") + ylab ("mC (%)") + theme_bw()
contextPlot
ggsave("C_par/C_par_context.pdf", contextPlot, units = "in", height = 2.5, width = 2)

CG_whole_genome_results <- CG_whole_genome_results[order(CG_whole_genome_results$context)]
CG_whole_genome_results$angle_start <- seq(0, 360-360/(4^(width_of_context-1)), 360/(4^(width_of_context-1)))
CG_whole_genome_results$angle_end <- seq(360/(4^(width_of_context-1)), 360, 360/(4^(width_of_context-1)))
CG_whole_genome_results[order(CG_whole_genome_results$mCG, decreasing = TRUE),]

biggestMCG <- max(CG_whole_genome_results$mCG)
ggplot(data = CG_whole_genome_results, mapping = aes(x = context, y = mCG)) +
  ylim(-10,biggestMCG) +
  # Bars around the circle in polar coordinates
  geom_col() +
  coord_polar() +  # Use polar coordinates
  theme_bw()  # Remove background and axes for a clean look

ggplot(data = CG_whole_genome_results, mapping = aes(x = plus1, color = minus1, y = mCG)) +
  geom_boxplot() + theme_bw() + xlab("Context") + ylab("mC (%)")

ggplot(data = CG_whole_genome_results, mapping = aes(x = plus1, color = plus2, y = mCG)) +
  geom_boxplot() + theme_bw() + xlab("Context") + ylab("mC (%)")



##### Find genome composition: #####

#Get gene GenomicRanges object from GFF3 file:
txdb <- makeTxDbFromGFF("C_par/Cyapar1_GeneCatalog_genes_20200807.gff3", format="gff3")
CParGenes <- genes(txdb)
repeatGR <- loadRepeatAnnotation(REPEAT_GTF_PATH = "C_par/Cyapar1_AssemblyScaffolds.fasta.out.gtf.kimura")

only_promoters <- GenomicFeatures::promoters(txdb, upstream = 500, 0, columns = "gene_id")
only_introns <- GenomicFeatures::intronicParts(txdb)
only_exons <- GenomicFeatures::exons(txdb)

#Trim promoter granges down to bits that don't overlap with introns or exons or repeats
trimPromoters <- function(promoterGranges, only_genes){
  #Remove strand information
  strand(promoterGranges) <- "*"
  #Simplify granges to describe non-overlapping regions
  promoterGranges <- trim(GenomicRanges::reduce(promoterGranges))
  
  #Do the same for only_genes
  strand(only_genes) <- "*"
  only_genes <- trim(GenomicRanges::reduce(only_genes))
  
  return(GenomicRanges::setdiff(promoterGranges, only_genes))
}
trimmedPromoters <- trimPromoters(only_promoters, only_introns)
trimmedPromoters <- trimPromoters(trimmedPromoters, only_exons)
trimmedPromoters <- trimPromoters(trimmedPromoters, repeatGR)

#Trim introns down to bits that aren't exons or repeats:
trimmedIntrons <- trimPromoters(only_introns, only_exons)
trimmedIntrons <- trimPromoters(trimmedIntrons, repeatGR)

#Trim exons down to bits that aren't repeats:
trimmedExons <- trimPromoters(only_exons, repeatGR)

#Create genome Granges, then trim down to bits that are not promoters, exons, introns, or repeats
genome_space <- GRanges(seqnames = names(SLENGTHS),
                        ranges = IRanges(start = 1, end = SLENGTHS))[-c(705, 706)] # Example chromosome lengths

nonRepetitiveGenomeSpace <- trimPromoters(genome_space, trimmedPromoters)
nonRepetitiveGenomeSpace <- trimPromoters(nonRepetitiveGenomeSpace, trimmedExons)
nonRepetitiveGenomeSpace <- trimPromoters(nonRepetitiveGenomeSpace, trimmedIntrons)
nonRepetitiveGenomeSpace <- trimPromoters(nonRepetitiveGenomeSpace, repeatGR)

#Get mCG stats for a list of grange objects:

GR_list <- list(trimmedPromoters, trimmedExons, trimmedIntrons, repeatGR, nonRepetitiveGenomeSpace)
GR_names <- c("Promoters", "Exons", "Introns", "Repeats", "Intergenic (non-repetitive)")

data <- lapply(seq_along(GR_list), function(x){
  strand(GR_list[[x]]) <- "*"
  GR_region <- GenomicRanges::reduce(GR_list[[x]])
  BS_region <- subsetByOverlaps(CG_BS_CG, GR_region)
  fractionDt <- data.table(Feature = GR_names[x],
                           mcg = sum(assays(BS_region)$M, na.rm = TRUE)/sum(assays(BS_region)$Cov, na.rm = TRUE))
  return(fractionDt)
})

totalFractionData <- rbindlist(data)


totalFractionData$Feature <- factor(totalFractionData$Feature, levels=c("Promoters", "Exons", "Introns", "Repeats", "Intergenic (non-repetitive)"))

totalFractionPlot <- ggplot(data = totalFractionData, mapping = aes (x = Feature, y = (100*mcg), fill = Feature)) +
  geom_col() + ylab("m (%)") + xlab("Genomic feature") +
  geom_hline(yintercept=CG_rate, linetype="dashed", color = "black") + 
  geom_hline(yintercept=falsePositiveRate, linetype="dashed", color = "red") + 
  coord_cartesian(ylim=c(0, 50)) + ylab("mCG (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(values = c("#209410", "#97b831", "#f0db65", "#b30900", "#8097ff"))
totalFractionPlot
ggsave("C_par/C_par_mainFeatures_mCG.pdf", totalFractionPlot, height = 4, width = 5, units = "in")

genomicSize <- function(GR){
  strand(GR) <- "*"
  GR <- GenomicRanges::reduce(GR)
  GR <- trim(GR)
  return(sum(width(GR)))
}

genomeComposition <- data.frame(Feature = c("Promoters", "Exons", "Introns", "Repeats", "Intergenic (non-repetitive)"),
                                Size = c(genomicSize(trimmedPromoters),
                                         genomicSize(trimmedExons),
                                         genomicSize(trimmedIntrons),
                                         genomicSize(repeatGR),
                                         genomicSize(nonRepetitiveGenomeSpace)),
                                Genome = "Genome")

genomeComposition$Feature <- factor(totalFractionData$Feature, levels=c("Intergenic (non-repetitive)", "Repeats", "Introns", "Exons", "Promoters"))
genomeCompositionPlot <- ggplot(genomeComposition, mapping = aes(y = (Size/1000000), fill = Feature, x = Genome)) +
  geom_col() + ylab("Genome size (mb)") + coord_flip() + theme_bw() + 
  scale_fill_manual(values = c("#8097ff", "#b30900", "#f0db65", "#97b831", "#209410"))
genomeCompositionPlot
ggsave("C_par/C_par_mainFeatures_composition.pdf", genomeCompositionPlot, height = 1.5, width = 6, units = "in")

##### Genes analysis (including LGT): #####

#Get gene GenomicRanges object from GFF3 file:
txdb <- makeTxDbFromGFF("C_par/Cyapar1_GeneCatalog_genes_20200807.gff3", format="gff3")
CParGenes <- GenomicFeatures::genes(txdb)
CParGenes$gene_id <- sapply(strsplit(CParGenes$gene_id, "\\|"), function(x) x[3])

mcStats <- get_mc_stats_from_granges(BS = CG_BS_CG, GR = CParGenes)
CHG_mcStats <- get_mc_stats_from_granges(BS = CG_BS_CHG, GR = CParGenes)
CHH_mcStats <- get_mc_stats_from_granges(BS = CG_BS_CHH, GR = CParGenes)

#Get methylation levels per gene:
geneMethData <- data.table(gene_id = CParGenes$tx_name, 
                           mCG = as.vector(mcStats[[4]]),
                           met_CG = as.vector(mcStats[[2]]),
                           cov_CG = as.vector(mcStats[[3]]),
                           CG_sites = as.vector(mcStats[[1]]),
                           mCHG = as.vector(CHG_mcStats[[4]]),
                           met_CHG = as.vector(CHG_mcStats[[2]]),
                           cov_CHG = as.vector(CHG_mcStats[[3]]),
                           CHG_sites = as.vector(CHG_mcStats[[1]]),
                           mCHH = as.vector(CHH_mcStats[[4]]),
                           met_CHH = as.vector(CHH_mcStats[[2]]),
                           cov_CHH = as.vector(CHH_mcStats[[3]]),
                           CHH_sites = as.vector(CHH_mcStats[[1]]))
write.table(geneMethData, file = "C_par/C_par_gene_meth.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

#Load in taxonomy data for C. par genes:
taxDataPath <- "C_par/Cyapar1_with_taxonomy.txt"
taxData <- fread(taxDataPath, select = c(1, 11, 18, 19, 20, 21, 22, 23, 24))
taxData <- taxData[taxData$V23 != "Cyanophora"]
taxData$V1 <- sub("mRNA_", "", taxData$V1)
taxData <- taxData[taxData$V18 %in% c("Eukaryota", "Archaea", "Bacteria", "Viruses"),]
taxData$V18 <- ifelse(taxData$V18 %in% c("Archaea", "Bacteria"), "Prokaryota", taxData$V18)
gene_taxonomies <- data.table(gene_id = NULL, domains = NULL)
for (gene in unique(taxData$V1)) {
  subsetTaxData <- taxData[taxData$V1 == gene]
  domain <- paste(sort(unique(subsetTaxData$V18), decreasing = TRUE), collapse = ", ")
  newRow <- data.frame(gene_id = gene, domains = domain)
  gene_taxonomies <- rbind(gene_taxonomies, newRow)
}

write.table(gene_taxonomies, file = "C_par/C_par_gene_taxonomic_categories.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
categorySummary <- gene_taxonomies %>% group_by(domains) %>% summarise(C_par = n())
write.table(categorySummary, file = "C_par/C_par_categorySummary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#Add expression data to table:
expressionData <- fread("C_par/Cyapar_EBI.genes.results", select = c(1,6))
expressionData$gene_id <- sub("gene_", "", expressionData$gene_id)
geneMethData <- merge(geneMethData, expressionData, by = "gene_id", all = TRUE)

#See distribution of methylation across genes:
geneMethData
ggplot(geneMethData, mapping = aes(x = mCG*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + geom_vline(xintercept = 50, color = "red", linetype = "dashed")

#Describe genes as methylated or unmethylated using a cut-off, e.g. 50%
geneMethData$methylationStatus <- ifelse(geneMethData$mCG >= 0.5, "Methylated", "Unmethylated")
geneMethData$methylationStatus <- factor(geneMethData$methylationStatus, levels = c("Unmethylated", "Methylated"))

#Visualise relationship between mCG and TPM:
ggplot(data = geneMethData[geneMethData$TPM != 0,], mapping = aes(x = mCG, y = log2(TPM))) +
  geom_point(alpha = 0.1) + theme_bw() + xlab("mCG (%)") + ylab("log2(TPM)")

expressedGenes <- geneMethData[geneMethData$TPM != 0,]
expressedGenes$decile <- ntile(expressedGenes$TPM, 10)
nonExpressedGenes <- geneMethData[geneMethData$TPM == 0,]
nonExpressedGenes$decile <- 0
geneMethData <- rbind(expressedGenes, nonExpressedGenes)
geneMethData$decile <- as.character(geneMethData$decile)
geneMethData$decile <- factor(geneMethData$decile, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

#Export gene beds, split by decile:
for (i in c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) {
  decileTable <- geneMethData[geneMethData$decile == i,]
  decile_granges <- CParGenes[CParGenes$gene_id %in% decileTable$gene_id]
  export.bed(decile_granges, paste0("C_par/", i, "_decile_genes.bed"))
}

mCG_tpm_decile <- ggplot(geneMethData, mapping = aes(x = decile, y = mCG * 100, group = decile)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("TPM decile") + ylab("mCG (%)")

ggsave("C_par/C_par_mCG_tpm_decile.pdf", mCG_tpm_decile, width = 3, height = 3, units = "in")

#Add LGT description to methylation data.
gene_taxonomies <- merge(gene_taxonomies, geneMethData, "gene_id")
taxData <- merge(taxData, geneMethData, by.x = "V1", by.y = "gene_id")

#Check distribution of methylation over genes for which we could assign domain of origin
ggplot(gene_taxonomies, mapping = aes(x = mCG*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + geom_vline(xintercept = 50, color = "red", linetype = "dashed")
## Result: Enriched in unmethylated genes

totalGeneTaxonomyFrequency <- gene_taxonomies %>% group_by(methylationStatus, domains) %>% summarise(frequency = n()) %>%
  mutate(frequency_proportion = frequency / sum(frequency)) %>% ungroup() %>% as.data.frame()
totalGeneTaxonomyFrequency <- totalGeneTaxonomyFrequency[!is.na(totalGeneTaxonomyFrequency$methylationStatus),]
# LGT_mCG_col <- gene_taxonomies %>% group_by(domains) %>% summarise(total_mCG = sum(met)/sum(cov), count = n()) %>% ggplot(mapping = aes(x = domains, y = total_mCG*100, fill = domains)) +
#   geom_col() + 
#   geom_hline(yintercept = falsePositiveRate*100, linetype = "dashed", color = "red") +
#   geom_hline(yintercept = globalLevel*100, linetype = "dashed", color = "black") +
#   theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
#   ylab("mCG (%)") + xlab("Domain of origin") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# LGT_mCG_col
# ggsave("figures/acanthamoeba/initial_emseq/LGT/LGT_mCG_col.pdf", LGT_mCG_col, height = 10, width = 13, units = "cm")

count_data <- gene_taxonomies %>%
  group_by(domains) %>%
  summarise(count = sum(!is.na(mCG)))

LGT_mCG_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = mCG*100, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("mCG (%)") + xlab("Domain of origin") + 
  # coord_cartesian(ylim=c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  
  geom_text(data = count_data, aes(x = domains, y = 100, label = count), vjust = -0.5)
LGT_mCG_bars

count_data <- gene_taxonomies %>%
  group_by(domains) %>%
  summarise(count = sum(!is.na(mCHG)))

LGT_mCHG_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = mCHG*100, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("mCHG (%)") + xlab("Domain of origin") + 
  coord_cartesian(ylim=c(0, 5.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  
  geom_text(data = count_data, aes(x = domains, y = 5, label = count), vjust = -0.5)
LGT_mCHG_bars

count_data <- gene_taxonomies %>%
  group_by(domains) %>%
  summarise(count = sum(!is.na(mCHH)))

LGT_mCHH_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = mCHH*100, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("mCHH (%)") + xlab("Domain of origin") + 
  coord_cartesian(ylim=c(0, 4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  
  geom_text(data = count_data, aes(x = domains, y = 3.8, label = count), vjust = -0.5)
LGT_mCHH_bars

ggsave("C_par/C_par_LGT_bars_CG.pdf", LGT_mCG_bars, height = 10, width = 13, units = "cm")
ggsave("C_par/C_par_LGT_bars_CHG.pdf", LGT_mCHG_bars, height = 10, width = 13, units = "cm")
ggsave("C_par/C_par_LGT_bars_CHH.pdf", LGT_mCHH_bars, height = 10, width = 13, units = "cm")

count_data <- gene_taxonomies %>%
  group_by(domains) %>%
  summarise(count = sum(!is.na(TPM)))

LGT_TPM_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = TPM, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("TPM") + xlab("Domain of origin") + coord_cartesian(ylim=c(0, 60)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  
  geom_text(data = count_data, aes(x = domains, y = 52, label = count), vjust = -0.5)
LGT_TPM_bars
ggsave("C_par/Cpar_TPM_bars.pdf", LGT_mCG_bars, height = 10, width = 13, units = "cm")

ggplot(gene_taxonomies, mapping = aes(x = mCG)) +
  geom_histogram(bins = 50, fill = "lightgrey", color = "black") +
  theme_bw()

ggplot(gene_taxonomies, mapping = aes(x = mCG, fill = domains)) +
  geom_histogram(bins = 50, color = "black") +
  theme_bw()

sum(totalGeneTaxonomyFrequency[totalGeneTaxonomyFrequency$methylationStatus == "Methylated",]$frequency)
sum(totalGeneTaxonomyFrequency[totalGeneTaxonomyFrequency$methylationStatus == "Unmethylated",]$frequency)

count_data <- totalGeneTaxonomyFrequency %>%
  group_by(methylationStatus) %>%
  summarise(count = sum(frequency))

LGT <- ggplot() +
  geom_col(totalGeneTaxonomyFrequency, mapping = aes(x = methylationStatus, y = frequency_proportion*100, fill = domains)) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("Proportion (%)") + xlab("Methylation status (>50% mCG)") +
  
  geom_text(data = count_data, aes(x = methylationStatus, y = 100, label = count), vjust = -0.5)
LGT
ggsave("C_par/C_par_LGT.pdf", LGT, height = 15.14, width = 13.69, units = "cm")

# bestMethylatedViralHits <- taxData[taxData$V18 == "Viruses",] %>%
#   group_by(V1) %>%
#   slice_min(order_by = V11, with_ties = FALSE) %>%
#   ungroup() %>%
#   filter(!is.na(methylationStatus) & methylationStatus == "Methylated")
# 
# bestMethylatedViralHits %>% ggplot(mapping = aes(x = V19)) +
#   geom_bar() + theme_bw() + xlab("Viral phylum of methylated genes") + ylab("Count")
# 
# bestMethylatedViralHits %>% filter(V19 == "Nucleocytoviricota") %>% ggplot(mapping = aes(x = V22)) +
#   geom_bar() + theme_bw() + xlab("Viral order of Nucleocytoviricota methylated genes") + ylab("Count")
# 
# bestMethylatedViralHits %>% filter(V19 == "Nucleocytoviricota") %>% ggplot(mapping = aes(x = V23)) +
#   geom_bar() + theme_bw() + xlab("Viral genus of Nucleocytoviricota methylated genes") + ylab("Count")



#Calculate enrichment of groups within 
data_grouped <- totalGeneTaxonomyFrequency %>% 
  group_by(domains) %>% 
  summarise(
    # Total_Frequency = sum(frequency),
    Frequency_50pc_mCG = sum(frequency[methylationStatus == "Methylated"]),
    Frequency_Unmethylated = sum(frequency[methylationStatus == "Unmethylated"]),
    Proportion_50pc_mCG = sum(frequency_proportion[methylationStatus == "Methylated"]),
    Proportion_Unmethylated = sum(frequency_proportion[methylationStatus == "Unmethylated"])
  )
total_50pc_mCG <- sum(data_grouped$Frequency_50pc_mCG)
total_Unmethylated <- sum(data_grouped$Frequency_Unmethylated)

data_grouped <- data_grouped %>%
  rowwise() %>%
  mutate(p_value = fisher.test(matrix(c(Frequency_50pc_mCG,
                                        total_50pc_mCG - Frequency_50pc_mCG,
                                        Frequency_Unmethylated,
                                        total_Unmethylated - Frequency_Unmethylated),
                                      nrow = 2), alternative = "two.sided")$p.value) %>%
  ungroup()
data_grouped
write.table(data_grouped, file = "C_par/C_par_methylated_category_enrichment.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

##### Repeat analysis #####
repeatGR <- loadRepeatAnnotation(REPEAT_GTF_PATH = "C_par/Cyapar1_AssemblyScaffolds.fasta.out.gtf.kimura")

repeatFamilyContributions <- data.frame()
for (family in unique(repeatGR$family)) {
  familyGranges <- repeatGR[repeatGR$family == family]
  export.bed(familyGranges, paste0("C_par/C_par_", family, "_repeats.bed"))
  repeatFamilyContributions <- rbind(repeatFamilyContributions, data.frame(family = family, bases = sum(width(repeatGR[repeatGR$family == family]))))
}
genomeSize <- sum(SLENGTHS[1:(length(SLENGTHS)-2)])

repeatFamilyContributions$contribution <- 100*repeatFamilyContributions$bases/genomeSize
repeatFamilyContributions <- repeatFamilyContributions[repeatFamilyContributions$family %in% c("LINE", "DNA", "LTR", "Unknown"),]
repeatFamilyContributionsPlot <- ggplot(repeatFamilyContributions, mapping = aes(x = family, y = contribution)) +
  geom_col(fill = "black") + ylab("Genome %") + xlab("Repeat sub-class") + theme_bw()
repeatFamilyContributionsPlot
ggsave("C_par/C_par_repeatFamilyContribution.pdf", plot = repeatFamilyContributionsPlot, width = 6, height = 6, units = "cm")

repeatMeth <- get_mc_stats_from_granges(CG_BS_CG, repeatGR)
repeatMethDF <- data.frame(repeat_id = repeatGR$gene_id,
                           family = repeatGR$family,
                           kimura = repeatGR$kimura,
                           mCG = as.vector(repeatMeth[[4]]),
                           M = as.vector(repeatMeth[[2]]),
                           Cov = as.vector(repeatMeth[[3]]))

repeatMethDF <- repeatMethDF[repeatMethDF$family %in% c("LINE", "LTR", "DNA"),]

kimura_histogram_split_by_family <- ggplot(data=repeatMethDF, mapping = aes(x = kimura)) +
  geom_histogram(binwidth=1, fill = "black") +
  facet_wrap(family ~ ., scales = "free_y") +
  theme_bw() +
  xlab("Kimura score") +
  ylab("Count")
kimura_histogram_split_by_family
ggsave(filename = "C_par/C_par_kimura_histogram_split_by_family.pdf", plot = kimura_histogram_split_by_family, width = 5, height = 2, units = "in")

repeatSubclassMethylation <- repeatMethDF %>% group_by(family) %>% 
  summarise(mCG = 100*(sum(M, na.rm=TRUE)/sum(Cov, na.rm=TRUE)),
            n = n()) %>% 
  ggplot(mapping = aes(x = family, y = mCG)) +
  geom_col(fill = "black") + theme_bw() + xlab("Repeat sub-class") + 
  ylab("mCG (%)") + 
  geom_hline(yintercept = globalMethylation*100, color = "red", linetype = "dashed") +
  geom_text(aes(x = family, y = mCG + 3, label = n))
repeatSubclassMethylation
ggsave(filename = "C_par/C_par_repeatSubclassMethylation.pdf", plot = repeatSubclassMethylation, width = 2, height = 2, units = "in")

repeatMethDF <- repeatMethDF[!is.na(repeatMethDF$kimura),]

repeatMethDF$kimura_range <- cut_kimura_into_ranges(repeatMethDF$kimura)
repeatSubclassKimuramCG <- repeatMethDF %>% group_by(family, kimura_range) %>% 
  summarise(mCG = 100*(sum(M, na.rm=TRUE)/sum(Cov, na.rm=TRUE)),
            n = n()) %>% 
  ggplot(mapping = aes(x = kimura_range, y = mCG)) +
  geom_col(fill = "black") +
  facet_wrap(. ~ family) + theme_bw() + xlab("Repeat sub-class") + 
  ylab("mCG (%)") + 
  geom_hline(yintercept = globalMethylation*100, color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(aes(x = kimura_range, y = mCG + 4, label = n))
repeatSubclassKimuramCG
ggsave(filename = "C_par/C_par_repeatSubclassKimuramCG.pdf", plot = repeatSubclassKimuramCG, width = 8, height = 4, units = "in")

kimura_range_vs_methylation_split_by_family <- ggplot(data=repeatMethDF, mapping=aes(x=kimura_range, y=mCG)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_violin(alpha=0) +
  facet_wrap(. ~ family) +
  theme_bw() +
  xlab("Kimura score") +
  ylab("Methylation (%)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
kimura_range_vs_methylation_split_by_family
ggsave(filename = "C_par/C_par_repeatSubclassKimura_box.pdf", plot = kimura_range_vs_methylation_split_by_family, width = 8, height = 4, units = "in")

########## Create ideogram ########

sliding_window <- function(SLENGTHS,
                           BS,
                           WINDOW_SIZE){
  ORG_SLENGTHS <- SLENGTHS[!names(SLENGTHS) %in% c("NC_001416.1", "M77789.2")]
  #Create GRanges object describing each organism scaffold
  contig_gr <- GRanges(seqnames = names(ORG_SLENGTHS), ranges = IRanges(start = 1, end = ORG_SLENGTHS), seqlengths = ORG_SLENGTHS)
  
  #Create GRanges object describing each scaffold, split into windows of a set width
  sliding_windows_150bp_GR <- slidingWindows(contig_gr, width = WINDOW_SIZE, step = WINDOW_SIZE) %>% unlist()
  
  #Get methylation statistics based on this sliding windows GRanges object
  Sliding_windows_mC <- get_mc_stats_from_granges(BS, sliding_windows_150bp_GR)
  
  colData <- DataFrame(Treatment=c("Static"),
                       row.names=c("Static"))
  SE <- SummarizedExperiment(assays=list(mCG_met = as.matrix(Sliding_windows_mC[[2]]),
                                         mCG_cov = as.matrix(Sliding_windows_mC[[3]]),
                                         mCG_fraction=as.matrix(Sliding_windows_mC[[4]])),
                             rowRanges=sliding_windows_150bp_GR, colData=colData)
  rowRanges(SE)$C_position_per_gene <- Sliding_windows_mC[[1]]
  return(SE)
}

only_genes <- CParGenes
strand(only_genes) <- "*"

#Load in taxonomy data for A. cas genes:

taxData <- taxData[taxData$V18 %in% c("Viruses"),]

viralGenes <- only_genes[only_genes$gene_id %in% taxData$V1]

virus <- data.table(Type = "Viral",
                    Shape = "box",
                    Chr=as.vector(seqnames(viralGenes)),
                    Start=start(viralGenes),
                    End=end(viralGenes),
                    color = "000000")

# virus$Chr <- sub("^QPMI01_", "", virus$Chr)

karyotype<-data.table(Chr=names(SLENGTHS),Start=0,End=SLENGTHS)
# karyotype$Chr <- sub("^QPMI01_", "", karyotype$Chr)
mostViralChrs <- (virus %>% group_by(Chr) %>% summarise(n = n()) %>% arrange(desc(n)))[1:12,]$Chr
# largestChrs <- (karyotype %>% arrange(desc(End)))[1:24,]$Chr

karyotype <- karyotype[karyotype$Chr %in% mostViralChrs,]

methylation_windows <- sliding_window(SLENGTHS[names(SLENGTHS) %in% karyotype$Chr], CG_BS_CG[seqnames(CG_BS_CG) %in% karyotype$Chr,], 500)

addDensity <- function(SE, gr){
  densities <- c()
  SE_gr <- granges(SE)
  windowWidths <- width(SE_gr)
  iterations <- length(SE_gr)
  for (i in seq_len(iterations)) {
    regionsInRange <- GenomicRanges::intersect(SE_gr[i], gr)
    overlapWithGR <- sum(width(regionsInRange))
    densities <- c(densities, overlapWithGR)
  }
  return(densities/windowWidths)
}

rowRanges(methylation_windows)$geneDensity <- addDensity(methylation_windows, only_genes)

methylation <- data.table(Chr=as.vector(seqnames(methylation_windows)),
                          Start=start(methylation_windows),
                          End=end(methylation_windows),
                          Value=as.vector(assays(methylation_windows)$mCG_fraction)*100,
                          Color="000000")

geneDensity <- data.table(Chr=as.vector(seqnames(methylation_windows)),
                          Start=start(methylation_windows),
                          End=end(methylation_windows),
                          Value=rowRanges(methylation_windows)$geneDensity*100,
                          Color="000000")

# methylation <- methylation[!is.na(methylation$Value),]

# methylation$Chr <- sub("^scaffold_", "", methylation$Chr)
# geneDensity$Chr <- sub("^scaffold_", "", geneDensity$Chr)
# repeatDensity$Chr <- sub("^scaffold_", "", repeatDensity$Chr)
# karyotype$Chr <- sub("^scaffold_", "", karyotype$Chr)

# ideogram(karyotype, 
#          overlaid = geneDensity, 
#          label = repeatDensity,
#          label_type = "heatmap",
#          output = file.path(OUTPUT_DIR, paste("ideogram_gene_repeat.svg", sep="")), 
#          colorset1 = c("#f7f7f7", "#e34a33"), 
#          colorset2 = c("#f7f7f7", "#2c7fb8"))
# 
# ideogram(karyotype, 
#          overlaid = geneDensity, 
#          label = methylation,
#          label_type = "line",
#          output = file.path(OUTPUT_DIR, paste("ideogram_geneDensity_mCG.svg", sep="")), 
#          colorset1 = c("#f7f7f7", "#e34a33"))
# 
# ideogram(karyotype, 
#          overlaid = repeatDensity, 
#          label = methylation,
#          label_type = "line",
#          output = file.path(OUTPUT_DIR, paste("ideogram_repeatDensity_mCG.svg", sep="")), 
#          colorset1 = c("#f7f7f7", "#e34a33"))
# 
# 
# ideogram(karyotype, 
#          overlaid=repeatDensity, 
#          label = virus,
#          label_type = "marker",
#          output = file.path(OUTPUT_DIR, "ideogram_virus_repeatDensity.svg"), 
#          colorset1 = c("#ffffff", "#b88dff", "#3700ff"))

virus <- virus[virus$Chr %in% karyotype$Chr,]
virus[1,]$Shape <- "circle"
virus[2,]$Shape <- "triangle"
virus[1,]$color <- "080808"
virus[2,]$color <- "090909"
virus[1,]$Type <- "test"
virus[2,]$Type <- "fdasfd"


ideogram(karyotype, 
         overlaid=methylation, 
         label = virus,
         label_type = "marker",
         output = file.path("./", "ideogram_virus_mCG.svg"), 
         colorset1 = c("#ffffff", "#ff9f9b", "#ff0041"))

ideogram(karyotype, 
         overlaid=geneDensity, 
         label = virus,
         label_type = "marker",
         output = file.path(OUTPUT_DIR, "ideogram_virus_geneDensity.svg"), 
         colorset1 = c("#ffffff", "#e7da91", "#c2b800"))

ideogram(karyotype = karyotype, 
         overlaid = repeatDensity, 
         label = methylation, 
         label_type = "polygon", 
         colorset1 = c("#ffffff", "#b88dff", "#3700ff"),
         output = file.path(OUTPUT_DIR, "ideogram_repeatDensity_methLine.svg"))

# convertSVG(file.path(OUTPUT_DIR, paste("ideogram_virus_m", title, ".svg", sep="")), file = file.path(OUTPUT_DIR, paste("ideogram_virus_m", title, ".png", sep="")), device = "png")
