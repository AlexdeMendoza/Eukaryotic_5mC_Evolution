source("loadFunctions.R")

CGMAP_PATH = "A_cas/EMSeq_A_cas_merged.CGmap.gz"
ORGANISM_FASTA_PATH = "A_cas/Neff_assembly.fa"
width_of_context = 5
FASTA_INDEX_PATH = "A_cas/Neff_assembly_pUC19_lambda.fasta.fai"

##### Load methylation data: #####

CG_LOCI <- build_cg_loci_from_fasta_mod(ORGANISM_FASTA_PATH, width_of_context, positiveCntrlPath = "pUC19.fasta", negativeCntrlPath = "lambda.fasta")
SLENGTHS <- index2slengths(FASTA_INDEX_PATH)

methylationData <- loadMethylationData()

CGMAP_SUMMARY <- methylationData[[1]]
CG_METH_ARRAY <- methylationData[[2]]

CG_BS <- build_bsseq(SLENGTHS, CGMAP_SUMMARY, CG_METH_ARRAY)
A_cas_BS <- CG_BS[!seqnames(CG_BS) %in% c("M77789.2", "NC_001416.1")]
CG_BS_CG <- A_cas_BS[rowRanges(A_cas_BS)$dinucleotide == "CG",]
export_cg_map(CG_BS_CG, "A_cas/A_cas_CG.CGmap")

singleBaseMethDistribution <- as.data.frame(assays(CG_BS_CG[assays(CG_BS_CG)$Cov >= 10 & assays(CG_BS_CG)$Fraction >= 1/24,])$Fraction)
singleBaseMethDistribution <- ggplot(singleBaseMethDistribution, mapping = aes(x = V1)) + 
  geom_histogram(aes(y = after_stat(density)), color = "black", fill ="grey", breaks = seq(0, 1, by = 1/24)) + 
  theme_classic() +
  xlab("mCG / CG") + ylab("Density")
singleBaseMethDistribution
ggsave("A_cas/A_cas_cytosine_meth_ratios.pdf", singleBaseMethDistribution, width = 3, height = 3, units = "in")

##### Run kmer analysis #####

CG_kmer_tree <- generate_kmer_tree("c", width_of_context)

CG_whole_genome_results <- extract_kmer_stats(kmer_tree=CG_kmer_tree, BS=, width_of_context)
CG_whole_genome_results$mCG <- CG_whole_genome_results$mCG*100
CG_whole_genome_results <- CG_whole_genome_results[order(CG_whole_genome_results$mCG, decreasing = TRUE),]

CG_whole_genome_results$plus1 <- substr(CG_whole_genome_results$context, 4, 4)
CG_whole_genome_results$plus2 <- substr(CG_whole_genome_results$context, 5, 5)
CG_whole_genome_results$minus1 <- substr(CG_whole_genome_results$context, 2, 2)
CG_whole_genome_results$minus2 <- substr(CG_whole_genome_results$context, 1, 1)

dinucleotide <- CG_whole_genome_results %>% group_by(plus1) %>% summarise(weighted_mCG = sum(mCG * n) / sum(n))

ggplot(dinucleotide, mapping = aes(x = plus1, y = weighted_mCG)) +
  geom_col() + 
  # geom_hline(yintercept = falsePositiveRate, color = "red", linetype = "dashed") + 
  theme_bw() + xlab("CX") + ylab("mC (%)")

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

##### Genes analysis #####

#Get gene GenomicRanges object from GFF3 file:
GENE_GTF_PATH <- "A_cas/Neff_annotations.gtf"
txdb <- makeTxDbFromGFF(GENE_GTF_PATH, format="gtf")

genes <- GenomicFeatures::genes(txdb)

mcStats <- get_mc_stats_from_granges(BS = CG_BS_CG, GR = genes)

#Get methylation levels per gene:
geneMethData <- data.table(gene_id = genes$gene_id, mCG = as.vector(mcStats[[4]]), 
                           met = as.vector(mcStats[[2]]), 
                           cov = as.vector(mcStats[[3]]), 
                           CG_sites = as.vector(mcStats[[1]]), 
                           gene_length = width(genes))
write.table(geneMethData, file = "A_cas/A_cas_gene_meth.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

#Load in taxonomy data for genes:
taxDataPath <- "A_cas/A_cas_Neff_tax.tsv"
taxData <- fread(taxDataPath, select = c(1, 11, 18, 19, 20, 21, 22, 23, 24))
taxData$V1 <- gsub(".mRNA.1", "", taxData$V1)
taxData <- taxData[taxData$V23 != "Acanthamoeba"]
taxData <- taxData[taxData$V18 %in% c("Eukaryota", "Archaea", "Bacteria", "Viruses"),]
taxData$V18 <- ifelse(taxData$V18 %in% c("Archaea", "Bacteria"), "Prokaryota", taxData$V18)
gene_taxonomies <- data.table(gene_id = NULL, domains = NULL)
for (gene in unique(taxData$V1)) {
  subsetTaxData <- taxData[taxData$V1 == gene]
  domain <- paste(sort(unique(subsetTaxData$V18), decreasing = TRUE), collapse = ", ")
  newRow <- data.frame(gene_id = gene, domains = domain)
  gene_taxonomies <- rbind(gene_taxonomies, newRow)
}

write.table(gene_taxonomies, file = "A_cas/A_cas_gene_taxonomic_categories.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
categorySummary <- gene_taxonomies %>% group_by(domains) %>% summarise(A_cas = n())
write.table(categorySummary, file = "A_cas/A_cas_categorySummary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#Add expression data to table:
expressionData <- fread("A_cas/Acas_EBI.genes.results", select = c(1,6))
geneMethData <- merge(geneMethData, expressionData, by = "gene_id", all = TRUE)

#See distribution of methylation across genes:
cutOff <- 0.2
ggplot(geneMethData, mapping = aes(x = mCG*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + geom_vline(xintercept = cutOff*100, color = "red", linetype = "dashed")

#Describe genes as methylated or unmethylated using a cut-off, e.g. 0.5%
geneMethData$methylationStatus <- ifelse(geneMethData$mCG >= cutOff, "Methylated", "Unmethylated")
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
  decile_granges <- genes[genes$gene_id %in% decileTable$gene_id]
  export.bed(decile_granges, paste0("A_cas/A_cas_", i, "_decile_genes.bed"))
}

mCG_tpm_decile <- ggplot(geneMethData, mapping = aes(x = decile, y = mCG * 100, group = decile)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("TPM decile") + ylab("mCG (%)")
mCG_tpm_decile
ggsave("A_cas/A_cas_mCG_tpm_decile.pdf", mCG_tpm_decile, width = 3, height = 3, units = "in")

#Add LGT description to methylation data.
gene_taxonomies <- merge(gene_taxonomies, geneMethData, "gene_id")
taxData <- merge(taxData, geneMethData, by.x = "V1", by.y = "gene_id")
categorySummaryMethylation <- gene_taxonomies %>% group_by(domains, methylationStatus) %>% summarise(A_cas = n())
write.table(categorySummaryMethylation, file = "A_cas/A_cas_categorySummaryMethylation.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#Check distribution of methylation over genes for which we could assign domain of origin
ggplot(gene_taxonomies, mapping = aes(x = mCG*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + geom_vline(xintercept = cutOff*100, color = "red", linetype = "dashed")

totalGeneTaxonomyFrequency <- gene_taxonomies %>% group_by(methylationStatus, domains) %>% summarise(frequency = n()) %>%
  mutate(frequency_proportion = frequency / sum(frequency)) %>% ungroup() %>% as.data.frame()
totalGeneTaxonomyFrequency <- totalGeneTaxonomyFrequency[!is.na(totalGeneTaxonomyFrequency$methylationStatus),]

falsePositiveRate <- sum(assays(CG_BS_CG[seqnames(CG_BS_CG) == "NC_001416.1"])$M, na.rm = TRUE) / sum(assays(CG_BS_CG[seqnames(CG_BS_CG) == "NC_001416.1"])$Cov, na.rm = TRUE)
globalLevel <- sum(assays(CG_BS_CG[!seqnames(CG_BS_CG) %in% c("NC_001416.1", "M77789.2")])$M, na.rm = TRUE) / sum(assays(CG_BS_CG[!seqnames(CG_BS_CG) %in% c("NC_001416.1", "M77789.2")])$Cov, na.rm = TRUE)
LGT_mCG_col <- gene_taxonomies %>% group_by(domains) %>% summarise(total_mCG = sum(met, na.rm = TRUE)/sum(cov, na.rm = TRUE), count = n()) %>% ggplot(mapping = aes(x = domains, y = total_mCG*100, fill = domains)) +
  geom_col() +
  geom_hline(yintercept = falsePositiveRate*100, linetype = "dashed", color = "red") +
  geom_hline(yintercept = globalLevel*100, linetype = "dashed", color = "black") +
  theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("mCG (%)") + xlab("Domain of origin") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
LGT_mCG_col
ggsave("A_cas/A_casLGT_mCG_col.pdf", LGT_mCG_col, height = 10, width = 13, units = "cm")

count_data <- gene_taxonomies %>%
  group_by(domains) %>%
  summarise(count = sum(!is.na(mCG)))

LGT_mCG_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = mCG*100, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("mCG (%)") + xlab("Domain of origin") + coord_cartesian(ylim=c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  
  geom_text(data = count_data, aes(x = domains, y = 100, label = count), vjust = -0.5)
LGT_mCG_bars
ggsave("A_cas/A_cas_LGT_bars.pdf", LGT_mCG_bars, height = 10, width = 13, units = "cm")

count_data <- gene_taxonomies %>%
  group_by(domains) %>%
  summarise(count = sum(!is.na(TPM)))

LGT_TPM_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = TPM, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("TPM") + xlab("Domain of origin") + coord_cartesian(ylim=c(0, 120)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  
  geom_text(data = count_data, aes(x = domains, y = 100, label = count), vjust = -0.5)
LGT_TPM_bars
ggsave("A_cas/A_cas_TPM_boxs.pdf", LGT_TPM_bars, height = 10, width = 13, units = "cm")

count_data <- totalGeneTaxonomyFrequency %>%
  group_by(methylationStatus) %>%
  summarise(count = sum(frequency))

LGT <- ggplot() +
  geom_col(data = totalGeneTaxonomyFrequency, mapping = aes(x = methylationStatus, y = frequency_proportion*100, fill = domains)) + 
  theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("Proportion (%)") + xlab(paste0("Methylation status (>", cutOff*100, "% mCG)")) +
  geom_text(data = count_data, aes(x = methylationStatus, y = 100, label = count), vjust = -0.5)
LGT
ggsave("A_cas/A_cas_LGT.pdf", LGT, height = 15.14, width = 13.69, units = "cm")


#Calculate enrichment of groups within 
data_grouped <- totalGeneTaxonomyFrequency %>% 
  group_by(domains) %>% 
  summarise(
    # Total_Frequency = sum(frequency),
    Frequency_methylated = sum(frequency[methylationStatus == "Methylated"]),
    Frequency_unmethylated = sum(frequency[methylationStatus == "Unmethylated"]),
    Proportion_methylated = sum(frequency_proportion[methylationStatus == "Methylated"]),
    Proportion_unmethylated = sum(frequency_proportion[methylationStatus == "Unmethylated"])
  )
total_methylated <- sum(data_grouped$Frequency_methylated)
total_unmethylated <- sum(data_grouped$Frequency_unmethylated)

data_grouped <- data_grouped %>%
  rowwise() %>%
  mutate(p_value = fisher.test(matrix(c(Frequency_methylated,
                                        total_methylated - Frequency_methylated,
                                        Frequency_unmethylated,
                                        total_unmethylated - Frequency_unmethylated),
                                      nrow = 2), alternative = "two.sided")$p.value) %>%
  ungroup()
data_grouped
write.table(data_grouped, file = "A_cas/A_cas_methylated_category_enrichment.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

##### Repeat analysis #####
repeatGR <- loadRepeatAnnotation(REPEAT_GTF_PATH = "A_cas/Neff_assembly.fa.out.gtf.kimura")
repeatGR <- repeatGR[width(repeatGR) >= 500]
repeatFamilyContributions <- data.frame()
for (family in unique(repeatGR$family)) {
  familyGranges <- repeatGR[repeatGR$family == family]
  export.bed(familyGranges, paste0("A_cas/A_cas_", family, "_repeats.bed"))
  repeatFamilyContributions <- rbind(repeatFamilyContributions, data.frame(family = family, bases = sum(width(repeatGR[repeatGR$family == family]))))
}
genomeSize <- sum(SLENGTHS[1:(length(SLENGTHS)-2)])

repeatFamilyContributions$contribution <- 100*repeatFamilyContributions$bases/genomeSize
repeatFamilyContributions <- repeatFamilyContributions[repeatFamilyContributions$family %in% c("LINE", "DNA", "LTR", "Unknown"),]
repeatFamilyContributionsPlot <- ggplot(repeatFamilyContributions, mapping = aes(x = family, y = contribution)) +
  geom_col(fill = "black") + ylab("Genome %") + xlab("Repeat sub-class") + theme_bw()
repeatFamilyContributionsPlot
ggsave("A_cas/A_cas_repeatFamilyContribution.pdf", plot = repeatFamilyContributionsPlot, width = 6, height = 6, units = "cm")

repeatMeth <- get_mc_stats_from_granges(CG_BS_CG, repeatGR)
repeatMethDF <- data.frame(repeat_id = repeatGR$gene_id,
                           family = repeatGR$family,
                           kimura = repeatGR$kimura,
                           mCG = as.vector(repeatMeth[[4]]),
                           M = as.vector(repeatMeth[[2]]),
                           Cov = as.vector(repeatMeth[[3]]))

repeatMethDF_filter <- repeatMethDF[repeatMethDF$family %in% c("LINE", "LTR", "DNA"),]

kimura_histogram_split_by_family <- ggplot(data=repeatMethDF_filter, mapping = aes(x = kimura)) +
  geom_histogram(binwidth=1, fill = "black") +
  facet_wrap(family ~ ., scales = "free_y") +
  theme_bw() +
  xlab("Kimura score") +
  ylab("Count")
kimura_histogram_split_by_family
ggsave(filename = "A_cas/A_cas_kimura_histogram_split_by_family.pdf", plot = kimura_histogram_split_by_family, width = 5, height = 2, units = "in")

repeatSubclassMethylation <- repeatMethDF_filter %>% group_by(family) %>% 
  summarise(mCG = 100*(sum(M, na.rm=TRUE)/sum(Cov, na.rm=TRUE)),
            n = n()) %>% 
  ggplot(mapping = aes(x = family, y = mCG)) +
  geom_col(fill = "black") + theme_bw() + xlab("Repeat sub-class") + 
  ylab("mCG (%)") + 
  geom_hline(yintercept = globalLevel*100, color = "red", linetype = "dashed") +
  geom_text(aes(x = family, y = mCG + 0.4, label = n))
repeatSubclassMethylation
ggsave(filename = "A_cas/A_cas_repeatSubclassMethylation.pdf", plot = repeatSubclassMethylation, width = 2, height = 2, units = "in")

repeatMethDF_filter <- repeatMethDF_filter[!is.na(repeatMethDF_filter$kimura),]

repeatMethDF_filter$kimura_range <- cut_kimura_into_ranges(repeatMethDF_filter$kimura)
repeatSubclassKimuramCG <- repeatMethDF_filter %>% group_by(family, kimura_range) %>% 
  summarise(mCG = 100*(sum(M, na.rm=TRUE)/sum(Cov, na.rm=TRUE)),
            n = n()) %>% 
  ggplot(mapping = aes(x = kimura_range, y = mCG)) +
  geom_col(fill = "black") +
  facet_wrap(. ~ family) + theme_bw() + xlab("Repeat sub-class") + 
  ylab("mCG (%)") + 
  geom_hline(yintercept = globalLevel*100, color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(aes(x = kimura_range, y = mCG + 0.8, label = n))
repeatSubclassKimuramCG
ggsave(filename = "A_cas/A_cas_repeatSubclassKimuramCG.pdf", plot = repeatSubclassKimuramCG, width = 8, height = 4, units = "in")

kimura_range_vs_methylation_split_by_family <- ggplot(data=repeatMethDF_filter, mapping=aes(x=kimura_range, y=mCG)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_violin(alpha=0) +
  facet_wrap(. ~ family) +
  theme_bw() +
  xlab("Kimura score") +
  ylab("Methylation (%)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
kimura_range_vs_methylation_split_by_family
ggsave(filename = "A_cas/A_cas_repeatSubclassKimura_box.pdf", plot = kimura_range_vs_methylation_split_by_family, width = 8, height = 4, units = "in")

##### This section is to filter the repeat gtf for mapping expression data: ######

repeat_gtf <- fread("A_cas/Neff_assembly.fa.out.gtf.kimura")
  
repeat_gr <- GRanges(seqnames = repeat_gtf$V1, ranges=IRanges(start=repeat_gtf$V4, end = repeat_gtf$V5), strand = repeat_gtf$V7, id = repeat_gtf$V9)

# 1) removing the repeats that are completely within genes
## From which genes do I want to eliminate repeats? The ones I'm SURE are genes. Therefore I'm only going to use genes that are not completely within repeats 
flattened_repeats <- GenomicRanges::reduce(repeat_gr, ignore.strand=TRUE)
genes_gtf <- fread("A_cas/Neff_annotations.gtf")
genes_gr <- GRanges(seqnames = genes_gtf$V1, ranges=IRanges(start=genes_gtf$V4, end = genes_gtf$V5), strand = genes_gtf$V7, id = genes_gtf$V9)
genes_not_entirely_covered_by_repeats <- subsetByOverlaps(genes_gr, flattened_repeats, type=c("within"), invert=TRUE, ignore.strand=TRUE)

##Now I want to exclude repeats that are totally covered by this subset of gene:
repeat_gr <- subsetByOverlaps(repeat_gr, genes_not_entirely_covered_by_repeats, type=c("within"), invert=TRUE, ignore.strand=TRUE)

## Finally, exclude those that are < 200bp:
repeat_gr <- repeat_gr[width(repeat_gr) >= 200]
repeat_gtf <- repeat_gtf[repeat_gtf$V9 %in% repeat_gr$id]
repeat_gtf <- repeat_gtf[order(repeat_gtf[,c('V1','V4','V5','V7','V3','V8')]),]
repeat_gtf <- repeat_gtf[!duplicated(repeat_gtf[,c('V1','V4','V5','V3')]),]
write.table(repeat_gtf, file = "A_cas/Neff_assembly.fa.out.gtf.kimura.morethan200.notingenes", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


##### Compare conditions ######

# Methylation level comparisons #####

drugTreatmentSamples <- c(rep("DMSO", 3), rep("Aza", 3), rep("Decitabine", 2), rep("Zebularine", 3))
drugTreatmentCGmaps <- c("A_cas/DMSO1_2_EKDL210007760-1a-7UDI280-AK890_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/DMSO2_2_EKDL210007760-1a-AK903-GH07_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/DMSO3_2_EKDL210007760-1a-AK907-AK1864_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/s5Aza1_2_EKDL210007760-1a-7UDI228-7UDI200_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/s5Aza2_2_EKDL210007760-1a-AK889-7UDI229_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/s5Aza3_2_EKDL210007760-1a-7UDI237-AK1867_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/Decitabine1_2_EKDL210007760-1a-AK1853-7UDI283_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/Decitabine2_2_EKDL210007760-1a-AK850-AK2430_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/Zebularine1_2_EKDL210007760-1a-AK845-GD06_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/Zebularine2_2_EKDL210007760-1a-AK1868-AK1852_HLV2TDSX2_L2_1.fq.gz.CGmap.gz",
                         "A_cas/Zebularine3_2_EKDL210007760-1a-AK869-AK847_HLV2TDSX2_L2_1.fq.gz.CGmap.gz")

drugTreatmentBS <- lapply(drugTreatmentCGmaps, function(x){
  drugTreatmentMethylationData <- loadMethylationData(x)
  CGMAP_SUMMARY <- drugTreatmentMethylationData[[1]]
  CG_METH_ARRAY <- drugTreatmentMethylationData[[2]]
  CG_BS <- build_bsseq(SLENGTHS, CGMAP_SUMMARY, CG_METH_ARRAY)
  CG_BS_CG <- CG_BS[rowRanges(CG_BS)$dinucleotide == "CG",]
  return(CG_BS_CG)
})

pUC19_rates <- vapply(drugTreatmentBS, FUN = function(x){
  pUC19 <- x[seqnames(x) == "M77789.2"]
  return(sum(assays(pUC19)$M, na.rm = TRUE)/sum(assays(pUC19)$Cov, na.rm = TRUE))
}, FUN.VALUE = numeric(1))

lambda_rates <- vapply(drugTreatmentBS, FUN = function(x){
  lambda <- x[seqnames(x) == "NC_001416.1"]
  return(sum(assays(lambda)$M, na.rm = TRUE)/sum(assays(lambda)$Cov, na.rm = TRUE))
}, FUN.VALUE = numeric(1))

global_rates <- vapply(drugTreatmentBS, FUN = function(x){
  A_cas <- x[!seqnames(x) %in% c("NC_001416.1", "M77789.2")]
  return(sum(assays(A_cas)$M, na.rm = TRUE)/sum(assays(A_cas)$Cov, na.rm = TRUE))
}, FUN.VALUE = numeric(1))

drugTreatmentResults <- data.frame(treatment = drugTreatmentSamples,
                                   truePositive = pUC19_rates,
                                   falseNegative = lambda_rates,
                                   uncorrectedLevels = global_rates,
                                   alexCorrection = global_rates - lambda_rates,
                                   lukeCorrection =  (global_rates-lambda_rates) / (pUC19_rates-lambda_rates))

drugTreatmentResults$treatment <- factor(drugTreatmentResults$treatment, levels = c("DMSO", "Aza", "Decitabine", "Zebularine"))

globalLevelsScatter <- ggplot(drugTreatmentResults, mapping = aes(x = treatment, y = uncorrectedLevels*100)) +
  geom_point(size = 1, position = position_dodge2(width=0.5)) + ylim(0, 2) + theme_bw() + 
  ylab("Uncorrected mCG (%)") + xlab("Replicate") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
globalLevelsScatter

ggplot(drugTreatmentResults, mapping = aes(x = treatment, y = lukeCorrection*100)) +
  geom_point(size = 1, position = position_dodge2(width=0.5)) + ylim(0, 2) + theme_bw() + 
  ylab("Uncorrected mCG (%)") + xlab("Replicate") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

t.test(drugTreatmentResults[drugTreatmentResults$treatment == "DMSO",]$lukeCorrection, drugTreatmentResults[drugTreatmentResults$treatment == "Aza",]$lukeCorrection)
t.test(drugTreatmentResults[drugTreatmentResults$treatment == "DMSO",]$lukeCorrection, drugTreatmentResults[drugTreatmentResults$treatment == "Decitabine",]$lukeCorrection)

t.test(drugTreatmentResults[drugTreatmentResults$treatment == "DMSO",]$alexCorrection, drugTreatmentResults[drugTreatmentResults$treatment == "Aza",]$alexCorrection)
t.test(drugTreatmentResults[drugTreatmentResults$treatment == "DMSO",]$alexCorrection, drugTreatmentResults[drugTreatmentResults$treatment == "Decitabine",]$alexCorrection)

wilcox.test(drugTreatmentResults[drugTreatmentResults$treatment == "DMSO",]$lukeCorrection, drugTreatmentResults[drugTreatmentResults$treatment == "Aza",]$lukeCorrection)
wilcox.test(drugTreatmentResults[drugTreatmentResults$treatment == "DMSO",]$lukeCorrection, drugTreatmentResults[drugTreatmentResults$treatment == "Decitabine",]$lukeCorrection)

wilcox.test(drugTreatmentResults[drugTreatmentResults$treatment == "DMSO",]$alexCorrection, drugTreatmentResults[drugTreatmentResults$treatment == "Aza",]$alexCorrection)
wilcox.test(drugTreatmentResults[drugTreatmentResults$treatment == "DMSO",]$alexCorrection, drugTreatmentResults[drugTreatmentResults$treatment == "Decitabine",]$alexCorrection)

##### DE analysis of genes: #####

# Load expression data, and exclude genes: #

geneMethData <- geneMethData[!is.na(geneMethData$mCG),]
methylatedIDs <- geneMethData[geneMethData$mCG >= cutOff,]$gene_id

EXPRESSION_DATA_PATH = c("A_cas/DMSO1.bam.cntTable.cleaned",
                         "A_cas/DMSO2.bam.cntTable.cleaned",
                         "A_cas/DMSO3.bam.cntTable.cleaned",
                         "A_cas/Dec1.bam.cntTable.cleaned",
                         "A_cas/Dec2.bam.cntTable.cleaned",
                         "A_cas/Dec3.bam.cntTable.cleaned",
                         "A_cas/Aza1.bam.cntTable.cleaned",
                         "A_cas/Aza2.bam.cntTable.cleaned",
                         "A_cas/Aza3.bam.cntTable.cleaned")
SAMPLES <- c("DMSO1", "DMSO2", "DMSO3", "Dec1", "Dec2", "Dec3", "Aza1", "Aza2", "Aza3")
expressionData <- load_expression_data(EXPRESSION_DATA_PATH)
geneExpressionData <- expressionData[expressionData[,ncol(expressionData)] %in% genes$gene_id,]
methylatedGeneExpressionData <- geneExpressionData[geneExpressionData$gene_id %in% methylatedIDs,]

# Make Upset showing the number of repeats expressed with tpm >=1 in different samples: #
geneExpressionUpsetData <- geneExpressionData[,1:ncol(geneExpressionData)-1]
row.names(geneExpressionUpsetData) <- geneExpressionData[,ncol(geneExpressionData)]
geneExpressionUpsetData <- geneExpressionUpsetData>0
geneExpressionUpsetData[geneExpressionUpsetData == TRUE] <- 1
geneExpressionUpsetData[geneExpressionUpsetData == FALSE] <- 0
colnames(geneExpressionUpsetData) <- SAMPLES
geneExpressionUpsetData <- as.data.frame(geneExpressionUpsetData)
upset(geneExpressionUpsetData, sets = colnames(geneExpressionUpsetData), nintersects = NA)
reducedUpsetData <- data.frame(row.names = row.names(geneExpressionUpsetData),
                               DMSO = ifelse(rowSums(geneExpressionUpsetData[,1:3])>=2, 1, 0),
                               Aza = ifelse(rowSums(geneExpressionUpsetData[,4:6])>=2, 1, 0),
                               Dec = ifelse(rowSums(geneExpressionUpsetData[,7:9])>=2, 1, 0))
treatmentUpset <- upset(reducedUpsetData, sets = colnames(reducedUpsetData), order.by = c("degree"))
pdf(file="A_cas/treatmentUpset.pdf",
    width = 5, height = 5) # or other device
treatmentUpset
dev.off()

# Make Upset showing the number of repeats expressed with tpm >=1 in different samples, but looking only at methylated genes: #
methylatedGeneExpressionUpsetData <- methylatedGeneExpressionData[,1:ncol(methylatedGeneExpressionData)-1]
row.names(methylatedGeneExpressionUpsetData) <- methylatedGeneExpressionData[,ncol(methylatedGeneExpressionData)]
methylatedGeneExpressionUpsetData <- methylatedGeneExpressionUpsetData>0
methylatedGeneExpressionUpsetData[methylatedGeneExpressionUpsetData == TRUE] <- 1
methylatedGeneExpressionUpsetData[methylatedGeneExpressionUpsetData == FALSE] <- 0
colnames(methylatedGeneExpressionUpsetData) <- SAMPLES
methylatedGeneExpressionUpsetData <- as.data.frame(methylatedGeneExpressionUpsetData)
upset(methylatedGeneExpressionUpsetData, sets = colnames(methylatedGeneExpressionUpsetData), nintersects = NA)
reducedUpsetData <- data.frame(row.names = row.names(methylatedGeneExpressionUpsetData),
                               DMSO = ifelse(rowSums(methylatedGeneExpressionUpsetData[,1:3])>=2, 1, 0),
                               Aza = ifelse(rowSums(methylatedGeneExpressionUpsetData[,4:6])>=2, 1, 0),
                               Dec = ifelse(rowSums(methylatedGeneExpressionUpsetData[,7:9])>=2, 1, 0))
treatmentUpset <- upset(reducedUpsetData, sets = colnames(reducedUpsetData), order.by = c("degree"))
pdf(file="A_cas/methylatedTreatmentUpset.pdf",
    width = 5, height = 5) # or other device
treatmentUpset
dev.off()

# Get IDs for genes which have hits against viruses: #
viralTaxData <- taxData[taxData$V18 %in% c("Viruses"),]
viralGenesIDs <- unique(viralTaxData$V1)

# Get IDs for genes which do not have orthologues in C3 strain of A. cas (indicating LGT): #
genesInC3 <- fread("A_cas/Neff_in_C3.tsv", header = TRUE, col.names = "gene_id")$gene_id
genesInC3 <- sub("\\.mRNA\\.[0-9]+$", "", genesInC3)
genesNotInC3 <- genes[!genes$gene_id %in% genesInC3]$gene_id

#Get IDs for genes within GEVEs:
GEVEs <- import.bed("A_cas/GEVEs.formatted.bed")
GEVE_gene_IDs <- subsetByOverlaps(genes, GEVEs)$gene_id

##### Run DEseq and make graphs of fold changes etc #####

pfam_numbers <- fread("A_cas/Neff_annotations.pep.noStop.1iso.PfamA.domtblout.4cols", sep = "\t", fill = TRUE, col.names = c("Pfam","PFID","gene_id","evalue")) %>% .[.$evalue < 0.001,]
pfam_numbers$gene_id <- sapply(pfam_numbers$gene_id, function(x){
  return(strsplit(x, ".mRNA")[[1]][1])
})

geneID2GO <- readMappings(file = "A_cas/acanthamoeba.eggnog.GOs", sep = " ", IDsep = ",")
names(geneID2GO) <- sub("\\.mRNA\\.[0-9]+$", "", names(geneID2GO))
geneUniverse <- names(geneID2GO)

bothTreatmentsDE <- data.table()
# For each of Aza and Dec, call DE genes and carry out eggnog enrichment:
controlSamples <- c(1, 2, 3)
#treatmentDrug <- "Aza"
for (treatmentDrug in c("Dec", "Aza")) {
  if (treatmentDrug == "Dec") {
    treatmentSamples <- c(4, 5, 6)
  }else{
    treatmentSamples <- c(7, 8, 9)
  }
  geneDEseqSummary <- summariseDESEQ(expressionData, controlSamples, treatmentSamples, SAMPLES) %>% 
    mutate(Treatment = treatmentDrug) %>%
    filter(!is.na(padj),
           !is.na(log2FoldChange),
           transcript_id %in% genes$gene_id) %>% arrange(padj) #%>% head(10)
  write.table(geneDEseqSummary, paste0("A_cas/A_cas_genes_DEseq_", treatmentDrug, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  #Do PFAM / EggNog enrichment
  upregulated_ids <- (geneDEseqSummary %>% filter(padj < 0.05, log2FoldChange > 0))$transcript_id
  upregulated_pfams <- pfamEnrichments(ids = upregulated_ids) %>% arrange(fishp)
  # View(upregulated_pfams[upregulated_pfams$fishp < 0.05,])
  upregulated_pfams$Pfam <- factor(upregulated_pfams$Pfam, levels= rev(c(upregulated_pfams$Pfam)))
  upregulated_pfams <- upregulated_pfams[1:60,]
  gg_upregulated_pfam <- ggplot(upregulated_pfams, aes(x = Pfam, y = -log10(fishp))) + 
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab("Pfam domain") +
    ylab("-log10(p-value)") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() + coord_flip()
  gg_upregulated_pfam
  ggsave(gg_upregulated_pfam, 
         filename = paste0("A_cas/A_cas_PFAMs_upregulated_", treatmentDrug, ".pdf"),
         height = 12, width = 7, units = "cm")
  
  upregulated_GOs <- get_GOs(geneNames = upregulated_ids, GOlist = geneID2GO)
  # View(upregulated_GOs[upregulated_GOs$classicFisher < 0.05 & upregulated_GOs$Ontology == "BP",])
  # View(upregulated_GOs[upregulated_GOs$classicFisher < 0.05 & upregulated_GOs$Ontology == "MF",])
  upregulated_GOs <- upregulated_GOs %>%
    group_by(Ontology) %>%
    slice_min(order_by = classicFisher, n = 30) %>%
    ungroup()
  gg_upregulated_GOs <- plot_GOs(df = upregulated_GOs, name = "Upregulated genes")
  ggsave(paste0("A_cas/A_cas_eggnogUpregulatedGO_", treatmentDrug, ".pdf"),
         plot = gg_upregulated_GOs, height = 12, width = 15)
  
  geneDEseqSummary <- merge(geneDEseqSummary, geneMethData, by.x = "transcript_id", by.y = "gene_id")
  geneDEseqSummary$NeffExclusive <- geneDEseqSummary$transcript_id %in% genesNotInC3
  geneDEseqSummary$Viral <- geneDEseqSummary$transcript_id %in% viralGenesIDs
  geneDEseqSummary$GEVE <- geneDEseqSummary$transcript_id %in% GEVE_gene_IDs
  bothTreatmentsDE <- rbind(bothTreatmentsDE, geneDEseqSummary)
  
  #mCG vs FC, colored by significance:
  mCG_expression_change <- ggplot(geneDEseqSummary, mapping = aes(x = mCG*100, y = log2FoldChange, color = padj < 0.05)) +
    geom_point(size = 1) + xlab("Untreated mCG (%)") + ylab("log2 Fold Change") + theme_bw() + geom_hline(yintercept = 0)
  mCG_expression_change
  ggsave(mCG_expression_change, 
         filename = paste0("A_cas/A_cas_scatter_", treatmentDrug, ".pdf"),
         height = 8, width = 11, units = "cm")
  
  #mCG vs FC, colored by presence in C3:
  mCG_expression_change <- ggplot(geneDEseqSummary, mapping = aes(x = mCG*100, y = log2FoldChange, color = NeffExclusive)) +
    geom_point(size = 1) + xlab("Untreated mCG (%)") + ylab("log2 Fold Change") + theme_bw() + geom_hline(yintercept = 0)
  mCG_expression_change
  ggsave(mCG_expression_change, 
         filename = paste0("A_cas/A_cas_C3presence_scatter_", treatmentDrug, ".pdf"),
         height = 8, width = 11, units = "cm")
  
  #mCG vs FC, colored by presence in GEVEs:
  mCG_expression_change <- ggplot(geneDEseqSummary, mapping = aes(x = mCG*100, y = log2FoldChange, color = GEVE)) +
    geom_point(size = 1) + xlab("Untreated mCG (%)") + ylab("log2 Fold Change") + theme_bw() + geom_hline(yintercept = 0)
  mCG_expression_change
  
  ggsave(mCG_expression_change, 
         filename = paste0("A_cas/A_cas_GEVE_genes_scatter_", treatmentDrug, ".pdf"),
         height = 8, width = 11, units = "cm")
  
  #mCG vs FC, colored by viral or not:
  mCG_expression_change <- ggplot(geneDEseqSummary, mapping = aes(x = mCG*100, y = log2FoldChange, color = Viral)) +
    geom_point(size = 1) + xlab("Untreated mCG (%)") + ylab("log2 Fold Change") + theme_bw() + geom_hline(yintercept = 0)
  mCG_expression_change
  ggsave(mCG_expression_change, 
         filename = paste0("A_cas/A_cas_viral_scatter_", treatmentDrug, ".pdf"),
         height = 8, width = 11, units = "cm")
  
  #Proportion of upregulated/downregulated, split by methylation status:
  geneDEseqSummaryBarData <- geneDEseqSummary %>%
    filter(padj < 0.05) %>%
    mutate(Direction = ifelse(log2FoldChange < 0, "Down", ifelse(log2FoldChange > 0, "Up", "Neutral")),
           ThresholdmCG = ifelse(mCG < 0.2, "< 20%", ">= 20%"),
           Direction = factor(Direction, levels = c("Down", "Up")),
           ThresholdmCG = factor(ThresholdmCG, levels = c("< 20%", ">= 20%"))) %>%
    group_by(Direction, ThresholdmCG) %>%
    summarise(count = n(), .groups = 'drop') %>%
    complete(Direction, ThresholdmCG, fill = list(count = 0))
  geneDEseqSummaryBar <- ggplot(geneDEseqSummaryBarData, mapping = aes(x = ThresholdmCG, fill = Direction, y = count)) +
    geom_col(position = "dodge") + xlab("Untreated mCG (%)") + ylab("log2 Fold Change") + theme_bw()
  geneDEseqSummaryBar
  ggsave(geneDEseqSummaryBar, 
         filename = paste0("A_cas/A_cas_bar_", treatmentDrug, ".pdf"),
         height = 8, width = 8, units = "cm")
  
  #Proportion of upregulated/downregulated, split by presence in C3:
  geneDEseqSummaryBarData <- geneDEseqSummary %>%
    filter(padj < 0.05) %>%
    mutate(Direction = ifelse(log2FoldChange < 0, "Down", ifelse(log2FoldChange > 0, "Up", "Neutral")),
           Direction = factor(Direction, levels = c("Down", "Up"))) %>%
    group_by(Direction, NeffExclusive) %>%
    summarise(count = n(), .groups = 'drop') %>%
    complete(Direction, NeffExclusive, fill = list(count = 0))
  geneDEseqSummaryBar <- ggplot(geneDEseqSummaryBarData, mapping = aes(x = NeffExclusive, fill = Direction, y = count)) +
    geom_col(position = "dodge") + xlab("Genes exclusive to Neff (not in C3)") + ylab("log2 Fold Change") + theme_bw()
  geneDEseqSummaryBar
  ggsave(geneDEseqSummaryBar, 
         filename = paste0("A_cas/A_cas_C3_bar_", treatmentDrug, ".pdf"),
         height = 8, width = 8, units = "cm")
  
  #Proportion of upregulated/downregulated, split by presence in C3:
  geneDEseqSummaryBarData <- geneDEseqSummary %>%
    filter(padj < 0.05) %>%
    mutate(Direction = ifelse(log2FoldChange < 0, "Down", ifelse(log2FoldChange > 0, "Up", "Neutral")),
           Direction = factor(Direction, levels = c("Down", "Up"))) %>%
    group_by(Direction, GEVE) %>%
    summarise(count = n(), .groups = 'drop') %>%
    complete(Direction, GEVE, fill = list(count = 0))
  geneDEseqSummaryBar <- ggplot(geneDEseqSummaryBarData, mapping = aes(x = GEVE, fill = Direction, y = count)) +
    geom_col(position = "dodge") + xlab("Genes in GEVE") + ylab("log2 Fold Change") + theme_bw()
  geneDEseqSummaryBar
  ggsave(geneDEseqSummaryBar, 
         filename = paste0("A_cas/A_cas_GEVE_genes_bar_", treatmentDrug, ".pdf"),
         height = 8, width = 8, units = "cm")
  
  #Proportion of upregulated/downregulated, split by whether gene is viral or not:
  geneDEseqSummaryBarData <- geneDEseqSummary %>%
    filter(padj < 0.05) %>%
    mutate(Direction = ifelse(log2FoldChange < 0, "Down", ifelse(log2FoldChange > 0, "Up", "Neutral")),
           Direction = factor(Direction, levels = c("Down", "Up"))) %>%
    group_by(Direction, Viral) %>%
    summarise(count = n(), .groups = 'drop') %>%
    complete(Direction, Viral, fill = list(count = 0))
  geneDEseqSummaryBar <- ggplot(geneDEseqSummaryBarData, mapping = aes(x = Viral, fill = Direction, y = count)) +
    geom_col(position = "dodge") + xlab("Viral genes") + ylab("log2 Fold Change") + theme_bw()
  geneDEseqSummaryBar
  ggsave(geneDEseqSummaryBar, 
         filename = paste0("A_cas/A_cas_viral_bar_", treatmentDrug, ".pdf"),
         height = 8, width = 8, units = "cm")
}

##### Compare treatments #####

bothTreatmentsDF <- as.data.frame(bothTreatmentsDE[,c(1,3,7,10,11)])

compareTreatments <- pivot_wider(bothTreatmentsDF, 
                                 names_from = Treatment, 
                                 values_from = c(log2FoldChange, padj),
                                 names_sep = "_")

bothTreatmentsDF$Direction <- ifelse(bothTreatmentsDF$log2FoldChange > 0, "Up", "Down")
bothTreatmentsDF$methylationStatus <- ifelse(bothTreatmentsDF$mCG >= 0.2, "Methylated", "Unmethylated")
bothTreatmentsDF$Treatment <- ifelse(bothTreatmentsDF$Treatment == "Aza", "Aza", "Dec")

significantGenesFC <- bothTreatmentsDF %>% filter(padj < 0.05) %>% group_by(Treatment, Direction, methylationStatus) %>% summarise(count = n(), .groups = 'drop') %>%
  complete(Direction, Treatment, methylationStatus, fill = list(count = 0)) %>% ggplot(mapping = aes(x = methylationStatus, y = count, fill = Direction)) +
  geom_col(position = "dodge") + theme_bw() + ylab("DE Genes") + facet_grid(cols = vars(Treatment)) +xlab("Methylation status")
significantGenesFC
ggsave("A_cas/significantGenes_UpDown.pdf", significantGenesFC,
       height = 7, width = 13, units = "cm")

compareTreatments$NeffExclusive <- compareTreatments$transcript_id %in% genesNotInC3
compareTreatments$Viral <- compareTreatments$transcript_id %in% viralGenesIDs
compareTreatments$Methylated <- compareTreatments$mCG >= 0.2
compareTreatments$GEVE <- compareTreatments$transcript_id %in% GEVE_gene_IDs

non_NAs <- compareTreatments[!is.na(compareTreatments$log2FoldChange_Dec) & !is.na(compareTreatments$log2FoldChange_Aza),]
cor(non_NAs$log2FoldChange_Dec, non_NAs$log2FoldChange_Aza, method = "pearson")

non_NAs$alpha = ifelse(pmin(non_NAs$padj_Dec, non_NAs$padj_Aza)<0.05, 1, 0.01)

#Compare Dec vs Aza change, highlighting methylated genes:
geneDecAzaMethylation <- ggplot() +
  geom_point(data = non_NAs[!non_NAs$Methylated,], mapping = aes(x = log2FoldChange_Aza,
                                                                                     y = log2FoldChange_Dec,
                                                                                     alpha = alpha),
             color = "darkblue", size = 0.5, show.legend = FALSE) + 
  geom_point(data = non_NAs[non_NAs$Methylated,], mapping = aes(x = log2FoldChange_Aza,
                                                                                    y = log2FoldChange_Dec,
                                                                                    alpha = alpha),
             color = "red", size = 0.5, show.legend = FALSE) +
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw()
geneDecAzaMethylation

ggsave("A_cas/A_cas_Aza_vs_Dec_methylation_FC_sig_eitherTrt.pdf",
       geneDecAzaMethylation,
       height = 9, width = 9, units = "cm")

#Compare Dec vs Aza change, highlighting methylated (>20% mCG) genes:
geneDecAzaMethylation <- ggplot() +
  geom_point(data = compareTreatments[!compareTreatments$Methylated,], mapping = aes(x = log2FoldChange_Aza,
                                                                                     y = log2FoldChange_Dec),
             color = "darkblue", size = 0.5) + 
  geom_point(data = compareTreatments[compareTreatments$Methylated,], mapping = aes(x = log2FoldChange_Aza,
                                                                                    y = log2FoldChange_Dec),
             color = "red", size = 0.5) +
  geom_hline(yintercept = mean(compareTreatments[compareTreatments$Methylated,]$log2FoldChange_Dec, na.rm = TRUE), color = "red", linetype = "dashed") + 
  geom_vline(xintercept = mean(compareTreatments[compareTreatments$Methylated,]$log2FoldChange_Aza, na.rm = TRUE), color = "red", linetype = "dashed") + 
  geom_hline(yintercept = mean(compareTreatments[!compareTreatments$Methylated,]$log2FoldChange_Dec, na.rm = TRUE), color = "darkblue", linetype = "dashed") + 
  geom_vline(xintercept = mean(compareTreatments[!compareTreatments$Methylated,]$log2FoldChange_Aza, na.rm = TRUE), color = "darkblue", linetype = "dashed") + 
  xlab("log2FC 5-azacytidine | 20% mCG genes highlighted") + ylab("log2FC Decitabine") + theme_bw()
geneDecAzaMethylation
ggsave("A_cas/Aza_vs_Dec_methylation.pdf",
       geneDecAzaMethylation,
       height = 9, width = 9, units = "cm")

compareTreatments[compareTreatments$log2FoldChange_Dec > 4 & compareTreatments$log2FoldChange_Aza > 4 & compareTreatments$Methylated,]

##### PFAM / GO enrichment of methylated / viral genes #####
#Get methylated geneIDs:
methylatedGenes <- geneMethData[geneMethData$methylationStatus == "Methylated",]
methylatedGeneIDs <- methylatedGenes$gene_id

methylated_pfams <- pfamEnrichments(ids = methylatedGeneIDs) %>% arrange(fishp)
methylated_pfams[methylated_pfams$fishp < 0.05,]
methylated_pfams$Pfam <- factor(methylated_pfams$Pfam, levels= rev(c(methylated_pfams$Pfam)))
methylated_pfams <- methylated_pfams[1:60,]
gg_methylated_pfam <- ggplot(methylated_pfams, aes(x = Pfam, y = -log10(fishp))) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Pfam domain") +
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() + coord_flip()
gg_methylated_pfam
ggsave(gg_methylated_pfam, 
       filename = paste0("A_cas/A_casPFAMs_methylated_genes.pdf"),
       height = 12, width = 7, units = "cm")

methylated_GOs <- get_GOs(geneNames = methylatedGeneIDs, GOlist = geneID2GO)
# View(methylated_GOs[methylated_GOs$classicFisher < 0.05 & methylated_GOs$Ontology == "BP",])
# View(methylated_GOs[methylated_GOs$classicFisher < 0.05 & methylated_GOs$Ontology == "MF",])
methylated_GOs <- methylated_GOs %>%
  group_by(Ontology) %>%
  slice_min(order_by = classicFisher, n = 30) %>%
  ungroup()
gg_methylated_GOs <- plot_GOs(df = methylated_GOs, name = "methylated genes")
gg_methylated_GOs
ggsave("A_cas/A_cas_eggnog_methylated_genes_GO.pdf",
       plot = gg_methylated_GOs, height = 12, width = 15)

## Do for viral genes:
viral_pfams <- pfamEnrichments(ids = viralGenesIDs) %>% arrange(fishp)
# View(viral_pfams[viral_pfams$fishp < 0.05,])
viral_pfams$Pfam <- factor(viral_pfams$Pfam, levels= rev(c(viral_pfams$Pfam)))
viral_pfams <- viral_pfams[1:60,]
gg_viral_pfam <- ggplot(viral_pfams, aes(x = Pfam, y = -log10(fishp))) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Pfam domain") +
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() + coord_flip()
gg_viral_pfam
ggsave(gg_viral_pfam, 
       filename = paste0("A_cas/A_cas_PFAMs_viral_genes.pdf"),
       height = 12, width = 7, units = "cm")

viral_GOs <- get_GOs(geneNames = viralGenesIDs, GOlist = geneID2GO)
# View(viral_GOs[viral_GOs$classicFisher < 0.05 & viral_GOs$Ontology == "BP",])
# View(viral_GOs[viral_GOs$classicFisher < 0.05 & viral_GOs$Ontology == "MF",])
viral_GOs <- viral_GOs %>%
  group_by(Ontology) %>%
  slice_min(order_by = classicFisher, n = 30) %>%
  ungroup()
gg_viral_GOs <- plot_GOs(df = viral_GOs, name = "viral genes")
gg_viral_GOs
ggsave(paste0("A_cas/A_cas_eggnog_viral_genes_GO.pdf"),
       plot = gg_viral_GOs, height = 12, width = 15)

#Repeat for genes in GEVEs:

GEVE_pfams <- pfamEnrichments(ids = GEVE_gene_IDs) %>% arrange(fishp)
# GEVE_pfams[GEVE_pfams$fishp < 0.05,]
GEVE_pfams$Pfam <- factor(GEVE_pfams$Pfam, levels= rev(c(GEVE_pfams$Pfam)))
GEVE_pfams <- GEVE_pfams[1:60,]
gg_GEVE_pfam <- ggplot(GEVE_pfams, aes(x = Pfam, y = -log10(fishp))) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Pfam domain") +
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() + coord_flip()
gg_GEVE_pfam
ggsave(gg_GEVE_pfam, 
       filename = paste0("A_cas/A_cas_PFAMs_GEVE_genes.pdf"),
       height = 12, width = 7, units = "cm")

GEVE_GOs <- get_GOs(geneNames = GEVE_gene_IDs, GOlist = geneID2GO)
# GEVE_GOs[GEVE_GOs$classicFisher < 0.05 & GEVE_GOs$Ontology == "BP",]
# GEVE_GOs[GEVE_GOs$classicFisher < 0.05 & GEVE_GOs$Ontology == "MF",]
GEVE_GOs <- GEVE_GOs %>%
  group_by(Ontology) %>%
  slice_min(order_by = classicFisher, n = 30) %>%
  ungroup()
gg_GEVE_GOs <- plot_GOs(df = GEVE_GOs, name = "GEVE genes")
gg_GEVE_GOs
ggsave(paste0("A_cas/A_cas_eggnog_GEVE_genes_GO.pdf"),
       plot = gg_GEVE_GOs, height = 12, width = 15)

##### DE analysis of repeats between drug treatments: #####

methylatedRepeatIDs <- repeatMethDF[repeatMethDF$mCG >= cutOff,]$repeat_id

#Load expression data, and exclude genes:

repeatExpressionData <- expressionData[!expressionData[,ncol(expressionData)] %in% genes$gene_id,]
methylatedRepeatExpressionData <- repeatExpressionData[repeatExpressionData$gene_id %in% methylatedRepeatIDs,]

#Make Upset showing the number of repeats expressed with tpm >=1 in different samples:
repeatExpressionUpsetData <- repeatExpressionData[,1:ncol(repeatExpressionData)-1]
row.names(repeatExpressionUpsetData) <- repeatExpressionData[,ncol(repeatExpressionData)]
repeatExpressionUpsetData <- repeatExpressionUpsetData>0
repeatExpressionUpsetData[repeatExpressionUpsetData == TRUE] <- 1
repeatExpressionUpsetData[repeatExpressionUpsetData == FALSE] <- 0
colnames(repeatExpressionUpsetData) <- SAMPLES
repeatExpressionUpsetData <- as.data.frame(repeatExpressionUpsetData)
upset(repeatExpressionUpsetData, sets = colnames(repeatExpressionUpsetData))
reducedUpsetData <- data.frame(row.names = row.names(repeatExpressionUpsetData),
                               DMSO = ifelse(rowSums(repeatExpressionUpsetData[,1:3])>=3, 1, 0),
                               Aza = ifelse(rowSums(repeatExpressionUpsetData[,4:6])>=3, 1, 0),
                               Dec = ifelse(rowSums(repeatExpressionUpsetData[,7:9])>=3, 1, 0))
treatmentUpset <- upset(reducedUpsetData, sets = colnames(reducedUpsetData), order.by = c("degree"))
pdf(file="A_cas/A_cas_repeats_treatmentUpset.pdf",
    width = 5, height = 5) # or other device
treatmentUpset
dev.off()

#Make Upset showing the number of repeats expressed with tpm >=1 in different samples, using methylated TEs:
methylatedRepeatExpressionUpsetData <- methylatedRepeatExpressionData[,1:ncol(methylatedRepeatExpressionData)-1]
row.names(methylatedRepeatExpressionUpsetData) <- methylatedRepeatExpressionData[,ncol(methylatedRepeatExpressionData)]
methylatedRepeatExpressionUpsetData <- methylatedRepeatExpressionUpsetData>0
methylatedRepeatExpressionUpsetData[methylatedRepeatExpressionUpsetData == TRUE] <- 1
methylatedRepeatExpressionUpsetData[methylatedRepeatExpressionUpsetData == FALSE] <- 0
colnames(methylatedRepeatExpressionUpsetData) <- SAMPLES
methylatedRepeatExpressionUpsetData <- as.data.frame(methylatedRepeatExpressionUpsetData)
upset(methylatedRepeatExpressionUpsetData, sets = colnames(methylatedRepeatExpressionUpsetData))
dev.off()
methylatedReducedUpsetData <- data.frame(row.names = row.names(methylatedRepeatExpressionUpsetData),
                                         DMSO = ifelse(rowSums(methylatedRepeatExpressionUpsetData[,1:3])>=2, 1, 0),
                                         Aza = ifelse(rowSums(methylatedRepeatExpressionUpsetData[,4:6])>=2, 1, 0),
                                         Dec = ifelse(rowSums(methylatedRepeatExpressionUpsetData[,7:9])>=2, 1, 0))
methylatedTreatmentUpset <- upset(methylatedReducedUpsetData, sets = colnames(methylatedReducedUpsetData), order.by = c("degree"))
pdf(file="A_cas/A_cas_repeats_methylatedTreatmentUpset.pdf",
    width = 5, height = 5) # or other device
methylatedTreatmentUpset
dev.off()

# Empty data.frame to store info for comparing treatments
bothTreatmentsDE <- data.table()

# For each of Aza and Dec, call DE repeats and merge with methylation stats from the repeat SummarizedExperiment object:
# treatmentDrug = "Aza"
controlSamples <- c(1, 2, 3)
for (treatmentDrug in c("Dec", "Aza")) {
  if (treatmentDrug == "Dec") {
    treatmentSamples <- c(4, 5, 6)
  }else{
    treatmentSamples <- c(7, 8, 9)
  }
  
  DEseqSummary <- summariseDESEQ(expressionData, controlSamples, treatmentSamples, SAMPLES) %>% 
    # filter(padj < 0.05) %>%
    mutate(Treatment = treatmentDrug,
           Repeat = !transcript_id %in% genes$gene_id,
           Direction = ifelse(log2FoldChange<0, "Down", "Up"))
  bothTreatmentsDE <- rbind(bothTreatmentsDE, DEseqSummary)
}
#####

write.table(bothTreatmentsDE, "A_cas/A_cas_DE_table.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
bothTreatmentsDE <- bothTreatmentsDE[bothTreatmentsDE$Repeat,]
bothTreatmentsDF <- as.data.frame(bothTreatmentsDE[,c(1,3,7,10,12)])
compareTreatments <- pivot_wider(bothTreatmentsDF, 
                                 names_from = Treatment, 
                                 values_from = c(log2FoldChange, padj, Direction),
                                 names_sep = "_")
bothTreatmentsDF$Methylated <- bothTreatmentsDF$transcript_id %in% methylatedRepeatIDs
significantRepeatsFC <- bothTreatmentsDF %>% filter(padj < 0.05) %>% group_by(Treatment, Direction) %>% summarise(count = n(), .groups = 'drop') %>%
  complete(Direction, Treatment, fill = list(count = 0)) %>% ggplot(mapping = aes(x = Treatment, y = count, fill = Direction)) +
  geom_col(position = "dodge") + theme_bw() + ylab("DE TEs")
significantRepeatsFC

allRepeatsFC <- bothTreatmentsDF %>% filter(!is.na(Direction)) %>% group_by(Treatment, Direction) %>% summarise(count = n(), .groups = 'drop') %>%
  complete(Direction, Treatment, fill = list(count = 0)) %>% ggplot(mapping = aes(x = Treatment, y = count, fill = Direction)) +
  geom_col(position = "dodge") + theme_bw() + ylab("DE TEs")
allRepeatsFC
ggsave("A_cas/A_cas_allRepeatsFC.pdf", allRepeatsFC,
       height = 5, width = 6, units = "cm")
ggsave("A_cas/A_cas_significantRepeatsFC.pdf", significantRepeatsFC,
       height = 5, width = 6, units = "cm")

ggplot(compareTreatments, mapping = aes(x = log2FoldChange_Aza, y = log2FoldChange_Dec)) +
  geom_point(size = 0.5) + ylim(-8,8) + xlim(-8,8) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw()

ggplot(compareTreatments, mapping = aes(x = log2FoldChange_Aza, y = log2FoldChange_Dec)) +
  geom_point(size = 0.5) + ylim(-8,8) + xlim(-8,8) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw() + stat_cor(method="pearson")

compareTreatmentsWithMeth <- merge(compareTreatments, repeatMethDF, by.x = "transcript_id", by.y = "repeat_id")
compareTreatmentsWithMeth$Methylated <- compareTreatmentsWithMeth$mCG >= cutOff
compareTreatmentsWithMeth$Young <- compareTreatmentsWithMeth$kimura <= 15
ggplot(compareTreatmentsWithMeth, mapping = aes(x = log2FoldChange_Aza, y = log2FoldChange_Dec, colour = mCG)) +
  geom_point() + ylim(-3.5,3.5) + xlim(-3.5,3.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw()

repeatDecAzaMethylation <- ggplot() +
  geom_point(data = compareTreatmentsWithMeth[!compareTreatmentsWithMeth$Methylated,], mapping = aes(x = log2FoldChange_Aza,
                                                                                                     y = log2FoldChange_Dec),
             color = "darkblue", size = 0.5) + 
  geom_point(data = compareTreatmentsWithMeth[compareTreatmentsWithMeth$Methylated,], mapping = aes(x = log2FoldChange_Aza,
                                                                                                    y = log2FoldChange_Dec),
             color = "red", size = 0.5) +
  geom_hline(yintercept = mean(compareTreatmentsWithMeth[compareTreatmentsWithMeth$Methylated,]$log2FoldChange_Dec, na.rm = TRUE), color = "red", linetype = "dashed") + 
  geom_vline(xintercept = mean(compareTreatmentsWithMeth[compareTreatmentsWithMeth$Methylated,]$log2FoldChange_Aza, na.rm = TRUE), color = "red", linetype = "dashed") + 
  geom_hline(yintercept = mean(compareTreatmentsWithMeth[!compareTreatmentsWithMeth$Methylated,]$log2FoldChange_Dec, na.rm = TRUE), color = "darkblue", linetype = "dashed") + 
  geom_vline(xintercept = mean(compareTreatmentsWithMeth[!compareTreatmentsWithMeth$Methylated,]$log2FoldChange_Aza, na.rm = TRUE), color = "darkblue", linetype = "dashed") + 
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw()
repeatDecAzaMethylation
ggsave("A_cas/A_cas_repeats_Aza_vs_Dec_methylation.pdf",
       repeatDecAzaMethylation,
       height = 9, width = 9, units = "cm")

repeatDecAzaKimura <- ggplot() +
  geom_point(data = compareTreatmentsWithMeth[!compareTreatmentsWithMeth$Young,], mapping = aes(x = log2FoldChange_Aza,
                                                                                                y = log2FoldChange_Dec),
             color = "darkblue", size = 0.5) + 
  geom_point(data = compareTreatmentsWithMeth[compareTreatmentsWithMeth$Young,], mapping = aes(x = log2FoldChange_Aza,
                                                                                               y = log2FoldChange_Dec),
             color = "red", size = 0.5) +
  geom_hline(yintercept = mean(compareTreatmentsWithMeth[compareTreatmentsWithMeth$Young,]$log2FoldChange_Dec, na.rm = TRUE), color = "red", linetype = "dashed") + 
  geom_vline(xintercept = mean(compareTreatmentsWithMeth[compareTreatmentsWithMeth$Young,]$log2FoldChange_Aza, na.rm = TRUE), color = "red", linetype = "dashed") + 
  geom_hline(yintercept = mean(compareTreatmentsWithMeth[!compareTreatmentsWithMeth$Young,]$log2FoldChange_Dec, na.rm = TRUE), color = "darkblue", linetype = "dashed") + 
  geom_vline(xintercept = mean(compareTreatmentsWithMeth[!compareTreatmentsWithMeth$Young,]$log2FoldChange_Aza, na.rm = TRUE), color = "darkblue", linetype = "dashed") + 
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw()
repeatDecAzaKimura
ggsave("A_cas/A_cas_repeats_Aza_vs_Dec_kimura.pdf",
       repeatDecAzaKimura,
       height = 9, width = 9, units = "cm")

ggplot(compareTreatmentsWithMeth, mapping = aes(x = log2FoldChange_Aza, y = log2FoldChange_Dec, colour = kimura)) +
  geom_point() + ylim(-3.5,3.5) + xlim(-3.5,3.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw()

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
