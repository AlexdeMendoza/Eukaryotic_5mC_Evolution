source("loadFunctions.R")

CGMAP_PATH = "N_gruberi/EMSeq_N_gru_merged.CGmap.gz"
ORGANISM_FASTA_PATH = "N_gruberi/Ngru_medaka_RagTag.fasta"
width_of_context = 5
FASTA_INDEX_PATH = "N_gruberi/Ngru_medaka_RagTag_pUC19_lambda.fasta.fai"

##### Load methylation data: #####

CG_LOCI <- build_cg_loci_from_fasta_mod(ORGANISM_FASTA_PATH, width_of_context, positiveCntrlPath = "pUC19.fasta", negativeCntrlPath = "lambda.fasta")
SLENGTHS <- index2slengths(FASTA_INDEX_PATH)

methylationData <- loadMethylationData()

CGMAP_SUMMARY <- methylationData[[1]]
CG_METH_ARRAY <- methylationData[[2]]

CG_BS <- build_bsseq(SLENGTHS, CGMAP_SUMMARY, CG_METH_ARRAY)
N_gru_BS <- CG_BS[!seqnames(CG_BS) %in% c("M77789.2", "NC_001416.1")]
lambda_BS <- CG_BS[seqnames(CG_BS) == "NC_001416.1"]
globalLevel <- sum(assays(N_gru_BS)$M, na.rm = TRUE) / sum(assays(N_gru_BS)$Cov, na.rm = TRUE)
falsePositiveRate <- sum(assays(lambda_BS)$M, na.rm = TRUE) / sum(assays(lambda_BS)$Cov, na.rm = TRUE)

CG_BS_CG <- N_gru_BS[rowRanges(N_gru_BS)$dinucleotide == "CG",]
singleBaseMethDistribution <- as.data.frame(assays(CG_BS_CG[assays(CG_BS_CG)$Cov >= 10 & assays(CG_BS_CG)$Fraction >= 1/24,])$Fraction)
singleBaseMethDistribution <- ggplot(singleBaseMethDistribution, mapping = aes(x = V1)) + 
  geom_histogram(aes(y = after_stat(density)), color = "black", fill ="grey", breaks = seq(0, 1, by = 1/24)) + 
  theme_classic() +
  xlab("mCG / CG") + ylab("Density")
singleBaseMethDistribution
ggsave("N_gruberi/N_gru_cytosine_meth_ratios.pdf", singleBaseMethDistribution, width = 3, height = 3, units = "in")

singleBaseMethDistribution <- as.data.frame(assays(N_gru_BS[assays(N_gru_BS)$Cov >= 10 & assays(N_gru_BS)$Fraction >= 1/24,])$Fraction)
singleBaseMethDistribution <- ggplot(singleBaseMethDistribution, mapping = aes(x = V1)) + 
  geom_histogram(aes(y = after_stat(density)), color = "black", fill ="grey", breaks = seq(0, 1, by = 1/24)) + 
  theme_classic() +
  xlab("mC / C") + ylab("Density")
singleBaseMethDistribution
ggsave("N_gruberi/N_gru_cytosine_meth_ratios_CN.pdf", singleBaseMethDistribution, width = 3, height = 3, units = "in")

##### Run kmer analysis #####

CG_kmer_tree <- generate_kmer_tree("c", width_of_context)

CG_whole_genome_results <- extract_kmer_stats(kmer_tree=CG_kmer_tree, BS=CG_BS, width_of_context)
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

##### Gene analysis: #####

#Get gene GenomicRanges object from GFF3 file:
GENE_GTF_PATH <- "N_gruberi/NgruPASA.agatPlusLiftoff.gtf"
txdb <- makeTxDbFromGFF(GENE_GTF_PATH, format="gtf")
transcripts <- transcripts(txdb, columns = c("gene_id", "tx_name"))

mcStats <- get_mc_stats_from_granges(BS = CG_BS, GR = transcripts)

#Get methylation levels per gene:
geneMethData <- data.table(transcript_id = transcripts$tx_name, 
                           gene_id = unlist(transcripts$gene_id), 
                           mC = as.vector(mcStats[[4]]),
                           met = as.vector(mcStats[[2]]),
                           cov = as.vector(mcStats[[3]]),
                           C_sites = as.vector(mcStats[[1]]))

write.table(geneMethData, file = "N_gruberi/N_gruberi_gene_meth.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

#Load in taxonomy data for genes:
taxDataPath <- "N_gruberi/N_gru_tax.tsv"
taxData <- fread(taxDataPath, select = c(1, 11, 18, 19, 20, 21, 22, 23, 24))
taxData <- taxData[taxData$V23 != "Naegleria"]
taxData <- taxData[taxData$V18 %in% c("Eukaryota", "Archaea", "Bacteria", "Viruses"),]
taxData$V18 <- ifelse(taxData$V18 %in% c("Archaea", "Bacteria"), "Prokaryota", taxData$V18)
gene_taxonomies <- data.table(gene_id = NULL, domains = NULL)
for (gene in unique(taxData$V1)) {
  subsetTaxData <- taxData[taxData$V1 == gene]
  domain <- paste(sort(unique(subsetTaxData$V18), decreasing = TRUE), collapse = ", ")
  newRow <- data.frame(gene_id = gene, domains = domain)
  gene_taxonomies <- rbind(gene_taxonomies, newRow)
}

write.table(gene_taxonomies, file = "N_gruberi/N_gruberi_gene_taxonomic_categories.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
categorySummary <- gene_taxonomies %>% group_by(domains) %>% summarise(N_gru = n())
write.table(categorySummary, file = "N_gruberi/N_gruberi_categorySummary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#Add expression data to table:
expressionData <- fread("N_gruberi/Ngru_EBI.genes.results", select = c(1,6))
geneMethData <- merge(geneMethData, expressionData, by = "gene_id", all = TRUE)

#See distribution of methylation across genes:
cutOff <- 0.005
ggplot(geneMethData, mapping = aes(x = mC*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + geom_vline(xintercept = cutOff*100, color = "red", linetype = "dashed") + scale_y_log10()

#Describe genes as methylated or unmethylated using a cut-off, e.g. 0.5%
geneMethData$methylationStatus <- ifelse(geneMethData$mC >= cutOff, "Methylated", "Unmethylated")
geneMethData$methylationStatus <- factor(geneMethData$methylationStatus, levels = c("Unmethylated", "Methylated"))

#Visualise relationship between mCG and TPM:
ggplot(data = geneMethData[geneMethData$TPM != 0,], mapping = aes(x = mC, y = log2(TPM))) +
  geom_point(alpha = 0.1) + theme_bw() + xlab("mCG (%)") + ylab("log2(TPM)")

expressedGenes <- geneMethData[geneMethData$TPM != 0,]
expressedGenes$decile <- ntile(expressedGenes$TPM, 10)
nonExpressedGenes <- geneMethData[geneMethData$TPM == 0,]
nonExpressedGenes$decile <- "0"
geneMethData <- rbind(expressedGenes, nonExpressedGenes)
geneMethData$decile <- as.character(geneMethData$decile)
geneMethData$decile <- factor(geneMethData$decile, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

#Export gene beds, split by decile:
for (i in c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) {
  decileTable <- geneMethData[geneMethData$decile == i,]
  decile_granges <- transcripts[unlist(transcripts$gene_id) %in% decileTable$gene_id]
  export.bed(decile_granges, paste0("N_gruberi/N_gru_", i, "_decile_genes.bed"))
}

mCG_tpm_decile <- ggplot(geneMethData, mapping = aes(x = decile, y = mC * 100, group = decile)) +
  geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim=c(0, 0.3)) + theme_bw() + xlab("TPM decile") + ylab("mCG (%)")
mCG_tpm_decile
ggsave("N_gruberi/N_gru_mCG_tpm_decile.pdf", mCG_tpm_decile, width = 3, height = 3, units = "in")

#Add LGT description to methylation data.
gene_taxonomies <- merge(gene_taxonomies, geneMethData, by.x = "gene_id", by.y = "transcript_id")
taxData <- merge(taxData, geneMethData, by.x = "V1", by.y = "transcript_id")
categorySummaryMethylation <- gene_taxonomies %>% group_by(domains, methylationStatus) %>% summarise(N_gru = n())
write.table(categorySummaryMethylation, file = "N_gruberi/N_gru_categorySummaryMethylation.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#Check distribution of methylation over genes for which we could assign domain of origin
ggplot(gene_taxonomies, mapping = aes(x = mC*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + geom_vline(xintercept = cutOff*100, color = "red", linetype = "dashed")

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


LGT_mC_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = mC*100, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("mC (%)") + xlab("Domain of origin") + coord_cartesian(ylim=c(0, 0.3)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
LGT_mC_bars
ggsave("N_gruberi/N_gruberi_LGT_bars.pdf", LGT_mC_bars, height = 10, width = 13, units = "cm")

LGT_TPM_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = TPM, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("TPM") + xlab("Domain of origin") + coord_cartesian(ylim=c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
LGT_TPM_bars
ggsave("N_gruberi/N_gru_TPM_bars.pdf", LGT_TPM_bars, height = 10, width = 13, units = "cm")

ggplot(gene_taxonomies, mapping = aes(x = mC)) +
  geom_histogram(bins = 50, fill = "lightgrey", color = "black") +
  theme_bw()

ggplot(gene_taxonomies, mapping = aes(x = mC, fill = domains)) +
  geom_histogram(bins = 50, color = "black") +
  theme_bw()

count_data <- totalGeneTaxonomyFrequency %>%
  group_by(methylationStatus) %>%
  summarise(count = sum(frequency))

LGT <- ggplot() +
  geom_col(data = totalGeneTaxonomyFrequency, mapping = aes(x = methylationStatus, y = frequency_proportion*100, fill = domains)) + 
  theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("Proportion (%)") + xlab(paste0("Methylation status (>", cutOff*100, "% mC)")) +
  geom_text(data = count_data, aes(x = methylationStatus, y = 100, label = count), vjust = -0.5)
LGT
ggsave("N_gruberi/N_gruberi_LGT.pdf", LGT, height = 15.14, width = 13.69, units = "cm")


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
write.table(data_grouped, file = "N_gruberi/N_gru_methylated_category_enrichment.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

##### Repeat analysis #####
repeatGR <- loadRepeatAnnotation(REPEAT_GTF_PATH = "N_gruberi/Ngru_medaka_RagTag.fasta.out.noSatellite.gtf.kimura")
repeatFamilyContributions <- data.frame()
for (family in unique(repeatGR$family)) {
  familyGranges <- repeatGR[repeatGR$family == family]
  export.bed(familyGranges, paste0("N_gruberi/N_gru_", family, "_repeats.bed"))
  repeatFamilyContributions <- rbind(repeatFamilyContributions, data.frame(family = family, bases = sum(width(repeatGR[repeatGR$family == family]))))
}
genomeSize <- sum(SLENGTHS[1:(length(SLENGTHS)-2)])

repeatFamilyContributions$contribution <- 100*repeatFamilyContributions$bases/genomeSize
repeatFamilyContributions <- repeatFamilyContributions[repeatFamilyContributions$family %in% c("LINE", "DNA", "LTR", "Unknown"),]
repeatFamilyContributionsPlot <- ggplot(repeatFamilyContributions, mapping = aes(x = family, y = contribution)) +
  geom_col(fill = "black") + ylab("Genome %") + xlab("Repeat sub-class") + theme_bw()
repeatFamilyContributionsPlot
ggsave("N_gruberi/N_gru_repeatFamilyContribution.pdf", plot = repeatFamilyContributionsPlot, width = 6, height = 6, units = "cm")

repeatMeth <- get_mc_stats_from_granges(CG_BS, repeatGR)
repeatMethDF <- data.frame(repeat_id = repeatGR$gene_id,
                           family = repeatGR$family,
                           kimura = repeatGR$kimura,
                           mC = as.vector(repeatMeth[[4]]),
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
ggsave(filename = "N_gruberi/N_gru_kimura_histogram_split_by_family.pdf", plot = kimura_histogram_split_by_family, width = 5, height = 2, units = "in")

repeatSubclassMethylation <- repeatMethDF_filter %>% group_by(family) %>% 
  summarise(mCG = 100*(sum(M, na.rm=TRUE)/sum(Cov, na.rm=TRUE)),
            n = n()) %>% 
  ggplot(mapping = aes(x = family, y = mCG)) +
  geom_col(fill = "black") + theme_bw() + xlab("Repeat sub-class") + 
  ylab("mCG (%)") + 
  geom_hline(yintercept = globalLevel*100, color = "red", linetype = "dashed") +
  geom_text(aes(x = family, y = mCG + 0.1, label = n))
repeatSubclassMethylation
ggsave(filename = "N_gruberi/N_gru_repeatSubclassMethylation.pdf", plot = repeatSubclassMethylation, width = 2, height = 2, units = "in")

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
  geom_text(aes(x = kimura_range, y = mCG + 0.7, label = n))
repeatSubclassKimuramCG
ggsave(filename = "N_gruberi/N_gru_repeatSubclassKimuramCG.pdf", plot = repeatSubclassKimuramCG, width = 8, height = 4, units = "in")

kimura_range_vs_methylation_split_by_family <- ggplot(data=repeatMethDF_filter, mapping=aes(x=kimura_range, y=mC)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_violin(alpha=0) +
  facet_wrap(. ~ family) +
  theme_bw() +
  xlab("Kimura score") +
  ylab("Methylation (%)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
kimura_range_vs_methylation_split_by_family
ggsave(filename = "N_gruberi/N_gru_repeatSubclassKimura_box.pdf", plot = kimura_range_vs_methylation_split_by_family, width = 8, height = 4, units = "in")

kimura_range_vs_methylation <- ggplot(data=repeatMethDF_filter, mapping=aes(x=kimura_range, y=mC*100)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  xlab("Kimura score") +
  ylab("Methylation (%)") +
  # coord_cartesian(ylim=c(0, 4)) + 
  geom_hline(yintercept = 100*globalLevel, color = "black", linetype = "dashed") + 
  geom_hline(yintercept = 100*falsePositiveRate, color = "red", linetype = "dashed")
kimura_range_vs_methylation
ggsave(filename = "N_gruberi/N_gru_kimura_range_vs_methylation.pdf", plot = kimura_range_vs_methylation, width = 6, height = 6, units = "cm")

repeatKimuramC <- repeatMethDF_filter %>% group_by(kimura_range) %>% 
  summarise(mC = 100*(sum(M, na.rm=TRUE)/sum(Cov, na.rm=TRUE)),
            n= n()) %>% 
  ggplot(mapping = aes(x = kimura_range, y = mC)) +
  geom_col() + theme_bw() + xlab("Kimura score") + 
  ylab("mC (%)") + 
  geom_hline(yintercept = globalLevel*100, color = "black", linetype = "dashed") + 
  geom_hline(yintercept = falsePositiveRate*100, color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(aes(x = kimura_range, y = mC + 0.2, label = n))
repeatKimuramC
ggsave(filename = "N_gruberi/N_gru_repeatKimuramC.pdf", plot = repeatKimuramC, width = 6, height = 6, units = "cm")

##### This section is to filter the repeat gtf for mapping expression data: ######

repeat_gtf <- fread("N_gruberi/Ngru_medaka_RagTag.fasta.out.noSatellite.gtf.kimura")

repeat_gr <- GRanges(seqnames = repeat_gtf$V1, ranges=IRanges(start=repeat_gtf$V4, end = repeat_gtf$V5), strand = repeat_gtf$V7, id = repeat_gtf$V9)

# 1) removing the repeats that are completely within genes
## From which genes do I want to eliminate repeats? The ones I'm SURE are genes. Therefore I'm only going to use genes that are not completely within repeats 
flattened_repeats <- GenomicRanges::reduce(repeat_gr, ignore.strand=TRUE)
genes_gtf <- fread("N_gruberi/NgruPASA.agatPlusLiftoff.gtf")
genes_gr <- GRanges(seqnames = genes_gtf$V1, ranges=IRanges(start=genes_gtf$V4, end = genes_gtf$V5), strand = genes_gtf$V7, id = genes_gtf$V9)
genes_not_entirely_covered_by_repeats <- subsetByOverlaps(genes_gr, flattened_repeats, type=c("within"), invert=TRUE, ignore.strand=TRUE)

##Now I want to exclude repeats that are totally covered by this subset of gene:
repeat_gr <- subsetByOverlaps(repeat_gr, genes_not_entirely_covered_by_repeats, type=c("within"), invert=TRUE, ignore.strand=TRUE)

## Finally, exclude those that are < 200bp:
repeat_gr <- repeat_gr[width(repeat_gr) >= 200]
repeat_gtf <- repeat_gtf[repeat_gtf$V9 %in% repeat_gr$id]
repeat_gtf <- repeat_gtf[order(repeat_gtf[,c('V1','V4','V5','V7','V3','V8')]),]
repeat_gtf <- repeat_gtf[!duplicated(repeat_gtf[,c('V1','V4','V5','V3')]),]
write.table(repeat_gtf, file = "N_gruberi/Ngru_medaka_RagTag.fasta.out.noSatellite.gtf.kimura.morethan200.notingenes", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


##### DE analysis of genes: #####

pfam_numbers <- fread("N_gruberi/NgruPASA.agatPlusLiftoff.pep.PfamA.domtblout.4col", sep = "\t", fill = TRUE, col.names = c("Pfam","PFID","tx_name","evalue")) %>% .[.$evalue < 0.001,]

txByGene <- transcriptsBy(txdb, by = "gene")
longestTranscripts <- rbindlist(lapply(names(txByGene), function(gene) {
  transcripts <- txByGene[[gene]]
  return(data.table(tx_name = transcripts[width(transcripts)==max(width(transcripts))]$tx_name[1],
                    gene_id = gene))
}))
pfam_numbers  <- merge(pfam_numbers, longestTranscripts, by = "tx_name")

geneID2GO <- readMappings(file = "N_gruberi/naegleria.eggnog.GOs", sep = " ", IDsep = ",")
geneID2GO <- geneID2GO[names(geneID2GO) %in% longestTranscripts$tx_name]
names(geneID2GO) <- longestTranscripts$gene_id[match(names(geneID2GO), longestTranscripts$tx_name)]
geneUniverse <- names(geneID2GO)

#Load expression data, and exclude genes: #
gene_gr <- GenomicFeatures::genes(txdb)
mC_stats <- get_mc_stats_from_granges(CG_BS, gene_gr)
genes_mC <- data.table(gene_id = gene_gr$gene_id, 
                       mC = as.vector(mC_stats[[4]]),
                       cov = as.vector(mC_stats[[3]]))

genes_mC <- genes_mC[!is.na(genes_mC$mC),]

EXPRESSION_DATA_PATH = c("N_gruberi/DMSO1.bam.cntTable.cleaned",
                         "N_gruberi/DMSO2.bam.cntTable.cleaned",
                         "N_gruberi/DMSO3.bam.cntTable.cleaned",
                         "N_gruberi/Dec1.bam.cntTable.cleaned",
                         "N_gruberi/Dec2.bam.cntTable.cleaned",
                         "N_gruberi/Dec3.bam.cntTable.cleaned",
                         "N_gruberi/Aza1.bam.cntTable.cleaned",
                         "N_gruberi/Aza2.bam.cntTable.cleaned",
                         "N_gruberi/Aza3.bam.cntTable.cleaned")
SAMPLES <- c("DMSO1", "DMSO2", "DMSO3", "Dec1", "Dec2", "Dec3", "Aza1", "Aza2", "Aza3")
expressionData <- load_expression_data(EXPRESSION_DATA_PATH)
geneExpressionData <- expressionData[expressionData[,ncol(expressionData)] %in% gene_gr$gene_id,]

# Run DEseq and make graphs of fold changes etc #

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
           transcript_id %in% gene_gr$gene_id) %>% arrange(padj)
  write.table(geneDEseqSummary, paste0("N_gruberi/N_gru_genes_DEseq_", treatmentDrug, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
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
         filename = paste0("N_gruberi/N_gru_PFAMs_upregulated_", treatmentDrug, ".pdf"),
         height = 12, width = 7, units = "cm")
  
  upregulated_GOs <- get_GOs(geneNames = upregulated_ids, GOlist = geneID2GO)
  # View(upregulated_GOs[upregulated_GOs$classicFisher < 0.05 & upregulated_GOs$Ontology == "BP",])
  # View(upregulated_GOs[upregulated_GOs$classicFisher < 0.05 & upregulated_GOs$Ontology == "MF",])
  upregulated_GOs <- upregulated_GOs %>%
    group_by(Ontology) %>%
    slice_min(order_by = classicFisher, n = 30) %>%
    ungroup()
  gg_upregulated_GOs <- plot_GOs(df = upregulated_GOs, name = "Upregulated genes")
  ggsave(paste0("N_gruberi/N_gru_eggnogUpregulatedGO_", treatmentDrug, ".pdf"),
         plot = gg_upregulated_GOs, height = 12, width = 15)
  
  #Look at relationship between methylation and change in expression:
  geneDEseqSummary <- merge(geneDEseqSummary, genes_mC, by.x = "transcript_id", by.y = "gene_id")
  bothTreatmentsDE <- rbind(bothTreatmentsDE, geneDEseqSummary)
  
  #mCG vs FC, colored by significance:
  mCG_expression_change <- ggplot(geneDEseqSummary, mapping = aes(x = mC*100, y = log2FoldChange, color = padj < 0.05)) +
    geom_point(size = 1) + xlab("Untreated mC (%)") + ylab("log2 Fold Change") + theme_bw() + geom_hline(yintercept = 0)
  mCG_expression_change
  ggsave(mCG_expression_change, 
         filename = paste0("N_gruberi/N_gru_scatter_", treatmentDrug, ".pdf"),
         height = 8, width = 11, units = "cm")
  
  #Proportion of upregulated/downregulated, split by methylation status:
  geneDEseqSummaryBarData <- geneDEseqSummary %>%
    filter(padj < 0.05) %>%
    mutate(Direction = ifelse(log2FoldChange < 0, "Down", ifelse(log2FoldChange > 0, "Up", "Neutral")),
           ThresholdmC = ifelse(mC < cutOff, "Unmethylated", "Methylated"),
           Direction = factor(Direction, levels = c("Down", "Up")),
           ThresholdmC = factor(ThresholdmC, levels = c("Unmethylated", "Methylated"))) %>%
    group_by(Direction, ThresholdmC) %>%
    summarise(count = n(), .groups = 'drop') %>%
    complete(Direction, ThresholdmC, fill = list(count = 0))
  geneDEseqSummaryBar <- ggplot(geneDEseqSummaryBarData, mapping = aes(x = ThresholdmC, fill = Direction, y = count)) +
    geom_col(position = "dodge") + xlab("Untreated mC (%)") + ylab("log2 Fold Change") + theme_bw()
  geneDEseqSummaryBar
  ggsave(geneDEseqSummaryBar, 
         filename = paste0("N_gruberi/DE_genes_UpDown_bar_", treatmentDrug, ".pdf"),
         height = 8, width = 8, units = "cm")
}

# Compare aza vs decitabine

bothTreatmentsDF <- as.data.frame(bothTreatmentsDE[,c(1,3,7,10)])
compareTreatments <- pivot_wider(bothTreatmentsDF, 
                                 names_from = Treatment, 
                                 values_from = c(log2FoldChange, padj),
                                 names_sep = "_")
compareTreatments <- merge(compareTreatments, genes_mC, by.x = "transcript_id", by.y = "gene_id")

#By looking at genes in IGV, 5% (0.05) seems to be a good threshold for "real" methylated genes
compareTreatments$Methylated <- compareTreatments$mC >= cutOff

compareTreatments <- compareTreatments[!is.na(compareTreatments$log2FoldChange_Dec),]
compareTreatments <- compareTreatments[!is.na(compareTreatments$log2FoldChange_Aza),]
compareTreatments$alpha = ifelse(pmin(compareTreatments$padj_Dec, compareTreatments$padj_Aza)<0.05, 1, 0.01)

#Compare Dec vs Aza change, highlighting methylated genes:
geneDecAzaMethylation <- ggplot() +
  geom_point(data = compareTreatments[!compareTreatments$Methylated,], mapping = aes(x = log2FoldChange_Aza,
                                                                                     y = log2FoldChange_Dec,
                                                                                     alpha = alpha),
             color = "darkblue", size = 0.5, show.legend = FALSE) + 
  geom_point(data = compareTreatments[compareTreatments$Methylated,], mapping = aes(x = log2FoldChange_Aza,
                                                                                    y = log2FoldChange_Dec,
                                                                                    alpha = alpha),
             color = "red", size = 0.5, show.legend = FALSE) +
  # geom_hline(yintercept = mean(compareTreatments[compareTreatments$Methylated & compareTreatments$alpha == 1,]$log2FoldChange_Dec, na.rm = TRUE), color = "red", linetype = "dashed") + 
  # geom_vline(xintercept = mean(compareTreatments[compareTreatments$Methylated & compareTreatments$alpha == 1,]$log2FoldChange_Aza, na.rm = TRUE), color = "red", linetype = "dashed") + 
  # geom_hline(yintercept = mean(compareTreatments[!compareTreatments$Methylated & compareTreatments$alpha == 1,]$log2FoldChange_Dec, na.rm = TRUE), color = "darkblue", linetype = "dashed") + 
  # geom_vline(xintercept = mean(compareTreatments[!compareTreatments$Methylated & compareTreatments$alpha == 1,]$log2FoldChange_Aza, na.rm = TRUE), color = "darkblue", linetype = "dashed") + 
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw()
geneDecAzaMethylation
ggsave("N_gruberi/Aza_vs_Dec_methylation_FC_eitherTrt.pdf",
       geneDecAzaMethylation,
       height = 9, width = 9, units = "cm")


##### PFAM / GO enrichment of methylated / viral genes #####
#Get methylated geneIDs:
methylatedGenes <- genes_mC[genes_mC$mC > cutOff,]
methylatedGeneIDs <- methylatedGenes$gene_id

#Get viral gene IDs:
taxData <- fread("N_gruberi/N_gru_tax.tsv", select = c(1, 18, 23))
tx_gene <- mcols(transcripts(txdb, columns = c("gene_id", "tx_name")))
taxData <- merge(taxData, tx_gene, by.y = "tx_name", by.x = "V1")
viralTaxData <- taxData[taxData$V18=="Viruses",]
viralGeneIDs <- unique(unlist(viralTaxData$gene_id))

methylated_pfams <- pfamEnrichments(ids = methylatedGeneIDs) %>% arrange(fishp)
# View(methylated_pfams[methylated_pfams$fishp < 0.05,])
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
       filename = paste0("N_gruberi/PFAMs_methylated_genes.pdf"),
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
ggsave(paste0("N_gruberi/eggnog_methylated_genes_GO.pdf"),
       plot = gg_methylated_GOs, height = 12, width = 15)

## Do for viral genes:
viral_pfams <- pfamEnrichments(ids = viralGeneIDs) %>% arrange(fishp)
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
       filename = paste0("N_gruberi/N_gru_PFAMs_viral_genes.pdf"),
       height = 12, width = 7, units = "cm")

viral_GOs <- get_GOs(geneNames = viralGeneIDs, GOlist = geneID2GO)
# View(viral_GOs[viral_GOs$classicFisher < 0.05 & viral_GOs$Ontology == "BP",])
# View(viral_GOs[viral_GOs$classicFisher < 0.05 & viral_GOs$Ontology == "MF",])
viral_GOs <- viral_GOs %>%
  group_by(Ontology) %>%
  slice_min(order_by = classicFisher, n = 30) %>%
  ungroup()
gg_viral_GOs <- plot_GOs(df = viral_GOs, name = "viral genes")
gg_viral_GOs
ggsave(paste0("N_gruberi/N_gru_eggnog_viral_genes_GO.pdf"),
       plot = gg_viral_GOs, height = 12, width = 15)

##### DE analysis of repeats: #####

methylatedRepeatIDs <- repeatMethDF[repeatMethDF$mC >= cutOff,]$repeat_id

#Load expression data, and exclude genes:

repeatExpressionData <- expressionData[!expressionData[,ncol(expressionData)] %in% gene_gr$gene_id,]
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
dev.off()
reducedUpsetData <- data.frame(row.names = row.names(repeatExpressionUpsetData),
                               DMSO = ifelse(rowSums(repeatExpressionUpsetData[,1:3])>=3, 1, 0),
                               Aza = ifelse(rowSums(repeatExpressionUpsetData[,4:6])>=3, 1, 0),
                               Dec = ifelse(rowSums(repeatExpressionUpsetData[,7:9])>=3, 1, 0))
treatmentUpset <- upset(reducedUpsetData, sets = colnames(reducedUpsetData), order.by = c("degree"))
pdf(file="N_gruberi/N_gru_repeats_treatmentUpset.pdf",
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
pdf(file="N_gruberi/N_gru_repeats_methylatedTreatmentUpset.pdf",
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
           Repeat = !transcript_id %in% gene_gr$gene_id,
           Direction = ifelse(log2FoldChange<0, "Down", "Up"))
  bothTreatmentsDE <- rbind(bothTreatmentsDE, DEseqSummary)
}
#####

write.table(bothTreatmentsDE, "N_gruberi/N_gru_DE_table.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
bothTreatmentsDE <- bothTreatmentsDE[bothTreatmentsDE$Repeat,]
bothTreatmentsDF <- as.data.frame(bothTreatmentsDE[,c(1,3,7,10,12)])
compareTreatments <- pivot_wider(bothTreatmentsDF, 
                                 names_from = Treatment, 
                                 values_from = c(log2FoldChange, padj, Direction),
                                 names_sep = "_")

significantRepeatsFC <- bothTreatmentsDF %>% filter(padj < 0.05) %>% group_by(Treatment, Direction) %>% summarise(count = n(), .groups = 'drop') %>%
  complete(Direction, Treatment, fill = list(count = 0)) %>% ggplot(mapping = aes(x = Treatment, y = count, fill = Direction)) +
  geom_col(position = "dodge") + theme_bw() + ylab("DE TEs")
significantRepeatsFC
allRepeatsFC <- bothTreatmentsDF %>% filter(!is.na(Direction)) %>% group_by(Treatment, Direction) %>% summarise(count = n(), .groups = 'drop') %>%
  complete(Direction, Treatment, fill = list(count = 0)) %>% ggplot(mapping = aes(x = Treatment, y = count, fill = Direction)) +
  geom_col(position = "dodge") + theme_bw() + ylab("TEs")
allRepeatsFC
ggsave("N_gruberi/N_gru_allRepeatsFC.pdf", allRepeatsFC,
       height = 5, width = 6, units = "cm")
ggsave("N_gruberi/N_gru_significantRepeatsFC.pdf", significantRepeatsFC,
       height = 5, width = 6, units = "cm")

ggplot(compareTreatments, mapping = aes(x = log2FoldChange_Aza, y = log2FoldChange_Dec)) +
  geom_point(size = 0.5) + ylim(-8,8) + xlim(-8,8) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw()

ggplot(compareTreatments, mapping = aes(x = log2FoldChange_Aza, y = log2FoldChange_Dec)) +
  geom_point(size = 0.5) + ylim(-8,8) + xlim(-8,8) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlab("log2FC 5-azacytidine") + ylab("log2FC Decitabine") + theme_bw() + stat_cor(method="pearson")

compareTreatmentsWithMeth <- merge(compareTreatments, repeatMethDF, by.x = "transcript_id", by.y = "repeat_id")
compareTreatmentsWithMeth$Methylated <- compareTreatmentsWithMeth$mC >= cutOff
compareTreatmentsWithMeth$Young <- compareTreatmentsWithMeth$kimura <= 15
ggplot(compareTreatmentsWithMeth, mapping = aes(x = log2FoldChange_Aza, y = log2FoldChange_Dec, colour = mC)) +
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
ggsave("N_gruberi/N_gru_repeats_Aza_vs_Dec_methylation.pdf",
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
ggsave("N_gruberi/N_gru_repeats_Aza_vs_Dec_kimura.pdf",
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
