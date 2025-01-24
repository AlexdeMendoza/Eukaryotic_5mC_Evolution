##### Load libraries and functions #####

source("loadFunctions.R")

##### Load methylation data: #####

CGMAP_PATH = "P_patens/SRR18164483_1.fastq.CGmap.gz"
ORGANISM_FASTA_PATH = "P_patens/Physcomitrium_patens_formerly_Physcomitrella_patens.formatted.faa"
width_of_context = 5
FASTA_INDEX_PATH = "P_patens/Physcomitrium_patens_formerly_Physcomitrella_patens.formatted.faa.fai"

CG_LOCI <- build_cg_loci_from_fasta_mod(ORGANISM_FASTA_PATH, width_of_context)
SLENGTHS <- index2slengths(FASTA_INDEX_PATH)

methylationData <- loadMethylationData()

CGMAP_SUMMARY <- methylationData[[1]]
CG_METH_ARRAY <- methylationData[[2]]

CG_BS <- build_bsseq(SLENGTHS, CGMAP_SUMMARY, CG_METH_ARRAY)

CG_BS_CG <- CG_BS[rowRanges(CG_BS)$dinucleotide == "CG",]
CG_BS_CHG <- CG_BS[substr(rowRanges(CG_BS)$context, 4, 4) != "g" & substr(rowRanges(CG_BS)$context, 5, 5) == "g",]
CG_BS_CHH <- CG_BS[substr(rowRanges(CG_BS)$context, 4, 4) != "g" & substr(rowRanges(CG_BS)$context, 5, 5) != "g",]

export_cg_map(CG_BS_CG, "P_patens/P_patens_CG.CGmap")
export_cg_map(CG_BS_CHG, "P_patens/P_patens_CHG.CGmap")
export_cg_map(CG_BS_CHH, "P_patens/P_patens_CHH.CGmap")

##### Run kmer analysis #####

CG_kmer_tree <- generate_kmer_tree("c", width_of_context)

CG_whole_genome_results <- extract_kmer_stats(kmer_tree=CG_kmer_tree, BS=CG_BS, width_of_context) # CG_whole_genome_results <- readRDS(file.path(SCRATCH_DIR, "CG_whole_genome_results.rds"))

CG_whole_genome_results$mCG <- CG_whole_genome_results$mCG*100
CG_whole_genome_results <- CG_whole_genome_results[order(CG_whole_genome_results$mCG, decreasing = TRUE),]
CG_whole_genome_results
CG_whole_genome_results$plus1 <- substr(CG_whole_genome_results$context, 4, 4)
CG_whole_genome_results$plus2 <- substr(CG_whole_genome_results$context, 5, 5)
CG_whole_genome_results$minus1 <- substr(CG_whole_genome_results$context, 2, 2)
CG_whole_genome_results$minus2 <- substr(CG_whole_genome_results$context, 1, 1)

dinucleotide <- CG_whole_genome_results %>% group_by(plus1) %>% summarise(weighted_mCG = sum(mCG * n) / sum(n))

dinucleotide_plot <- ggplot(dinucleotide, mapping = aes(x = plus1, y = weighted_mCG)) +
  geom_col() + 
  # geom_hline(yintercept = falsePositiveRate, color = "red", linetype = "dashed") + 
  theme_bw() + xlab("CX") + ylab("mC (%)")
dinucleotide_plot

CG_global_level <- sum(assays(CG_BS_CG)$M, na.rm = TRUE) / sum(assays(CG_BS_CG)$Cov, na.rm = TRUE)
CHG_global_level <- sum(assays(CG_BS_CHG)$M, na.rm = TRUE) / sum(assays(CG_BS_CHG)$Cov, na.rm = TRUE)
CHH_global_level <- sum(assays(CG_BS_CHH)$M, na.rm = TRUE) / sum(assays(CG_BS_CHH)$Cov, na.rm = TRUE)

contextLevels <- data.frame(context = c("CG", "CHG", "CHH"),
                            mC = c(CG_global_level, CHG_global_level, CHH_global_level))

contextPlot <- ggplot(contextLevels, mapping = aes(x = context, y = mC*100)) +
  geom_col() + xlab("Cytosine context") + ylab ("mC (%)") + theme_bw()
contextPlot
ggsave("P_patens/P_patens_context.pdf", contextPlot, units = "in", height = 2.5, width = 2)

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

##### Load in methylation data for genes: #####

#Get gene GenomicRanges object from GFF3 file:
GENE_GTF_PATH <- "P_patens/P_patens_CoGe_V3.agat.gtf"
txdb <- makeTxDbFromGFF(GENE_GTF_PATH, format="gtf")
transcripts <- transcripts(txdb, columns = c("gene_id", "tx_name"))

transcripts$gene_id <- unlist(transcripts$gene_id)
#Keep just the longest transcript per gene:
longestTranscript <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
for (gene in unique(transcripts$gene_id)) {
  geneTranscript <- transcripts[transcripts$gene_id == gene]
  maxLength <- max(width(geneTranscript))
  longestTranscript <- c(longestTranscript, geneTranscript[width(geneTranscript) == maxLength][1])
}


GEVE_prodigal_genes <- import.bed(con = "P_patens/P_pat_GEVE_prodigal.bed")
GEVE_prodigal_genes$tx_name <- GEVE_prodigal_genes$name
longestTranscript <- c(longestTranscript, GEVE_prodigal_genes)
mcStats <- get_mc_stats_from_granges(BS = CG_BS_CG, GR = longestTranscript)
CHG_mcStats <- get_mc_stats_from_granges(BS = CG_BS_CHG, GR = longestTranscript)
CHH_mcStats <- get_mc_stats_from_granges(BS = CG_BS_CHH, GR = longestTranscript)

#Get methylation levels per gene:
geneMethData <- data.table(gene_id = longestTranscript$tx_name, 
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

write.table(geneMethData, file = "P_patens/P_patens_gene_meth.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

#Load in taxonomy data for genes:
taxDataPath <- "P_patens/P_patens_tax.tsv"
taxData <- fread(taxDataPath, select = c(1, 11, 18, 19, 20, 21, 22, 23, 24))
taxData <- taxData[taxData$V23 != "Physcomitrium"]
taxData <- taxData[taxData$V18 %in% c("Eukaryota", "Archaea", "Bacteria", "Viruses"),]
taxData$V18 <- ifelse(taxData$V18 %in% c("Archaea", "Bacteria"), "Prokaryota", taxData$V18)
taxData <- taxData[taxData$V1 %in% longestTranscript$tx_name,]

viralRecall_tax <- data.frame(V1 = GEVE_prodigal_genes$tx_name,
                              V18 = "viralRecall")
taxData <- rbind(taxData, viralRecall_tax, fill = TRUE)

gene_taxonomies <- data.table(gene_id = NULL, domains = NULL)
for (gene in unique(taxData$V1)) {
  subsetTaxData <- taxData[taxData$V1 == gene]
  domain <- paste(sort(unique(subsetTaxData$V18), decreasing = TRUE), collapse = ", ")
  newRow <- data.frame(gene_id = gene, domains = domain)
  gene_taxonomies <- rbind(gene_taxonomies, newRow)
}
write.table(gene_taxonomies, file = "P_patens/P_patens_gene_taxonomic_categories.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
categorySummary <- gene_taxonomies %>% group_by(domains) %>% summarise(P_pat = n())
write.table(categorySummary, file = "P_patens/P_pat_categorySummary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# #Add expression data to table:
# expressionData <- fread("../Cyapar_EBI.genes.results", select = c(1,6))
# expressionData$gene_id <- sub("gene_", "", expressionData$gene_id)
# geneMethData <- merge(geneMethData, expressionData, by = "gene_id")

#See distribution of methylation across genes:
cutOff <- 0.1
ggplot(geneMethData, mapping = aes(x = mCG*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + 
  geom_vline(xintercept = cutOff*100, color = "red", linetype = "dashed")
ggplot(geneMethData, mapping = aes(x = mCHG*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + 
  geom_vline(xintercept = cutOff*100, color = "red", linetype = "dashed")
ggplot(geneMethData, mapping = aes(x = mCHH*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + 
  geom_vline(xintercept = cutOff*100, color = "red", linetype = "dashed")

#Describe genes as methylated or unmethylated using a cut-off, e.g. 0.5%
geneMethData$CG_methylationStatus <- ifelse(geneMethData$mCG >= cutOff, "Methylated", "Unmethylated")
geneMethData$CHG_methylationStatus <- ifelse(geneMethData$mCHG >= cutOff, "Methylated", "Unmethylated")
geneMethData$CHH_methylationStatus <- ifelse(geneMethData$mCHH >= cutOff, "Methylated", "Unmethylated")
geneMethData$CG_methylationStatus <- factor(geneMethData$CG_methylationStatus, levels = c("Unmethylated", "Methylated"))
geneMethData$CHG_methylationStatus <- factor(geneMethData$CHG_methylationStatus, levels = c("Unmethylated", "Methylated"))
geneMethData$CHH_methylationStatus <- factor(geneMethData$CHH_methylationStatus, levels = c("Unmethylated", "Methylated"))

# #Visualise relationship between mCG and TPM:
# ggplot(data = geneMethData[geneMethData$TPM != 0,], mapping = aes(x = mCG, y = log2(TPM))) +
#   geom_point(alpha = 0.1) + theme_bw() + xlab("mCG (%)") + ylab("log2(TPM)")

#Add LGT description to methylation data.
gene_taxonomies <- merge(gene_taxonomies, geneMethData, "gene_id")
taxData <- merge(taxData, geneMethData, by.x = "V1", by.y = "gene_id")

#Check distribution of methylation over genes for which we could assign domain of origin
ggplot(gene_taxonomies, mapping = aes(x = mCG*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + geom_vline(xintercept = cutOff*100, color = "red", linetype = "dashed")

CG_totalGeneTaxonomyFrequency <- gene_taxonomies %>% group_by(CG_methylationStatus, domains) %>% summarise(frequency = n()) %>%
  mutate(frequency_proportion = frequency / sum(frequency)) %>% ungroup() %>% as.data.frame()
CG_totalGeneTaxonomyFrequency <- CG_totalGeneTaxonomyFrequency[!is.na(CG_totalGeneTaxonomyFrequency$CG_methylationStatus),]

CHG_totalGeneTaxonomyFrequency <- gene_taxonomies %>% group_by(CHG_methylationStatus, domains) %>% summarise(frequency = n()) %>%
  mutate(frequency_proportion = frequency / sum(frequency)) %>% ungroup() %>% as.data.frame()
CHG_totalGeneTaxonomyFrequency <- CHG_totalGeneTaxonomyFrequency[!is.na(CHG_totalGeneTaxonomyFrequency$CHG_methylationStatus),]

CHH_totalGeneTaxonomyFrequency <- gene_taxonomies %>% group_by(CHH_methylationStatus, domains) %>% summarise(frequency = n()) %>%
  mutate(frequency_proportion = frequency / sum(frequency)) %>% ungroup() %>% as.data.frame()
CHH_totalGeneTaxonomyFrequency <- CHH_totalGeneTaxonomyFrequency[!is.na(CHH_totalGeneTaxonomyFrequency$CHH_methylationStatus),]

# LGT_mCG_col <- gene_taxonomies %>% group_by(domains) %>% summarise(total_mCG = sum(met)/sum(cov), count = n()) %>% ggplot(mapping = aes(x = domains, y = total_mCG*100, fill = domains)) +
#   geom_col() +
#   geom_hline(yintercept = falsePositiveRate*100, linetype = "dashed", color = "red") +
#   geom_hline(yintercept = globalLevel*100, linetype = "dashed", color = "black") +
#   theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
#   ylab("mCG (%)") + xlab("Domain of origin") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# LGT_mCG_col
# ggsave("figures/acanthamoeba/initial_emseq/LGT/LGT_mCG_col.pdf", LGT_mCG_col, height = 10, width = 13, units = "cm")

gene_taxonomies$domains <- factor(gene_taxonomies$domains, levels = c("Eukaryota", "Prokaryota", "Prokaryota, Eukaryota", "Viruses", "Viruses, Eukaryota", "Viruses, Prokaryota", "viralRecall"))


count_data <- gene_taxonomies %>%
  group_by(domains) %>%
  summarise(count = sum(!is.na(mCG)))

LGT_mCG_bars <- ggplot(gene_taxonomies, aes(x = domains, y = mCG * 100, fill = domains)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#64a858")) +
  ylab("mCG (%)") +
  xlab("Domain of origin") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  
  geom_text(data = count_data, aes(x = domains, y = 100, label = count), vjust = -0.5)
LGT_mCG_bars

count_data <- gene_taxonomies %>%
  group_by(domains) %>%
  summarise(count = sum(!is.na(mCHG)))

LGT_mCHG_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = mCHG*100, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#64a858")) +
  ylab("mCHG (%)") + xlab("Domain of origin") + 
  # coord_cartesian(ylim=c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_text(data = count_data, aes(x = domains, y = 100, label = count), vjust = -0.5)
LGT_mCHG_bars

count_data <- gene_taxonomies %>%
  group_by(domains) %>%
  summarise(count = sum(!is.na(mCHH)))

LGT_mCHH_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = mCHH*100, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#64a858")) +
  ylab("mCHH (%)") + xlab("Domain of origin") + 
  coord_cartesian(ylim=c(0, 80)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_text(data = count_data, aes(x = domains, y = 75, label = count), vjust = -0.5)
LGT_mCHH_bars

ggsave("P_patens/P_patens_LGT_bars_CG.pdf", LGT_mCG_bars, height = 10, width = 13, units = "cm")
ggsave("P_patens/P_patens_LGT_bars_CHG.pdf", LGT_mCHG_bars, height = 10, width = 13, units = "cm")
ggsave("P_patens/P_patens_LGT_bars_CHH.pdf", LGT_mCHH_bars, height = 10, width = 13, units = "cm")

# LGT_TPM_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = TPM, fill = domains)) +
#   geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
#   ylab("TPM") + xlab("Domain of origin") + coord_cartesian(ylim=c(0, 60)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# LGT_TPM_bars
# ggsave("Cpar_TPM_bars.pdf", LGT_mCG_bars, height = 10, width = 13, units = "cm")

ggplot(gene_taxonomies, mapping = aes(x = mCG)) +
  geom_histogram(bins = 50, fill = "lightgrey", color = "black") +
  theme_bw()

ggplot(gene_taxonomies, mapping = aes(x = mCG, fill = domains)) +
  geom_histogram(bins = 50, color = "black") +
  theme_bw()

CG_totalGeneTaxonomyFrequency$domains <- factor(CG_totalGeneTaxonomyFrequency$domains, levels = c("Eukaryota", "Prokaryota", "Prokaryota, Eukaryota", "Viruses", "Viruses, Eukaryota", "Viruses, Prokaryota", "viralRecall"))
CHG_totalGeneTaxonomyFrequency$domains <- factor(CHG_totalGeneTaxonomyFrequency$domains, levels = c("Eukaryota", "Prokaryota", "Prokaryota, Eukaryota", "Viruses", "Viruses, Eukaryota", "Viruses, Prokaryota", "viralRecall"))
CHH_totalGeneTaxonomyFrequency$domains <- factor(CHH_totalGeneTaxonomyFrequency$domains, levels = c("Eukaryota", "Prokaryota", "Prokaryota, Eukaryota", "Viruses", "Viruses, Eukaryota", "Viruses, Prokaryota", "viralRecall"))

count_data <- CG_totalGeneTaxonomyFrequency %>%
  group_by(CG_methylationStatus) %>%
  summarise(count = sum(frequency))

CG_LGT <- ggplot() +
  geom_col(data = CG_totalGeneTaxonomyFrequency, mapping = aes(x = CG_methylationStatus, y = frequency_proportion*100, fill = domains)) + 
  theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#64a858")) +
  ylab("Proportion (%)") + xlab(paste0("Methylation status (>", cutOff*100, "% mCG)")) +
  geom_text(data = count_data, aes(x = CG_methylationStatus, y = 100, label = count), vjust = -0.5)
CG_LGT

count_data <- CHG_totalGeneTaxonomyFrequency %>%
  group_by(CHG_methylationStatus) %>%
  summarise(count = sum(frequency))

CHG_LGT <- ggplot() +
  geom_col(data = CHG_totalGeneTaxonomyFrequency, mapping = aes(x = CHG_methylationStatus, y = frequency_proportion*100, fill = domains)) + 
  theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#64a858")) +
  ylab("Proportion (%)") + xlab(paste0("Methylation status (>", cutOff*100, "% mCHG)")) +
  geom_text(data = count_data, aes(x = CHG_methylationStatus, y = 100, label = count), vjust = -0.5)
CHG_LGT

count_data <- CHH_totalGeneTaxonomyFrequency %>%
  group_by(CHH_methylationStatus) %>%
  summarise(count = sum(frequency))

CHH_LGT <- ggplot() +
  geom_col(data = CHH_totalGeneTaxonomyFrequency, mapping = aes(x = CHH_methylationStatus, y = frequency_proportion*100, fill = domains)) + 
  theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#64a858")) +
  ylab("Proportion (%)") + xlab(paste0("Methylation status (>", cutOff*100, "% mCHH)")) +
  geom_text(data = count_data, aes(x = CHH_methylationStatus, y = 100, label = count), vjust = -0.5)
CHH_LGT

ggsave("P_patens/P_patens_CG_LGT.pdf", CG_LGT, height = 15.14, width = 13.69, units = "cm")
ggsave("P_patens/P_patens_CHG_LGT.pdf", CHG_LGT, height = 15.14, width = 13.69, units = "cm")
ggsave("P_patens/P_patens_CHH_LGT.pdf", CHH_LGT, height = 15.14, width = 13.69, units = "cm")

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



#Calculate enrichment of domains within methylated genes
CG_data_grouped <- CG_totalGeneTaxonomyFrequency %>% 
  group_by(domains) %>% 
  summarise(
    # Total_Frequency = sum(frequency),
    Frequency_Methylated = sum(frequency[CG_methylationStatus == "Methylated"]),
    Frequency_Unmethylated = sum(frequency[CG_methylationStatus == "Unmethylated"]),
    Proportion_Methylated = sum(frequency_proportion[CG_methylationStatus == "Methylated"]),
    Proportion_Unmethylated = sum(frequency_proportion[CG_methylationStatus == "Unmethylated"])
  )
total_Methylated <- sum(CG_data_grouped$Frequency_Methylated)
total_Unmethylated <- sum(CG_data_grouped$Frequency_Unmethylated)

CG_data_grouped <- CG_data_grouped %>%
  rowwise() %>%
  mutate(p_value = fisher.test(matrix(c(Frequency_Methylated,
                                        total_Methylated - Frequency_Methylated,
                                        Frequency_Unmethylated,
                                        total_Unmethylated - Frequency_Unmethylated),
                                      nrow = 2), alternative = "two.sided")$p.value) %>%
  ungroup()
CG_data_grouped
write.table(CG_data_grouped, file = "P_patens/CG_methylated_category_enrichment.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

CHG_data_grouped <- CHG_totalGeneTaxonomyFrequency %>% 
  group_by(domains) %>% 
  summarise(
    # Total_Frequency = sum(frequency),
    Frequency_Methylated = sum(frequency[CHG_methylationStatus == "Methylated"]),
    Frequency_Unmethylated = sum(frequency[CHG_methylationStatus == "Unmethylated"]),
    Proportion_Methylated = sum(frequency_proportion[CHG_methylationStatus == "Methylated"]),
    Proportion_Unmethylated = sum(frequency_proportion[CHG_methylationStatus == "Unmethylated"])
  )
total_Methylated <- sum(CHG_data_grouped$Frequency_Methylated)
total_Unmethylated <- sum(CHG_data_grouped$Frequency_Unmethylated)

CHG_data_grouped <- CHG_data_grouped %>%
  rowwise() %>%
  mutate(p_value = fisher.test(matrix(c(Frequency_Methylated,
                                        total_Methylated - Frequency_Methylated,
                                        Frequency_Unmethylated,
                                        total_Unmethylated - Frequency_Unmethylated),
                                      nrow = 2), alternative = "two.sided")$p.value) %>%
  ungroup()
CHG_data_grouped
write.table(CHG_data_grouped, file = "P_patens/CHG_methylated_category_enrichment.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

data_grouped <- CHH_totalGeneTaxonomyFrequency %>% 
  group_by(domains) %>% 
  summarise(
    # Total_Frequency = sum(frequency),
    Frequency_Methylated = sum(frequency[CHH_methylationStatus == "Methylated"]),
    Frequency_Unmethylated = sum(frequency[CHH_methylationStatus == "Unmethylated"]),
    Proportion_Methylated = sum(frequency_proportion[CHH_methylationStatus == "Methylated"]),
    Proportion_Unmethylated = sum(frequency_proportion[CHH_methylationStatus == "Unmethylated"])
  )
total_Methylated <- sum(data_grouped$Frequency_Methylated)
total_Unmethylated <- sum(data_grouped$Frequency_Unmethylated)

CHH_data_grouped <- data_grouped %>%
  rowwise() %>%
  mutate(p_value = fisher.test(matrix(c(Frequency_Methylated,
                                        total_Methylated - Frequency_Methylated,
                                        Frequency_Unmethylated,
                                        total_Unmethylated - Frequency_Unmethylated),
                                      nrow = 2), alternative = "two.sided")$p.value) %>%
  ungroup()
CHH_data_grouped
write.table(CHH_data_grouped, file = "P_patens/CHH_methylated_category_enrichment.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

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
