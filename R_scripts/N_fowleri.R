source("loadFunctions.R")

CGMAP_PATH = "N_fowleri/Nfow.R9_5mC.methyl.context.bed.CGmap.gz"
ORGANISM_FASTA_PATH = "N_fowleri/Naegleria_fowleri_GCF_008403515.1.gDNA.fasta"
width_of_context = 5
FASTA_INDEX_PATH = "N_fowleri/Naegleria_fowleri_GCF_008403515.1.gDNA.fasta.fai"
taxDataPath <- "N_fowleri/N_fowleri_tax.txt"

##### Load methylation data: #####

CG_LOCI <- build_cg_loci_from_fasta_mod(ORGANISM_FASTA_PATH, width_of_context)
SLENGTHS <- index2slengths(FASTA_INDEX_PATH)

methylationData <- loadMethylationData()

CGMAP_SUMMARY <- methylationData[[1]]
CG_METH_ARRAY <- methylationData[[2]]

CG_BS <- build_bsseq(SLENGTHS, CGMAP_SUMMARY, CG_METH_ARRAY)
CG_BS_CG <- CG_BS[rowRanges(CG_BS)$dinucleotide == "CG",]

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

##### Load in methylation data for genes: #####

#Get gene GenomicRanges object from GFF3 file:
GENE_GTF_PATH <- "N_fowleri/Naegleria_fowleri_GCF_008403515.1.gtf"
txdb <- makeTxDbFromGFF(GENE_GTF_PATH, format="gtf")
transcripts <- transcripts(txdb)
transcripts$tx_name <- sub("rna-", "", transcripts$tx_name)

CG_BS_CG <- CG_BS[rowRanges(CG_BS)$dinucleotide == "CG",]
mcStats <- get_mc_stats_from_granges(BS = CG_BS_CG, GR = transcripts)

#Get methylation levels per gene:
geneMethData <- data.table(gene_id = transcripts$tx_name, mCG = mcStats[[4]])
colnames(geneMethData)[2] <- "mCG" 
write.table(geneMethData, file = "N_fowleri/N_fowleri_gene_meth.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

#Load in taxonomy data for genes:
taxData <- fread(taxDataPath, select = c(1, 11, 18, 19, 20, 21, 22, 23, 24))
taxData <- taxData[taxData$V23 != "Naegleria"]
taxData$V1 <- sub(".*\\|(.*)\\|.*", "\\1", taxData$V1)
taxData <- taxData[taxData$V18 %in% c("Eukaryota", "Archaea", "Bacteria", "Viruses"),]
taxData$V18 <- ifelse(taxData$V18 %in% c("Archaea", "Bacteria"), "Prokaryote", taxData$V18)
gene_taxonomies <- data.table(gene_id = NULL, domains = NULL)
for (gene in unique(taxData$V1)) {
  subsetTaxData <- taxData[taxData$V1 == gene]
  domain <- paste(sort(unique(subsetTaxData$V18), decreasing = TRUE), collapse = ", ")
  newRow <- data.frame(gene_id = gene, domains = domain)
  gene_taxonomies <- rbind(gene_taxonomies, newRow)
}

write.table(gene_taxonomies, file = "N_fowleri/N_fowleri_gene_taxonomic_categories.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# #Add expression data to table:
# expressionData <- fread("../Cyapar_EBI.genes.results", select = c(1,6))
# expressionData$gene_id <- sub("gene_", "", expressionData$gene_id)
# geneMethData <- merge(geneMethData, expressionData, by = "gene_id")

#See distribution of methylation across genes:
cutOff <- 0.005
ggplot(geneMethData, mapping = aes(x = mCG*100)) + 
  geom_histogram(color = "black", fill = "lightgrey", bins = 50) + 
  theme_bw() + xlab("mCG (%)") + ylab("Count") + geom_vline(xintercept = cutOff*100, color = "red", linetype = "dashed")

#Describe genes as methylated or unmethylated using a cut-off, e.g. 0.5%
geneMethData$methylationStatus <- ifelse(geneMethData$mCG >= cutOff, "Methylated", "Unmethylated")
geneMethData$methylationStatus <- factor(geneMethData$methylationStatus, levels = c("Unmethylated", "Methylated"))

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


LGT_mCG_bars <- ggplot(gene_taxonomies, mapping = aes(x = domains, y = mCG*100, fill = domains)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("mCG (%)") + xlab("Domain of origin") + coord_cartesian(ylim=c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
LGT_mCG_bars
ggsave("N_fowleri/N_fowl_LGT_bars.pdf", LGT_mCG_bars, height = 10, width = 13, units = "cm")

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

sum(totalGeneTaxonomyFrequency[totalGeneTaxonomyFrequency$methylationStatus == "Methylated",]$frequency)
sum(totalGeneTaxonomyFrequency[totalGeneTaxonomyFrequency$methylationStatus == "Unmethylated",]$frequency)
LGT <- ggplot(totalGeneTaxonomyFrequency, mapping = aes(x = methylationStatus, y = frequency_proportion*100, fill = domains)) +
  geom_col() + theme_bw() + scale_fill_manual(values = c("#E69F00", "#ab2320", "#ff6e6e", "#56B4E9", "#003f5c", "#2a6dad", "#8097ff")) +
  ylab("Proportion (%)") + xlab(paste0("Methylation status (>", cutOff*100, "% mCG)"))
LGT
ggsave("N_fowleri/N_fowl_LGT.pdf", LGT, height = 15.14, width = 13.69, units = "cm")

bestMethylatedViralHits <- taxData[taxData$V18 == "Viruses",] %>%
  group_by(V1) %>%
  slice_min(order_by = V11, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!is.na(methylationStatus) & methylationStatus == "Methylated")

bestMethylatedViralHits %>% ggplot(mapping = aes(x = V19)) +
  geom_bar() + theme_bw() + xlab("Viral phylum of methylated genes") + ylab("Count")

bestMethylatedViralHits %>% filter(V19 == "Nucleocytoviricota") %>% ggplot(mapping = aes(x = V22)) +
  geom_bar() + theme_bw() + xlab("Viral order of Nucleocytoviricota methylated genes") + ylab("Count")

bestMethylatedViralHits %>% filter(V19 == "Nucleocytoviricota") %>% ggplot(mapping = aes(x = V23)) +
  geom_bar() + theme_bw() + xlab("Viral genus of Nucleocytoviricota methylated genes") + ylab("Count")



#Calculate enrichment of groups within 
data_grouped <- totalGeneTaxonomyFrequency %>% 
  group_by(domains) %>% 
  summarise(
    # Total_Frequency = sum(frequency),
    Frequency_point5pc_mCG = sum(frequency[methylationStatus == "Methylated"]),
    Frequency_Unmethylated = sum(frequency[methylationStatus == "Unmethylated"]),
    Proportion_point5pc_mCG = sum(frequency_proportion[methylationStatus == "Methylated"]),
    Proportion_Unmethylated = sum(frequency_proportion[methylationStatus == "Unmethylated"])
  )
total_point5pc_mCG <- sum(data_grouped$Frequency_point5pc_mCG)
total_Unmethylated <- sum(data_grouped$Frequency_Unmethylated)

data_grouped <- data_grouped %>%
  rowwise() %>%
  mutate(p_value = fisher.test(matrix(c(Frequency_point5pc_mCG,
                                        total_point5pc_mCG - Frequency_point5pc_mCG,
                                        Frequency_Unmethylated,
                                        total_Unmethylated - Frequency_Unmethylated),
                                      nrow = 2), alternative = "two.sided")$p.value) %>%
  ungroup()
data_grouped

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
