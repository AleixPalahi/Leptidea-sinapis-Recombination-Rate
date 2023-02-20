# We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)

# We load the data.
gene_density_1Mb <- read_table2("cds_1Mb.txt", col_names = FALSE)
RR_1Mb <- read_table2("RR-1Mb-sorted.txt", col_names = TRUE)

scaff_29_RR <- read_table2("scaffold-29-1Mb.txt", col_names = FALSE)

windows_RR <- rbind (scaff_1_RR, scaff_10_RR, scaff_11_RR, scaff_12_RR, scaff_13_RR, scaff_14_RR, scaff_15_RR, scaff_16_RR, scaff_17_RR, scaff_18_RR, scaff_19_RR, scaff_2_RR, scaff_20_RR, scaff_21_RR, scaff_22_RR, scaff_23_RR, scaff_24_RR, scaff_25_RR, scaff_26_RR, scaff_27_RR, scaff_28_RR, scaff_29_RR, scaff_3_RR, scaff_4_RR, scaff_5_RR, scaff_6_RR, scaff_7_RR, scaff_8_RR, scaff_9_RR)

gene_density_rr_1Mb <- cbind(RR_1Mb$Chromosome, RR_1Mb$Start, RR_1Mb$End, RR_1Mb$RR, gene_density_1Mb$X5)

gene_density_rr_1Mb <- as.tibble(gene_density_rr_1Mb)

names(gene_density_rr_1Mb)[1] <- "Chromosome"
names(gene_density_rr_1Mb)[2] <- "Initial_position"
names(gene_density_rr_1Mb)[3] <- "Final_position"
names(gene_density_rr_1Mb)[4] <- "Rate"
names(gene_density_rr_1Mb)[5] <- "Gene_positions"

gene_density_rr_1Mb$Rate <- as.numeric(as.character(gene_density_rr_1Mb$Rate))
gene_density_rr_1Mb$Initial_position <- as.numeric(as.character(gene_density_rr_1Mb$Initial_position))
gene_density_rr_1Mb$Final_position <- as.numeric(as.character(gene_density_rr_1Mb$Final_position))
gene_density_rr_1Mb$Gene_positions <- as.numeric(as.character(gene_density_rr_1Mb$Gene_positions))


gene_density_rr_1Mb$Rate <- gene_density_rr_1Mb$Rate / 100
gene_density_rr_1Mb$cM <- gene_density_rr_1Mb$Rate * 100000000

gene_density_rr_1Mb[613, 3] = 24821562

gene_density_rr_1Mb$Gene_density <- gene_density_rr_1Mb$Gene_positions / (gene_density_rr_1Mb$Final_position - gene_density_rr_1Mb$Initial_position)

# We test for normality of the variables.
Normality_gene <- shapiro.test(gene_density_rr_1Mb$Gene_density)
Normality_gene # p = 2.2e-16 (non-normal distribution)
Normality_RR <- shapiro.test(gene_density_rr_1Mb$cM)
Normality_RR # p = 2.2e-16 (non-normal distribution)

gene_density_rr_1Mb_filtered <- gene_density_rr_1Mb %>% filter(Gene_density <= 0.21) %>% filter(cM <= 36)

Spearman_test <- cor.test(gene_density_rr_1Mb_filtered$cM, gene_density_rr_1Mb_filtered$Gene_density, method = 'spearman', alternative = "two.sided", exact = FALSE)
Spearman_test

max(gene_density_rr_1Mb_filtered$Gene_density)

# Now we can plot the results.
q <- ggplot(gene_density_rr_1Mb_filtered, aes(x=cM, y=Gene_density)) + theme_light() +
  geom_point(shape=20, size=2.5, color="darkgrey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black", fill ="lightgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold")) +
  geom_label(label="rho = -0.036, p = 0.381", x = 25, y = 0.20, size = 7) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("CDS (proportion of window)") 
q

png("gene-density-RR-1Mb.png", units="in", width=6, height=4, res=300)
p
dev.off()

# Now we can do the analyses at a chromosome level.
cds_chromosome <- gene_density_1Mb %>% group_by(X1) %>%
  summarize(Positions = sum(X5))

# We update the file with the chromosome-level data and upload it.
cds_chromosome <- read_excel("chromosome_stats.xlsx", sheet = "All_chromosomes")

cds_chromosome$`RR (cM/Mb)` <- as.numeric(as.character(cds_chromosome$`RR (cM/Mb)`))
cds_chromosome$`Size (bp)` <- as.numeric(as.character(cds_chromosome$`Size (bp)`))
cds_chromosome$Gene_positions <- as.numeric(as.character(cds_chromosome$Gene_positions))

cds_chromosome$Gene_density <- cds_chromosome$Gene_positions / cds_chromosome$`Size (bp)`
cds_chromosome$cM <- cds_chromosome$`RR (cM/Mb)` * 100000000

# We test for normality of the variables.
Normality_gene <- shapiro.test(cds_chromosome$Gene_density)
Normality_gene # p = 0.3825 (normal distribution)
Normality_RR <- shapiro.test(cds_chromosome$cM)
Normality_RR # p = 0.0315 (non-normal distribution)

Spearman_test <- cor.test(cds_chromosome$cM, cds_chromosome$Gene_density, method = 'pearson', alternative = "greater", exact = FALSE)
Spearman_test

# Now we can plot the results.
p <- ggplot(cds_chromosome, aes(x=cM, y=Gene_density)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") +
  geom_label(label="rho = 0.096\np = 0.3101", x = 14, y = 0.050) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("CDS density") 
p
