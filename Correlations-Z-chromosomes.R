# We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(psych)
library(reshape)

# We load the data.
info <- read.delim("features_1Mb.txt", header = TRUE, sep = "\t")

info$Length <- as.numeric(as.character(info$Length))
info$Rec_Rate <- as.numeric(as.character(info$Rec_Rate))
info$GC <- as.numeric(as.character(info$GC))
info$Gene_Density <- as.numeric(as.character(info$Gene_Density))
info$DNA <- as.numeric(as.character(info$DNA))
info$LINE <- as.numeric(as.character(info$LINE))
info$SINE <- as.numeric(as.character(info$SINE))
info$LTR <- as.numeric(as.character(info$LTR))

info$Gene_Density <- (info$Gene_Density) / (info$Length)
info$Rec_Rate <- (info$Rec_Rate) * 100000000
info$DNA <- (info$DNA) / (info$Length)
info$LINE <- (info$LINE) / (info$Length)
info$SINE <- (info$SINE) / (info$Length)
info$LTR <- (info$LTR) / (info$Length)

info <- info %>% filter(Rec_Rate < 40)
info <- info %>% filter(Stable == "N")

# We do the tests
Pearson_DNA <- cor.test(info$Rec_Rate, info$DNA, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_DNA

Pearson_SINE <- cor.test(info$Rec_Rate, info$SINE, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_SINE

Pearson_LINE <- cor.test(info$Rec_Rate, info$LINE, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_LINE

Pearson_LTR <- cor.test(info$Rec_Rate, info$LTR, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_LTR

Pearson_GC <- cor.test(info$Rec_Rate, info$GC, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_GC

Pearson_Gene_Density <- cor.test(info$Rec_Rate, info$Gene_Density, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_Gene_Density
