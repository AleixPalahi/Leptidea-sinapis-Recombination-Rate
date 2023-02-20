install.packages("psych")
install.packages("reshape")

# We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(psych)
library(reshape)

#We load all the data necessary for each 1Mb window.
gene_density_1Mb <- read_table2("cds_1Mb.txt", col_names = FALSE)
GC_1Mb <- read_table2("GC_Mb_windows.txt", col_names = TRUE)
DNA_windows <- read_table2("TE_DNA-recrate_1Mb_coverage.txt", col_names = FALSE)
LINE_windows <- read_table2("TE_LINE-recrate_1Mb_coverage.txt", col_names = FALSE)
SINE_windows <- read_table2("TE_SINE-recrate_1Mb_coverage.txt", col_names = FALSE)
LTR_windows <- read_table2("TE_LTR-recrate_1Mb_coverage.txt", col_names = FALSE)

#We have to sort the chromosome column, 
DNA_windows$X1 <- gsub("^.{0,13}", "", DNA_windows$X1)
DNA_windows$X1 <- gsub("^1$", "01", DNA_windows$X1)
DNA_windows$X1 <- gsub("^2$", "02", DNA_windows$X1)
DNA_windows$X1 <- gsub("^3$", "03", DNA_windows$X1)
DNA_windows$X1 <- gsub("^4$", "04", DNA_windows$X1)
DNA_windows$X1 <- gsub("^5$", "05", DNA_windows$X1)
DNA_windows$X1 <- gsub("^6$", "06", DNA_windows$X1)
DNA_windows$X1 <- gsub("^7$", "07", DNA_windows$X1)
DNA_windows$X1 <- gsub("^8$", "08", DNA_windows$X1)
DNA_windows$X1 <- gsub("^9$", "09", DNA_windows$X1)
DNA_windows <- DNA_windows[order(DNA_windows$X1,DNA_windows$X2),]

SINE_windows$X1 <- gsub("^.{0,13}", "", SINE_windows$X1)
SINE_windows$X1 <- gsub("^1$", "01", SINE_windows$X1)
SINE_windows$X1 <- gsub("^2$", "02", SINE_windows$X1)
SINE_windows$X1 <- gsub("^3$", "03", SINE_windows$X1)
SINE_windows$X1 <- gsub("^4$", "04", SINE_windows$X1)
SINE_windows$X1 <- gsub("^5$", "05", SINE_windows$X1)
SINE_windows$X1 <- gsub("^6$", "06", SINE_windows$X1)
SINE_windows$X1 <- gsub("^7$", "07", SINE_windows$X1)
SINE_windows$X1 <- gsub("^8$", "08", SINE_windows$X1)
SINE_windows$X1 <- gsub("^9$", "09", SINE_windows$X1)
SINE_windows <- SINE_windows[order(SINE_windows$X1,SINE_windows$X2),]

LINE_windows$X1 <- gsub("^.{0,13}", "", LINE_windows$X1)
LINE_windows$X1 <- gsub("^1$", "01", LINE_windows$X1)
LINE_windows$X1 <- gsub("^2$", "02", LINE_windows$X1)
LINE_windows$X1 <- gsub("^3$", "03", LINE_windows$X1)
LINE_windows$X1 <- gsub("^4$", "04", LINE_windows$X1)
LINE_windows$X1 <- gsub("^5$", "05", LINE_windows$X1)
LINE_windows$X1 <- gsub("^6$", "06", LINE_windows$X1)
LINE_windows$X1 <- gsub("^7$", "07", LINE_windows$X1)
LINE_windows$X1 <- gsub("^8$", "08", LINE_windows$X1)
LINE_windows$X1 <- gsub("^9$", "09", LINE_windows$X1)
LINE_windows <- LINE_windows[order(LINE_windows$X1,LINE_windows$X2),]

LTR_windows$X1 <- gsub("^.{0,13}", "", LTR_windows$X1)
LTR_windows$X1 <- gsub("^1$", "01", LTR_windows$X1)
LTR_windows$X1 <- gsub("^2$", "02", LTR_windows$X1)
LTR_windows$X1 <- gsub("^3$", "03", LTR_windows$X1)
LTR_windows$X1 <- gsub("^4$", "04", LTR_windows$X1)
LTR_windows$X1 <- gsub("^5$", "05", LTR_windows$X1)
LTR_windows$X1 <- gsub("^6$", "06", LTR_windows$X1)
LTR_windows$X1 <- gsub("^7$", "07", LTR_windows$X1)
LTR_windows$X1 <- gsub("^8$", "08", LTR_windows$X1)
LTR_windows$X1 <- gsub("^9$", "09", LTR_windows$X1)
LTR_windows <- LTR_windows[order(LTR_windows$X1,LTR_windows$X2),]

names(GC_1Mb)[6] <- "A"
names(GC_1Mb)[7] <- "C"
names(GC_1Mb)[8] <- "G"
names(GC_1Mb)[9] <- "T"
GC_1Mb$GC_Content <- (GC_1Mb$C + GC_1Mb$G) / (GC_1Mb$A + GC_1Mb$G + GC_1Mb$C + GC_1Mb$T)

GC21 <- GC_1Mb %>% filter(`#1_usercol` == "HiC_scaffold_21")
gene_density21  <- gene_density_1Mb %>% filter(X1 == "HiC_scaffold_21")
DNA21 <- DNA_windows %>% filter(X1 == 21)
LINE21 <- LINE_windows %>% filter(X1 == 21)
SINE21 <- SINE_windows %>% filter(X1 == 21)
LTR21 <- LTR_windows %>% filter(X1 == 21)

info21 <- cbind(SINE21$X1, GC21$`2_usercol`, GC21$`3_usercol`, GC21$`12_seq_len`, LINE21$X4, GC21$GC_Content, gene_density21$X5, DNA21$X6, LINE21$X6, SINE21$X6, LTR21$X6)
info21 <- data.frame(info21)



#We create the matrix with all the info for each 1Mb window.
info <- cbind(SINE_windows$X1, GC_1Mb$`2_usercol`, GC_1Mb$`3_usercol`, GC_1Mb$`12_seq_len`, LINE_windows$X4, GC_1Mb$GC_Content, gene_density_1Mb$X5, DNA_windows$X6, LINE_windows$X6, SINE_windows$X6, LTR_windows$X6)
info <- data.frame(info)












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

# We explore the data.
meltData <- melt(info)
p <- ggplot(meltData, aes(factor(variable), value))
p + geom_boxplot() + facet_wrap(~variable, scale="free")
p #Not many outliers, max 10-12 per variable (out of 610...)

info_no_outliers <- info %>% filter(GC < 0.34)
info_no_outliers <- info_no_outliers %>% filter(Gene_Density < 0.1) 
info_no_outliers <- info_no_outliers %>% filter(DNA < 0.125) 
info_no_outliers <- info_no_outliers %>% filter(DNA > 0.05)
info_no_outliers <- info_no_outliers %>% filter(LINE < 0.30) 
info_no_outliers <- info_no_outliers %>% filter(LINE >0.17) 
info_no_outliers <- info_no_outliers %>% filter(0.02 < SINE)
info_no_outliers <- info_no_outliers %>% filter(LTR < 0.125)



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

model <- lm(Rec_Rate ~ GC + Gene_Density + DNA + LINE + SINE + LTR, data = info)
summary(model)
# We eliminate variable with the highest p-value (LTR)
model <- lm(Rec_Rate ~ GC + Gene_Density + DNA + LINE + SINE, data = info)
summary(model)
# We eliminate variable with the highest p-value (DNA)
model <- lm(Rec_Rate ~ GC + Gene_Density + LINE + SINE, data = info)
summary(model)
# We eliminate variable with the highest p-value (Gene_Density)
model <- lm(Rec_Rate ~ GC + LINE + SINE, data = info)
summary(model)
# We eliminate variable with the highest p-value (LTR)
model <- lm(Rec_Rate ~ LINE + SINE, data = info)
summary(model)



singleDNA <- lm(Rec_Rate ~ DNA, data = info)
summary(singleDNA)

model <- lm(Rec_Rate ~ GC + Gene_Density + DNA + LINE + SINE + LTR, data = info_no_outliers)
summary(model)
# We eliminate variable with the highest p-value (LTR)
model <- lm(Rec_Rate ~ GC + Gene_Density + DNA + LINE + SINE, data = info_no_outliers)
summary(model)
# We eliminate variable with the highest p-value (GC)
model <- lm(Rec_Rate ~ Gene_Density + DNA + LINE + SINE, data = info_no_outliers)
summary(model)
# We eliminate variable with the highest p-value (Gene_Density)
model <- lm(Rec_Rate ~ DNA + LINE + SINE, data = info_no_outliers)
summary(model)
# We eliminate variable with the highest p-value (LTR)
model <- lm(Rec_Rate ~ LINE + SINE, data = info_no_outliers)
summary(model)



