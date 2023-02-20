install.packages("vegan")
# We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(vegan)
library(cowplot)

# We load the data.
GC_windows <- read_table2("GC_100kb.txt", col_names = TRUE)

# We rename some variables.
names(GC_windows)[1] <- "Chromosome"
names(GC_windows)[2] <- "Initial_Position"
names(GC_windows)[3] <- "Final_Position"
names(GC_windows)[12] <- "Length"

# And correct the calculation on GC percentage.
GC_windows$GC_Content <- (GC_windows$`7_num_C` + GC_windows$`8_num_G`) / (GC_windows$`6_num_A` + GC_windows$`8_num_G` + GC_windows$`7_num_C` + GC_windows$`9_num_T`)

# Now we load the RR data for each window in the genome.
RR <- read.delim("GC_RR_windows.txt")

# Modify the variables to have them the recombination rate in cM/Mb.
RR <- RR %>% mutate(cM = Rec_Rate / 0.00000001)
RR <- RR %>% mutate(GC = GC_content * 100)

GC_RR_windows <- cbind(RR$Chromosome, GC_windows$Initial_Position, GC_windows$Final_Position, RR$cM, GC_windows$GC_Content)

GC_RR <- data.frame(GC_RR_windows)

# Shapiro-Wilk tests for normal distribution must be limited to 5000 samples, so we need to randomly subsample 5000 rows of
# the approx. 6000 we have. We can do this with the package vegan. We create a subsample for GC and another for RR.
Subsample_GC <- GC_RR[sample(1:nrow(GC_RR), 5000, replace=FALSE),]
Subsample_RR <- GC_RR[sample(1:nrow(GC_RR), 5000, replace=FALSE),]

Subsample <- GC_RR[sample(1:nrow(GC_RR), 0, replace=FALSE),]

# We test for normality of the variables.
Normality_GC <- shapiro.test(Subsample_GC$GC_content)
Normality_GC # p = 2.2e-16 (non-normal distribution)
Normality_cM <- shapiro.test(Subsample_RR$cM)
Normality_cM # p = 2.2e-16 (non-normal distribution)

# Since both parameters show a non-normal distribution, we must use a Spearman correlation to assess the association between
# both.

# Now the correlation scores and significancy are computed.
Spearman_test <- cor.test(GC_RR$X5, GC_RR$X4, method = 'spearman', alternative = "less", exact = FALSE)
Spearman_test
# The association between GC content and the recombination rate is SIGNIFICANT (p-value = 6.4e-05) and NEGATIVE (rho = -0.0494)
# Since the p-value loses a bit of importance and significance when the 

# Now we can plot the results.
p <- ggplot(GC_RR, aes(x=X5, y=X4)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.2) + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") +
  xlab ("GC Content") + ylab ("Recombination rate (cM/Mb)") 
p




# The same analyses can be repeated for each chromosome.

# First, the chromosome-specific dataset must be obtained.
chr_22 <- GC_RR %>% filter(X1 == 2)

# We test for normality of the variables.
Normality_chr_GC <- shapiro.test(chr_29$GC)
Normality_chr_GC
Normality_chr_cM <- shapiro.test(chr_29$cM)
Normality_chr_cM

# Now the correlation scores and significancy are computed.
Spearman_test <- cor.test(chr_2$X5, chr_2$X4, method = 'spearman', alternative = "less", exact = FALSE)
Spearman_test

# Now we can plot the results.
p <- ggplot(GC_RR, aes(x=X5, y=X4)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.2) + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") +
  xlab ("GC Content") + ylab ("Recombination rate (cM/Mb)") + xlim(0.29,0.39) + ylim(0,45)
p

summary(lm(X5~X4, GC_RR))
ggplot(GC_RR[GC_RR$X1 == "2",], aes(x=X2, y=scale(X4)))+geom_line()+geom_line(aes( y=scale(X5) ), col = "Blue")




GC_RR_terminal <- read_table2("GC_RR_terminal.txt", col_names = TRUE)

GC_RR_terminal <- GC_RR_terminal %>% mutate(cM = Rec_Rate * 100000000)

# Now the correlation scores and significancy are computed.
Spearman_test <- cor.test(GC_RR_terminal$cM, GC_RR_terminal$GC_content, method = 'pearson', alternative = "less", exact = FALSE)
Spearman_test

p <- ggplot(GC_RR_terminal, aes(x=GC_content, y=cM)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.2) + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") +
  xlab ("GC Content") + ylab ("Recombination rate (cM/Mb)")
p



# We can repeat the association essay on a 1Mb scale.
GC_1Mb <- read_table("GC_Mb_windows.txt", col_names = TRUE)
RR_1Mb <- read_table("RR-1Mb-sorted.txt", col_names = TRUE)


# We rename some variables.
names(GC_1Mb)[1] <- "Chromosome"
names(GC_1Mb)[2] <- "Initial_Position"
names(GC_1Mb)[3] <- "Final_Position"
names(GC_1Mb)[12] <- "Length"

names(RR_1Mb)[1] <- "Chromosome"
names(RR_1Mb)
# And correct the calculation on GC percentage.
GC_1Mb$GC_1Mb <- (GC_1Mb$`7_num_C` + GC_1Mb$`8_num_G`) / (GC_1Mb$`6_num_A` + GC_1Mb$`8_num_G` + GC_1Mb$`7_num_C` + GC_1Mb$`9_num_T`)

GC_RR_1Mb <- cbind(RR_1Mb$Chromosome, RR_1Mb$Start, RR_1Mb$End, RR_1Mb$RR, GC_1Mb$GC_1Mb)
GC_RR_1Mb <- as_tibble(GC_RR_1Mb)

names(GC_RR_1Mb)[1] <- "Chromosome"
names(GC_RR_1Mb)[2] <- "Initial_Position"
names(GC_RR_1Mb)[3] <- "Final_Position"
names(GC_RR_1Mb)[4] <- "RR"
names(GC_RR_1Mb)[5] <- "GC"

GC_RR_1Mb$RR <- GC_RR_1Mb$RR / 100
GC_RR_1Mb$cM <- GC_RR_1Mb$RR * 100000000

max(GC_RR_1Mb$RR)

GC_RR_1Mb <- GC_RR_1Mb[-426,]

Spearman_test <- cor.test(GC_RR_1Mb$RR, GC_RR_1Mb$GC, method = 'spearman', alternative = "two.sided", exact = FALSE)
Spearman_test

p <- ggplot(GC_RR_1Mb, aes(x=cM, y=GC)) + theme_light() +
  geom_point(shape=20, size=2.5, color="darkgrey", alpha=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold")) +
  #stat_smooth(alpha=0.8) +
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black", fill="lightgrey") +
  geom_label(label="rho = -0.066, p = 0.100", x = 25, y = 0.358, size = 7) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("GC Content")
p


grid <- plot_grid(p, q, labels = "AUTO")

png("GC-gene-density-RR-1Mb.png", units="in", width=10, height=6, res=300)
grid
dev.off()

ggsave("GC-gene-density-RR-1Mb.png", grid, width=12, height=4, units="in", dpi=300)

p <- ggplot(GC_RR_1Mb, aes(x=Initial_Position, y=scale(RR))) + theme_light() +
  geom_line(size=1.5, color="grey", alpha=1) +
  geom_line(aes(y=scale(GC)), color="blue")+
  #stat_smooth() +
  #geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") +
  xlab ("RR") + ylab ("GC")
p


GC_1Mb$Unmasked <- (GC_1Mb$`6_num_A` + GC_1Mb$`8_num_G` + GC_1Mb$`7_num_C` + GC_1Mb$`9_num_T`)/1000000
GC_RR_1Mb <- cbind(RR_1Mb$Chromosome, RR_1Mb$Start, RR_1Mb$End, RR_1Mb$RR, GC_1Mb$GC_1Mb, GC_1Mb$Unmasked)
GC_RR_1Mb <- as_tibble(GC_RR_1Mb)

names(GC_RR_1Mb)[1] <- "Chromosome"
names(GC_RR_1Mb)[2] <- "Initial_Position"
names(GC_RR_1Mb)[3] <- "Final_Position"
names(GC_RR_1Mb)[4] <- "RR"
names(GC_RR_1Mb)[5] <- "GC"
names(GC_RR_1Mb)[6] <- "Unmasked"

GC_RR_1Mb$RR <- GC_RR_1Mb$RR / 100
GC_RR_1Mb$cM <- GC_RR_1Mb$RR * 100000000

max(GC_RR_1Mb$cM)
GC_RR_1Mb <- GC_RR_1Mb[-426,]

Spearman_test <- cor.test(GC_RR_1Mb$cM, GC_RR_1Mb$Unmasked, method = 'spearman', alternative = "greater", exact = FALSE)
Spearman_test

p <- ggplot(GC_RR_1Mb, aes(x=Unmasked, y=cM)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=1) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") +
  geom_label(label="rho = -0.0752, p = 0.0314", x = 0.53, y = 32) +
  xlab ("Prop. Unmasked sites") + ylab ("Recombination rate (cM/Mb)")
p
