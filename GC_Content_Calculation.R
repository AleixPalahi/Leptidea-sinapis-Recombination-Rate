#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(PopGenome)
library(readxl)

# First, we are going to focus on the GC content in the whole genome.

# We load the data.
GC_windows <- read_table2("GC_100kb.txt", col_names = TRUE)

# We rename some variables.
names(GC_windows)[1] <- "Chromosome"
names(GC_windows)[2] <- "Initial_Position"
names(GC_windows)[3] <- "Final_Position"
names(GC_windows)[12] <- "Length"

# And correct the calculation on GC percentage.
GC_windows$GC_Content <- (GC_windows$`7_num_C` + GC_windows$`8_num_G`) / (GC_windows$`6_num_A` + GC_windows$`8_num_G` + GC_windows$`7_num_C` + GC_windows$`9_num_T`)

# Now we calculate the weighted GC content taking into account the length of each window 
# (the terminal windows are smaller than the rest!!!!)
sum(GC_windows$Length) # 600915751
GC_windows$Weight <- (GC_windows$Length / 600915751) * GC_windows$GC_Content
sum(GC_windows$Weight) # The genome-wide GC content is 32.65%


# We can now calculate it for every chromosome.
chr_stats <- read_excel("chromosome_stats.xlsx", sheet = "All_chromosomes", col_names = TRUE)
autosomes_stats <- read_excel("chromosome_stats.xlsx", sheet = "Autosomes", col_names = TRUE)

names(chr_stats)[3] <- "Length"
names(chr_stats)[4] <- "RR"
names(chr_stats)[11] <- "GC"

chr_stats$RR <- chr_stats$RR * 100000000
autosomes_stats$RR <- autosomes_stats$RR * 100000000

chr_stats$Length <- chr_stats$Length / 1000000
autosomes_stats$Length <- autosomes_stats$Length / 1000000

# Now the check the normality of RR, chromosome length and GC content.
Normality_chr_length <- shapiro.test(as.numeric(chr_stats$`Size (bp)`))
Normality_chr_length
Normality_chr_GC <- shapiro.test(as.numeric(chr_stats$`GC Content`))
Normality_chr_GC
Normality_chr_RR <- shapiro.test(as.numeric(chr_stats$`RR (cM/Mb)`))
Normality_chr_RR

Spearman_test <- cor.test(chr_stats$RR, chr_stats$GC, method = 'spearman', alternative = "greater", exact = FALSE)
Spearman_test # Non-significant for all (rho=0.1207, p=0.2664), autosomes as well (rho=0.0933, p=0.3251)

Spearman_test <- cor.test(chr_stats$Length, chr_stats$GC, method = 'spearman', alternative = "less", exact = FALSE)
Spearman_test # Non-significant for all (rho=-0.2458, p=0.09934), autosomes too (rho=-0.3005, p=0.06789)

Spearman_test <- cor.test(chr_stats$Length, chr_stats$RR, method = 'spearman', alternative = "less", exact = FALSE)
Spearman_test # Significant with all chromosomes (rho=-0.3182, p=0.04625), autosomes too (rho=0.0354, p=0.03542)

a <- ggplot(chr_stats, aes(x=Length, y=RR)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey") + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") +
  geom_label(label="p = 0.0463", x = 30, y = 15) +
  xlab ("Chromosome Size (Mb)") + ylab ("GC Content")
a


# Now, we can calculate the GC content in the hotspots.
GC_Hotspots <- read_table2("GC_hotspots.txt", col_names = TRUE)

names(GC_Hotspots)[1] <- "Chromosome"
names(GC_Hotspots)[2] <- "Initial_Position"
names(GC_Hotspots)[3] <- "Final_Position"
names(GC_Hotspots)[12] <- "Length"

GC_Hotspots$GC_Content <- (GC_Hotspots$`7_num_C` + GC_Hotspots$`8_num_G`) / (GC_Hotspots$`6_num_A` + GC_Hotspots$`8_num_G` + GC_Hotspots$`7_num_C` + GC_Hotspots$`9_num_T`)

total_length <- sum(GC_Hotspots$Length)

GC_Hotspots$Weight <- (GC_Hotspots$Length / 5131007) * GC_Hotspots$GC_Content
GC_Hotspots <- na.omit(GC_Hotspots) # 35 hotspots lay in masked regions (exons or TEs)
sum(GC_Hotspots$Weight) # The GC content in hotspots is 32.43%


# Finally, we can work on the flanking regions of each hotspots (5 kb on each side).
GC_left_flanks <- read_table2("GC_left_flanks.txt", col_names = TRUE)
GC_right_flanks <- read_table2("GC_right_flanks.txt", col_names = TRUE)

names(GC_left_flanks)[1] <- "Chromosome"
names(GC_left_flanks)[2] <- "Initial_Position"
names(GC_left_flanks)[3] <- "Final_Position"
names(GC_left_flanks)[12] <- "Length"

names(GC_right_flanks)[1] <- "Chromosome"
names(GC_right_flanks)[2] <- "Initial_Position"
names(GC_right_flanks)[3] <- "Final_Position"
names(GC_right_flanks)[12] <- "Length"

GC_right_flanks$GC_Content <- (GC_right_flanks$`7_num_C` + GC_right_flanks$`8_num_G`) / (GC_right_flanks$`6_num_A` + GC_right_flanks$`8_num_G` + GC_right_flanks$`7_num_C` + GC_right_flanks$`9_num_T`)


GC_flanks <- rbind(GC_left_flanks, GC_right_flanks)

GC_flanks$GC_Content <- (GC_flanks$`7_num_C` + GC_flanks$`8_num_G`) / (GC_flanks$`6_num_A` + GC_flanks$`8_num_G` + GC_flanks$`7_num_C` + GC_flanks$`9_num_T`)

flanks_length <- sum(GC_flanks$Length)

GC_flanks$Weight <- (GC_flanks$Length / 31229301) * GC_flanks$GC_Content
GC_flanks <- na.omit(GC_flanks) # 4 flanking regions are completely masked.
sum(GC_flanks$Weight) # The GC content in hotspots is 32.45%

# CDS have 0.449 GC content, based in the older assembly.
