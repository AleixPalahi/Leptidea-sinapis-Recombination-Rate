# We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)

# We load the data.
codons <- read_table2("all-cds-pyrho.txt", col_names = FALSE)

names(codons)[1] <- "Chromosome"
names(codons)[2] <- "Initial_position"
names(codons)[3] <- "Final_position"
names(codons)[4] <- "Rate"
codons$cM <- as.numeric(codons$Rate * 100000000)

utr <- read_table2("all-utr.txt", col_names = FALSE)

names(utr)[1] <- "Chromosome"
names(utr)[2] <- "Initial_position"
names(utr)[3] <- "Final_position"
names(utr)[4] <- "Rate"
utr$cM <- as.numeric(utr$Rate * 100000000)

introns <- read_table2("all-introns-uniq.txt", col_names = FALSE)

names(introns)[1] <- "Chromosome"
names(introns)[2] <- "Initial_position"
names(introns)[3] <- "Final_position"
names(introns)[4] <- "Rate"
introns$cM <- as.numeric(introns$Rate * 100000000)

intergenic <- read_table2("all-intergenic.txt", col_names = FALSE)

names(intergenic)[1] <- "Chromosome"
names(intergenic)[2] <- "Initial_position"
names(intergenic)[3] <- "Final_position"
names(intergenic)[4] <- "Rate"
intergenic$cM <- as.numeric(intergenic$Rate * 100000000)

result1 <- wilcox.test(codons$cM, intergenic$cM, alternative = "two.sided")
print(result1$p.value)
result2 <- wilcox.test(utr$cM, intergenic$cM, alternative = "two.sided")
print(result2$p.value)
result3 <- wilcox.test(introns$cM, intergenic$cM, alternative = "two.sided")
print(result3$p.value)




LINE <- read_table2("all-LINE-pyrho.txt", col_names = FALSE)

names(LINE)[1] <- "Chromosome"
names(LINE)[2] <- "Initial_position"
names(LINE)[3] <- "Final_position"
names(LINE)[4] <- "Rate"
LINE$cM <- as.numeric(LINE$Rate * 100000000)

SINE <- read_table2("all-SINE-pyrho.txt", col_names = FALSE)

names(SINE)[1] <- "Chromosome"
names(SINE)[2] <- "Initial_position"
names(SINE)[3] <- "Final_position"
names(SINE)[4] <- "Rate"
SINE$cM <- as.numeric(SINE$Rate * 100000000)

LTR <- read_table2("all-LTR-pyrho.txt", col_names = FALSE)

names(LTR)[1] <- "Chromosome"
names(LTR)[2] <- "Initial_position"
names(LTR)[3] <- "Final_position"
names(LTR)[4] <- "Rate"
LTR$cM <- as.numeric(LTR$Rate * 100000000)

DNA <- read_table2("all-DNA-pyrho.txt", col_names = FALSE)

names(DNA)[1] <- "Chromosome"
names(DNA)[2] <- "Initial_position"
names(DNA)[3] <- "Final_position"
names(DNA)[4] <- "Rate"
DNA$cM <- as.numeric(DNA$Rate * 100000000)


result4 <- wilcox.test(LINE$cM, intergenic$cM, alternative = "two.sided")
print(result4)
result5 <- wilcox.test(SINE$cM, intergenic$cM, alternative = "greater")
print(result5)
result6 <- wilcox.test(LTR$cM, intergenic$cM, alternative = "two.sided")
print(result6)
result7 <- wilcox.test(DNA$cM, intergenic$cM, alternative = "two.sided")
print(result7)



