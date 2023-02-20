#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(PopGenome)

# Load the data and name the columns.
Recombination_rate_scaff_1 <- read_table2("scaff-1-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_1)[1] <- "Initial_Position"
names(Recombination_rate_scaff_1)[2] <- "Final_Position"
names(Recombination_rate_scaff_1)[3] <- "Rec_Rate"

Recombination_rate_scaff_1 <- Recombination_rate_scaff_1 %>% mutate(Mb = Initial_Position / 1000000)
Recombination_rate_scaff_1 <- Recombination_rate_scaff_1 %>% mutate(cM = Rec_Rate / 0.00000001)

SNP_Density_scaff_1 <- read_table2("scaff-1-100kb.snpden")
SNP_Density_scaff_1$BIN_MID <- SNP_Density_scaff_1$BIN_START + 50000
SNP_Density_scaff_1$Density_Normalized <- SNP_Density_scaff_1$`VARIANTS/KB`*20
SNP_Density_scaff_1$Mb <- SNP_Density_scaff_1$BIN_MID / 1000000

RR1 <- data.frame(Recombination_rate_scaff_1)
SNP1 <- data.frame(SNP_Density_scaff_1)

p <- ggplot() + theme_light () + geom_line(data=RR1, aes(x=Mb, y=cM), color="black") + geom_line(data=SNP1, aes(x=Mb, y=Density_Normalized), color="lightgrey")
p
