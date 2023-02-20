#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(PopGenome)
library("readxl")

Windows_4 <- read_table2("scaffold-29-100kb.txt", col_names = FALSE)

names(Windows_4)[1] <- "Initial_Position"
names(Windows_4)[2] <- "Final_Position"
names(Windows_4)[3] <- "Rec_Rate"

Windows_4 <- Windows_4 %>% mutate(Rec_Rate = Rec_Rate / 100)
Windows_4 <- Windows_4 %>% mutate(cM = Rec_Rate / 0.00000001)

Windows_4 <- Windows_4 %>% filter(cM <= 0.75)

Recombination_rate_scaff_4 <- read_table2("scaff-29-miss0.8.hdf-27.rmap", col_names = FALSE)

names(Recombination_rate_scaff_4)[1] <- "Initial_Position"
names(Recombination_rate_scaff_4)[2] <- "Final_Position"
names(Recombination_rate_scaff_4)[3] <- "Rec_Rate"

Recombination_rate_scaff_4 <- Recombination_rate_scaff_4 %>% mutate(cM = Rec_Rate / 0.00000001)

Hotspots_4 <- Recombination_rate_scaff_4 %>% filter(Initial_Position >= 11000000) %>% filter(Initial_Position <= 13000000) %>% filter(cM <= 0.75)


# We load the data of all the hotspots detected.
coldspots <- read_excel("Coldspots_Length.xlsx", sheet = "Coldspots_CORRECT")

# 

Coldspots_filtered <- coldspots %>% filter(Length >= 100000) # 70 coldspots are very large (longer than 100000 bp)
Coldspots_filtered <- coldspots %>% filter(Length <= 10000) # 554 coldspots are short (shorter than 10000 bp)

mean(coldspots$Length) # 29768 bp long
median(coldspots$Length) # 12480 bp long
quantile(coldspots$Length)
sum(coldspots$Length)

ggplot(coldspots, aes(x=Length)) + geom_histogram(binwidth = 2500) # It is clearly visible how hotspots tend to be of moderate RR, while just some exceed 200 cM/Mb


# Regarding the distribution of hotspots, we can check whether they accumulate in the A or Z. We start by loading the data.
coldspots_per_scaffold <- read_excel("Coldspots_Length.xlsx", sheet = "Coldspots_per_scaffold")

#Now we transform the length to Mb
coldspots_per_scaffold <- coldspots_per_scaffold %>% mutate(Mb = Length / 1000000)
coldspots_per_scaffold <- coldspots_per_scaffold %>% mutate(coldspot_density = Coldspots / Mb)

with(coldspots_per_scaffold, shapiro.test(coldspot_density[Type == "A"])) # p = 0.56e-06
with(coldspots_per_scaffold, shapiro.test(coldspot_density[Type == "Z"])) # p = 0.523
Variances_Homogeneity <- var.test(coldspot_density ~ Type, data = coldspots_per_scaffold)
Variances_Homogeneity # p = 0.887 -> Both groups are normally distributed, and the variances are homogeneous.
t_test <- t.test(coldspot_density ~ Type, data = coldspots_per_scaffold, var.equal = TRUE, alternative = "less")
t_test # Higher frequency in Z chromosomes, but it is non-significant (p = 0.231)

wilcox <- wilcox.test(coldspot_density ~ Type, data = coldspots_per_scaffold, alternative = "less")
wilcox



# How abundant are the coldspots in the subtelomeric regions?

coldspots <- read_excel("Coldspots_Length.xlsx", sheet = "Coldspots_CORRECT")
coldspots_central <- coldspots %>% filter(Position == "Central")
coldspot_terminal <- coldspots %>% filter(Position == "Telomere")

sum(coldspots_central$Length) # 26899442
sum(coldspot_terminal$Length) # 11293533

sum <- sum(coldspots_central$Length) + sum(coldspot_terminal$Length)

percentage_coldspots_central <- sum(coldspots_central$Length) / sum * 100 # 70.4%
percentage_coldspots_terminal <- sum(coldspot_terminal$Length) / sum * 100 # 29.6%

coldspots_central_long <- coldspots_central %>% filter(Length > 100000)
mean(coldspots_central$Length)
mean
