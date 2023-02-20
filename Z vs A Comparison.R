#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
# Load the data.
RR_length <- read_excel("scaff-length-RR.xlsx")
names(RR_length)[4] <- "Type"
names(RR_length)[5] <- "RR_Own_hdf"

# Modify the variables to have them in Mb and cM/Mb.
RR_length <- RR_length %>% mutate(Mb = Size / 1000000)
RR_length <- RR_length %>% mutate(cM = RR_Own_hdf / 0.00000001)

# We test for normality of the variables.
with(RR_length, shapiro.test(cM[Type == "A"])) # p = 0.07383
with(RR_length, shapiro.test(cM[Type == "Z"])) # p = 0.9299

# We test for homogeneity of variances.
Variances_Homogeneity <- var.test(cM ~ Type, data = RR_length) # p = 0.0847
Variances_Homogeneity
# Now we can compute the t-test to see if the means for Z and Autosomes are statistically different or not.
t_test <- t.test(cM ~ Type, data = RR_length, var.equal = TRUE, alternative = "less")
t_test # No significant (p = 0.6230) difference between Z and Autosomes although Z have a lower RR (7.65 in A - 7.03 in Z).

wilcox <- wilcox.test(cM ~ Type, data=RR_length)
wilcox # No significant again (p = 0.9195)

# Now we plot the results.
p <- ggplot(RR_length, aes(x=Type, y=cM, color = Type)) + geom_boxplot() + theme_light() +
  ylab("Rec. rate (cM/Mb)") + xlab("Chromosome Type")
p







