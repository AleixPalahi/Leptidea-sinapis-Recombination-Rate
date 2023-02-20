install.packages("readxl")
install.packages("ggpubr")
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
Normality_Mb <- shapiro.test(RR_length$Mb)
Normality_Mb # p = 0.4246 (normal distribution)
Normality_cM <- shapiro.test(RR_length$cM)
Normality_cM # p = 0.05921 (normal distribution)

# Now the correlation scores and significancy are computed.
Pearson_test <- cor.test(log(RR_length$Mb), log(RR_length$cM), method = 'pearson', alternative = "less")
Pearson_test
Spearman_test <- cor.test(RR_length$Mb, RR_length$cM, method = 'spearman', alternative = "less")
Spearman_test
# We see that the correlation between chromosome size and recombination rate is NOT SIGNIFICANT.


# Now we can plot the results.
colors <- c("orange", "blue")
p <- ggplot(RR_length, aes(x=Mb, y=cM)) + theme_light() +
  theme(legend.title = element_blank(), legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  geom_point(shape=20, size=4, aes(colour = factor(Type))) + 
  scale_color_manual(values=c("orange", "blue")) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="grey") +
  geom_label(x = 32.5, y=14.5, label = c("rho = -0.292\np = 0.062"), colour = "black") +
  xlab ("Chromosome length (Mb)") + ylab ("Recombination rate (cM/Mb)") 
p

png("RR-Length-chromosomes.png", units="in", width=6, height=4, res=300)
p
dev.off()


# We can subsample the autosomes, as a more clear relationship can be observed in them.
Autosomes <- subset(RR_length, Type =='A')

Normality_Mb <- shapiro.test(Autosomes$Mb)
Normality_Mb
Normality_cM <- shapiro.test(Autosomes$cM)
Normality_cM

# Now the correlation scores and significancy are computed.
Pearson_test <- cor.test(Autosomes$Mb, Autosomes$cM, method = 'pearson', alternative = "less")
Pearson_test
Spearman_test <- cor.test(Autosomes$Mb, Autosomes$cM, method = 'spearman', alternative = "less")
Spearman_test

# Now we can plot the results.
p <- ggplot(Autosomes, aes(x=Mb, y=cM, color=Type)) + theme_light() +
  geom_point(shape=20, size=2.5) + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="grey") +
  xlab ("Chromosome length (Mb)") + ylab ("Recombination rate (cM/Mb)")
p

ggsave




# We can also calculate the genome-wide RR.
genome_size <- 600915751 # Total size of the reference genome.
RR_length$Weighted_RR <- (RR_length$Size / genome_size) * RR_length$cM
RR <- sum(RR_length$Weighted_RR) # GwRR = 7.37598 cM/Mb






