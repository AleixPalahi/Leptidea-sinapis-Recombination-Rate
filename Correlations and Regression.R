# Load packages
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)

#Load data
info <- read.delim("features_1Mb.txt", header = TRUE, sep = "\t")

#Turn variables into numerical
info$Length <- as.numeric(as.character(info$Length))
info$Rec_Rate <- as.numeric(as.character(info$Rec_Rate))
info$GC <- as.numeric(as.character(info$GC))
info$Gene_Density <- as.numeric(as.character(info$Gene_Density))
info$DNA <- as.numeric(as.character(info$DNA))
info$LINE <- as.numeric(as.character(info$LINE))
info$SINE <- as.numeric(as.character(info$SINE))
info$LTR <- as.numeric(as.character(info$LTR))

#Transform them
info$Gene_Density <- (info$Gene_Density) / (info$Length)
info$Rec_Rate <- (info$Rec_Rate) * 100000000
info$DNA <- (info$DNA) / (info$Length)
info$LINE <- (info$LINE) / (info$Length)
info$SINE <- (info$SINE) / (info$Length)
info$LTR <- (info$LTR) / (info$Length)

#Filter the clear outlier
info <- info %>% filter(Rec_Rate < 40)

#CORRELATIONS
Pearson_DNA <- cor.test(info$Rec_Rate, info$DNA, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_DNA
Pearson_SINE <- cor.test(info$Rec_Rate, info$SINE, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_SINE
Pearson_LINE <- cor.test(info$Rec_Rate, info$LINE, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_LINE
Pearson_LTR <- cor.test(info$Rec_Rate, info$LTR, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_LTR

# Now we can plot the results.
p <- ggplot(info, aes(x=Rec_Rate, y=DNA)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = 0.08, p = 0.04", x = 26, y = 0.1385, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("DNA (fraction of window)") 
p

q <- ggplot(info, aes(x=Rec_Rate, y=LTR)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = -0.15, p = 2.00e-4", x = 26, y = 0.087, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("LTR (fraction of window)") 
q

r <- ggplot(info, aes(x=Rec_Rate, y=SINE)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = 0.29, p = 3.66e-13", x = 26, y = 0.0755, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("SINE (fraction of window)") 
r

s <- ggplot(info, aes(x=Rec_Rate, y=LINE)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = -0.19, p = 3.22e-6", x = 26, y = 0.34, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("LINE (fraction of window)") 
s

grid <- plot_grid(p, q, r, s, labels = "AUTO")
ggsave("TEclasses-RR-1Mb.png", plot = last_plot(), units= c("in"), width = 10, height = 6, dpi = 300)





#REGRESSION

# Explore the data.
meltData <- melt(info)
p <- ggplot(meltData, aes(factor(variable), value))
p + geom_boxplot() + facet_wrap(~variable, scale="free")
p

# Remove outliers.
info_no_outliers <- info %>% filter(GC < 0.34)
info_no_outliers <- info_no_outliers %>% filter(Gene_Density < 0.1) 
info_no_outliers <- info_no_outliers %>% filter(DNA < 0.125) 
info_no_outliers <- info_no_outliers %>% filter(DNA > 0.05)
info_no_outliers <- info_no_outliers %>% filter(LINE < 0.30) 
info_no_outliers <- info_no_outliers %>% filter(LINE >0.17) 
info_no_outliers <- info_no_outliers %>% filter(0.02 < SINE)
info_no_outliers <- info_no_outliers %>% filter(LTR < 0.125)

# Model including outliers.
model1 <- lm(Rec_Rate ~ GC + Gene_Density + DNA + LINE + SINE + LTR, data = info)
summary(model1)
# We eliminate variable with the highest p-value (LTR)
model1 <- lm(Rec_Rate ~ GC + Gene_Density + DNA + LINE + SINE, data = info)
summary(model1)
# We eliminate variable with the highest p-value (DNA)
model1 <- lm(Rec_Rate ~ GC + Gene_Density + LINE + SINE, data = info)
summary(model1)
# We eliminate variable with the highest p-value (Gene_Density)
model1 <- lm(Rec_Rate ~ GC + LINE + SINE, data = info)
summary(model1)
# We eliminate variable with the highest p-value (LTR)
model1 <- lm(Rec_Rate ~ LINE + SINE, data = info)
summary(model1)

# Model excluding outliers.
model2 <- lm(Rec_Rate ~ GC + Gene_Density + DNA + LINE + SINE + LTR, data = info_no_outliers)
summary(model2)
# We eliminate variable with the highest p-value (LTR)
model2 <- lm(Rec_Rate ~ GC + Gene_Density + DNA + LINE + SINE, data = info_no_outliers)
summary(model2)
# We eliminate variable with the highest p-value (GC)
model2 <- lm(Rec_Rate ~ Gene_Density + DNA + LINE + SINE, data = info_no_outliers)
summary(model2)
# We eliminate variable with the highest p-value (Gene_Density)
model2 <- lm(Rec_Rate ~ DNA + LINE + SINE, data = info_no_outliers)
summary(model2)
# We eliminate variable with the highest p-value (LTR)
model2 <- lm(Rec_Rate ~ LINE + SINE, data = info_no_outliers)
summary(model2)