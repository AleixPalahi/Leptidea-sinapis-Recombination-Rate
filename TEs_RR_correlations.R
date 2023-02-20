# We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(vegan)
library(cowplot)

# We load the data.
DNA_windows <- read_table2("TE_DNA-recrate_1Mb_coverage.txt", col_names = FALSE)

DNA_windows <- DNA_windows %>% mutate(cM = X4 / 0.00000001)

DNA_windows <- DNA_windows %>% filter(cM < 40)

# Now the correlation scores and significancy are computed.
Spearman_test <- cor.test(DNA_windows$cM, DNA_windows$X8, method = 'pearson', alternative = "greater", exact = FALSE)
Spearman_test

# Now we can plot the results.
p <- ggplot(DNA_windows, aes(x=cM, y=X8)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = 0.184, p = 2.32e-06", x = 26, y = 0.138, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("DNA (fraction of window)") 
p


# We load the data.
LTR_windows <- read_table2("TE_LTR-recrate_1Mb_coverage.txt", col_names = FALSE)

LTR_windows <- LTR_windows %>% mutate(cM = X4 / 0.00000001)

LTR_windows <- LTR_windows %>% filter(cM < 40)

# Now the correlation scores and significancy are computed.
Spearman_test <- cor.test(LTR_windows$cM, LTR_windows$X8, method = 'spearman', alternative = "greater", exact = FALSE)
Spearman_test

# Now we can plot the results.
q <- ggplot(LTR_windows, aes(x=cM, y=X8)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = -0.051, p = 0.11", x = 26, y = 0.058, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("LTR (fraction of window)") 
q





# We load the data.
SINE_windows <- read_table2("TE_SINE-recrate_1Mb_coverage.txt", col_names = FALSE)

SINE_windows <- SINE_windows %>% mutate(cM = X4 / 0.00000001)

SINE_windows <- SINE_windows %>% filter(cM < 40)

# Now the correlation scores and significancy are computed.
Spearman_test <- cor.test(SINE_windows$cM, SINE_windows$X8, method = 'spearman', alternative = "greater", exact = FALSE)
Spearman_test

# Now we can plot the results.
r <- ggplot(SINE_windows, aes(x=cM, y=X8)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = 0.307, p = 4.07e-15", x = 26, y = 0.076, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("SINE (fraction of window)") 
r




# We load the data.
LINE_windows <- read_table2("TE_LINE-recrate_1Mb_coverage.txt", col_names = FALSE)

LINE_windows <- LINE_windows %>% mutate(cM = X4 / 0.00000001)

LINE_windows <- LINE_windows %>% filter(cM < 40)

# Now the correlation scores and significancy are computed.
Spearman_test <- cor.test(LINE_windows$cM, LINE_windows$X8, method = 'spearman', alternative = "less", exact = FALSE)
Spearman_test

# Now we can plot the results.
s <- ggplot(LINE_windows, aes(x=cM, y=X8)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = -0.044, p = 0.14", x = 26, y = 0.3, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("LINE (fraction of window)") 
s

grid <- plot_grid(p, q, r, s, labels = "AUTO")
grid

png("TE-RR-1Mb.png", units="in", width=10, height=6, res=300)
grid
dev.off()







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
  geom_label(label="rho = 0.09, p = 0.03", x = 28.85, y = 0.1385, size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("DNA (fraction of window)") +
  scale_x_continuous(expand = c(0.01,0.01))
p

q <- ggplot(info, aes(x=Rec_Rate, y=LTR)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = -0.11, p = 8.50e-3", x = 27.7, y = 0.088, size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("LTR (fraction of window)") +
  scale_x_continuous(expand = c(0.01,0.01))
q

r <- ggplot(info, aes(x=Rec_Rate, y=SINE)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = 0.29, p = 3.31e-13", x = 27.65, y =0.0755, size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("SINE (fraction of window)") +
  scale_x_continuous(expand = c(0.01,0.01))
r

s <- ggplot(info, aes(x=Rec_Rate, y=LINE)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = -0.19, p = 3.01e-6", x = 27.8, y = 0.338, size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("LINE (fraction of window)") +
  scale_x_continuous(expand = c(0.01,0.01))
s

grid <- plot_grid(p, q, r, s, labels = "AUTO", label_size = 10)
ggsave("TEclasses-RR-1Mb.png", plot = last_plot(), units= c("in"), width = 10, height = 6, dpi = 300)


Pearson_gene <- cor.test(info$Rec_Rate, info$Gene_Density, method = 'spearman', alternative = "less", exact = FALSE)
Pearson_gene
Pearson_GC <- cor.test(info$Rec_Rate, info$GC, method = 'spearman', alternative = "two.sided", exact = FALSE)
Pearson_GC

t <- ggplot(info, aes(x=Rec_Rate, y=Gene_Density)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = -0.04, p = 0.17", x = 29.5, y = 0.215, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  scale_x_continuous(expand = c(0.01,0.01)) + scale_y_continuous(expand = c(0.01,0.01)) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("CDS (Fraction of window)") + ylim(0,0.22)
t

u <- ggplot(info, aes(x=Rec_Rate, y=GC)) + theme_light() +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = -0.07, p = 0.10", x = 29.5, y = 0.358, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) + 
  scale_x_continuous(expand = c(0.01,0.01)) + scale_y_continuous(expand = c(0.001,0.001)) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("GC Content")
u

grid <- plot_grid(u, t, labels = "AUTO", label_size = 10)
grid
ggsave("GC-Gene-1Mb.png", plot = last_plot(), units= c("in"), width = 10, height = 3, dpi = 300)
