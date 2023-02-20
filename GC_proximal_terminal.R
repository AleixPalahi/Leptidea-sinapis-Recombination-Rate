#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(PopGenome)
library(readxl)

# First, we are going to focus on the GC content in the whole genome.

# We load the data.
GC_terminal <- read_table2("GC_content_terminal_windows.txt", col_names = TRUE)
GC_proximal <- read_table2("GC_content_proximal_windows.txt", col_names = TRUE)

# And correct the calculation on GC percentage.
GC_terminal$GC_Content <- (GC_terminal$`7_num_C` + GC_terminal$`8_num_G`) / (GC_terminal$`6_num_A` + GC_terminal$`8_num_G` + GC_terminal$`7_num_C` + GC_terminal$`9_num_T`)
GC_proximal$GC_Content <- (GC_proximal$`7_num_C` + GC_proximal$`8_num_G`) / (GC_proximal$`6_num_A` + GC_proximal$`8_num_G` + GC_proximal$`7_num_C` + GC_proximal$`9_num_T`)

# We check for the significancy of the differences in GC content in the subtelomeres and the
# central positions of the genome.
wilcox <- wilcox.test(GC_terminal$GC_Content, GC_proximal$GC_Content)
wilcox
wilcox$p.value # Significant difference between the subtelomeric regions and the rest of 
               # the genome (p = 5.72645e-22)
mean(GC_terminal$GC_Content) # 33.46 %
mean(GC_proximal$GC_Content) # 32.61 %

# Now we plot the results.

# To do it, we must indicate in each case if a window is proximal or terminal, and merge both
# data frames.
GC_terminal$Type <- "Terminal"
GC_proximal$Type <- "Proximal"

GC_proximal_terminal <- rbind(GC_terminal, GC_proximal)
GC_terminal <- GC_terminal[-10,]

max(GC_proximal_terminal$GC_Content)
GC_proximal_terminal <- GC_proximal_terminal[-40,] 

p <- ggplot(GC_proximal_terminal, aes(x=Type, y=GC_Content, color = Type)) + geom_boxplot(show.legend = FALSE) + theme_light() +
  ylab("GC Content (%)") + xlab("Relative Position") +
  annotate("text", x=1, y=0.28, label = "32.61%") + annotate("text", x=2, y=0.28, label = "33.38%")
p




# We load the data.
RR_terminal <- read_table2("GC_RR_terminal.txt", col_names = TRUE)
RR_proximal <- read_table2("GC_RR_proximal.txt", col_names = TRUE)

RR_terminal$Type <- "Terminal"
RR_proximal$Type <- "Proximal"

RR_proximal_terminal <- rbind(RR_terminal, RR_proximal)

wilcox <- wilcox.test(Rec_Rate ~ Type, data = RR_proximal_terminal,
                   exact = FALSE)
wilcox
wilcox$p.value
