install.packages("ggplot2")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("Rmisc")
install.packages("PopGenome")
install.packages("data.table")

library(ggplot2)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(PopGenome)
library(data.table)

cds_positions <- read_table("all-cds-pyrho.txt", sep = t, col_names = FALSE)

names(cds_positions)[1] <- "Chromosome"
names(cds_positions)[2] <- "Initial_position"
names(cds_positions)[3] <- "Final_position"
names(cds_positions)[4] <- "Rate"

cds_positions_expanded <- setDT(cds_positions)[, .(Position = Initial_position:Final_position, Rate = Rate), .(id = 1:nrow(cds_positions))]
cds_positions_expanded <- cds_positions_expanded %>% mutate(cM = Rate * 100000000)
max(cds_positions_expanded$cM)
cds_positions_expanded <- cds_positions_expanded %>% filter(cM <= 708)

sum(cds_positions_expanded$Length)
mean(cds_positions_expanded$cM) #7.324 cM/Mb (5.514 cM/Mb)
median(cds_positions_expanded$cM) # 0.913 cM/Mb (0.759 cM/Mb)
sd(cds_positions_expanded$cM) #22.67 cM/Mb (117.200 cM/Mb)
CI(cds_positions_expanded$cM, ci=0.95) # 7.333676  -  7.315082 cM/Mb (5.525 - 5.502)

hist(cds_positions_expanded$cM, breaks =1000, xlim = c(0,80))

p <- cds_positions_expanded %>%
  filter(cM<74) %>%
  ggplot(aes(x=cM)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  annotate("text", label = "mean = 5.514 cM/Mb\nmedian = 0.759 cM/Mb", x = 30, y = 3e+06) +
  labs(y= "CDS Position count", x = "Recombination rate (cM/Mb)") + xlim(-2,40) +
  theme(
    plot.title = element_text(size=15)
  )
p




intergenic_positions <- read_table("all-intergenic.txt", col_names = FALSE)

names(intergenic_positions)[1] <- "Chromosome"
names(intergenic_positions)[2] <- "Initial_position"
names(intergenic_positions)[3] <- "Final_position"
names(intergenic_positions)[4] <- "Rate"

intergenic_positions$Length <- intergenic_positions$Final_position - intergenic_positions$Initial_position
max(intergenic_positions$Length)


intergenic_positions_expanded <- setDT(intergenic_positions)[, .(Position = Initial_position:Final_position, Rate = Rate), .(id = 1:nrow(intergenic_positions))]
intergenic_positions_expanded <- intergenic_positions_expanded %>% mutate(cM = Rate * 100000000)
max(intergenic_positions_expanded$cM)
intergenic_positions_expanded <- intergenic_positions_expanded %>% filter(cM <= 10000)

sum(intergenic_positions_expanded$Length)
mean(intergenic_positions_expanded$cM) #7.357 cM/Mb
median(intergenic_positions_expanded$cM) # 1.473 cM/Mb
sd(intergenic_positions_expanded$cM) #24.071 cM/Mb
CI(intergenic_positions_expanded$cM, ci=0.95) # 7.360 - 7.355 cM/Mb

hist(intergenic_positions_expanded$cM, breaks =1000, xlim = c(0,80))

p <- intergenic_positions_expanded %>%
  filter(cM<74) %>%
  ggplot(aes(x=cM)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  annotate("text", label = "mean = 7.522 cM/Mb\nmedian = 1.473 cM/Mb", x = 30, y = 4.5e+07) +
  labs(y= "CDS Position count", x = "Recombination rate (cM/Mb)") + xlim(-2,40) +
  theme(
    plot.title = element_text(size=15)
  )
p








introns_positions <- read_table("all-intergenic.txt", col_names = FALSE)

names(introns_positions)[1] <- "Chromosome"
names(introns_positions)[2] <- "Initial_position"
names(introns_positions)[3] <- "Final_position"
names(introns_positions)[4] <- "Rate"

introns_positions$Length <- introns_positions$Final_position - introns_positions$Initial_position
max(introns_positions$Length)
introns_positions <- introns_positions %>% filter(Length <= 5000)


introns_positions_expanded <- setDT(introns_positions)[, .(Position = Initial_position:Final_position, Rate = Rate), .(id = 1:nrow(introns_positions))]
introns_positions_expanded <- introns_positions_expanded %>% mutate(cM = Rate * 100000000)
max(introns_positions_expanded$cM)
introns_positions_expanded <- introns_positions_expanded %>% filter(cM <= 708)

sum(cds_positions_expanded$Length)
mean(introns_positions_expanded$cM) #7.520 cM/Mb
median(introns_positions_expanded$cM) # 1.312 cM/Mb
sd(introns_positions_expanded$cM) #19.405 cM/Mb
CI(introns_positions_expanded$cM, ci=0.95) # 7.523 - 7.517 cM/Mb

hist(introns_positions_expanded$cM, breaks =1000, xlim = c(0,80))

p <- introns_positions_expanded %>%
  filter(cM<74) %>%
  ggplot(aes(x=cM)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  annotate("text", label = "mean = 7.520 cM/Mb\nmedian = 1.312 cM/Mb", x = 30, y = 4.5e+07) +
  labs(y= "CDS Position count", x = "Recombination rate (cM/Mb)") + xlim(-2,40) +
  theme(
    plot.title = element_text(size=15)
  )
p









utr_positions <- read.delim("all-utr.txt", header = FALSE)

names(utr_positions)[1] <- "Chromosome"
names(utr_positions)[2] <- "Initial_position"
names(utr_positions)[3] <- "Final_position"
names(utr_positions)[4] <- "Rate"

utr_positions$Length <- utr_positions$Final_position - utr_positions$Initial_position

utr_positions_expanded <- setDT(utr_positions)[, .(Position = Initial_position:Final_position, Rate = Rate), .(id = 1:nrow(utr_positions))]
utr_positions_expanded <- utr_positions_expanded %>% mutate(cM = Rate * 100000000)
max(utr_positions_expanded$cM)
utr_positions_expanded <- utr_positions_expanded %>% filter(cM <= 708)

mean(utr_positions_expanded$cM) #6.264 cM/Mb
median(utr_positions_expanded$cM) # 1.095 cM/Mb
sd(utr_positions_expanded$cM) #18.309 cM/Mb
CI(utr_positions_expanded$cM, ci=0.95) # 6.325 - 6.203 cM/Mb

hist(introns_positions_expanded$cM, breaks =1000, xlim = c(0,80))

p <- utr_positions_expanded %>%
  filter(cM<74) %>%
  ggplot(aes(x=cM)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  annotate("text", label = "mean = 6.264 cM/Mb\nmedian = 1.095 cM/Mb", x = 30, y = 90000) +
  labs(y= "CDS Position count", x = "Recombination rate (cM/Mb)") + xlim(-2,40) +
  theme(
    plot.title = element_text(size=15)
  )
p




SINE_positions <- read.delim("all-SINE-pyrho.txt", header = FALSE)

names(SINE_positions)[1] <- "Chromosome"
names(SINE_positions)[2] <- "Initial_position"
names(SINE_positions)[3] <- "Final_position"
names(SINE_positions)[4] <- "Rate"

SINE_positions$Length <- SINE_positions$Final_position - SINE_positions$Initial_position

SINE_positions_expanded <- setDT(SINE_positions)[, .(Position = Initial_position:Final_position, Rate = Rate), .(id = 1:nrow(SINE_positions))]
SINE_positions_expanded <- SINE_positions_expanded %>% mutate(cM = Rate * 100000000)
max(SINE_positions_expanded$cM)
SINE_positions_expanded <- SINE_positions_expanded %>% filter(cM <= 708)

mean(SINE_positions_expanded$cM) #7.399 cM/Mb
median(SINE_positions_expanded$cM) # 1.630 cM/Mb
sd(SINE_positions_expanded$cM) #18.410 cM/Mb
CI(SINE_positions_expanded$cM, ci=0.95) # 7.411 - 7.388 cM/Mb

hist(SINE_positions_expanded$cM, breaks =1000, xlim = c(0,80))

p <- SINE_positions_expanded %>%
  filter(cM<74) %>%
  ggplot(aes(x=cM)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  annotate("text", label = "mean = 7.482 cM/Mb\nmedian = 1.630 cM/Mb", x = 30, y = 2500000) +
  labs(y= "SINE Position count", x = "Recombination rate (cM/Mb)") + xlim(-2,40) +
  theme(
    plot.title = element_text(size=15)
  )
p




LTR_positions <- read.delim("all-LTR-pyrho.txt", header = FALSE)

names(LTR_positions)[1] <- "Chromosome"
names(LTR_positions)[2] <- "Initial_position"
names(LTR_positions)[3] <- "Final_position"
names(LTR_positions)[4] <- "Rate"

LTR_positions$Length <- LTR_positions$Final_position - LTR_positions$Initial_position

LTR_positions_expanded <- setDT(LTR_positions)[, .(Position = Initial_position:Final_position, Rate = Rate), .(id = 1:nrow(LTR_positions))]
LTR_positions_expanded <- LTR_positions_expanded %>% mutate(cM = Rate * 100000000)
max(LTR_positions_expanded$cM)
LTR_positions_expanded <- LTR_positions_expanded %>% filter(cM <= 708)

mean(LTR_positions_expanded$cM) #7.753 cM/Mb
median(LTR_positions_expanded$cM) # 1.387 cM/Mb
sd(LTR_positions_expanded$cM) #23.244 cM/Mb
CI(LTR_positions_expanded$cM, ci=0.95) # 7.771 - 7.734 cM/Mb

hist(LTR_positions_expanded$cM, breaks =1000, xlim = c(0,80))

p <- LTR_positions_expanded %>%
  filter(cM<74) %>%
  ggplot(aes(x=cM)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  annotate("text", label = "mean = 7.767 cM/Mb\nmedian = 1.387 cM/Mb", x = 30, y = 1000000) +
  labs(y= "SINE Position count", x = "Recombination rate (cM/Mb)") + xlim(-2,40) +
  theme(
    plot.title = element_text(size=15)
  )
p




DNA_positions <- read.delim("all-DNA-pyrho.txt", header = FALSE)

names(DNA_positions)[1] <- "Chromosome"
names(DNA_positions)[2] <- "Initial_position"
names(DNA_positions)[3] <- "Final_position"
names(DNA_positions)[4] <- "Rate"

DNA_positions$Length <- DNA_positions$Final_position - DNA_positions$Initial_position

DNA_positions_expanded <- setDT(DNA_positions)[, .(Position = Initial_position:Final_position, Rate = Rate), .(id = 1:nrow(DNA_positions))]
DNA_positions_expanded <- DNA_positions_expanded %>% mutate(cM = Rate * 100000000)
max(DNA_positions_expanded$cM)
DNA_positions_expanded <- DNA_positions_expanded %>% filter(cM <= 708)

mean(DNA_positions_expanded$cM) #6.794 cM/Mb
median(DNA_positions_expanded$cM) # 1.347 cM/Mb
sd(DNA_positions_expanded$cM) #18.717 cM/Mb
CI(DNA_positions_expanded$cM, ci=0.95) # 6.802 - 6.787 cM/Mb

hist(DNA_positions_expanded$cM, breaks =1000, xlim = c(0,80))

p <- DNA_positions_expanded %>%
  filter(cM<74) %>%
  ggplot(aes(x=cM)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  annotate("text", label = "mean = 7.767 cM/Mb\nmedian = 1.387 cM/Mb", x = 30, y = 6000000) +
  labs(y= "SINE Position count", x = "Recombination rate (cM/Mb)") + xlim(-2,40) +
  theme(
    plot.title = element_text(size=15)
  )
p




LINE_positions <- read.delim("all-LINE-pyrho.txt", header = FALSE)

names(LINE_positions)[1] <- "Chromosome"
names(LINE_positions)[2] <- "Initial_position"
names(LINE_positions)[3] <- "Final_position"
names(LINE_positions)[4] <- "Rate"

DNA_positions$Length <- DNA_positions$Final_position - DNA_positions$Initial_position

LINE_positions_expanded <- setDT(LINE_positions)[, .(Position = Initial_position:Final_position, Rate = Rate), .(id = 1:nrow(LINE_positions))]
LINE_positions_expanded <- LINE_positions_expanded %>% mutate(cM = Rate * 100000000)
max(LINE_positions_expanded$cM)
LINE_positions_expanded <- LINE_positions_expanded %>% filter(cM <= 708)

mean(LINE_positions_expanded$cM) #6.562 cM/Mb
median(LINE_positions_expanded$cM) # 1.218 cM/Mb
sd(LINE_positions_expanded$cM) # 18,498 cM/Mb
CI(LINE_positions_expanded$cM, ci=0.95) # 6.568 - 6.557 cM/Mb

hist(LINE_positions_expanded$cM, breaks =1000, xlim = c(0,80))

p <- LINE_positions_expanded %>%
  filter(cM<74) %>%
  ggplot(aes(x=cM)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  annotate("text", label = "mean = 6.621 cM/Mb\nmedian = 1.218 cM/Mb", x = 30, y = 10000000) +
  labs(y= "SINE Position count", x = "Recombination rate (cM/Mb)") + xlim(-2,40) +
  theme(
    plot.title = element_text(size=15)
  )
p



RR_all_features <- read.delim("all_features.txt")

colors <- c("Lightblue"="Lightblue", "lightgrey"="lightgrey", "lightgreen"="lightgreen")
p <- ggplot(RR_all_features) +
  geom_bar( aes(x= reorder(Feature, Order), y=Mean), stat="identity", fill=RR_all_features$Type, alpha=0.5) +
  theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  geom_errorbar( aes(x=Feature, ymin=Lower_CI, ymax=Upper_CI), width=0.3, colour="orange", alpha=0.9, size=0.5) +
  geom_hline(yintercept=7.37, linetype="dashed", color = "black", size=0.5) +
  labs(y= "Recombination rate (cM/Mb)", x = "Feature")
  #scale_color_manual(values = colors)
p

p <- ggplot(RR_all_features) +
  geom_bar(aes(x= reorder(Feature, Order), y=Mean), stat="identity", alpha=0.5) +
  theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        axis.text=element_text(size=8), axis.title=element_text(size=8,face="bold")) +
  geom_errorbar( aes(x=Feature, ymin=Lower_CI, ymax=Upper_CI), width=0.3, colour="orange", alpha=0.9, size=0.5) +
  geom_hline(yintercept=7.37, linetype="dashed", color = "black", size=0.3) +
  labs(y= "Recombination rate (cM/Mb)", x = "Feature")
#scale_color_manual(values = colors)
p


png("RR-in-features.png", units="in", width=6, height=4, res=300)
p
dev.off()

