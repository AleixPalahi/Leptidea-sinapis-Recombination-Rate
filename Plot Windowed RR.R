#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(PopGenome)
library(ggpubr)
library(vegan)
library(cowplot)

# Load the data and name the columns.
Windows_scaff_1 <- read_table2("scaffold-17-100kb-own.txt", col_names = FALSE)
names(Windows_scaff_1)[1] <- "Window_Start"
names(Windows_scaff_1)[2] <- "Window_End"
names(Windows_scaff_1)[3] <- "Window_RR"

Windows_scaff_1<- Windows_scaff_1 %>% mutate(Window_Mid = Window_Start + 50000)
Windows_scaff_1 <- Windows_scaff_1 %>% mutate(Mb = Window_Mid / 1000000)
Windows_scaff_1 <- Windows_scaff_1 %>% mutate(cM = Window_RR * 100000000)

# Transform it to be able to plot it together.

p <- ggplot(Windows_scaff_1, aes(x=Mb, y=cM)) + theme_light()
p <- p + geom_line(color="grey") + scale_x_continuous(limits = c(0,11.565203), breaks = seq(0,11.565203, by = 5))
p <- p + scale_y_continuous(limits = c(0,30), breaks = c(0,5,10,15,20,25,30)) + ylab("Rec.Rate (cM/Mb)")
p


# We can generate all the plots first
names(scaff_5_RR)[1] <- "Window_Start"
names(scaff_5_RR)[2] <- "Window_End"
names(scaff_5_RR)[3] <- "Window_RR"

scaff_5_RR <- scaff_5_RR %>% mutate(Window_Mid = Window_Start + 50000)
scaff_5_RR <- scaff_5_RR %>% mutate(Mb = Window_Mid / 1000000)
scaff_5_RR <- scaff_5_RR %>% mutate(cM = Window_RR * 100000000)

plot1 <- ggplot(scaff_1_RR, aes(x=Mb, y=cM)) + theme_light() +
  geom_line(color="grey", size = 1) + scale_x_continuous(limits = c(0,34.455698), breaks = seq(0,35, by = 5)) +
  scale_y_continuous(limits = c(0,15), breaks = c(0,5,10,15,20,25,30)) + ylab("Rec.Rate (cM/Mb)") + ggtitle("Chromosome 1")
plot1


# We can plot all chromosomes together in a grid.

# We load the data.
rr_1Mb <- read.table("features_1Mb.txt", header = TRUE)

# Transform them properly.
rr_1Mb$Chr <- factor(rr_1Mb$Chr, levels = c("Z1", "02", "03", "04", "05", "Z2", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "Z3", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"))
rr_1Mb <- rr_1Mb %>% mutate(Mid = Start + 500000)
rr_1Mb <- rr_1Mb %>% mutate(Mb = Mid / 1000000)
rr_1Mb <- rr_1Mb %>% mutate(cM = Rec_Rate * 100000000)

# Create an object for all plot in one.
facet <- ggplot(rr_1Mb, aes(Mb, cM)) + theme_light() + geom_line(color="black", size = 1) + 
  ylab("Recombination rate (cM/Mb)") + xlab ("Position (Mb)")
facet

# Use vars() to supply the variable we want to separate for (chromosome in this case).
facet_grid <- facet + facet_wrap(vars(Chr), ncol = 5, nrow = 6, scales = "free") +
  theme(strip.background = element_rect(fill="lightgrey")) + scale_x_continuous(limits = c(0,34.455698)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  theme(strip.text = element_text(colour = 'black', size = 12, face="bold"))
facet_wrap()

facet_grid
# We save the final plot.
ggsave("rr-1Mb-all-chromosomes.png", facet_grid, width=10, height=12, units="in", dpi=300)



Recombination_rate_scaff_24 <- read_table("scaff-24-miss0.8.own-hdf.rmap", col_names = FALSE)
names(Recombination_rate_scaff_24)[1] <- "Initial_Position"
names(Recombination_rate_scaff_24)[2] <- "Final_Position"
names(Recombination_rate_scaff_24)[3] <- "Rec_Rate"

data24 = data.frame(Recombination_rate_scaff_24)
data24 <- data24 %>% mutate(cM = Rec_Rate / 0.00000001)
data24 <- data24 %>% mutate(Initial_Position = Initial_Position / 1000000)
data24 <- data24 %>% mutate(Final_Position = Final_Position / 1000000)

Windows_scaff_24 <- read_table("scaffold-24-100kb.txt", col_names = FALSE)
names(Windows_scaff_24)[1] <- "Window_Start"
names(Windows_scaff_24)[2] <- "Window_End"
names(Windows_scaff_24)[3] <- "Window_RR"

Windows_scaff_24 <- Windows_scaff_24 %>% mutate(Window_Mid = Window_Start + 50000)
Windows_scaff_24 <- Windows_scaff_24 %>% mutate(Mb = Window_Mid / 1000000)
Windows_scaff_24 <- Windows_scaff_24 %>% mutate(cM = X4 * 1000000)
Windows_scaff_24 <- Windows_scaff_24 %>% mutate(BG_cM = X5 / 100)
Windows_scaff_24 <- Windows_scaff_24 %>% mutate(BG_cM = BG_cM * 100000000)
Windows_scaff_24 <- Windows_scaff_24 %>% mutate(BG10_cM = BG_cM * 10)

Windows_scaff_24 = data.frame(Windows_scaff_24)





p <- ggplot() +
  geom_line(aes(x=Initial_Position, y=cM), data24, color = "lightgrey", alpha = 1) +
  geom_line(aes(x = Mb, y = BG10_cM), Windows_scaff_24, color= "orange", size = 1) +
  theme_light () + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  scale_x_continuous(expand = c(0,0), limits = c(0,13.684430), breaks = seq(0,13.684430, by = 1)) + 
  scale_y_continuous(expand = c(0.5,0)) +
  xlab("Position (Mb)") + ylab("Rec. Rate (cM/Mb)") +ylim(0,305) +
  geom_segment(aes(x = 0.648605, xend = 0.648605, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 1.977204, xend = 1.977204, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 2.251781, xend = 2.251781, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 2.436298, xend = 2.436298, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 2.458160, xend = 2.458160, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 2.520921, xend = 2.520921, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 2.537270, xend = 2.537270, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 2.784382, xend = 2.784382, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 2.789670, xend = 2.789670, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.029627, xend = 3.029627, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.208062, xend = 3.208062, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.220324, xend = 3.220324, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.233714, xend = 3.233714, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.242459, xend = 3.242459, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.252970, xend = 3.252970, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.416733, xend = 3.416733, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.534869, xend = 3.534869, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.544242, xend = 3.544242, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.591364, xend = 3.591364, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.695547, xend = 3.695547, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 3.802419, xend = 3.802419, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 4.058482, xend = 4.058482, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 4.338929, xend = 4.338929, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 5.670887, xend = 5.670887, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 5.689594, xend = 5.689594, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 5.824567, xend = 5.824567, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 6.213987, xend = 6.213987, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 7.053773, xend = 7.053773, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 7.183370, xend = 7.183370, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 7.764328, xend = 7.764328, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 8.457574, xend = 8.457574, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 8.928313, xend = 8.928313, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 8.932633, xend = 8.932633, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 8.955143, xend = 8.955143, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 8.970670, xend = 8.970670, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 8.982644, xend = 8.982644, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 8.983777, xend = 8.983777, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 9.143393, xend = 9.143393, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 9.152175, xend = 9.152175, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 9.200790, xend = 9.200790, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 9.567110, xend = 9.567110, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.019687, xend = 10.019687, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.098411, xend = 10.098411, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.384596, xend = 10.384596, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.445720, xend = 10.445720, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.488494, xend = 10.488494, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.518053, xend = 10.518053, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.564323, xend = 10.564323, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.859051, xend = 10.859051, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.863866, xend = 10.863866, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.943830, xend = 10.943830, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 10.961220, xend = 10.961220, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 11.199027, xend = 11.199027, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 11.200003, xend = 11.200003, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 11.336718, xend = 11.336718, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 12.374145, xend = 12.374145, y = -Inf, yend = 0.1), size = 0.3, color = "red") +
  geom_segment(aes(x = 13.087962, xend = 13.087962, y = -Inf, yend = 0.1), size = 0.3, color = "red")
p


grid <- plot_grid(facet_grid, p, labels = "AUTO", label_fontface = "bold", label_size = 10, nrow = 2, ncol = 1, rel_heights = c(3,1))
grid
ggsave("Figure1.png", plot = last_plot(), units= c("in"), width = 12, height = 15, dpi = 300)


