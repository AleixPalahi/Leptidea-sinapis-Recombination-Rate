# This script is intended to plot the recombination rate with SNP density, to see if the 
# regions of invariable RR are due to a lack of resolution or are indeed regions of constant RR

#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)

#Here is the code needed to load the data and name the different columns of the table.
Recombination_rate_scaff_29 <- read_table2("scaff-29-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_29)[1] <- "Initial_Position"
names(Recombination_rate_scaff_29)[2] <- "Final_Position"
names(Recombination_rate_scaff_29)[3] <- "Rec_Rate"

SNP_Density_scaff_29 <- read_table2("scaff-29-100kb.snpden")
names(SNP_Density_scaff_29)[2] <- "Position"
names(SNP_Density_scaff_29)[4] <- "Variants"


# We can transform them into data frames now.
data29 = data.frame(Recombination_rate_scaff_29)
snp_den_29 = data.frame(SNP_Density_scaff_29)

# We can convert the positions to Mb, by dividing their values by 10^6, and create a new column for such values.
data29 <- data29 %>% mutate(Mb = Initial_Position / 1000000)
data29 <- data29 %>% mutate(cM = Rec_Rate / 0.00000001)

snp_den_29 <- snp_den_29 %>% mutate(Position = Position / 1000000)

# We plot the results.
p <- ggplot(data29, aes(x=Mb, y=cM)) + theme_light ()
p <- p + ggtitle("RR and Marker Density - Scaffold 29") + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,11.165203), breaks = seq(0,11.165203, by = 5))
p <- p + scale_y_continuous(limits = c(0,100), breaks = c(0,10,25,50,100)) + ylab("Rec.Rate (cM/Mb)")
p

r <- ggplot(snp_den_29, aes(x=Position, y=Variants)) 
r <- r + theme_light () +
  geom_line(color="blue") + scale_x_continuous(limits = c(0,11.165203), breaks = seq(0,11.165203, by = 5))
r <- r + scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10)) + xlab("Position (Mb)") + ylab("Variant Density (SNPs/kb)")
r

# We now unify all plots and represent them in a single figure.
grid.arrange(p, r, ncol = 1, nrow = 2)








rr_1Mb <- read_table("rr-1Mb-corrected.txt", col_names = FALSE)
density_1Mb <- read_table("SNPden-1Mb-corrected.txt", col_names = FALSE)

rr_density_1Mb <- cbind(rr_1Mb$X1, rr_1Mb$X2, rr_1Mb$X3, rr_1Mb$X4, density_1Mb$X3)
rr_density_1Mb <- data.frame(rr_density_1Mb)

names(rr_density_1Mb)[1] <- "Chromosome"
names(rr_density_1Mb)[2] <- "Initial_Position"
names(rr_density_1Mb)[3] <- "Final_Position"
names(rr_density_1Mb)[4] <- "RR"
names(rr_density_1Mb)[5] <- "Markers"

rr_density_1Mb <- rr_density_1Mb %>% mutate(Length = Final_Position - Initial_Position)
rr_density_1Mb <- rr_density_1Mb %>% mutate(cM = RR * 100000000)
rr_density_1Mb <- rr_density_1Mb %>% mutate(Density = Markers / Length)

rr_density_1Mb_filtered <- rr_density_1Mb %>% filter(cM < 30)

rr_density_no_begin <- rr_density_1Mb %>% filter(Initial_Position > 500000)
rr_density_no_subtelomeric <- rr_density_no_begin[-c(34, 66, 96, 123,150, 176, 201, 226, 250, 273, 295, 317, 338, 358, 378, 398, 417, 435, 452, 468, 484, 499, 513, 526, 539, 551, 563, 574, 585), ] 
rr_density_no_subtelomeric <- rr_density_no_subtelomeric %>% filter(Density > 0.001)
rr_density_no_subtelomeric <- rr_density_no_subtelomeric %>% filter(cM < 40)

spearman <- cor.test(rr_density_no_subtelomeric$Density, rr_density_no_subtelomeric$cM, method = "spearman", alternative = "two.sided", exact = FALSE)
spearman

p <- ggplot(rr_density_no_subtelomeric, aes(x=cM, y=Density)) + theme_light() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  geom_point(shape=20, size=1.5, color="grey", alpha=0.9) + 
  geom_smooth(method=lm, se=TRUE, linetype="dashed", color="black") +
  geom_label(label="rho = 0.121, p = 0.004", x = 29.5, y = 0.026, size = 4) +
  xlab ("Recombination rate (cM/Mb)") + ylab ("Number of markers (proportion of window)") +
  scale_x_continuous(expand = c(0.01,0.01))
p

png("RR-marker-density-filtered.png", units="in", width=6, height=4, res=300)
p
dev.off()

ggsave("RR-marker-density-filtered.png", p, width=6, height=4, units="in", dpi=300)



