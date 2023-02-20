install.packages("dplyr")
library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(vegan)
library(cowplot)

coldspots <- read_excel("Coldspots_Length.xlsx", sheet = "Coldspots_CORRECT")

chr <- read_excel("Scaff-Length-RR.xlsx")
chr$center <- chr$Size / 2

coldspots_1 <- filter(coldspots, Scaffold == "HiC_scaffold_1")
coldspots_2 <- filter(coldspots, Scaffold == "HiC_scaffold_2")
coldspots_3 <- filter(coldspots, Scaffold == "HiC_scaffold_3")
coldspots_4 <- filter(coldspots, Scaffold == "HiC_scaffold_4")
coldspots_5 <- filter(coldspots, Scaffold == "HiC_scaffold_5")
coldspots_6 <- filter(coldspots, Scaffold == "HiC_scaffold_6")
coldspots_7 <- filter(coldspots, Scaffold == "HiC_scaffold_7")
coldspots_8 <- filter(coldspots, Scaffold == "HiC_scaffold_8")
coldspots_9 <- filter(coldspots, Scaffold == "HiC_scaffold_9")
coldspots_10 <- filter(coldspots, Scaffold == "HiC_scaffold_10")
coldspots_11 <- filter(coldspots, Scaffold == "HiC_scaffold_11")
coldspots_12 <- filter(coldspots, Scaffold == "HiC_scaffold_12")
coldspots_13 <- filter(coldspots, Scaffold == "HiC_scaffold_13")
coldspots_14 <- filter(coldspots, Scaffold == "HiC_scaffold_14")
coldspots_15 <- filter(coldspots, Scaffold == "HiC_scaffold_15")
coldspots_16 <- filter(coldspots, Scaffold == "HiC_scaffold_16")
coldspots_17 <- filter(coldspots, Scaffold == "HiC_scaffold_17")
coldspots_18 <- filter(coldspots, Scaffold == "HiC_scaffold_18")
coldspots_19 <- filter(coldspots, Scaffold == "HiC_scaffold_19")
coldspots_20 <- filter(coldspots, Scaffold == "HiC_scaffold_20")
coldspots_21 <- filter(coldspots, Scaffold == "HiC_scaffold_21")
coldspots_22 <- filter(coldspots, Scaffold == "HiC_scaffold_22")
coldspots_23 <- filter(coldspots, Scaffold == "HiC_scaffold_23")
coldspots_24 <- filter(coldspots, Scaffold == "HiC_scaffold_24")
coldspots_25 <- filter(coldspots, Scaffold == "HiC_scaffold_25")
coldspots_26 <- filter(coldspots, Scaffold == "HiC_scaffold_26")
coldspots_27 <- filter(coldspots, Scaffold == "HiC_scaffold_27")
coldspots_28 <- filter(coldspots, Scaffold == "HiC_scaffold_28")
coldspots_29 <- filter(coldspots, Scaffold == "HiC_scaffold_29")

coldspots_1$Distance_center <- (abs(coldspots_1$Initial_pos - 17227849) / 17227849 /2)
coldspots_2$Distance_center <- (abs(coldspots_2$Initial_pos - 16422112) / 16422112 /2)
coldspots_3$Distance_center <- (abs(coldspots_3$Initial_pos - 15272082) / 15272082 /2)
coldspots_4$Distance_center <- (abs(coldspots_4$Initial_pos - 13711016) / 13711016 /2)
coldspots_5$Distance_center <- (abs(coldspots_5$Initial_pos - 13611650) / 13611650 /2)
coldspots_6$Distance_center <- (abs(coldspots_6$Initial_pos - 13296701) / 13296701 /2)
coldspots_7$Distance_center <- (abs(coldspots_7$Initial_pos - 12709236) / 12709236 /2)
coldspots_8$Distance_center <- (abs(coldspots_8$Initial_pos - 12635893) / 12635893 /2)
coldspots_9$Distance_center <- (abs(coldspots_9$Initial_pos - 12410781) / 12410781 /2)
coldspots_10$Distance_center <- (abs(coldspots_10$Initial_pos - 11848176) / 11848176 /2)
coldspots_11$Distance_center <- (abs(coldspots_11$Initial_pos - 11299846) / 11299846 /2)
coldspots_12$Distance_center <- (abs(coldspots_12$Initial_pos - 11154678) / 11154678 /2)
coldspots_13$Distance_center <- (abs(coldspots_13$Initial_pos - 10918758) / 10918758 /2)
coldspots_14$Distance_center <- (abs(coldspots_14$Initial_pos - 10501169) / 10501169 /2)
coldspots_15$Distance_center <- (abs(coldspots_15$Initial_pos - 10385810) / 10385810 /2)
coldspots_16$Distance_center <- (abs(coldspots_16$Initial_pos - 10177264) / 10177264 /2)
coldspots_17$Distance_center <- (abs(coldspots_17$Initial_pos - 9703488) / 9703488 /2)
coldspots_18$Distance_center <- (abs(coldspots_18$Initial_pos - 9257705) / 9257705 /2)
coldspots_19$Distance_center <- (abs(coldspots_19$Initial_pos - 8886558) / 8886558 /2)
coldspots_20$Distance_center <- (abs(coldspots_20$Initial_pos - 8153522) / 8153522 /2)
coldspots_21$Distance_center <- (abs(coldspots_21$Initial_pos - 8122023) / 8122023 /2)
coldspots_22$Distance_center <- (abs(coldspots_22$Initial_pos - 7892680) / 7892680 /2)
coldspots_23$Distance_center <- (abs(coldspots_23$Initial_pos - 7412015) / 7412015 /2)
coldspots_24$Distance_center <- (abs(coldspots_24$Initial_pos - 6842215) / 6842215 /2)
coldspots_25$Distance_center <- (abs(coldspots_25$Initial_pos - 6735720) / 6735720 /2)
coldspots_26$Distance_center <- (abs(coldspots_26$Initial_pos - 6268922) / 6268922 /2)
coldspots_27$Distance_center <- (abs(coldspots_27$Initial_pos - 6188244) / 6188244 /2)
coldspots_28$Distance_center <- (abs(coldspots_28$Initial_pos - 5829164) / 5829164 /2)
coldspots_29$Distance_center <- (abs(coldspots_29$Initial_pos - 5582602) / 5582602 /2)

coldspots_bra <- rbind(coldspots_1, coldspots_2, coldspots_3, coldspots_4, coldspots_5, coldspots_6, coldspots_7, coldspots_8, coldspots_9, coldspots_10, coldspots_11, coldspots_12, coldspots_13, coldspots_14, coldspots_15, coldspots_16, coldspots_17, coldspots_18, coldspots_19, coldspots_20, coldspots_21, coldspots_22, coldspots_23, coldspots_24, coldspots_25, coldspots_26, coldspots_27, coldspots_28, coldspots_29)

coldspots_bra$percentil <- substr(coldspots_bra$Distance_center, 3, 4)
coldspots_bra$percentil <- as.numeric(coldspots_bra$percentil)
coldspots_bra$percentil <- coldspots_bra$percentil + 1
summary_cold <- coldspots_bra %>% count(percentil)
summary_cold

summary_cold$percentil <- as.numeric(summary_cold$percentil)
summary_cold$n <- as.numeric(summary_cold$n)
cor.test(summary_cold$percentil, summary_cold$n, alternative = "two.sided")

q <- ggplot(summary_cold, aes(x=percentil, y=n)) + theme_light() + 
  geom_area(color = "blue", fill = "lightblue", alpha = 0.2) +
  geom_vline(xintercept = 45, linetype = "dashed", size = 1) +
  ylab("Number of Coldspots") + xlab("Relative position (Percentile)") + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
q

central <- data.frame(rep(c("C", "T"), times= c(36,5)))
type <- cbind(summary, central)
colnames(type)[3] <- "Type"

wilcox.test(n ~ Type, data = type, exact = FALSE) # W = 3, p = 5.62e-4





###### Figure for Coldspot and Hotspot distribution (must run both scripts) ########

grid <- plot_grid(p, q, labels = "AUTO", label_size = 10, nrow = 2)
grid
ggsave("coldspot-hotspot-percentage.png", plot = last_plot(), units= c("in"), width = 10, height = 8, dpi = 300)



