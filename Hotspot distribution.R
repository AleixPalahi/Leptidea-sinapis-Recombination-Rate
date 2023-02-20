install.packages("dplyr")
library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)

hotspots <- read_excel("Hotspots_Length.xlsx", sheet = "HOTSPOTS_CORRECT")
Hotspots_filtered <- hotspots %>% filter(Length >= 750) %>% filter(RR >= 25) %>% filter(Length <= 10000)

chr <- read_excel("Scaff-Length-RR.xlsx")
chr$center <- chr$Size / 2

hotspots_1 <- filter(Hotspots_filtered, Scaffold == 1)
hotspots_2 <- filter(Hotspots_filtered, Scaffold == 2)
hotspots_3 <- filter(Hotspots_filtered, Scaffold == 3)
hotspots_4 <- filter(Hotspots_filtered, Scaffold == 4)
hotspots_5 <- filter(Hotspots_filtered, Scaffold == 5)
hotspots_6 <- filter(Hotspots_filtered, Scaffold == 6)
hotspots_7 <- filter(Hotspots_filtered, Scaffold == 7)
hotspots_8 <- filter(Hotspots_filtered, Scaffold == 8)
hotspots_9 <- filter(Hotspots_filtered, Scaffold == 9)
hotspots_10 <- filter(Hotspots_filtered, Scaffold == 10)
hotspots_11 <- filter(Hotspots_filtered, Scaffold == 11)
hotspots_12 <- filter(Hotspots_filtered, Scaffold == 12)
hotspots_13 <- filter(Hotspots_filtered, Scaffold == 13)
hotspots_14 <- filter(Hotspots_filtered, Scaffold == 14)
hotspots_15 <- filter(Hotspots_filtered, Scaffold == 15)
hotspots_16 <- filter(Hotspots_filtered, Scaffold == 16)
hotspots_17 <- filter(Hotspots_filtered, Scaffold == 17)
hotspots_18 <- filter(Hotspots_filtered, Scaffold == 18)
hotspots_19 <- filter(Hotspots_filtered, Scaffold == 19)
hotspots_20 <- filter(Hotspots_filtered, Scaffold == 20)
hotspots_21 <- filter(Hotspots_filtered, Scaffold == 21)
hotspots_22 <- filter(Hotspots_filtered, Scaffold == 22)
hotspots_23 <- filter(Hotspots_filtered, Scaffold == 23)
hotspots_24 <- filter(Hotspots_filtered, Scaffold == 24)
hotspots_25 <- filter(Hotspots_filtered, Scaffold == 25)
hotspots_26 <- filter(Hotspots_filtered, Scaffold == 26)
hotspots_27 <- filter(Hotspots_filtered, Scaffold == 27)
hotspots_28 <- filter(Hotspots_filtered, Scaffold == 28)
hotspots_29 <- filter(Hotspots_filtered, Scaffold == 29)

hotspots_1$Distance_center <- (abs(hotspots_1$Initial_pos - 17227849) / 17227849 /2)
hotspots_2$Distance_center <- (abs(hotspots_2$Initial_pos - 16422112) / 16422112 /2)
hotspots_3$Distance_center <- (abs(hotspots_3$Initial_pos - 15272082) / 15272082 /2)
hotspots_4$Distance_center <- (abs(hotspots_4$Initial_pos - 13711016) / 13711016 /2)
hotspots_5$Distance_center <- (abs(hotspots_5$Initial_pos - 13611650) / 13611650 /2)
hotspots_6$Distance_center <- (abs(hotspots_6$Initial_pos - 13296701) / 13296701 /2)
hotspots_7$Distance_center <- (abs(hotspots_7$Initial_pos - 12709236) / 12709236 /2)
hotspots_8$Distance_center <- (abs(hotspots_8$Initial_pos - 12635893) / 12635893 /2)
hotspots_9$Distance_center <- (abs(hotspots_9$Initial_pos - 12410781) / 12410781 /2)
hotspots_10$Distance_center <- (abs(hotspots_10$Initial_pos - 11848176) / 11848176 /2)
hotspots_11$Distance_center <- (abs(hotspots_11$Initial_pos - 11299846) / 11299846 /2)
hotspots_12$Distance_center <- (abs(hotspots_12$Initial_pos - 11154678) / 11154678 /2)
hotspots_13$Distance_center <- (abs(hotspots_13$Initial_pos - 10918758) / 10918758 /2)
hotspots_14$Distance_center <- (abs(hotspots_14$Initial_pos - 10501169) / 10501169 /2)
hotspots_15$Distance_center <- (abs(hotspots_15$Initial_pos - 10385810) / 10385810 /2)
hotspots_16$Distance_center <- (abs(hotspots_16$Initial_pos - 10177264) / 10177264 /2)
hotspots_17$Distance_center <- (abs(hotspots_17$Initial_pos - 9703488) / 9703488 /2)
hotspots_18$Distance_center <- (abs(hotspots_18$Initial_pos - 9257705) / 9257705 /2)
hotspots_19$Distance_center <- (abs(hotspots_19$Initial_pos - 8886558) / 8886558 /2)
hotspots_20$Distance_center <- (abs(hotspots_20$Initial_pos - 8153522) / 8153522 /2)
hotspots_21$Distance_center <- (abs(hotspots_21$Initial_pos - 8122023) / 8122023 /2)
hotspots_22$Distance_center <- (abs(hotspots_22$Initial_pos - 7892680) / 7892680 /2)
hotspots_23$Distance_center <- (abs(hotspots_23$Initial_pos - 7412015) / 7412015 /2)
hotspots_24$Distance_center <- (abs(hotspots_24$Initial_pos - 6842215) / 6842215 /2)
hotspots_25$Distance_center <- (abs(hotspots_25$Initial_pos - 6735720) / 6735720 /2)
hotspots_26$Distance_center <- (abs(hotspots_26$Initial_pos - 6268922) / 6268922 /2)
hotspots_27$Distance_center <- (abs(hotspots_27$Initial_pos - 6188244) / 6188244 /2)
hotspots_28$Distance_center <- (abs(hotspots_28$Initial_pos - 5829164) / 5829164 /2)
hotspots_29$Distance_center <- (abs(hotspots_29$Initial_pos - 5582602) / 5582602 /2)

hotspots_bra <- rbind(hotspots_1, hotspots_2, hotspots_3, hotspots_4, hotspots_5, hotspots_6, hotspots_7, hotspots_8, hotspots_9, hotspots_10, hotspots_11, hotspots_12, hotspots_13, hotspots_14, hotspots_15, hotspots_16, hotspots_17, hotspots_18, hotspots_19, hotspots_20, hotspots_21, hotspots_22, hotspots_23, hotspots_24, hotspots_25, hotspots_26, hotspots_27, hotspots_28, hotspots_29)

hotspots_bra$percentil <- substr(hotspots_bra$Distance_center, 3, 4)
hotspots_bra$percentil <- as.numeric(hotspots_bra$percentil)
hotspots_bra$percentil <- hotspots_bra$percentil + 1
summary_hot <- hotspots_bra %>% count(percentil)
summary_hot

summary_hot$percentil <- as.numeric(summary_hot$percentil)
summary_hot$n <- as.numeric(summary_hot$n)
cor.test(summary_hot$percentil, summary_hot$n, alternative = "two.sided")

p <- ggplot(summary_hot, aes(x=percentil, y=n)) + theme_light() + 
  geom_area(color = "darkorange", fill = "orange", alpha = 0.2) + 
  geom_vline(xintercept = 45, linetype = "dashed", size = 1) +
  ylab("Number of Hotspots") + xlab("Relative position (Percentile)") + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
p

central <- data.frame(rep(c("C", "T"), times= c(45,5)))
type <- cbind(summary, central)
colnames(type)[3] <- "Type"

wilcox.test(n ~ Type, data = type, exact = FALSE) # W = 220.5, p = 4.94e-4
