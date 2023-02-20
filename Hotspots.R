#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(PopGenome)
library("readxl")

Recombination_rate_scaff_4 <- read_table("scaff-17-miss0.8.own-hdf.rmap", col_names = FALSE)

names(Recombination_rate_scaff_4)[1] <- "Initial_Position"
names(Recombination_rate_scaff_4)[2] <- "Final_Position"
names(Recombination_rate_scaff_4)[3] <- "Rec_Rate"

Recombination_rate_scaff_4 <- Recombination_rate_scaff_4 %>% mutate(cM = Rec_Rate / 0.00000001)
Recombination_rate_scaff_4 <- Recombination_rate_scaff_4 %>% mutate(Length = Final_Position - Initial_Position)

Hotspots_4 <- Recombination_rate_scaff_4 %>% filter(Initial_Position >= 5000000) %>% filter(Initial_Position <= 5100000) %>% filter(cM >= 210)
Hotspots_4 <- Hotspots_4 %>% filter(cM <= 1000)
Hotspots_4 <- Hotspots_4 %>% filter(Initial_Position >= 5000000) %>% filter(Initial_Position <= 5100000) %>% filter(Length >= 750) %>% filter(Length <= 10000)



# We load the data of all the hotspots detected.
hotspots <- read_excel("Hotspots_Length.xlsx", sheet = "HOTSPOTS_CORRECT")

# We can plot the potential association between the hotspot length and the RR, to see if maybe the shortest hotspots have significantly high or low RR.
ggplot(hotspots, aes(x=Length, y=RR)) + theme_light() + geom_point(shape=20, size=2.5) + xlab ("Length (bp)") + ylab ("Recombination rate (cM/Mb)")

# We see that some hotspots have abnormally high RR, and they are all of very short length. We can set the threshold at 750 bp then.

Hotspots_filtered <- hotspots %>% filter(Length >= 750) %>% filter(RR >= 25) %>% filter(Length <= 10000) # 3124 hotspots remain.

ggplot(Hotspots_filtered, aes(x=Length, y=RR)) + theme_light() + geom_point(shape=20, size=2.5) + xlab ("Length (bp)") + ylab ("Recombination rate (cM/Mb)")

# There seems to be a tendency for longer hotspots to have lower RR. We can investigate this.
shapiro.test(Hotspots_filtered$Length)
shapiro.test(Hotspots_filtered$RR) # Both variables are not normal, so we must use non-parametric correlations.
cor.test(Hotspots_filtered$Length, Hotspots_filtered$RR, method = 'spearman', alternative = "less", exact = FALSE)
cor.test(Hotspots_filtered$Length, Hotspots_filtered$RR, method = 'kendall', alternative = "less", exact = FALSE)
# Both non-parametric methods show a significant negative correlation between hotspots length and RR. This is confirmed by the fact that most outliers for RR are gone when we apply the 750 bp threshold.

mean(Hotspots_filtered$Length) # 1656 bp long
sd(Hotspots_filtered$Length) # 1193 bp
median(Hotspots_filtered$Length) # 1242 bp long
quantile(Hotspots_filtered$Length)

mean(Hotspots_filtered$RR) # 94.139 cM/Mb
sd(Hotspots_filtered$RR) # 62.5 cM/Mb
median(Hotspots_filtered$RR) # 77.121 cM/Mb
quantile(Hotspots_filtered$RR)
more_than_400 <- Hotspots_filtered %>% filter(RR >= 400)
max(Hotspots_filtered$RR)
sum(Hotspots_filtered$Length)
ggplot(Hotspots_filtered, aes(x=RR)) + geom_histogram(binwidth = 25) # It is clearly visible how hotspots tend to be of moderate RR, while just some exceed 200 cM/Mb


# Regarding the distribution of hotspots, we can check whether they accumulate in the A or Z. We start by loading the data.
hotspots_per_scaffold <- read_excel("Hotspots_Length.xlsx", sheet = "Hotspots_per_scaffold")

#Now we transform the length to Mb
hotspots_per_scaffold <- hotspots_per_scaffold %>% mutate(Mb = Length / 1000000)
hotspots_per_scaffold <- hotspots_per_scaffold %>% mutate(Hotspot_density = CORRECT_filtered / Mb)


with(hotspots_per_scaffold, shapiro.test(Hotspot_density[Type == "A"])) # p = 0.172
with(hotspots_per_scaffold, shapiro.test(Hotspot_density[Type == "Z"])) # p = 0.120
Variances_Homogeneity <- var.test(Hotspot_density ~ Type, data = hotspots_per_scaffold)
Variances_Homogeneity # p = 0.882 -> Both groups are normally distributed, and the variances are homogeneous.
t_test <- t.test(Hotspot_density ~ Type, data = hotspots_per_scaffold, var.equal = TRUE, alternative = "less")
t_test # Higher frequency in Z chromosomes, but it is non-significant (p = 0.199)



write.table(Hotspots_filtered, file = "hotspots_filtered.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)


# How much of the recombination takes place within the defined hotspots?
Hotspots_filtered <- read_table("hotspots_filtered.txt", col_names = FALSE)

names(Hotspots_filtered)[1] <- "Scaffold"
names(Hotspots_filtered)[2] <- "Initial_Position"
names(Hotspots_filtered)[3] <- "Final_Position"
names(Hotspots_filtered)[4] <- "Length"
names(Hotspots_filtered)[5] <- "RR"

Hotspots_filtered$Mb <- Hotspots_filtered$Length / 1000000
Hotspots_filtered$Weight <- Hotspots_filtered$RR * Hotspots_filtered$Mb
sum(Hotspots_filtered$Weight) # 471.71 cM in the hotspots, out of 4422 cM in total
sum(Hotspots_filtered$Weight) / 4422 # 10.67% of recombination events within the hot-spots.








# Check weird peak in chr 17.
Recombination_rate_scaff_17 <- read_table("scaff-17-miss0.8.own-hdf.rmap", col_names = FALSE)

names(Recombination_rate_scaff_17)[1] <- "Initial_Position"
names(Recombination_rate_scaff_17)[2] <- "Final_Position"
names(Recombination_rate_scaff_17)[3] <- "Rec_Rate"

Recombination_rate_scaff_17$cM <- Recombination_rate_scaff_17$Rec_Rate *100000000
Recombination_rate_scaff_17$Length <- Recombination_rate_scaff_17$Final_Position - Recombination_rate_scaff_17$Initial_Position

Scaff_17 <- Recombination_rate_scaff_17 %>% filter(11000000 < Initial_Position) %>% filter(Initial_Position < 12000000)

library(PopGenome)
library(data.table)

Scaff_17$Initial_Position <- Scaff_17$Initial_Position + 1
scaff_17_positions <- setDT(Scaff_17)[, .(Position = Initial_Position:Final_Position, Rate = Rec_Rate), .(id = 1:nrow(Scaff_17))]
scaff_17_positions <- subset(scaff_17_positions, select = -c(id))

write.table(scaff_17_positions, file = "scaffold-17-positions-11-12Mb.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# Delete the last 5 digits of the initial position
scaff_17_positions$Window_mid = substr(scaff_17_positions$Position,1,nchar(scaff_17_positions$Position)-6)
# Correct for blank spaces at the beginning, and add 0 to them.
scaff_17_positions$Window_mid <- sub("^$", 0, scaff_17_positions$Window)
# Now multiply all values by 100000.
scaff_17_positions <- transform(scaff_17_positions, Window_mid=as.numeric(as.character(Window_mid)))
scaff_17_positions$Window_mid<- scaff_17_positions$Window_mid*1000000
# Add 50000 to all values, to make plotting more accurate.
scaff_17_positions$Window_mid<- scaff_17_positions$Window_mid + 500000
# Now we just need to do the sum of all the RR_Weight values in the intervals assigned.
means_scaff_1_1Mb <- scaff_17_positions %>% 
  group_by(Window_mid) %>% 
  mutate(mean = mean(Rate))




ggplot(Recombination_rate_scaff_17, aes(x = Initial_Position, y = cM)) + geom_line()



Scaff_17 <- Scaff_17 %>% filter(10000000 < Initial_Position) %>% filter(Initial_Position < 12000000) %>% filter(Length > 500) %>% filter(cM > 73)



