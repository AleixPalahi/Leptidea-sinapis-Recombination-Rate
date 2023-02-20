# The purpose of this script is to evaluate whether the reduction in recombination rate in the more distal positions of the 
# chromosome is in fact statistically significant.
install.packages("stringr")
#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(stringr)

# We load the data.
aleix <- function (name,file) {

name <- read_table2(file, col_names = FALSE)
#names(name)[1] <- "Window_Start"
#names(name)[2] <- "Window_End"
#names(name)[3] <- "Window_RR"

# We can extract each the rows for the first and last 5 lines of each chromosome, and then bind them into a single table.
chr_L <- head(name, 5)
chr_R <- tail(name, 5)
chr_1_ends <- rbind(chr_L, chr_R)

# We bind the 10 more distal windows from each chromosome into a single table.
#chr_ends <- rbind(chr_1_ends, chr_2_ends, chr_3_ends, chr_4_ends, chr_5_ends, chr_6_ends, chr_7_ends, chr_8_ends, chr_9_ends, chr_10_ends, chr_11_ends, chr_12_ends, chr_13_ends, chr_14_ends, chr_15_ends, chr_16_ends, chr_17_ends, chr_18_ends, chr_19_ends, chr_20_ends, chr_21_ends, chr_22_ends, chr_23_ends, chr_24_ends, chr_25_ends, chr_26_ends, chr_27_ends, chr_28_ends, chr_29_ends)

# We now have to separate all the windows that are not in the distal positions.
chr_mid <-  slice(name, 6:(n() - 5))
}
# And bind the proximal windows of all chromosomes.
#chr_mid <- rbind(chr_1_mid, chr_2_mid, chr_3_mid, chr_4_mid, chr_5_mid, chr_6_mid, chr_7_mid, chr_8_mid, chr_9_mid, chr_10_mid, chr_11_mid, chr_12_mid, chr_13_mid, chr_14_mid, chr_15_mid, chr_16_mid, chr_17_mid, chr_18_mid, chr_19_mid, chr_20_mid, chr_21_mid, chr_22_mid, chr_23_mid, chr_24_mid, chr_25_mid, chr_26_mid, chr_27_mid, chr_28_mid, chr_29_mid)

# In order to be able to distinguish the types of windows, we can add a column to the different tables with either "End" 
# or "Mid". 

#chr_ends$Type <- c("End")
#chr_mid$Type <- c("Mid")

# And now we bind everything into a single table.
#chr_mid_ends <- rbind(chr_mid, chr_ends)



# We test for normality of the variables.
#with(chr_mid_ends, shapiro.test(cM[Type == "Mid"])) 
#with(chr_mid_ends, shapiro.test(cM[Type == "End"])) # Normality is not observed in the RR distributions.

# We test for homogeneity of variances.
#Variances_Homogeneity <- var.test(cM ~ Type, data = chr_mid_ends) # p-value = 2.847e-13
#Variances_Homogeneity

# We proceed with the Wilcoxon test to check whether the differences are significant.
wilcox_L <- wilcox.test(chr_mid$Window_RR, chr_L$Window_RR, alternative="greater")
wilcox_L # Significant differences (p = 8.630715e-61).
wilcox_R <- wilcox.test(chr_mid$Window_RR, chr_R$Window_RR, alternative="greater")
wilcox_R # Significant differences (p = 8.630715e-61).

wilcox$p.value
return (c(wilcox_L$p.value, wilcox_R$p.value))


#The mean RR in the chromosome ends can also be obtained.
mean(chr_ends$cM) # Mean RR in ends = 2.46 cM/Mb.




chr23 <- read_table2("scaffold-23-100kb-corrected.txt", col_names = FALSE)
names(chr23)[1] <- "Window_Start"
names(chr23)[2] <- "Window_End"
names(chr23)[3] <- "Window_RR"
chr23$Window_RR <- as.numeric(chr23$Window_RR)

# We can extract each the rows for the first and last 5 lines of each chromosome, and then bind them into a single table.
chr_L <- head(chr23, 5)
chr_R <- tail(chr23, 5)
#chr_1_ends <- rbind(chr_L, chr_R)

# We bind the 10 more distal windows from each chromosome into a single table.
#chr_ends <- rbind(chr_1_ends, chr_2_ends, chr_3_ends, chr_4_ends, chr_5_ends, chr_6_ends, chr_7_ends, chr_8_ends, chr_9_ends, chr_10_ends, chr_11_ends, chr_12_ends, chr_13_ends, chr_14_ends, chr_15_ends, chr_16_ends, chr_17_ends, chr_18_ends, chr_19_ends, chr_20_ends, chr_21_ends, chr_22_ends, chr_23_ends, chr_24_ends, chr_25_ends, chr_26_ends, chr_27_ends, chr_28_ends, chr_29_ends)

# We now have to separate all the windows that are not in the distal positions.
chr_mid <-  slice(chr23, 6:(n() - 5))

wilcox_L <- wilcox.test(chr_mid$Window_RR, chr_L$Window_RR, alternative="greater")
wilcox_L # Significant differences (p = 8.630715e-61).
wilcox_R <- wilcox.test(chr_mid$Window_RR, chr_R$Window_RR, alternative="greater")
wilcox_R # Significant differences (p = 8.630715e-61).

wilcox$p.value





chr_1 <- read_table2("scaffold-1-100kb.txt", col_names = FALSE)

# We can extract each the rows for the first and last 5 lines of each chromosome, and then bind them into a single table.
chr_L <- head(chr_1, 5)
chr_R <- tail(chr_1, 5)
chr_1_ends <- rbind(chr_L, chr_R)
chr_ends <- rbind(chr_1_ends, chr_2_ends, chr_3_ends, chr_4_ends, chr_5_ends, chr_6_ends, chr_7_ends, chr_8_ends, chr_9_ends, chr_10_ends, chr_11_ends, chr_12_ends, chr_13_ends, chr_14_ends, chr_15_ends, chr_16_ends, chr_17_ends, chr_18_ends, chr_19_ends, chr_20_ends, chr_21_ends, chr_22_ends, chr_23_ends, chr_24_ends, chr_25_ends, chr_26_ends, chr_27_ends, chr_28_ends, chr_29_ends)

names(name)[1] <- "Window_Start"
names(name)[2] <- "Window_End"
names(name)[4] <- "Window_RR"

# We now have to separate all the windows that are not in the distal positions.
chr_1_mid <-  slice(chr1, 6:(n() - 5))

bu()