#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(PopGenome)
library(data.table)

# Load the data and name the columns.
Recombination_rate_scaff_1 <- read_table2("scaff-29-miss0.8.own-hdf.rmap", col_names = FALSE)
names(Recombination_rate_scaff_1)[1] <- "Initial_Position"
names(Recombination_rate_scaff_1)[2] <- "Final_Position"
names(Recombination_rate_scaff_1)[3] <- "Rec_Rate"

Recombination_rate_scaff_1$Initial_Position <- Recombination_rate_scaff_1$Initial_Position + 1

scaff_1_positions <- setDT(Recombination_rate_scaff_1)[, .(Position = Initial_Position:Final_Position, Rate = Rec_Rate), .(id = 1:nrow(Recombination_rate_scaff_1))]

scaff_1_positions <- subset(scaff_1_positions, select = -c(id))

write.table(scaff_1_positions, file = "scaffold-17-positions.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# Now we create the windows (100 kb long)

# Delete the last 5 digits of the initial position
scaff_1_positions$Window_mid = substr(scaff_1_positions$Position,1,nchar(scaff_1_positions$Position)-5)
# Correct for blank spaces at the beginning, and add 0 to them.
scaff_1_positions$Window_mid <- sub("^$", 0, scaff_1_positions$Window)
# Now multiply all values by 100000.
scaff_1_positions <- transform(scaff_1_positions, Window_mid=as.numeric(as.character(Window_mid)))
scaff_1_positions$Window_mid<- scaff_1_positions$Window_mid*100000
# Add 50000 to all values, to make plotting more accurate.
scaff_1_positions$Window_mid<- scaff_1_positions$Window_mid + 50000
# Now we just need to do the sum of all the RR_Weight values in the intervals assigned.
means_scaff_1_100kb <- scaff_1_positions %>% 
  group_by(Window_mid) %>% 
  mutate(mean = mean(Rate))

#To save the results, we extract some columns, and eliminate all rows with duplicate values (as they belong
# to the same interval and are redundant).

Window_mid <- means_scaff_1_100kb %>% pull(Window_mid)
Window_mid <- unique(Window_mid)  

Window_RR <- means_scaff_1_100kb %>% pull(mean)
Window_RR <- unique(Window_RR) 
Window_RR
Window_cM <- means_scaff_1_100kb %>% pull(cM)
Window_cM <- unique(Window_cM) 

# We recalculate the window start and end positions from the mid point.
Window_Start <- (Window_mid - 50000)
Window_End <- (Window_mid + 50000)

# We create a new data frame for all these extracted information.
scaffold_16_100kb = data.frame(Window_Start, Window_End, Window_RR)

# We output a txt file with the tab-separated columns with the information.
write.table(scaffold_16_100kb, file = "scaffold-17-100kb.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
