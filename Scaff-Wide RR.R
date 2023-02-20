#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(PopGenome)

# Load the data and name the columns.
Recombination_rate_scaff_3 <- read_table2("scaff-17-miss0.8.own-hdf-2.rmap", col_names = FALSE)
names(Recombination_rate_scaff_3)[1] <- "Initial_Position"
names(Recombination_rate_scaff_3)[2] <- "Final_Position"
names(Recombination_rate_scaff_3)[3] <- "Rec_Rate"

# We transform it into a data frame to manage the high volumes of data better.
scaff1_windows = data.frame(Recombination_rate_scaff_3)

# Create weighted value of RR for each interval (between markers).
scaff1_windows$Interval_length<- scaff1_windows$Final_Position - scaff1_windows$Initial_Position
scaff1_windows$Weight<- scaff1_windows$Interval_length / 19406976

# Multiply the RR per the weight, and then do the mean of this value.
scaff1_windows$RR_Weight<- scaff1_windows$Rec_Rate*scaff1_windows$Weight

# Now we just need to do the sum of all the RR_Weight values in the intervals assigned.
means_scaff1 <- sum(scaff1_windows$RR_Weight)
means_scaff1


