# FIRST, we can check for consistency among replicates of the optimize pyrho runs for the same scaffold.

# We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)

# We load the data and name the columns for each replicate (.rmap file)
Recombination_rate_scaff_5 <- read_table2("scaff-5-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_5)[1] <- "Initial_Position"
names(Recombination_rate_scaff_5)[2] <- "Final_Position"
names(Recombination_rate_scaff_5)[3] <- "Rec_Rate"

Recombination_rate_scaff_5_v2 <- read_table2("scaff-5-miss0.97-v2.rmap", col_names = FALSE)
names(Recombination_rate_scaff_5_v2)[1] <- "Initial_Position_2"
names(Recombination_rate_scaff_5_v2)[2] <- "Final_Position_2"
names(Recombination_rate_scaff_5_v2)[3] <- "Rec_Rate_2"

Recombination_rate_scaff_5_v3 <- read_table2("scaff-5-miss0.97-v3.rmap", col_names = FALSE)
names(Recombination_rate_scaff_5_v3)[1] <- "Initial_Position_3"
names(Recombination_rate_scaff_5_v3)[2] <- "Final_Position_3"
names(Recombination_rate_scaff_5_v3)[3] <- "Rec_Rate_3"

# We can transform them into data frames now.
data5 = data.frame(Recombination_rate_scaff_5)
data5_2 = data.frame(Recombination_rate_scaff_5_v2)
data5_3 = data.frame(Recombination_rate_scaff_5_v3)

# We can convert the positions to Mb, by dividing their values by 10^6, and create a new column for such values.
data5 <- data5 %>% mutate(Mb = Initial_Position / 1000000)
data5_2 <- data5_2 %>% mutate(Mb = Initial_Position_2 / 1000000)
data5_3 <- data5_3 %>% mutate(Mb = Initial_Position_3 / 1000000)

data5 <- data5 %>% mutate(cM = Rec_Rate / 0.00000001)
data5_2 <- data5_2 %>% mutate(cM = Rec_Rate_2 / 0.00000001)
data5_3 <- data5_3 %>% mutate(cM = Rec_Rate_3 / 0.00000001)

# Now it's time to plot them individually. We name the plots for each replicate p, q and r. We add different elements to each one so that when we unite them they are aesthetically cohesive 
p <- ggplot(data5, aes(x=Mb, y=cM)) 
p <- p + theme_light () + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,27.223299), breaks = seq(0,27.223299, by = 5))
p <- p + scale_y_continuous(limits = c(0,100), breaks = c(0,10,25,50,100)) + ylab("Rec. rate (cM/Mb)")
p

q <- ggplot(data5_2, aes(x=Mb, y=cM)) 
q <- q + theme_light () + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  geom_line(color="red") + scale_x_continuous(limits = c(0,27.223299), breaks = seq(0,27.223299, by = 5))
q <- q + scale_y_continuous(limits = c(0,100), breaks = c(0,10,25,50,100)) + ylab("Rec. rate (cM/Mb)")
q

r <- ggplot(data5_3, aes(x=Mb, y=cM)) 
r <- r + theme_light () +
  geom_line(color="blue") + scale_x_continuous(limits = c(0,27.223299), breaks = seq(0,27.223299, by = 5))
r <- r + scale_y_continuous(limits = c(0,100), breaks = c(0,10,25,50,100)) + xlab("Position (Mb)") + ylab("Rec. rate (cM/Mb)")
r

# We now unify all plots and represent them in a single figure.
grid.arrange(p, q, r, ncol = 1, nrow = 3)
