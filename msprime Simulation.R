#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)

#Here is the code needed to load the data and name the different columns of the table.
simulation <- read_table("sim-long-trial.rmap", col_names = FALSE)
names(simulation)[1] <- "Initial_Position"
names(simulation)[2] <- "Final_Position"
names(simulation)[3] <- "Rec_Rate"

#We transform it into a data frame to manage the high volumes of data better.
simulation.df = data.frame(simulation)

# We can convert the positions to Mb, by dividing their values by 10^6, and create a new column for such values.
simulation.df <- simulation.df %>% mutate(Mb = Initial_Position / 1000000)
simulation.df <- simulation.df %>% mutate(cM = Rec_Rate / 0.00000001)

hotspots <- simulation.df %>% filter(cM < 30) %>% filter(cM > 20)
# aside from unions between different independent simulation runs, the maximum RR is 20 cM/Mb (3-fold increase over GwRR).

mean(simulation.df$cM)
max(simulation.df$cM)

#Now we plot the data.
p <- ggplot(simulation.df, aes(x=Mb, y=cM)) +
  theme_light () + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +
  geom_line(color="grey", lwd = 0.2) + scale_x_continuous(limits = c(0,10), breaks = seq(0,10, by = 1)) +
  scale_y_continuous(limits = c(0,21), breaks = c(0,5,10,15,20), expand = c(0.02,0.02)) + scale_x_continuous(expand = c(0.03,0.03)) + xlab("Position (Mb)") + ylab("Recombination Rate (cM/Mb)")
p

#We generate the plot at high resolution.
png("msprime-simulation.png", units="in", width=6, height=4, res=300)
p
dev.off()
