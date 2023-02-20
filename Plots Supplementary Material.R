# We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)

# We load the data and name the columns for each replicate (.rmap file)
Recombination_rate_scaff_29 <- read_table2("scaffold-1-100kb.txt", col_names = FALSE)
names(Recombination_rate_scaff_29)[1] <- "Initial_Position"
names(Recombination_rate_scaff_29)[2] <- "Final_Position"
names(Recombination_rate_scaff_29)[3] <- "Rec_Rate"

# We can transform them into data frames now.
data29 = data.frame(Recombination_rate_scaff_29)

# We can convert the positions to Mb, by dividing their values by 10^6, and create a new column for such values.
data29 <- data29 %>% mutate(Mb = Initial_Position / 1000000)
data29 <- data29 %>% mutate(cM = Rec_Rate / 0.00000001)
data29 <- data29 %>% mutate(Mb = Mb + 0.05)


# Now it's time to plot them individually. We name the plots for each replicate p, q and r. We add different elements to each one so that when we unite them they are aesthetically cohesive 
p <- ggplot(data29, aes(x=Mb, y=cM)) 
p <- p + theme_light () +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,34.555698), breaks = seq(0,34.555698, by = 5))
p <- p + scale_y_continuous(limits = c(0,20), breaks = c(0,5,10,15,20,25)) + ylab("Rec. rate (cM/Mb)") + xlab("Position(Mb)")
p

