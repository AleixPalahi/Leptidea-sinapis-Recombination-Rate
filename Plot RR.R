#We load the packages necessary to plot the results.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)

#Here is the code needed to load the data and name the different columns of the table.
Recombination_rate_scaff_3 <- read_table2("sim-long-trial.rmap", col_names = FALSE)
names(Recombination_rate_scaff_3)[1] <- "Initial_Position"
names(Recombination_rate_scaff_3)[2] <- "Final_Position"
names(Recombination_rate_scaff_3)[3] <- "Rec_Rate"

#We transform it into a data frame to manage the high volumes of data better.
data1 = data.frame(Recombination_rate_scaff_3)

# We can convert the positions to Mb, by dividing their values by 10^6, and create a new column for such values.
data1 <- data1 %>% mutate(Mb = Initial_Position / 1000000)
data1 <- data1 %>% mutate(cM = Rec_Rate / 0.00000001)

#Now we plot the data.
p <- ggplot(data1, aes(x=Mb, y=cM)) 
p <- p + theme_light () + ggtitle("msprime simulation constant RR") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,10), breaks = seq(0,10, by = 1))
p <- p + scale_y_continuous(limits = c(0,7), breaks = c(0,2,4,6)) + xlab("Position (kb)") + ylab("Recombination Rate (cM/Mb)")
p


Recombination_rate_scaff_2 <- read_table2("scaff-2-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_2)[1] <- "Initial_Position"
names(Recombination_rate_scaff_2)[2] <- "Final_Position"
names(Recombination_rate_scaff_2)[3] <- "Rec_Rate"

data2 = data.frame(Recombination_rate_scaff_2)

p <- ggplot(data2, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 2") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,32844224), breaks = seq(0,32844224, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_3 <- read_table2("scaff-3-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_3)[1] <- "Initial_Position"
names(Recombination_rate_scaff_3)[2] <- "Final_Position"
names(Recombination_rate_scaff_3)[3] <- "Rec_Rate"

data3 = data.frame(Recombination_rate_scaff_3)

p <- ggplot(data3, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 3") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,30544165), breaks = seq(0,30544165, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_4 <- read_table2("scaff-4-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_4)[1] <- "Initial_Position"
names(Recombination_rate_scaff_4)[2] <- "Final_Position"
names(Recombination_rate_scaff_4)[3] <- "Rec_Rate"

data4 = data.frame(Recombination_rate_scaff_4)

p <- ggplot(data4, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 4") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,27422032), breaks = seq(0,27422032, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_5 <- read_table2("scaff-5-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_5)[1] <- "Initial_Position"
names(Recombination_rate_scaff_5)[2] <- "Final_Position"
names(Recombination_rate_scaff_5)[3] <- "Rec_Rate"

data5 = data.frame(Recombination_rate_scaff_5)

p <- ggplot(data5, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 5") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,27223299), breaks = seq(0,27223299, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_6 <- read_table2("scaff-6-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_6)[1] <- "Initial_Position"
names(Recombination_rate_scaff_6)[2] <- "Final_Position"
names(Recombination_rate_scaff_6)[3] <- "Rec_Rate"

data6 = data.frame(Recombination_rate_scaff_6)

p <- ggplot(data6, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 6") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,26593402), breaks = seq(0,26593402, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_7 <- read_table2("scaff-7-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_7)[1] <- "Initial_Position"
names(Recombination_rate_scaff_7)[2] <- "Final_Position"
names(Recombination_rate_scaff_7)[3] <- "Rec_Rate"

data7 = data.frame(Recombination_rate_scaff_7)

p <- ggplot(data7, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 7") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,25418471), breaks = seq(0,25418471, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_8 <- read_table2("scaff-8-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_8)[1] <- "Initial_Position"
names(Recombination_rate_scaff_8)[2] <- "Final_Position"
names(Recombination_rate_scaff_8)[3] <- "Rec_Rate"

data8 = data.frame(Recombination_rate_scaff_8)

p <- ggplot(data8, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 8") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,25271786), breaks = seq(0,25271786, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_9 <- read_table2("scaff-9-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_9)[1] <- "Initial_Position"
names(Recombination_rate_scaff_9)[2] <- "Final_Position"
names(Recombination_rate_scaff_9)[3] <- "Rec_Rate"

data9 = data.frame(Recombination_rate_scaff_9)

p <- ggplot(data9, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 9") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,24821562), breaks = seq(0,24821562, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_10 <- read_table2("scaff-10-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_10)[1] <- "Initial_Position"
names(Recombination_rate_scaff_10)[2] <- "Final_Position"
names(Recombination_rate_scaff_10)[3] <- "Rec_Rate"

data10 = data.frame(Recombination_rate_scaff_10)

p <- ggplot(data10, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 10") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,23696351), breaks = seq(0,23696351, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_11 <- read_table2("scaff-11-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_11)[1] <- "Initial_Position"
names(Recombination_rate_scaff_11)[2] <- "Final_Position"
names(Recombination_rate_scaff_11)[3] <- "Rec_Rate"

data11 = data.frame(Recombination_rate_scaff_11)

p <- ggplot(data11, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 11") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,22599691), breaks = seq(0,22599691, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_12 <- read_table2("scaff-12-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_12)[1] <- "Initial_Position"
names(Recombination_rate_scaff_12)[2] <- "Final_Position"
names(Recombination_rate_scaff_12)[3] <- "Rec_Rate"

data12 = data.frame(Recombination_rate_scaff_12)

p <- ggplot(data12, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 12") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,22309355), breaks = seq(0,22309355, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_13 <- read_table2("scaff-13-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_13)[1] <- "Initial_Position"
names(Recombination_rate_scaff_13)[2] <- "Final_Position"
names(Recombination_rate_scaff_13)[3] <- "Rec_Rate"

data13 = data.frame(Recombination_rate_scaff_13)

p <- ggplot(data13, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 13") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,21837515), breaks = seq(0,21837515, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_14 <- read_table2("scaff-14-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_14)[1] <- "Initial_Position"
names(Recombination_rate_scaff_14)[2] <- "Final_Position"
names(Recombination_rate_scaff_14)[3] <- "Rec_Rate"

data14 = data.frame(Recombination_rate_scaff_14)

p <- ggplot(data14, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 14") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,21002338), breaks = seq(0,21002338, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_15 <- read_table2("scaff-15-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_15)[1] <- "Initial_Position"
names(Recombination_rate_scaff_15)[2] <- "Final_Position"
names(Recombination_rate_scaff_15)[3] <- "Rec_Rate"

data15 = data.frame(Recombination_rate_scaff_15)

p <- ggplot(data15, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 15") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,20771620), breaks = seq(0,20771620, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_16 <- read_table2("scaff-16-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_16)[1] <- "Initial_Position"
names(Recombination_rate_scaff_16)[2] <- "Final_Position"
names(Recombination_rate_scaff_16)[3] <- "Rec_Rate"

data16 = data.frame(Recombination_rate_scaff_16)

p <- ggplot(data16, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 16") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,20354528), breaks = seq(0,20354528, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_17 <- read_table2("scaff-17-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_17)[1] <- "Initial_Position"
names(Recombination_rate_scaff_17)[2] <- "Final_Position"
names(Recombination_rate_scaff_17)[3] <- "Rec_Rate"

data17 = data.frame(Recombination_rate_scaff_17)

p <- ggplot(data17, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 17") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,19406976), breaks = seq(0,19406976, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_18 <- read_table2("scaff-18-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_18)[1] <- "Initial_Position"
names(Recombination_rate_scaff_18)[2] <- "Final_Position"
names(Recombination_rate_scaff_18)[3] <- "Rec_Rate"

data18 = data.frame(Recombination_rate_scaff_18)

p <- ggplot(data18, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 18") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,18515410), breaks = seq(0,18515410, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p


Recombination_rate_scaff_19 <- read_table2("scaff-19-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_19)[1] <- "Initial_Position"
names(Recombination_rate_scaff_19)[2] <- "Final_Position"
names(Recombination_rate_scaff_19)[3] <- "Rec_Rate"

data19 = data.frame(Recombination_rate_scaff_19)

p <- ggplot(data19, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 19") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,17773115), breaks = seq(0,17773115, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_20 <- read_table2("scaff-20-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_20)[1] <- "Initial_Position"
names(Recombination_rate_scaff_20)[2] <- "Final_Position"
names(Recombination_rate_scaff_20)[3] <- "Rec_Rate"

data20 = data.frame(Recombination_rate_scaff_20)

p <- ggplot(data20, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 20") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,16307044), breaks = seq(0,16307044, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_21 <- read_table2("scaff-21-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_21)[1] <- "Initial_Position"
names(Recombination_rate_scaff_21)[2] <- "Final_Position"
names(Recombination_rate_scaff_21)[3] <- "Rec_Rate"

data21 = data.frame(Recombination_rate_scaff_21)

p <- ggplot(data21, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 21") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,16244046), breaks = seq(0,16244046, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_22 <- read_table2("scaff-22-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_22)[1] <- "Initial_Position"
names(Recombination_rate_scaff_22)[2] <- "Final_Position"
names(Recombination_rate_scaff_22)[3] <- "Rec_Rate"

data22 = data.frame(Recombination_rate_scaff_22)

p <- ggplot(data22, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 22") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,15785361), breaks = seq(0,15785361, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_23 <- read_table2("scaff-23-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_23)[1] <- "Initial_Position"
names(Recombination_rate_scaff_23)[2] <- "Final_Position"
names(Recombination_rate_scaff_23)[3] <- "Rec_Rate"

data23 = data.frame(Recombination_rate_scaff_23)

p <- ggplot(data23, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 23") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,14824030), breaks = seq(0,14824030, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_24 <- read_table2("scaff-24-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_24)[1] <- "Initial_Position"
names(Recombination_rate_scaff_24)[2] <- "Final_Position"
names(Recombination_rate_scaff_24)[3] <- "Rec_Rate"

data24 = data.frame(Recombination_rate_scaff_24)

p <- ggplot(data24, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 24") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,13684430), breaks = seq(0,13684430, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_25 <- read_table2("scaff-25-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_25)[1] <- "Initial_Position"
names(Recombination_rate_scaff_25)[2] <- "Final_Position"
names(Recombination_rate_scaff_25)[3] <- "Rec_Rate"

data25 = data.frame(Recombination_rate_scaff_25)

p <- ggplot(data25, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 25") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,13471441), breaks = seq(0,13471441, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_26 <- read_table2("scaff-26-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_26)[1] <- "Initial_Position"
names(Recombination_rate_scaff_26)[2] <- "Final_Position"
names(Recombination_rate_scaff_26)[3] <- "Rec_Rate"

data26 = data.frame(Recombination_rate_scaff_26)

p <- ggplot(data26, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 26") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,12537843), breaks = seq(0,12537843, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_27 <- read_table2("scaff-27-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_27)[1] <- "Initial_Position"
names(Recombination_rate_scaff_27)[2] <- "Final_Position"
names(Recombination_rate_scaff_27)[3] <- "Rec_Rate"

data27 = data.frame(Recombination_rate_scaff_27)

p <- ggplot(data27, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 27") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,12376488), breaks = seq(0,12376488, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_28 <- read_table2("scaff-28-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_28)[1] <- "Initial_Position"
names(Recombination_rate_scaff_28)[2] <- "Final_Position"
names(Recombination_rate_scaff_28)[3] <- "Rec_Rate"

data28 = data.frame(Recombination_rate_scaff_28)

p <- ggplot(data28, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 28") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,11658327), breaks = seq(0,11658327, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p



Recombination_rate_scaff_29 <- read_table2("scaff-29-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_29)[1] <- "Initial_Position"
names(Recombination_rate_scaff_29)[2] <- "Final_Position"
names(Recombination_rate_scaff_29)[3] <- "Rec_Rate"

data29 = data.frame(Recombination_rate_scaff_29)

p <- ggplot(data29, aes(x=Initial_Position, y=Rec_Rate)) 
p <- p + theme_light () + ggtitle("Recombination Rate in Chromosome 29") +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,11165203), breaks = seq(0,11165203, by = 5000000))
p <- p + scale_y_log10() + xlab("Position (Kb)") + ylab("Recombination Rate (CO/bp)")
p










#Here is the code needed to load the data and name the different columns of the table.
Recombination_rate_scaff_1 <- read_table2("scaff-1-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_1)[1] <- "Initial_Position"
names(Recombination_rate_scaff_1)[2] <- "Final_Position"
names(Recombination_rate_scaff_1)[3] <- "Rec_Rate"

Recombination_rate_scaff_22 <- read_table2("scaff-22-miss0.97.rmap", col_names = FALSE)
names(Recombination_rate_scaff_22)[1] <- "Initial_Position"
names(Recombination_rate_scaff_22)[2] <- "Final_Position"
names(Recombination_rate_scaff_22)[3] <- "Rec_Rate"

#We transform it into a data frame to manage the high volumes of data better.
data1 = data.frame(Recombination_rate_scaff_1)
data22 = data.frame(Recombination_rate_scaff_22)

# We can convert the positions to Mb, by dividing their values by 10^6, and create a new column for such values.
data1 <- data1 %>% mutate(Mb = Initial_Position / 1000000)
data1 <- data1 %>% mutate(cM = Rec_Rate / 0.00000001)


data22 <- data22 %>% mutate(Mb = Initial_Position / 1000000)
data22 <- data22 %>% mutate(cM = Rec_Rate / 0.00000001)

#Now we plot the data.
p <- ggplot(data1, aes(x=Mb, y=cM)) 
p <- p + theme_light ()  + theme(axis.title.x = element_blank()) +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,34.455698), breaks = seq(0,34.455698, by = 5))
p <- p + scale_y_continuous(limits = c(0,100), breaks = c(0,10,25,50,100)) + ylab("Recombination Rate (cM/Mb)")
p

q <- ggplot(data22, aes(x=Mb, y=cM)) 
q <- q + theme_light () +
  geom_line(color="grey") + scale_x_continuous(limits = c(0,15.785361), breaks = seq(0,15.785361, by = 5))
q <- q + scale_y_continuous(limits = c(0,100), breaks = c(0,10,25,50,100)) + xlab("Position (Mb)") + ylab("Recombination Rate (cM/Mb)")
q

# We now unify all plots and represent them in a single figure.
grid.arrange(p, q, ncol = 1, nrow = 2)

