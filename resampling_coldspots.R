## 10k resampling of means
library("gridExtra")
library("ggplot2")
library(vegan)
library(cowplot)
## LINEs
genomedata <- read.table("TE_LINE-29768_coverage.txt", sep="\t", header=FALSE)
coldspotdata <- read.table("TE_LINE-coldspots_coverage.txt", sep="\t", header=FALSE)
meancoldspots <- mean(coldspotdata$V7)
meangenomewide <- mean(genomedata$V7)

## Sample 10k means from 1283 windows
results<-vector('list',10000)
for(i in 1:10000){
  x <- sample(genomedata$V7, 1283, replace = F)
  results[[i]]<-mean(x)  
}

## Calculate p-value
pvalue <- sum(unlist(results) < meancoldspots) /10000
pvaluetwotailed <- pvalue*2

## Plot results
df <- data.frame(difs = unlist(results))
charttitle <- paste("LINE density:", "p = 0.002" )

LINE <- ggplot(df, aes(x=difs)) + theme_light() +
  geom_histogram(color="white",fill="grey") +
  geom_vline(color="red",xintercept = meancoldspots) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(size=12, face="bold")) +
  labs(title = charttitle, tag = "D", x = "LINE (fraction of window)", y = "Count")
LINE

## SINEs
genomedata <- read.table("TE_SINE-29768_coverage.txt", sep="\t", header=FALSE)
coldspotdata <- read.table("TE_SINE-coldspots_coverage.txt", sep="\t", header=FALSE)
meancoldspots <- mean(coldspotdata$V7)
meangenomewide <- mean(genomedata$V7)

## Sample 10k means from 1283 windows
results<-vector('list',10000)
for(i in 1:10000){
  x <- sample(genomedata$V7, 1283, replace = F)
  results[[i]]<-mean(x)  
}

## Calculate p-value
pvalue <- sum(unlist(results) < meancoldspots) /10000
pvaluetwotailed <- pvalue*2

## Plot results
df <- data.frame(difs = unlist(results))
charttitle <- paste("SINE density:", "p < 0.001" )

SINE <- ggplot(df, aes(x=difs)) + theme_light() +
  geom_histogram(color="white",fill="grey") + 
  geom_vline(color="red",xintercept = meancoldspots) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(size=12, face="bold")) +
  labs(title = charttitle, tag = "C", x = "SINE (fraction of window)", y = "Count")
SINE


## LTRs
genomedata <- read.table("TE_LTR-29768_coverage.txt", sep="\t", header=FALSE)
coldspotdata <- read.table("TE_LTR-coldspots_coverage.txt", sep="\t", header=FALSE)
meancoldspots <- mean(coldspotdata$V7)
meangenomewide <- mean(genomedata$V7)

## Sample 10k means from 1283 windows
results<-vector('list',10000)
for(i in 1:10000){
  x <- sample(genomedata$V7, 1283, replace = F)
  results[[i]]<-mean(x)  
}

## Calculate p-value
pvalue <- sum(unlist(results) < meancoldspots) /10000
pvaluetwotailed <- pvalue*2

## Plot results
df <- data.frame(difs = unlist(results))
charttitle <- paste("LTR density:", "p = 0.007" )

LTR <- ggplot(df, aes(x=difs)) + theme_light() +
  geom_histogram(color="white",fill="grey") +
  geom_vline(color="red",xintercept = meancoldspots) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(size=12, face="bold")) +
  labs(title = charttitle, tag = "B", x = "LTR (fraction of window)", y = "Count")
LTR

## DNAs
genomedata <- read.table("TE_DNA-29768_coverage.txt", sep="\t", header=FALSE)
coldspotdata <- read.table("TE_DNA-coldspots_coverage.txt", sep="\t", header=FALSE)
meancoldspots <- mean(coldspotdata$V7)
meangenomewide <- mean(genomedata$V7)

## Sample 10k means from 1283 windows
results<-vector('list',10000)
for(i in 1:10000){
  x <- sample(genomedata$V7, 1283, replace = F)
  results[[i]]<-mean(x)  
}

## Calculate p-value
pvalue <- sum(unlist(results) < meancoldspots) /10000
pvaluetwotailed <- pvalue*2

## Plot results
df <- data.frame(difs = unlist(results))
charttitle <- paste("DNA density:", "p < 0.001")

DNA <- ggplot(df, aes(x=difs)) + theme_light() +
  geom_histogram(color="white",fill="grey") +
  geom_vline(color="red",xintercept = meancoldspots) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(size=12, face="bold")) +
  labs(title = charttitle, tag = "A", x = "DNA (fraction of window)", y = "Count")
DNA

# Create the final plot
p <- grid.arrange(DNA, LTR, SINE, LINE, ncol=2, nrow =2)
p

ggsave("Resampling-Coldspots.png", plot = p, units= c("in"), width = 10, height = 5, dpi = 300)

