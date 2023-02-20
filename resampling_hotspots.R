## 10k resampling of means
library("gridExtra")
library("ggplot2")

## LINEs
genomedata <- read.table("TE_LINE-1657_coverage.txt", sep="\t", header=FALSE)
hotspotdata <- read.table("TE_LINE-hotspots_coverage.txt", sep="\t", header=FALSE)
meanhotspots <- mean(hotspotdata$V7)
meangenomewide <- mean(genomedata$V7)

## Sample 10k means from 3124 windows
results<-vector('list',10000)
for(i in 1:10000){
  x <- sample(genomedata$V7, 3124, replace = F)
  results[[i]]<-mean(x)  
}

## Calculate p-value
pvalue <- sum(unlist(results) < meanhotspots) /10000
pvaluetwotailed <- pvalue*2

## Plot results
df <- data.frame(difs = unlist(results))
charttitle <- paste("LINE density:", "p < 0.001")

LINE <- ggplot(df, aes(x=difs)) + theme_light() +
  geom_histogram(color="white",fill="grey") +
  geom_vline(color="red",xintercept = meanhotspots) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(size=12, face="bold")) +
  labs(title = charttitle, tag = "D", x = "LINE (fraction of window)", y = "Count")
LINE

## SINEs
genomedata <- read.table("TE_SINE-1657_coverage.txt", sep="\t", header=FALSE)
hotspotdata <- read.table("TE_SINE-hotspots_coverage.txt", sep="\t", header=FALSE)
meanhotspots <- mean(hotspotdata$V7)
meangenomewide <- mean(genomedata$V7)

## Sample 10k means from 3124 windows
results<-vector('list',10000)
for(i in 1:10000){
  x <- sample(genomedata$V7, 3124, replace = F)
  results[[i]]<-mean(x)  
}

## Calculate p-value
pvalue <- sum(unlist(results) > meanhotspots) /10000
pvaluetwotailed <- pvalue*2

## Plot results
df <- data.frame(difs = unlist(results))
charttitle <- paste("SINE density:", "p < 0.001")

SINE <- ggplot(df, aes(x=difs)) + theme_light() +
  geom_histogram(color="white",fill="grey") +
  geom_vline(color="red",xintercept = meanhotspots) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(size=12, face="bold")) +
  labs(title = charttitle, tag = "C", x = "SINE (franction of window)", y = "Count")
SINE

## LTRs
genomedata <- read.table("TE_LTR-1657_coverage.txt", sep="\t", header=FALSE)
hotspotdata <- read.table("TE_LTR-hotspots_coverage.txt", sep="\t", header=FALSE)
meanhotspots <- mean(hotspotdata$V7)
meangenomewide <- mean(genomedata$V7)

## Sample 10k means from 3124 windows
results<-vector('list',10000)
for(i in 1:10000){
  x <- sample(genomedata$V7, 3124, replace = F)
  results[[i]]<-mean(x)  
}

## Calculate p-value
pvalue <- sum(unlist(results) < meanhotspots) /10000
pvaluetwotailed <- pvalue*2

## Plot results
df <- data.frame(difs = unlist(results))
charttitle <- paste("LTR density:", "p < 0.001")

LTR <- ggplot(df, aes(x=difs)) + theme_light() +
  geom_histogram(color="white",fill="grey") +
  geom_vline(color="red",xintercept = meanhotspots) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(size=12, face="bold")) +
  labs(title = charttitle, tag = "B", x = "LTR (fraction of window)", y = "Count")
LTR

## DNAs
genomedata <- read.table("TE_DNA-1657_coverage.txt", sep="\t", header=FALSE)
hotspotdata <- read.table("TE_DNA-hotspots_coverage.txt", sep="\t", header=FALSE)
meanhotspots <- mean(hotspotdata$V7)
meangenomewide <- mean(genomedata$V7)

## Sample 10k means from 3124 windows
results<-vector('list',10000)
for(i in 1:10000){
  x <- sample(genomedata$V7, 3124, replace = F)
  results[[i]]<-mean(x)  
}

## Calculate p-value
pvalue <- sum(unlist(results) < meanhotspots) /10000
pvaluetwotailed <- pvalue*2

## Plot results
df <- data.frame(difs = unlist(results))
charttitle <- paste("DNA density:", "p = 0.509" )

DNA <- ggplot(df, aes(x=difs)) + theme_light() +
  geom_histogram(color="white",fill="grey") +
  geom_vline(color="red",xintercept = meanhotspots) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(size=12, face="bold")) +
  labs(title = charttitle, tag = "A", x = "DNA (fraction of window)", y = "Count")
DNA


p <- grid.arrange(DNA, LTR, SINE, LINE, ncol=2, nrow =2)
p

ggsave("Resampling-Hotspots.png", plot = p, units= c("in"), width = 10, height = 5, dpi = 300)

