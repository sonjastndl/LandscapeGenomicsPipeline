
library(tidyverse)
library(gridExtra)
library(ggplot2)

args <- commandArgs(TRUE)

## set working directory
setwd(args[1])



## read AlleleFrequency Matrix
DATA=read.table(args[2],
  header=T)
AF=DATA[,3:ncol(DATA)]

## read Metadata
meta=read.csv(args[3],
  header=T)

## make sure that order of AF data and Meta data match
meta<-meta[match(colnames(AF),meta$sampleId),]

## make suuuuper simple function that fits a linear regression model and returns the p-values
lm.func <- function(x) {
      summary(lm(unlist(x)~y))[[4]][,4][2]
}

## fit linear models for both biovariables to AF data in each row (i.e. SNP) and append corresponding p-value to new column in raw Data table
for ( i in c("bio1","bio12","lat")){
  y=meta[[i]]
  p.val<-apply(AF,1,lm.func)
  ID=paste0(i,".pval")
  DATA[[ID]]<-p.val
}

## multiple testing problem??
y=runif(ncol(AF))
D=data.frame("X.1"=runif(nrow(AF)))
for (i in seq(2,ncol(AF),1)){
  D[[paste0("X.",i)]]<-runif(nrow(AF))
}

Bonf=0.05/nrow(DATA)
Test.p<-apply(D,1,lm.func)
Test.p[Test.p<Bonf]

pdf("results/Multtest_control.pdf",
  width=15,
  height=5)
hist(Test.p,breaks=100)
dev.off()

## Boferroni-correctd p-value threshold
Bonf=-log10(0.05/nrow(DATA))

## plot with ggplot
PLOT<-c("bio1.pval","bio12.pval","lat.pval") %>%
  map(function(z){
    PLOT.df<-data.frame(Pos = DATA$Pos/1000000, P.val = -log10(DATA[[z]]))
    pl <- ggplot(PLOT.df, aes(x = Pos, y = P.val)) +
      geom_point(alpha=0.3) +
      xlab("Genomic Position on 3R [Mbp]") +
      ylab("-log10(P.value)")+
      geom_hline(yintercept=Bonf,
        col="blue",
        lty=2)+
      geom_hline(yintercept=-log10(0.05),
        col="red",
        lty=2)+
      ggtitle(z)+
      theme_bw()
    return(pl)
  }) %>%
  arrangeGrob(grobs = ., nrow=3) %>% #add ncol/nrow argument here
grid.arrange()

ggsave("results/3R_Pvalues.png",
  PLOT,
  width=15,
  height=7)


