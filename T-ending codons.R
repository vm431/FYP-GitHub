library(datasets)
install.packages("pacman")
pacman::p_load(pacman, caret, lars, tidyverse, dplyr, GGally, ggplot2, ggthemes,ggvis, httr, plotly, shiny, tidyr)

codons = read.csv("C:\\Users\\Valee\\Desktop\\FYP-GitHub\\results\\data1.csv", header=TRUE, sep = ",")
colours = read.csv("C:\\Users\\Valee\\Desktop\\FYP\\Results\\colours.csv", header=TRUE, sep = ",")

nullx = 0
nully = 0
x = codons$IntergenicGC

# 4 synonymous codons
ACT = codons$ACT
TCT = codons$TCT
CCT = codons$CCT
CGT = codons$CGT
GCT = codons$GCT
GGT = codons$GGT
CTT = codons$CTT
GTT = codons$GTT
lineACT = abline(lm(ACT ~ x), col = "#FB9902", lwd = 1.5)
lineTCT = abline(lm(TCT ~ x), col = "#B2D732", lwd = 1.5)
lineCTT = abline(lm(CTT ~ x), col = "#347C98", lwd = 1.5)
lineCCT = abline(lm(CCT ~ x), col = "#0247FE", lwd = 1.5) #stable line
lineCGT = abline(lm(CGT ~ x), col = "#4424D6", lwd = 1.5) 
lineGTT = abline(lm(GTT ~ x), col = "#C21460", lwd = 1.5)
lineGCT = abline(lm(GCT ~ x), col = "#E67CCF", lwd = 1.5)
lineGGT = abline(lm(GGT ~ x), col = "#FF3366", lwd = 1.5)

# 2 synonymous codons
AAT = codons$AAT
AGT = codons$AGT
TAT = codons$TAT
TTT = codons$TTT
CAT = codons$CAT
GAT = codons$GAT
lineAAT = abline(lm(AAT ~ x), col = "#FE2712", lwd = 1.5)
lineAGT = abline(lm(AGT ~ x), col = "#e06377", lwd = 1.5)
lineTAT = abline(lm(TAT ~ x), col = "#d9ad7c", lwd = 1.5)
lineTTT = abline(lm(TTT ~ x), col = "#FEFE33", lwd = 1.5)
lineTGT = abline(lm(TGT ~ x), col = "#a2836e", lwd = 1.5) #lowest line
lineCAT = abline(lm(CAT ~ x), col = "#66B032", lwd = 1.5)
lineGAT = abline(lm(GAT ~ x), col = "#b8a9c9", lwd = 1.5)




plot(nullx, nully,
     main = "T-ending 4 synonymous codons usage over intergenic GC",
     xlab = "Intergenic GC content (%)",
     ylab = "T-ending 4 synonymous codons usage (%)",
     cex = 0.2,
     pch = 19,
     xlim = c(0,100),
     ylim = c(0,100))
legend(75, 90, legend=c(colours$codonT4),
       lty = 1,
       fill = colours$colour2, 
)

plot(nullx, nully,
     main = "T-ending 2 synonymous codons usage over intergenic GC",
     xlab = "Intergenic GC content (%)",
     ylab = "T-ending 2 synonymous codons usage (%)",
     cex = 0.2,
     pch = 19,
     xlim = c(0,100),
     ylim = c(0,100))
legend(75, 90, legend=c(colours$codonT2),
       lty = 1,
       fill = colours$colour1, 
)
