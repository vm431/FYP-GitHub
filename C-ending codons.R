library(datasets)
install.packages("pacman")
pacman::p_load(pacman, caret, lars, tidyverse, dplyr, GGally, ggplot2, ggthemes,ggvis, httr, plotly, shiny, tidyr)

codons = read.csv("C:\\Users\\Valee\\Desktop\\FYP-GitHub\\results\\data1.csv", header=TRUE, sep = ",")
colours = read.csv("C:\\Users\\Valee\\Desktop\\FYP\\Results\\colours.csv", header=TRUE, sep = ",")

nullx = 0
nully = 0
x = codons$IntergenicGC

# 4 synonymous codons

plot(nullx, nully,
     main = "C-ending 4 synonymous codons usage over intergenic GC",
     xlab = "Intergenic GC content (%)",
     ylab = "C-ending 4 synonymous codons usage (%)",
     cex = 0.2,
     pch = 19,
     xlim = c(0,100),
     ylim = c(0,100))
legend(0, 95, legend=c(colours$codonC4),
       lty = 1,
       fill = colours$colour2, 
)

ACC = codons$ACC
TCC = codons$TCC
CCC = codons$CCC
CGC = codons$CGC
GCC = codons$GCC
GGC = codons$GGC
CTC = codons$CTC
GTC = codons$GTC
lineACC = abline(lm(ACC ~ x), col = "#FB9902", lwd = 1.5)
lineTCC = abline(lm(TCC ~ x), col = "#B2D732", lwd = 1.5) #low line
lineCTC = abline(lm(CTC ~ x), col = "#347C98", lwd = 1.5)
lineCCC = abline(lm(CCC ~ x), col = "#0247FE", lwd = 1.5)
lineCGC = abline(lm(CGC ~ x), col = "#4424D6", lwd = 1.5) 
lineGTC = abline(lm(GTC ~ x), col = "#C21460", lwd = 1.5)
lineGCC = abline(lm(GCC ~ x), col = "#E67CCF", lwd = 1.5) #most used
lineGGC = abline(lm(GGC ~ x), col = "#FF3366", lwd = 1.5)


# 2 synonymous codons


plot(nullx, nully,
     main = "C-ending 2 synonymous codons usage over intergenic GC",
     xlab = "Intergenic GC content (%)",
     ylab = "C-ending 2 synonymous codons usage (%)",
     cex = 0.2,
     pch = 19,
     xlim = c(0,100),
     ylim = c(0,100))
legend(0, 95, legend=c(colours$codonC2),
       lty = 1,
       fill = colours$colour1, 
)

AAC = codons$AAC
AGC = codons$AGC
TAC = codons$TAC
TTC = codons$TTC
TGC = codons$TGC
CAC = codons$CAC
GAC = codons$GAC
lineAAC = abline(lm(AAC ~ x), col = "#FE2712", lwd = 1.5) #stable line
lineAGC = abline(lm(AGC ~ x), col = "#e06377", lwd = 1.5) #low line
lineTAC = abline(lm(TAC ~ x), col = "#d9ad7c", lwd = 1.5) #low line
lineTTC = abline(lm(TTC ~ x), col = "#FEFE33", lwd = 1.5)
lineTGC = abline(lm(TGC ~ x), col = "#a2836e", lwd = 1.5) #lowest line
lineGAC = abline(lm(GAC ~ x), col = "#b8a9c9", lwd = 1.5)



