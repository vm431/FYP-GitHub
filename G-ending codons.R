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
     main = "G-ending 4 synonymous codons usage over intergenic GC",
     xlab = "Intergenic GC content (%)",
     ylab = "G-ending 4 synonymous codons usage (%)",
     cex = 0.2,
     pch = 19,
     xlim = c(0,85),
     ylim = c(0,9))
legend(0, 8.9, legend=c(colours$codonG4),
       lty = 1,
       fill = colours$colour3, 
)

ACG = codons$ACG
TCG = codons$TCG
CCG = codons$CCG
CGG = codons$CGG
GCG = codons$GCG
GGG = codons$GGG
CTG = codons$CTG
GTG = codons$GTG
lineACG = abline(lm(ACG ~ x), col = "#FB9902", lwd = 1.5)
lineTCG = abline(lm(TCG ~ x), col = "#B2D732", lwd = 1.5)
lineCTG = abline(lm(CTG ~ x), col = "#347C98", lwd = 1.5)
lineCCG = abline(lm(CCG ~ x), col = "#0247FE", lwd = 1.5) 
lineCGG = abline(lm(CGG ~ x), col = "#4424D6", lwd = 1.5)
lineGTG = abline(lm(GTG ~ x), col = "#C21460", lwd = 1.5)
lineGCG = abline(lm(GCG ~ x), col = "#E67CCF", lwd = 1.5)
lineGGG = abline(lm(GGG ~ x), col = "#d9ad7c", lwd = 1.5)

# 2 synonymous codons

plot(nullx, nully,
     main = "G-ending 2 synonymous codons usage over intergenic GC",
     xlab = "Intergenic GC content (%)",
     ylab = "G-ending 2 synonymus codons usage (%)",
     cex = 0.2,
     pch = 19,
     xlim = c(0,85),
     ylim = c(0,9))
legend(0, 8.9, legend=c(colours$codonG2),
       lty = 1,
       fill = colours$colour, 
)

AAG = codons$AAG
AGG = codons$AGG
TTG = codons$TTG
CAG = codons$CAG
GAG = codons$GAG
lineAAG = abline(lm(AAG ~ x), col = "#FE2712", lwd = 1.5)
lineAGG = abline(lm(AGG ~ x), col = "#e06377", lwd = 1.5) #lowest line
lineTTG = abline(lm(TTG ~ x), col = "#FEFE33", lwd = 1.5) #goes down?
lineCAG = abline(lm(CAG ~ x), col = "#66B032", lwd = 1.5)
lineGAG = abline(lm(GAG ~ x), col = "#b8a9c9", lwd = 1.5)

