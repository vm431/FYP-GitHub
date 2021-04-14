library(datasets)
install.packages("pacman")
pacman::p_load(pacman, caret, lars, tidyverse, dplyr, GGally, ggplot2, ggthemes,ggvis, httr, plotly, shiny, tidyr)

codons = read.csv("C:\\Users\\Valee\\Desktop\\FYP-GitHub\\Results\\data1.csv", header=TRUE, sep = ",")
colours = read.csv("C:\\Users\\Valee\\Desktop\\FYP\\Results\\colours.csv", header=TRUE, sep = ",")

nullx = 0
nully = 0
x = codons$intergenicGC

# 4 synonymous codons

plot(nullx, nully,
     main = "A-ending 4 synonymous codons usage over intergenic GC",
     xlab = "Intergenic GC content (%)",
     ylab = "A-ending 4 synonymous codons usage (%)",
     cex = 0.2,
     pch = 19,
     xlim = c(0,85),
     ylim = c(0,10.5))
legend(70, 10, legend=c(colours$codonA4),
       lty = 1,
       fill = colours$colour3, 
)

ACA = codons$ACA
TCA = codons$TCA
CCA = codons$CCA
CGA = codons$CGA
GCA = codons$GCA
GGA = codons$GGA
CTA = codons$CTA
GTA = codons$GTA
lineACA = abline(lm(ACA ~ x), col = "#FB9902",lwd = 1.5)
lineTCA = abline(lm(TCA ~ x), col = "#B2D732",lwd = 1.5)
lineCTA = abline(lm(CTA ~ x), col = "#347C98",lwd = 1.5)
lineCCA = abline(lm(CCA ~ x), col = "#0247FE",lwd = 1.5)
lineCGA = abline(lm(CGA ~ x), col = "#4424D6",lwd = 1.5) #lowest line
lineGTA = abline(lm(GTA ~ x), col = "#C21460",lwd = 1.5)
lineGCA = abline(lm(GCA ~ x), col = "#E67CCF",lwd = 1.5)
lineGGA = abline(lm(GGA ~ x), col = "#FF3366",lwd = 1.5)

# 2 synonymous codons

plot(nullx, nully,
     main = "A-ending 2 synonymous codons usage over intergenic GC",
     xlab = "Intergenic GC content (%)",
     ylab = "A-ending 2 synonymous codons usage (%)",
     cex = 0.2,
     pch = 19,
     xlim = c(0,85),
     ylim = c(0,10.5))
legend(70, 10, legend=c(colours$codonA2),
       lty = 1,
       fill = colours$colour, 
)

AAA = codons$AAA
AGA = codons$AGA
TTA = codons$TTA
CAA = codons$CAA
GAA = codons$GAA
lineAAA = abline(lm(AAA ~ x), col = "#FE2712",lwd = 1.5)
lineAGA = abline(lm(AGA ~ x), col = "#e06377",lwd = 1.5)
lineTTA = abline(lm(TTA ~ x), col = "#FEFE33",lwd = 1.5)
lineCAA = abline(lm(CAA ~ x), col = "#66B032",lwd = 1.5)
lineGAA = abline(lm(GAA ~ x), col = "#b8a9c9",lwd = 1.5)

