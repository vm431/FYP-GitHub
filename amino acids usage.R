library(datasets)
install.packages("pacman")
pacman::p_load(pacman, caret, lars, tidyverse, dplyr, GGally, ggplot2, ggthemes,ggvis, httr, plotly, shiny, tidyr, ggrepel, dplyr)


GCdata = read.csv("C:\\Users\\Valee\\Desktop\\FYP\\Results\\GCdata.csv", header=TRUE, sep = ",")

#GC over intergenomic GC

sp <- plot(GCdata$IntergenomicGC, GCdata$GCcontent,
      main = "Genomic over Intergenomic GC content",
      xlab = "Intergenomic GC content",
      ylab = "Genomic GC content",
      col = "#7C3194",
      cex = 0.8,
      xlim = c(0,90),
      ylim = c(0,90),
      pch = 19,
      abline(lm(GCdata$GCcontent ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(55, 20, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))

# GC3 over intergenomic GC content

regline <- lm(GCdata$GC3content ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$GC3content,
     main = "Genomic GC3 over Intergenomic GC content",
     xlab = "Intergenomic GC content",
     ylab = "Genomic GC3 content",
     col = "#7C3194",
     xlim = c(0,90),
     ylim = c(0,100),
     pch = 19,
     cex = 0.8,
     abline(lm(GCdata$GC3content ~ GCdata$IntergenomicGC), col = "#431A50" )
     )
legend(55, 20, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))

#amino acid usage over intergenomic GC

regline <- lm(GCdata$GlyPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$GlyPer,
       main = "Glycine usage over Intergenomic GC content",
       xlab = "Intergenomic GC content",
       ylab = "Glycine usage",
       col = "#7C3194", 
       pch = 19,
       xlim = c(0,90),
       ylim = c(0,15),
       cex = 0.8,
       abline(lm(GCdata$GlyPer ~ GCdata$IntergenomicGC)))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                               paste("r2 =", rsq)))


regline <- lm(GCdata$PhePer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$PhePer,
           main = "Phenylalanineine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Phenylalanine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$PhePer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))



regline <- lm(GCdata$IlePer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$IlePer,
           main = "Isoleucine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Isoleucine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$IlePer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))



regline <- lm(GCdata$MetPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$MetPer,
           main = "Methionine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Methionine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$MetPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$ValPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$ValPer,
           main = "Valine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Valine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$ValPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$ProPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$ProPer,
           main = "Proline usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Proline usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$ProPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$ThrPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$ThrPer,
           main = "Threonine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Threonine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$ThrPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$AlaPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$AlaPer,
           main = "Alanine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Alanine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$AlaPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40,3, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$TyrPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$TyrPer,
           main = "Tyrosine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Tyrosine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$TyrPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$HisPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$HisPer,
           main = "Histidine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Histidine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$HisPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$GlnPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$GlnPer,
           main = "Glutamine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Glutamine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$GlnPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$AsnPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$AsnPer,
           main = "Asparagine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Asparagine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$AsnPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$LysPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$LysPer,
           main = "Lysine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Lysine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$LysPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$AspPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$AspPer,
           main = "Aspartic acid usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Aspartic acid usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$AspPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$GluPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$GluPer,
           main = "Glutamic acid usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Glutamic acid usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$GluPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$CysPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$CysPer,
           main = "Cysteine usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Cysteine usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$CysPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$TrpPer ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$TrpPer,
           main = "Tryptophan usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Tryptophan usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(GCdata$TrpPer ~ GCdata$IntergenomicGC), col = "#431A50" ))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))



LSRdata = read.csv("C:\\Users\\Valee\\Desktop\\FYP\\Results\\LSRdata.csv", header=TRUE, sep = ",")

regline <- lm(LSRdata$Ser2Per ~ LSRdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(LSRdata$IntergenomicGC, LSRdata$Ser2Per,
           main = "Serine 2 usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Serine 2 usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(LSRdata$Ser2Per ~ LSRdata$IntergenomicGC)))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(LSRdata$Arg2Per ~ LSRdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(LSRdata$IntergenomicGC, LSRdata$Arg2Per,
           main = "Arginine 2 usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Arginine 2 usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(LSRdata$Arg2Per ~ LSRdata$IntergenomicGC)))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))



regline <- lm(LSRdata$Arg4Per ~ LSRdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(LSRdata$IntergenomicGC, LSRdata$Arg4Per,
           main = "Arginine 4 usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Arginine 4 usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(LSRdata$Arg4Per ~ LSRdata$IntergenomicGC)))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(LSRdata$Ser4Per ~ LSRdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(LSRdata$IntergenomicGC, LSRdata$Ser4Per,
           main = "Serine 4 usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Serine 4 usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(LSRdata$Ser4Per ~ LSRdata$IntergenomicGC)))
legend(40, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(LSRdata$Leu4Per ~ LSRdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(LSRdata$IntergenomicGC, LSRdata$Leu4Per,
           main = "Leucine 4 usage over Intergenomic GC content",
           xlab = "Intergenomic GC content",
           ylab = "Leucine 4 usage",
           col = "#7C3194", 
           pch = 19,
           xlim = c(0,90),
           ylim = c(0,15),
           cex = 0.8,
           abline(lm(LSRdata$Leu4Per ~ LSRdata$IntergenomicGC)))
legend(5, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))


#amino acids slope over codon GC content                  

aaslope = read.csv("C:\\Users\\Valee\\Desktop\\FYP\\results\\aaslopefinal.csv", header=TRUE, sep = ",")
aminoacids = aaslope$AminoAcids
colour = aaslope$Colours
GC2 = aaslope$GC2codon
slope = aaslope$Slope

plot(GC2, slope, 
     main = "Relation between amino acids over intergenic GC slope and codon GC content at the first 2 positions",
     xlab = "Codon GC content at the first 2 positions",
     ylab = "Amino Acid usage over intergenic GC slope",
     abline(lm(aaslope$Slope ~ aaslope$GC2codon), col = "#431A50"))


sp <- ggplot(data = aaslope, aes(x = GC2, y = slope, color =colour),
             main = "Relation between amino acids over intergenic GC slope and codon GC content at the first 2 positions",
             xlab = "Codon GC content at the first 2 positions",
             ylab = "Amino Acid usage over intergenic GC slope",
             col = "#7C3194", 
             pch = 19,
             xlim = c(0,1),
             ylim = c(-0.2,0.2),
             cex = 0.8,
             abline(lm(aaslope$Slope ~ aaslope$GC2codon), col = "#431A50"))
            
legend(0.05, 0.2, legend=c(paste("Equation: y =", round(a, digits=3), ifelse(sign(b) == 1, "+", "-"), round(abs(b), digits=3), "· x"), 
                        paste("r2 =", rsq)))
sp+
   geom_label_repel(aes(label = aminoacids, col =colour), nudge_y = 0.005, nudge_x = 0.05)


regline <- lm(aaslope$Slope ~ aaslope$GC2codon)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)
regline
a

rlang::last_error()
