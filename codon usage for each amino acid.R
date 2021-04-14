library(datasets)
install.packages("pacman")
pacman::p_load(pacman, caret, lars, tidyverse, dplyr, GGally, ggplot2, ggthemes,ggvis, httr, plotly, shiny, tidyr)


GCdata = read.csv("C:\\Users\\Valee\\Desktop\\FYP\\Results\\GCdata.csv", header=TRUE, sep = ",")


regline <- lm(GCdata$PheGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$PheGC,
           main = "Phenylalanine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Phenylalanine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$PheGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$ThrGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$ThrGC,
           main = "Threonine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Threonine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$ThrGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$AlaGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$AlaGC,
           main = "Alanine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Alanine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$AlaGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 20, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$TyrGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$TyrGC,
           main = "Tyrosine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Tyrosine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$TyrGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$HisGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$HisGC,
           main = "Histidine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Histidine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$HisGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$GlnGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)
summary(regline)

sp <- plot(GCdata$IntergenomicGC, GCdata$GlnGC,
           main = "Glutamine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Glutamine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(-14.4,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$GlnGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))



regline <- lm(GCdata$AsnGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)
summary(regline)

sp <- plot(GCdata$IntergenomicGC, GCdata$AsnGC,
           main = "Asparagine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Asparagine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$AsnGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$LysGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$LysGC,
           main = "Lysine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Lysine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$LysGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$AspGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$AspGC,
           main = "Aspartic Acid GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Aspartic acid GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$AspGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))



regline <- lm(GCdata$GluGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$GluGC,
           main = "Glutamic Acid GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Glutamic acid GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$GluGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$CysGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$CysGC,
           main = "Cysteine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Cysteine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$CysGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$ValGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$ValGC,
           main = "Valine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Valine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$ValGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))



regline <- lm(GCdata$GlyGC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$GlyGC,
           main = "Glycine GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Glycine GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$GlyGC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$Leu2GC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$Leu2GC,
           main = "Leucine 2 GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Leucine 2 GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$Leu2GC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$Leu4GC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$Leu4GC,
           main = "Leucine 4 GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Leucine 4 GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$Leu4GC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$Ser2GC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$Ser2GC,
           main = "Serine 2 GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Serine 2 GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$Ser2GC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))



regline <- lm(GCdata$Ser4GC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$Ser4GC,
           main = "Serine 4 GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Serine 4 GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$Ser4GC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$Arg4GC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$Arg4GC,
           main = "Arginine 4 GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Arginine 4 GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$Arg4GC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))


regline <- lm(GCdata$Arg2GC ~ GCdata$IntergenomicGC)
a <- regline$coefficients[1]
b <- regline$coefficients[2]
rsq = format(summary(regline)$r.squared, digits = 3)

sp <- plot(GCdata$IntergenomicGC, GCdata$Arg2GC,
           main = "Arginine 2 GC3 codon percentage over intergenomic GC",
           xlab = "Intergenomic GC content",
           ylab = "Arginine 2 GC3 codon percentage",
           col = "#D577EE",
           xlim = c(0,90),
           ylim = c(0,100),
           pch = 19,
           cex = 0.8,
           abline(lm(GCdata$Arg2GC ~ GCdata$IntergenomicGC), col = "#431A50" )
)
legend(50, 15, legend=c(paste("Equation: y =", round(a, digits=0), ifelse(sign(b) == 1, "+", "-"), round(b, digits=3), "· x"), 
                        paste("r2 =", rsq)))
