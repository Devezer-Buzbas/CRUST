################
##
## @description Create plots using the Summary data from ABM simulations set
##              11,000 iterations and 1,000 iterations set as burn-in.
##
## @param None
##
## @return None
##
## @lastChange 2018-04-22
##
## @changes
##
################
# begin all
#-------------------------------------------------------------
#
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggthemes)
library(Hmisc)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra)


#############
## PATHS
#############
# Edit baseDir to the base directory of CRUST
baseDir = "."
inputDir = paste0(baseDir, "/data/summary")
outputDir = paste0(baseDir, "/data/plot")


#############
## LOAD SUMMARY DATA FILES
#############
# load summaryAIC 10,000
summaryAIC = read.csv(paste0(inputDir, "/summaryAIC10000.csv"),
    sep=";")

# load summarySC 10,000
summarySC = read.csv(paste0(inputDir, "/summarySC10000.csv"),
    sep=";")

# bind AIC and SC summary files and prepare for analysis
summary = rbind(summaryAIC, summarySC)
method = c(rep("AIC", nrow(summaryAIC)), rep("SC", nrow(summarySC)))
summary = cbind(summary, method)

#---------------------------------------------------
# begin Rearrange Summary Data

# Create categorical variable for scientist types
y = c()
for(i in 1:nrow(summary)){
  if (summary$nRey[i] == 300){
    y[i] = "Rey"
  } else if (summary$nTess[i] == 300){
    y[i] = "Tess"
  } else if (summary$nMave[i] == 300){
    y[i] = "Mave"
  } else if (summary$nBo[i] == 300){
    y[i] = "Bo"
  } else y[i] = "All"
}

summary$Scientist = y
#---------------------------------------------------

# Create ordinal variable for true model complexity

w = c()
for(i in 1:nrow(summary)){
  if (summary$tModel[i]=='Y ~ X1 + X2'){
    w[i]=1
  } else if (summary$tModel[i]=='Y ~ X1 + X2 + X3 + X1:X2'){
    w[i]=2
  } else if (summary$tModel[i]=='Y ~ X1 + X2 + X3 + X1:X2 + X1:X3 + X2:X3'){
    w[i]=3
  }
}

summary$TModel = w
#---------------------------------------------------
# clean up
summary = summary[c("replica", "TModel", "sigmaIndex", "Scientist",
        "selTrueModel", "departTrueModel", "replicated","firstTGM",
        "replicatedTGM", "replicatedNTGM", "method")]

names(summary) = c("replica","TrueModel", "Error", "Population", "PTMGM",
    "Stickiness", "ReproducibilityRate", "FTTMGM", "ReproducibilityRateTMGM",
    "ReproducibilityRateTMNGM", "Method")

 summary$TrueModel = as.factor(summary$TrueModel)
 summary$Error = as.factor(summary$Error)
 summary$Population = as.factor(summary$Population)
 summary$Method = as.factor(summary$Method)
 
# end Rearrange Summary Data
#--------------------------------------
# begin plots
 
dodge = position_dodge(width = 0.6)
mycol = c("#797878" ,"#FF5733" ,"#5DA9FE", "#FAEC3F" ,"#8F40FA")
#--------------------------------------

# Violin Plots
#--------------------------------------
# Figure 8C

p1 <- ggplot(summary, aes(x = Population, y = ReproducibilityRate,
            color = Population)) +
    geom_violin(position = dodge, width = 1) +
    scale_y_continuous(limits = c(-0.01, 1.1), breaks=c(0,1)) +
    scale_color_manual(values = c("#797878" ,"#FF5733" ,"#5DA9FE",
            "#FAEC3F" ,"#8F40FA"))
p11 <- p1 + stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),
    geom="pointrange", color="gray31")  

p2 <- ggplot(summary, aes(x = Population, y = ReproducibilityRateTMGM,
        color = Population))+
    geom_violin(position = dodge, width = 1) +
    scale_y_continuous(limits = c(-0.01, 1.1), breaks = c(0,1)) +
    scale_color_manual(values = c("#797878" ,"#FF5733" , "#5DA9FE",
            "#FAEC3F" ,"#8F40FA"))
p12 <- p2 + stat_summary(fun.data = mean_sdl, fun.args = list(mult=1),
    geom="pointrange", color="gray31")

p3 <- ggplot(summary, aes(x = Population, y = ReproducibilityRateTMNGM,
        color = Population)) +
    geom_violin(position = dodge, width = 1) +
    scale_y_continuous(limits = c(-0.01, 1.1), breaks = c(0, 1)) +
    scale_color_manual(values = c("#797878" ,"#FF5733" ,"#5DA9FE",
            "#FAEC3F" ,"#8F40FA"))
p13 <- p3 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31")

p = grid.arrange(p11, p12, p13, ncol = 3, nrow = 1)

ggsave(filename = paste0(outputDir, "/Fig8C.pdf"), plot = p,
    width = 15, height = 8, units = "in")

#--------------------------------------
# Figure 10B

p4 <- ggplot(summary, aes(x = Population, y = PTMGM, color = Population)) +
    geom_violin(position = dodge, width = 1) +
    scale_y_continuous(limits = c(-0.01, 1.1), breaks = c(0, 1)) +
    scale_color_manual(values = c("#797878" ,"#FF5733" ,"#5DA9FE",
            "#FAEC3F" ,"#8F40FA"))
p4 <- p4 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31")

ggsave(filename = paste0(outputDir, "/Fig10B.pdf"), plot = p4,
    width = 7.5, height = 4, units = "in")

#--------------------------------------
# Figure 10C

p5 <- ggplot(summary, aes(x = Population, y = Stickiness, color = Population)) +
    geom_violin(position = dodge, width = 1) +
    scale_y_continuous(limits = c(-0.01, 1.1), breaks = c(0, 1)) +
    scale_color_manual(values = c("#797878", "#FF5733", "#5DA9FE",
            "#FAEC3F" ,"#8F40FA"))
p5 <- p5 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31")

ggsave(filename = paste0(outputDir, "/Fig10C.pdf"), plot = p5,
    width = 7.5, height = 4, units = "in")

#--------------------------------------
# Figure 11

p2 <- ggplot(summary, aes(x = Error, y = PTMGM, fill = Method)) +
    geom_violin(position = dodge, width = 1, color='black') +
    scale_fill_brewer(palette="Spectral")
MethodByError <- p2 + stat_summary(fun.data = mean_sdl,
    fun.args = list(mult = 1), geom = "pointrange", color = "gray31",
    position = dodge)

p3 <- ggplot(summary, aes(x = TrueModel, y = PTMGM, fill = Method)) +
    geom_violin(position = dodge, width = 1, color = 'black') +
    scale_fill_brewer(palette = "Spectral")
MethodByTrueModel <- p3 + stat_summary(fun.data = mean_sdl,
    fun.args = list(mult = 1), geom = "pointrange", color = "gray31",
    position = dodge)

p <- grid.arrange(MethodByTrueModel, MethodByError, ncol = 1, nrow = 2)

ggsave(filename = paste0(outputDir, "/Fig11.pdf"), plot = p,
    width = 15, height = 8, units = "in")

#--------------------------------------
# Figure 12B, 12C, 12D

r1 <- ggplot(summary, aes(x = TrueModel, y = PTMGM, color = TrueModel)) +
    geom_violin(position = dodge, width = 1) +
    scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 1))
r1 <- r1 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31")

ggsave(filename = paste0(outputDir, "/Fig12B.pdf"), plot = r1,
    width = 7.5, height = 4, units = "in")

r2 <- ggplot(summary, aes(x = TrueModel, y = Stickiness, color = TrueModel)) +
    geom_violin(position = dodge, width = 1) +
    scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 1))
r2 <- r2 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31")  

ggsave(filename = paste0(outputDir, "/Fig12C.pdf"), plot = r2,
    width = 7.5, height = 4, units = "in")

r3 <- ggplot(summary, aes(x = TrueModel, y = ReproducibilityRateTMGM,
        color=TrueModel)) +
    geom_violin(position = dodge, width = 1) +
    scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 1))
r3 <- r3 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31")

ggsave(filename = paste0(outputDir, "/Fig12D.pdf"), plot = r3,
    width = 7.5, height = 4, units = "in")

#--------------------------------------
#SCATTER PLOT WITH MARGINAL GROUP DENSITIES
#--------------------------------------
# Figure 8A
#--------------------------------------
sp1 <- ggscatter(summary, x = "PTMGM", y = "ReproducibilityRate",
        color = "Population", palette = get_palette(mycol,5),
        size = 3, alpha = 0.6) +
    theme(legend.position = "none") +
    border()

xplot <- ggdensity(summary, "PTMGM", fill = "Population",
                   palette = get_palette(mycol, 5)) +
               theme(legend.position = "none")
yplot <- ggdensity(summary, "ReproducibilityRate", fill = "Population",
                   palette = get_palette(mycol, 5)) +
               theme(legend.position = "none") +
               rotate()
p <- plot_grid(xplot, NULL, sp1, yplot, ncol = 2, align = "hv",
    rel_widths = c(2, 1), rel_heights = c(1, 2))

ggsave(filename = paste0(outputDir, "/Fig8A.pdf"), plot = p,
    width = 7.5, height = 4, units = "in")

# end plots

#--------------------------------------
# begin summary statistics
#--------------------------------------
# SUMMARY STATISTICS PER POPULATION
#--------------------------------------
print(c("SUMMARY STATISTICS PER POPULATION (Rey, Tess, Mave, Bo, All)"))
print(c("PTMGM"))
y = matrix(0, nrow = nrow(summary) / 5, ncol = 5)
ind = (summary$Population == "Rey")
y[,1] = summary$PTMGM[ind]
ind = (summary$Population == "Tess")
y[,2] = summary$PTMGM[ind]
ind = (summary$Population == "Mave")
y[,3] = summary$PTMGM[ind]
ind = (summary$Population == "Bo")
y[,4] = summary$PTMGM[ind]
ind = (summary$Population == "All")
y[,5] = summary$PTMGM[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------
print(c("Stickiness"))
y = matrix(0, nrow = nrow(summary) / 5, ncol = 5)
ind = (summary$Population == "Rey")
y[,1] = summary$Stickiness[ind]
ind = (summary$Population == "Tess")
y[,2] = summary$Stickiness[ind]
ind = (summary$Population == "Mave")
y[,3] = summary$Stickiness[ind]
ind = (summary$Population == "Bo")
y[,4] = summary$Stickiness[ind]
ind = (summary$Population == "All")
y[,5] = summary$Stickiness[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------
print(c("ReproducibilityRate"))
y = matrix(0, nrow = nrow(summary) / 5, ncol = 5)
ind = (summary$Population == "Rey")
y[,1] = summary$ReproducibilityRate[ind]
ind = (summary$Population == "Tess")
y[,2] = summary$ReproducibilityRate[ind]
ind = (summary$Population == "Mave")
y[,3] = summary$ReproducibilityRate[ind]
ind = (summary$Population == "Bo")
y[,4] = summary$ReproducibilityRate[ind]
ind = (summary$Population == "All")
y[,5] = summary$ReproducibilityRate[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------
print(c("ReproducibilityRateTMGM"))
y = matrix(0, nrow = nrow(summary) / 5, ncol = 5)
ind = (summary$Population == "Rey")
y[,1] = summary$ReproducibilityRateTMGM[ind]
ind = (summary$Population == "Tess")
y[,2] = summary$ReproducibilityRateTMGM[ind]
ind = (summary$Population == "Mave")
y[,3] = summary$ReproducibilityRateTMGM[ind]
ind = (summary$Population == "Bo")
y[,4] = summary$ReproducibilityRateTMGM[ind]
ind = (summary$Population == "All")
y[,5] = summary$ReproducibilityRateTMGM[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------
# SUMMARY STATISTICS PER TRUE MODEL
#--------------------------------------
print(c("SUMMARY STATISTICS PER TRUE MODEL (Y = X1 + X2, Y = X1 + X2 + X3 + X1:X2, Y = X1 + X2 + X3 + X1:X2 + X1:X3 + X2:X3)"))
print(c("PTMGM"))
y = matrix(0, nrow = nrow(summary) / 3, ncol = 3)
ind = (summary$TrueModel == "1")
y[,1] = summary$PTMGM[ind]
ind = (summary$TrueModel == "2")
y[,2] = summary$PTMGM[ind]
ind = (summary$TrueModel == "3")
y[,3] = summary$PTMGM[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------
print(c("Stickiness"))
y = matrix(0, nrow = nrow(summary) / 3, ncol = 3)
ind = (summary$TrueModel == "1")
y[,1] = summary$Stickiness[ind]
ind = (summary$TrueModel == "2")
y[,2] = summary$Stickiness[ind]
ind = (summary$TrueModel == "3")
y[,3] = summary$Stickiness[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------
print(c("ReproducibilityRate"))
y = matrix(0, nrow = nrow(summary) / 3, ncol = 3)
ind = (summary$TrueModel == "1")
y[,1] = summary$ReproducibilityRate[ind]
ind = (summary$TrueModel == "2")
y[,2] = summary$ReproducibilityRate[ind]
ind = (summary$TrueModel == "3")
y[,3] = summary$ReproducibilityRate[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------
print(c("ReproducibilityRateTMGM"))
y = matrix(0, nrow = nrow(summary) / 3, ncol = 3)
ind = (summary$TrueModel == "1")
y[,1] = summary$ReproducibilityRateTMGM[ind]
ind = (summary$TrueModel == "2")
y[,2] = summary$ReproducibilityRateTMGM[ind]
ind = (summary$TrueModel == "3")
y[,3] = summary$ReproducibilityRateTMGM[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)

#--------------------------------------
# SUMMARY STATISTICS PER MODEL SEL STAT
#--------------------------------------
print(c("SUMMARY STATISTICS PER MODEL SEL STAT (AIC, SC)"))
print(c("PTMGM"))
y = matrix(0, nrow = nrow(summary) / 2, ncol = 2)
ind = (summary$Method == "AIC")
y[,1] = summary$PTMGM[ind]
ind = (summary$Method == "SC")
y[,2] = summary$PTMGM[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------

# end summary statistics

#--------------------------------------
# end all