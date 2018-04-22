################
##
## @description Create plots using the Summary data from ABM simulations set
##              11,000 iterations.
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
# load summaryAIC 11,000
summaryAIC = read.csv(paste0(inputDir, "/summaryAIC11000.csv"),
    sep=";")

# load summarySC 11,000
summarySC = read.csv(paste0(inputDir, "/summarySC11000.csv"),
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
  if (summary$tModel[i] == 'Y ~ X1 + X2'){
    w[i] = 1
  } else if (summary$tModel[i] == 'Y ~ X1 + X2 + X3 + X1:X2'){
    w[i] = 2
  } else if (summary$tModel[i] == 'Y ~ X1 + X2 + X3 + X1:X2 + X1:X3 + X2:X3'){
    w[i] = 3
  }
}

summary$TModel = w
#---------------------------------------------------
# clean up

summary = summary[c("replica", "TModel", "sigmaIndex",
                    "Scientist", "selTrueModel", "departTrueModel",
                    "replicated","firstTGM", "replicatedTGM",
                    "replicatedNTGM", "method")]

names(summary) <- c("replica","TrueModel", "Error", "Population",
    "PTMGM", "Stickiness", "ReproducibilityRate", "FTTMGM",
    "ReproducibilityRateTMGM", "ReproducibilityRateTMNGM", "Method")

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
# Figure 9
p4 <- ggplot(summary, aes(x = Population, y = FTTMGM, fill = TrueModel)) +
    geom_violin(position = dodge, width = 1)
p14 <- p4 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31", position = dodge) +
    theme(legend.position = "none")

p5 <- ggplot(summary, aes(x = Population, y = FTTMGM, fill = Error)) +
    geom_violin(position = dodge, width = 1)
p15 <- p5 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31", position = dodge) +
    theme(legend.position = "none")

p <- grid.arrange(p14, p15, ncol = 2, nrow = 1)

ggsave(filename = paste0(outputDir, "/Fig9.pdf"), plot = p,
    width = 15, height = 8, units = "in")

#--------------------------------------
# Figure 10A
p1 <- ggplot(summary, aes(x = Population, y = FTTMGM, color = Population)) +
    geom_violin(position = dodge, width = 1) +
    scale_color_manual(values = c("#797878" ,"#FF5733" ,"#5DA9FE",
            "#FAEC3F" ,"#8F40FA"))
p1 <- p1 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31")

ggsave(filename = paste0(outputDir, "/Fig10A.pdf"), plot = p1,
    width = 7.5, height = 4, units = "in")

#--------------------------------------
# Figure 12A
r2 <- ggplot(summary, aes(x = TrueModel, y = FTTMGM, color = TrueModel)) +
    geom_violin(position = dodge, width = 1)
r2 <- r2 + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange", color = "gray31")

ggsave(filename = paste0(outputDir, "/Fig12A.pdf"), plot = r2,
    width = 7.5, height = 4, units = "in")

#--------------------------------------
# SCATTER PLOT WITH MARGINAL GROUP DENSITIES
#--------------------------------------
# Figure 8B
#--------------------------------------
mycol <- c("#797878" ,"#FF5733" ,"#5DA9FE", "#FAEC3F" ,"#8F40FA")
sp2 <- ggscatter(summary, x = "FTTMGM", y = "ReproducibilityRate",
        color = "Population", palette = get_palette(mycol, 5),
        size = 3, alpha = 0.6) +
    border()
xplot <- ggdensity(summary, "PTMGM", fill = "Population",
    palette = get_palette(mycol, 5))
yplot <- ggdensity(summary, "ReproducibilityRate", fill = "Population",
        palette = get_palette(mycol, 5)) +
    rotate()
p <- plot_grid(xplot, NULL, sp2, yplot, ncol = 2, align = "hv",
    rel_widths = c(2, 1), rel_heights = c(1, 2))

ggsave(filename = paste0(outputDir, "/Fig8B.pdf"), plot = p,
    width = 7.5, height = 4, units = "in")
#--------------------------------------
# end plots

#--------------------------------------
# begin summary statistics

#--------------------------------------
#SUMMARY STATISTICS PER POPULATION
#--------------------------------------
print(c("SUMMARY STATISTICS PER POPULATION (Rey, Tess, Mave, Bo, All)"))
print(c("FTTMGM"))
y = matrix(0, nrow = nrow(summary) / 5, ncol = 5)
ind = (summary$Population == "Rey")
y[,1] = summary$FTTMGM[ind]
ind = (summary$Population == "Tess")
y[,2] = summary$FTTMGM[ind]
ind = (summary$Population == "Mave")
y[,3] = summary$FTTMGM[ind]
ind = (summary$Population == "Bo")
y[,4] = summary$FTTMGM[ind]
ind = (summary$Population == "All")
y[,5] = summary$FTTMGM[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------
# SUMMARY STATISTICS PER TRUE MODEL
#--------------------------------------
print(c("SUMMARY STATISTICS PER TRUE MODEL (Y = X1 + X2, Y = X1 + X2 + X3 + X1:X2, Y = X1 + X2 + X3 + X1:X2 + X1:X3 + X2:X3)"))
print(c("FTTMGM"))
y = matrix(0, nrow = nrow(summary) / 3, ncol = 3)
ind = (summary$TrueModel == "1")
y[,1] = summary$FTTMGM[ind]
ind = (summary$TrueModel == "2")
y[,2] = summary$FTTMGM[ind]
ind = (summary$TrueModel == "3")
y[,3] = summary$FTTMGM[ind]

apply(y, 2, median, na.rm = T)
apply(y, 2, IQR, na.rm = T)
#--------------------------------------

# end summary statistics

#--------------------------------------
# end all