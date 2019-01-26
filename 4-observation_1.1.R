###########################################################
# Observation
# 
# First created: 6/20/2018
#
# Contacts: Xin Jing <Xin.Jing@uvm.edu>
###########################################################

rm(list = ls())

# load library
library(TeachingDemos)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(cowplot)

char2seed("observations")

#################################################
# cleanSlopeDF 
# to clean the dataframe of slope estimate
# input: res.df, a matrix of slope stimate
# output: a neat dataframe of slope estimate
#------------------------------------------------
cleanSlopeDF <- function(x, ecosystem.types) {
  x <- data.frame(x)
  names(x) <- c("coef", "lci", "uci", "nfuncs", "ncombns")
  x$sig <- ifelse(x$lci > 0, "pos.sig", 
                  ifelse(x$uci < 0, "neg.sig", "neutral"))
  x$sig <- factor(x$sig, levels = c("pos.sig", "neutral", "neg.sig"))
  x$std <- rep(c("Averaging approach", "Summing approach", "Scaled by [0, 1]"),
               dim(x)[1] / 3)
  x$std <- factor(x$std,
                  levels = c("Averaging approach", "Summing approach", "Scaled by [0, 1]"),
                  labels = c("Averaging approach", "Summing approach", "Scaled by [0, 1]"))
  x$ecotypes <- rep(ecosystem.types, dim(x)[1])
  x$ecotypes <- factor(x$ecotypes)
  x <- data.frame(x)
}
#################################################
# getProbability
# calculate the probability of biodiversity with positive, neutral and negative effects on EMF
# input: res.df, a matrix of slope stimate
# output: get the probability
#------------------------------------------------
getProbability <- function(res.df, ecosystem.types) {
  df.nfuncs.sig <- ddply(res.df, .(std, nfuncs, sig), nrow)
  df.tot.obs <- ddply(df.nfuncs.sig, .(std, nfuncs), function(x) sum(x$V1))
  df.nfuncs.sig <- merge(df.nfuncs.sig, df.tot.obs, by = c("std", "nfuncs"))
  df.nfuncs.sig$prop <- with(df.nfuncs.sig, V1.x / V1.y)
  df.nfuncs.sig$ecotypes <- rep(ecosystem.types, dim(df.nfuncs.sig)[1])
  df.nfuncs.sig$ecotypes <- factor(df.nfuncs.sig$ecotypes)
  return(df.nfuncs.sig)
}
#################################################
# corrDescrip
# descriptive correlation matrics for ecosystem functions
# input: dat, a dataframe of ecosystem functions
#        func.names, names of ecosystem functions
# output: descriptive correlation matrics
#------------------------------------------------
corrDescrip <- function(dat, func.names) {
  cor.mat <- cor(dat[, func.names], use = "complete.obs")
  cat("Number of functions\n", length(func.names))
  cat("\nAverage correlation coefficients between ecosystem functions\n", mean(cor.mat[lower.tri(cor.mat)]))
  cat("\nRange correlation coefficients between ecosystem functions\n", range(cor.mat[lower.tri(cor.mat)]))
}

###########################################################
# Global drylands 
# Data from Maestre et al. 2012
# Full cite: Maestre, Fernando T., et al. "Plant species richness and ecosystem multifunctionality in global drylands." Science 335.6065 (2012): 214-218.

# load data
maestre.dat <- read.csv("./data/Maestre_Global_drylands_final_236_28_5_2018.csv")

# Z-scoring transformation for each function
vars4func.maestre <- c("TON", "BGL", "FOS", "P.HCL", "ORC", "AMO", "NIT",
               "AMI", "PRO", "PHE", "ARO", "HEX", "PEN", "NTR")
maestre.dat[, vars4func.maestre] <- apply(maestre.dat[, vars4func.maestre], 2, 
                                  function(x) {(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)})

# PCA analysis for climate
vars4clim <- c("AMT", "MDR", "TSE", "MAWM", "MTWEQ", "MTDQ",
               "RAI", "RASE", "RADQ", "RACQ")
clim.pca <- prcomp(maestre.dat[, vars4clim], scale = TRUE)
maestre.dat <- cbind(maestre.dat, clim.pca$x[, 1:4])

# Inspect the relationship between the number of functions and slope estimate
# using averaging, summing and scaled by [0, 1] approaches

res.df.drylands <- NULL
for (i in 1:length(vars4func.maestre)) {
  N <- dim(combn(vars4func.maestre, i))[2]  # using all combinations of functions
  if (N > 500) J <- 1:500 else J <- 1:N
  temp <- maestre.dat[, vars4func.maestre]
  for (j in J) {
    if (i == 1) {
      dat <- temp[, j]
      emf <- dat
    } else {
      dat <- sample(temp[vars4func.maestre], i)
      emf <- apply(dat, 1, function(x) {
        ret <- sum(x, na.rm = TRUE) / i
      })
    }
    
    for (l in c("av", "sm", "sc")) {
      if (l == "av") {
        emf <- emf
      } else if (l == "sm") {
        emf <- emf * i
      } else {
        emf <- (emf - min(emf, na.rm = TRUE)) / (max(emf, na.rm = TRUE) - min(emf, na.rm = TRUE))
      }
      df <- cbind(maestre.dat, emf)
      res <- matrix(NA, nrow = 1, ncol = 5)
      mod <- lm(emf ~ log(SR), data = df)
      res[, 1] <- coef(summary(mod))[2, 1]
      res[, 2] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
      res[, 3] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
      res[, 4] <- i
      res[, 5] <- j
      res.df.drylands <- rbind(res.df.drylands, res)
    }
  }
}

# inspcet BEMF-number of functions relationships
global.drylands <- cleanSlopeDF(res.df.drylands, ecosystem.types = "Global drylands")
# the proportion of functions with a significantly positive/neutral relationship with diversity
global.drylands.prob <- getProbability(global.drylands, ecosystem.types = "Global drylands")

###########################################################
# Tibetan grasslands
# Data from Jing et al. 2015
# Full cite: Jing, Xin, et al. "The links between ecosystem multifunctionality and above-and belowground biodiversity are mediated by climate." Nature communications 6 (2015).

# import data
tibet.dat <- read.csv("./data/jing_etal_Tibet.csv")

# data transformation
vars4func.tibet <- c("AGB", "PTN", "PTP", "BGB",
                    "SSOC", "SSTN", "TAN", "SSTP")
tibet.dat[, vars4func.tibet] <- apply(tibet.dat[, vars4func.tibet], 2, function(x) {
  x <- (x - mean(x, na.rm = T)) / sd(x, na.rm = TRUE)
})
vars4abiotic.tibet <- c("SM", "pH", "MAT")
tibet.dat[, vars4abiotic.tibet] <- apply(tibet.dat[, vars4abiotic.tibet], 2, function(x) {
  x <- (x - mean(x, na.rm = T)) / sd(x, na.rm = TRUE)
})

# for loops to calculate standardized coefficients
res.df.tibet <- NULL
for (i in 1:8) {
  N <- dim(combn(vars4func.tibet, i))[2]
  for (j in 1:N) {
    temp <- tibet.dat[vars4func.tibet]
    dat <- temp[, combn(vars4func.tibet, i)[, j]]
    dat <- data.frame(dat)
    emf <- apply(dat, 1, function(x) {
      ret <- sum(x) / i
    })
    for (l in c("av", "sm", "sc")) {
      if (l == "av") {
        emf <- emf
      } else if (l == "sm") {
        emf <- emf * i
      } else {
        emf <- (emf - min(emf)) / (max(emf) - min(emf))
      }
      df <- cbind(tibet.dat, emf)
      res <- matrix(NA, nrow = 1, ncol = 5)
      mod <- lm(emf ~ log(SR), data = df)
      res[, 1] <- coef(summary(mod))[2, 1]
      res[, 2] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
      res[, 3] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
      res[, 4] <- i
      res[, 5] <- j
      # cat(i, j, sep = "\n")
      res.df.tibet <- rbind(res.df.tibet, res)
    }
  }
}

# inspcet BEMF-number of functions relationships
tibetan.grasslands <- cleanSlopeDF(res.df.tibet, ecosystem.types = "Tibetan grasslands")
# the proportion of functions with a significantly positive/neutral relationship with diversity
tibetan.grasslands.prob <- getProbability(tibetan.grasslands, ecosystem.types = "Tibetan grasslands")

###########################################################
# BioDepth
# Data from Spehn et al. 2015
# Full cite: Spehn, E. M., et al. "Ecosystem effects of biodiversity manipulations in European grasslands." Ecological monographs 75.1 (2005): 37-63.
# import data
biodepth.dat <- read.csv("./data/all_biodepth.csv")
site.inf <- read.csv("./data/all_biodepth_site_vars.csv")

# data cleaning
biodepth.dat <- left_join(biodepth.dat, site.inf, by = "location")
biodepth.dat <- biodepth.dat %>% 
  mutate(fb = factor(NewBlock)) %>% 
  select(location, fb, latitude, longitude, altitude, 
         MAP, GST, pH, Diversity, biomassY3,
         root3, cotton3, wood3, N.Soil, N.g.m2, light3)
biodepth.dat[, -c(1:5, 9)] <- data.frame(apply(biodepth.dat[, -c(1:5, 9)], 2, function(x) {
  x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))
vars4func.biodepth <- c("biomassY3", "root3", 
               "N.Soil", "cotton3",
               "N.g.m2", "light3")

res.df.biodepth <- NULL
for (i in 1:6) {
  N <- dim(combn(vars4func.biodepth, i))[2]
  for (j in 1:N) {
    dat <- biodepth.dat[, combn(vars4func.biodepth, i)[, j]]
    dat <- data.frame(dat)
    emf <- apply(dat, 1, function(x) {
      ret <- sum(x, na.rm = TRUE) / i
    })
    for (l in c("av", "sm", "sc")) {
      if (l == "av") {
        emf <- emf
      } else if (l == "sm") {
        emf <- emf * i
      } else {
        emf <- (emf - min(emf, na.rm = TRUE)) / (max(emf, na.rm = TRUE) - min(emf, na.rm = TRUE))
      }
      df2 <- cbind(biodepth.dat, emf)
      res <- matrix(NA, nrow = 1, ncol = 5)
      mod <- lm(emf ~ log(Diversity),
                data = df2)
      res[, 1] <- coef(summary(mod))[2, 1]
      res[, 2] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
      res[, 3] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
      res[, 4] <- i
      res[, 5] <- j
      # cat(i, j, sep = "\n")
      res.df.biodepth <- rbind(res.df.biodepth, res)
    }
  }
}
# inspcet BEMF-number of functions relationships
biodepth.grasslands <- cleanSlopeDF(res.df.biodepth, ecosystem.types = "European grasslands")
# the proportion of functions with a significantly positive/neutral relationship with diversity
biodepth.grasslands.prob <- getProbability(biodepth.grasslands, ecosystem.types = "European grasslands")

###########################################################
# FunDivEUROPE
# Data from Ratcliffe et al. 2017
# Full cite: Ratcliffe, Sophia, et al. "Biodiversity and ecosystem functioning relations in European forests depend on environmental context." Ecology letters 20.11 (2017): 1414-1426.
# import data
fundiv.dat <- read.csv("./data/Ratcliffe_et_al_ELE_12849_ecosystem_function_variables.csv")
site.inf <- read.csv("./data/Ratcliffe_et_al_ELE_12849_regional_context_variables.csv")

# data cleaning
fundiv.dat <- left_join(fundiv.dat, site.inf, by = "region")
fundiv.dat[, 7:32] <- data.frame(apply(fundiv.dat[, 7:32], 2, function(x) {
  x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))
vars4func.fundiv <- names(fundiv.dat)[7:32]

# linear models
res.df.fundiv <- NULL
for (i in 1:26) {
  N <- dim(combn(vars4func.fundiv, i))[2]
  if (N > 500) J <- 1:500 else J <- 1:N
  k <- 0
  res1 <- numeric(length(J))
  res2 <- numeric(length(J))
  res3 <- numeric(length(J))
  res4 <- numeric(length(J))
  res5 <- numeric(length(J))

  for (j in J) {
    # make an empty data for each loop
    dat <- NULL
    emf <- NULL
    df2 <- NULL
    mod <- NULL

    # calculate EMF
    if (i == 1) {
      temp <- fundiv.dat[vars4func.fundiv]
      dat <- temp[, j]
      emf <- dat
    } else {
      dat <- sample(fundiv.dat[vars4func.fundiv], i)
      emf <- apply(dat, 1, function(x) {
        ret <- sum(x, na.rm = TRUE) / i
      })
    }
    
    for (l in c("av", "sm", "sc")) {
      if (l == "av") {
        emf <- emf
      } else if (l == "sm") {
        emf <- emf * i
      } else {
        emf <- (emf - min(emf, na.rm = TRUE)) / (max(emf, na.rm = TRUE) - min(emf, na.rm = TRUE))
      }
      df2 <- cbind(fundiv.dat, emf)
      
      # fit the model
      mod <- lm(emf ~ log(target_species_richness),
                data = df2)
      k <- k + 1
      res1[k] <- coef(summary(mod))[2, 1]
      res2[k] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
      res3[k] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
      res4[k] <- i
      res5[k] <- j
      # cat(i, j, sep = "\n")
    }
    res <- matrix(NA, nrow = 1, ncol = 5)
    res <- cbind(res1, res2, res3, res4, res5)
  }
  res.df.fundiv <- rbind(res.df.fundiv, res)
}

# inspcet BEMF-number of functions relationships
european.forests <- cleanSlopeDF(res.df.fundiv, ecosystem.types = "European forests")
# the proportion of functions with a significantly positive/neutral relationship with diversity
european.forests.prob <- getProbability(european.forests, ecosystem.types = "European forests")

###########################################################
# Jena grassland biodiversity experiment
# Data from Meyer et al. 2018
# Full cite: Meyer, Sebastian T., et al. "Biodiversity–multifunctionality relationships depend on identity and number of measured functions." Nature ecology & evolution 2.1 (2018): 44.
# import data
jena.dat <- read.csv("./data/Meyer_etal_NEE_Jena_exp.csv")
jena.funcs <- read.csv("./data/Meyer_etal_NEE_Jena_exp_funcs.csv")

# data cleaning
jena.dat <- left_join(jena.dat, jena.funcs, by = "plotcode")

# z-scoring transformation
jena.dat[, 15:96] <- data.frame(apply(jena.dat[, 15:96], 2, function(x) {
  x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))
vars4func.jena <- names(jena.dat)[15:96]

# Inspect the diversity effect and number of ecosystem functions
res.df.jena <- NULL
for (i in 1:82) {
  ifelse(i == 1, J <- 1:82, 
         ifelse(i == 81, J <- 1:82, 
                ifelse(i == 82, J <- 1, J <- 1:500)))
  k <- 0
  res1 <- res2 <- res3 <- res4 <- res5 <- res6 <- numeric(length(J)*3)
  
  for (j in J) {
    # make an empty data for each loop
    dat <- NULL
    emf <- NULL
    df2 <- NULL
    mod <- NULL
    if (i == 1) {
      temp <- jena.dat[vars4func.jena]
      dat <- temp[, j]
      emf <- dat
    } else {
      if (i == 81) {
        temp <- jena.dat[vars4func.jena]
        dat <- temp[, -j]
        emf <- apply(dat, 1, function(x) {
          ret <- sum(x, na.rm = TRUE) / 81  # average
        })
      } else {
        dat <- sample(jena.dat[vars4func.jena], i)
        emf <- apply(dat, 1, function(x) {
          ret <- sum(x, na.rm = TRUE) / i
        })
      }
    }
    for (l in c("av", "sm", "sc")) {
      if (l == "av") {
        emf <- emf
      } else if (l == "sm") {
        emf <- emf * i
      } else {
        emf <- (emf - min(emf)) / (max(emf) - min(emf))
      }
      df2 <- cbind(jena.dat, emf)
      
      # fit the model
      mod <- lm(emf ~ log(sowndiv), 
                data = df2)
      k <- k + 1
      res1[k] <- coef(summary(mod))[2, 1]
      res2[k] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
      res3[k] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
      res4[k] <- i
      res5[k] <- j
      # cat(i, j, l, sep = "\n")
    }
    res <- matrix(NA, ncol = 5, nrow = length(J)*3)
    res <- cbind(res1, res2, res3, res4, res5)
  }
  res.df.jena <- rbind(res.df.jena, res)
}
# inspcet BEMF-number of functions relationships
jena.grassland <- cleanSlopeDF(res.df.jena, ecosystem.types = "Jena grassland")
# the proportion of functions with a significantly positive/neutral relationship with diversity
jena.grassland.prob <- getProbability(jena.grassland, ecosystem.types = "Jena grassland")

###########################################################
# Plot
res.df.all <- rbind(jena.grassland, european.forests, 
                    biodepth.grasslands, tibetan.grasslands, global.drylands)
res.df.all$std <- factor(res.df.all$std, 
                         levels = c("Averaging approach", "Summing approach", "Scaled by [0, 1]"),
                         labels = c("Averaging approach", "Summing approach", "Scaling approach"))

ggplot(res.df.all, aes(x = nfuncs, y = coef)) +
  geom_jitter(width = 0.2, size = 0.8, 
              aes(color = sig, shape = sig)) +
  geom_smooth(method = "lm", lwd = 0.3, color = "black", se = FALSE) +
  scale_color_manual(values = c("#fb6a4a", "gray", "#4393c3")) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, col = "grey") +
  labs(x = "Number of functions considered", 
       y = "Biodiversity effect on EMF") +
  facet_grid(std ~ ecotypes, scales = "free") +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12.5),
        legend.position = 'none')
ggsave("./outputs/observations.pdf", width = 10, height = 5.6)
# ggsave("./outputs/observations.pdf", width = 11.55, height = 6.29)

# probablity plotting
res.df.all.prob <- rbind(jena.grassland.prob, european.forests.prob, 
                         biodepth.grasslands.prob, tibetan.grasslands.prob, 
                         global.drylands.prob)
ggplot(res.df.all.prob, aes(x = nfuncs, y = prop, color = sig)) +
  geom_point(aes(shape = sig), size = 1.5) +
  geom_smooth(method = "glm",
              # method.args = list(family = "binomial"),
              se = FALSE, size = 0.25) +
  scale_color_manual(values = c("#fb6a4a", "gray", "#4393c3")) +
  labs(x = "Number of functions considered", 
       y = "Probability of positive, neutral or negative\nbiodiversity effects on EMF") +
  facet_grid(std ~ ecotypes, scales = "free") +
  theme_bw(base_size = 16.5) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')
ggsave("./outputs/observations_probability.pdf", width = 11.65, height = 6.29)

###########################################################
# inspect biodiversity and EMF relationships
res.df.all$ecotypes <- factor(res.df.all$ecotypes,
                              levels = c("Jena grassland", "European forests",
                                         "European grasslands", "Tibetan grasslands",
                                         "Global drylands"))
# individual functions
p1 <- res.df.all %>% 
  filter(nfuncs < 2) %>% 
  filter(std == "Scaling approach") %>% 
  ggplot(aes(x = ecotypes, y = coef, fill = ecotypes, alpha = 0.66)) +
  geom_boxplot(coef = 1e30) +
  geom_jitter(width = 0.2, size = 1.5, shape = 1) +
  # scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")) +
  labs(x = "", y = "Biodiversity effect on\n individal ecosystem functions") +
  theme_classic(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# multifunctionality
p2 <- res.df.all %>% 
  filter(nfuncs != 1) %>%
  filter(std == "Scaling approach") %>% 
  ggplot(aes(coef, fill = ecotypes)) +
  geom_density(alpha = 0.66) +
  labs(x = "Biodiversity effect on EMF", y = "Density\n") +
  # scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")) +
  theme_classic(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.20, 0.65),
        legend.title = element_blank())

pdf("./outputs/boxplot_density.pdf", width = 6.8, height = 8.5)
plot_grid(p1, p2, nrow = 2, labels = c("a)", "b)"),
          rel_widths = c(1, 1))
dev.off()

###########################################################
# average correlation coefficients between ecosystem functions

# Jena grassland
corrDescrip(jena.dat, vars4func.jena)
# European forests
corrDescrip(fundiv.dat, vars4func.fundiv)
# European grasslands
corrDescrip(biodepth.dat, vars4func.biodepth)
# Tibetan grasslands
corrDescrip(tibet.dat, vars4func.tibet)
# Global drylands
corrDescrip(maestre.dat, vars4func.maestre)

###########################################################
#                   End of Script                         #
###########################################################