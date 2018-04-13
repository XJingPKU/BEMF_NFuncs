###########################################################
# The relationship between BEMF relationship and the number of functions
# 
# Data from the Jena project 
# First created: 11/20/2017
#
# Contacts: Xin Jing <Xin.Jing@uvm.edu>
###########################################################
options(useFancyQuotes = FALSE)
rm(list = ls())

# load library
library(plyr)
library(dplyr)
library(ggplot2)
library(TeachingDemos)

# create a seed for the random number generator
char2seed("randomsubset")

# import data
jena.exp <- read.csv("./data/Meyer_etal_NEE_Jena_exp.csv")
jena.funcs <- read.csv("./data/Meyer_etal_NEE_Jena_exp_funcs.csv")

# data cleaning
jena <- left_join(jena.exp, jena.funcs, by = "plotcode")

df <- jena
df[, 15:96] <- data.frame(apply(df[, 15:96], 2, function(x) {
  # x <- log(x + 1)
  x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

###########################################################
# using random subset

vars4func <- names(df)[15:96]
emf.dat <- df[vars4func]
names(emf.dat) <- paste("EF", 1:82, sep = "")
vars4func <- names(emf.dat)

###########################################################
# FUNCTION: getMatric
# calculate metric from random subset
# input: block, sowndiv, emf.dat and random subset (10 or 47)
# output: regression slope
#----------------------------------------------------------
getMetric <- function(block, sowndiv, emf.dat, N.random) {
  if (N.random <= 10) {
    dat <- emf.dat[, sample(seq_along(names(emf.dat)), N.random)]
    res.df <- NULL
    for (i in 1:N.random) {
      N <- dim(combn(N.random, i))[2]
      for (j in 1:N) {
        temp <- dat[, combn(N.random, i)[, j]]
        temp <- data.frame(temp)
        emf <- apply(temp, 1, function(x) {
          ret <- sum(x, na.rm = TRUE) / i
        })
        df2 <- cbind(block, sowndiv, emf)
        df2 <- data.frame(df2)
        df2$block <- factor(df2$block)
        mod <- lm(emf ~ block + log(sowndiv), data = df2)
        res <- matrix(NA, nrow = 1, ncol = 5)
        res[, 1] <- coef(summary(mod))[5, 1]
        res[, 2] <- coef(summary(mod))[5, 1] - 1.96 * coef(summary(mod))[5, 2]
        res[, 3] <- coef(summary(mod))[5, 1] + 1.96 * coef(summary(mod))[5, 2]
        res[, 4] <- i
        res[, 5] <- j
        cat(i, j, sep = "\n")
        res.df <- rbind(res.df, res)
      }
    }
    res.df <- data.frame(res.df)
    names(res.df) <- c("coef", "lci", "uci", "nfuncs", "ncombns")
    res.df$sig <- ifelse((res.df$lci >0), "pos.sig", 
                         ifelse(res.df$uci < 0, "neg.sig", "neutral"))
    res.df$sig <- factor(res.df$sig,
                         levels = c("pos.sig", "neutral", "neg.sig"))
    # mod.av <- lm(coef ~ nfuncs, data = res.df)
    # slope.av <- summary(mod.av)$coefficients[2, 1]
    if (length(res.df[res.df$sig == "pos.sig", ]$coef) <= 5) {
      slope.ps <- NA
    } else {
      mod.ps <- lm(coef ~ nfuncs, data = res.df[res.df$sig == "pos.sig", ])
      slope.ps <- summary(mod.ps)$coefficients[2, 1]
    }
    # return(list(av = slope.av, ps = slope.ps))
    return(list(ps = slope.ps))
  } else {
    dat <- emf.dat[, sample(seq_along(names(emf.dat)), N.random)]
    res.df <- NULL
    for (i in 1:N.random) {
      ifelse(i == 1, J <- 1:N.random, 
             ifelse(i == (N.random -1), J <- 1:N.random, 
                    ifelse(i == N.random, J <- 1, J <- 1:200)))
      k <- 0
      res1 <- numeric(length(J))
      res2 <- numeric(length(J))
      res3 <- numeric(length(J))
      res4 <- numeric(length(J))
      res5 <- numeric(length(J))
      for (j in J) {
        if (i == 1) {
          temp <- dat[, j]
          emf <- temp
        } else {
          if (i == (N.random - 1)) {
            temp <- dat[, -j]
            emf <- apply(temp, 1, function (x) {
              ret <- sum(x, na.rm = TRUE) / (N.random -1)
            })
          } else {
            temp <- sample(dat, i)
            emf <- apply(temp, 1, function(x) {
              ret <- sum(x, na.rm = TRUE) / i
            })
          }
        }
        df2 <- cbind(block, sowndiv, emf)
        df2 <- data.frame(df2)
        df2$block <- factor(df2$block)
        mod <- lm(emf ~ block + log(sowndiv), data = df2)
        k <- k + 1
        res1[k] <- coef(summary(mod))[5, 1]
        res2[k] <- coef(summary(mod))[5, 1] - 1.96 * coef(summary(mod))[5, 2]
        # res3[k] <- coef(summary(mod))[5, 1] + 1.96 * coef(summary(mod))[5, 2]
        res3[k] <- NA
        res4[k] <- i
        res5[k] <- j
        cat(i, j, sep = "\n")
      }
      res <- matrix(NA, nrow = 1, ncol = 5)
      res <- cbind(res1, res2, res3, res4, res5)
      res.df <- rbind(res.df, res)
    }
    res.df <- data.frame(res.df)
    names(res.df) <- c("coef", "lci", "uci", "nfuncs", "ncombns")
    res.df$sig <- ifelse((res.df$lci >0), "pos.sig", "no.pos.sig")
    res.df$sig <- factor(res.df$sig)
    # mod.av <- lm(coef ~ nfuncs, data = res.df)
    # slope.av <- summary(mod.av)$coefficients[2, 1]
    if (length(res.df[res.df$sig == "pos.sig", ]$coef) <= 5) {
      slope.ps <- NA
    } else {
      mod.ps <- lm(coef ~ nfuncs, data = res.df[res.df$sig == "pos.sig", ])
      slope.ps <- summary(mod.ps)$coefficients[2, 1]
    }
    # return(list(av = slope.av, ps = slope.ps))
    return(list(ps = slope.ps))
  }
}

###########################################################
# 5, 8, 10, 15, 20, 25, 30 random functions
getMetricRand <- function(nSim, N.random) {
  nSim <- nSim
  xSim <- rep(NA, nSim)
  for (i in seq_len(nSim)) {
    xSim[i] <- getMetric(block = df$block,
                          sowndiv = df$sowndiv,
                          emf.dat = emf.dat,
                          N.random = N.random)[1]
  }
  xSim <- unlist(xSim)
  xSim <- cbind(xSim, rep(N.random, nSim))
  colnames(xSim) <- c("xSim", "no.funcs")
  return(xSim)
}

sim5 <- getMetricRand(nSim = 200, N.random = 5)
sim8 <- getMetricRand(nSim = 200, N.random = 8)
sim10 <- getMetricRand(nSim = 200, N.random = 10)
sim15 <- getMetricRand(nSim = 200, N.random = 15)
sim20 <- getMetricRand(nSim = 200, N.random = 20)
sim25 <- getMetricRand(nSim = 200, N.random = 25)
sim30 <- getMetricRand(nSim = 200, N.random = 30)
sim35 <- getMetricRand(nSim = 200, N.random = 35)
sim40 <- getMetricRand(nSim = 200, N.random = 40)
sim45 <- getMetricRand(nSim = 200, N.random = 45)
sim54 <- getMetricRand(nSim = 200, N.random = 54)

df.sim <- rbind(sim5, sim8, sim10, sim15, sim20,
                sim25, sim30, sim35, sim40, sim45, sim54)
df.sim <- data.frame(df.sim)
# df.sim$xSim <- as.numeric(as.character(df.sim$sim10))
write.csv(df.sim, "./outputs/Jena_random_subset_nofuncs.csv")

# load the final data
df.sim2 <- read.csv("./outputs/Jena_random_subset_nofuncs_V1.csv")
# df.sim2$no.funcs <- factor(df.sim2$no.funcs)
ggplot(df.sim2, aes(x = no.funcs, y = xSim, group = no.funcs)) +
  geom_boxplot(outlier.shape = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_x_continuous(breaks = c(seq(0, 60, 10))) +
  ylim(c(-0.18, 0.01)) +
  labs(x = "Number of functions", y = "Slope estimate") +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank())
ggsave("outputs/Jena_randomization_test.pdf", height = 4.5, width = 6)

###########################################################