###########################################################
# The relationship between BEMF relationship and the number of functions
# 
# Data are from published and publicly available papers
# First created: 10/25/2017
#
# Contacts: Xin Jing <Xin.Jing@uvm.edu>
###########################################################
options(useFancyQuotes = FALSE)
rm(list = ls())

# load library
library(lme4)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)

###########################################################
# Berdugo et al. 2017 
# Full cite: Berdugo, Miguel, et al. "Plant spatial patterns identify alternative ecosystem multifunctionality states in global drylands." Nature Ecology & Evolution 1 (2017): 0003.
berdugo.dat <- read.csv("./data/Berdugo_etal_NEE_drylands.csv")
vars4func.berd <- c("NIT", "AMO", "OC", "TN", "AVP", "AMI", "PRO", "PEN",
                    "HEX", "ARO", "PHE", "PNT", "BGA", "PA", "IP", "TP")
berdugo.dat[, vars4func.berd] <- apply(berdugo.dat[, vars4func.berd], 2, 
                                       function(x) {(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)})

res.df.berd <- NULL
for (i in 1:16) {
  N <- dim(combn(vars4func.berd, i))[2]
  for (j in 1:N) {
    dat <- berdugo.dat[, combn(vars4func.berd, i)[, j]]
    dat <- data.frame(dat)
    emf <- apply(dat, 1, function(x) {
      ret <- sum(x) / i
    })
    df <- cbind(berdugo.dat, emf)
    res <- matrix(NA, nrow = 1, ncol = 5)
    mod <- lmer(emf ~ log(SR) + (1|Country), REML = TRUE, 
                control = lmerControl(optCtrl = list(maxfun = 2e5)),
                data = df)
    res[, 1] <- coef(summary(mod))[2, 1]
    res[, 2] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
    res[, 3] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
    res[, 4] <- i
    res[, 5] <- j
    cat(i, j, sep = "\n")
    res.df.berd <- rbind(res.df.berd, res)
  }
}
write.csv(res.df.berd, "./outputs/res.df.berd.csv")

###########################################################
# Gross et al. 2017
# Full cite: Gross, Nicolas, et al. "Functional trait diversity maximizes ecosystem multifunctionality." Nature ecology & evolution 1.5 (2017): 0132.

# import data
gross.dat <- read.csv("./data/Gross_etal_NEE_drylands.csv")

# data transformation
vars4func.gross <- c("BGL", "FOS", "AMP", "NTR", "I.NDVI")
gross.dat[, vars4func.gross] <- apply(gross.dat[, vars4func.gross], 2, function(x) {
  # ifelse(any(min(x, na.rm = T) < 0), log((x - min(x, na.rm = T)) + 1), log(x))
  x <- (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
})

# for loops to calculate standardized coefficients
res.df.gross <- NULL
for (i in 1:5) {
  N <- dim(combn(vars4func.gross, i))[2]
  for (j in 1:N) {
    dat <- gross.dat[, combn(vars4func.gross, i)[, j]]
    dat <- data.frame(dat)
    emf <- apply(dat, 1, function(x) {
      ret <- sum(x) / i
    })
    df <- cbind(gross.dat, emf)
    res <- matrix(NA, nrow = 1, ncol = 5)
    mod <- lmer(emf ~ log(SR) + (1|COU), REML = TRUE, 
                control = lmerControl(optCtrl = list(maxfun = 2e6)),
                data = df)
    res[, 1] <- coef(summary(mod))[2, 1]
    res[, 2] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
    res[, 3] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
    res[, 4] <- i
    res[, 5] <- j
    cat(i, j, sep = "\n")
    res.df.gross <- rbind(res.df.gross, res)
  }
}
write.csv(res.df.gross, "./outputs/res.df.gross.csv")

###########################################################
# Jing et al. 2015
# Full cite: Jing, Xin, et al. "The links between ecosystem multifunctionality and above-and belowground biodiversity are mediated by climate." Nature communications 6 (2015).

# import data
jing.dat <- read.csv("./data/jing_etal_Tibet.csv")

# data transformation
vars4func.jing <- c("AGB", "PTN", "PTP", "BGB",
                    "SSOC", "SSTN", "TAN", "SSTP")
# jing.dat$SR <- with(jing.dat, (log(SR) - mean(log(SR), na.rm = T)) / sd(log(SR), na.rm = T))
jing.dat[, vars4func.jing] <- apply(jing.dat[, vars4func.jing], 2, function(x) {
  # ifelse(any(min(x, na.rm = T) < 0), log((x - min(x, na.rm = T)) + 1), log(x))
  x <- (x - mean(x, na.rm = T)) / sd(x, na.rm = TRUE)
})

# for loops to calculate standardized coefficients
res.df.jing <- NULL
for (i in 1:8) {
  N <- dim(combn(vars4func.jing, i))[2]
  for (j in 1:N) {
    dat <- jing.dat[, combn(vars4func.jing, i)[, j]]
    dat <- data.frame(dat)
    emf <- apply(dat, 1, function(x) {
      ret <- sum(x) / i
    })
    df <- cbind(jing.dat, emf)
    res <- matrix(NA, nrow = 1, ncol = 5)
    mod <- lmer(emf ~ log(SR) + (1|Site), REML = TRUE, 
                control = lmerControl(optCtrl = list(maxfun = 2e6)),
                data = df)
    res[, 1] <- coef(summary(mod))[2, 1]
    res[, 2] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
    res[, 3] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
    res[, 4] <- i
    res[, 5] <- j
    cat(i, j, sep = "\n")
    res.df.jing <- rbind(res.df.jing, res)
  }
}
write.csv(res.df.jing, "./outputs/res.df.jing.csv")

###########################################################
# Quero et al. 2013
# Full cite: Quero, José L., et al. "On the importance of shrub encroachment by sprouters, climate, species richness and anthropic factors for ecosystem multifunctionality in semi-arid Mediterranean ecosystems." Ecosystems 16.7 (2013): 1248-1261.

# import data
quero.dat <- read.csv("./data/Quero_etal_Ecosystems_drylands.csv")

# data transformation
vars4func.quero <- c("ORC", "BGL", "HEX", "PEN", "TON",
                     "ATN", "AMI", "PRO", "AVP", "FOS")
# quero.dat$SR <- with(quero.dat, (log(SR) - mean(log(SR), na.rm = T)) / sd(log(SR), na.rm = T))
quero.dat[, vars4func.quero] <- apply(quero.dat[, vars4func.quero], 2, function(x) {
  # ifelse(any(min(x, na.rm = T) < 0), log((x - min(x, na.rm = T)) + 1), log(x))
  x <- (x - mean(x, na.rm = T)) / sd(x, na.rm = TRUE)
})

# for loops to calculate standardized coefficients
res.df.quero <- NULL
for (i in 1:10) {
  N <- dim(combn(vars4func.quero, i))[2]
  for (j in 1:N) {
    dat <- quero.dat[, combn(vars4func.quero, i)[, j]]
    dat <- data.frame(dat)
    emf <- apply(dat, 1, function(x) {
      ret <- sum(x) / i
    })
    df <- cbind(quero.dat, emf)
    res <- matrix(NA, nrow = 1, ncol = 5)
    mod <- lmer(emf ~ log(SR) + (1|VEG), REML = TRUE, 
                control = lmerControl(optCtrl = list(maxfun = 2e6)),
                data = df)
    res[, 1] <- coef(summary(mod))[2, 1]
    res[, 2] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
    res[, 3] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
    res[, 4] <- i
    res[, 5] <- j
    cat(i, j, sep = "\n")
    res.df.quero <- rbind(res.df.quero, res)
  }
}
write.csv(res.df.quero, "./outputs/res.df.quero.csv")

###########################################################
# to test whether BEMF relationship is related to the number of functions

# import data
df.berd <- read.csv("./outputs/res.df.berd.csv")
df.gross <- read.csv("./outputs/res.df.gross.csv")
df.jing <- read.csv("./outputs/res.df.jing.csv")
df.quero <- read.csv("./outputs/res.df.quero.csv")

# add short cites to each data set and significant signs
addCiteSigSigns <- function(x, case) {
  if (case == 1) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("Berdugo et al. 2017", dim(x)[1])
    x$cite <- as.factor(x$cite)
    x$sig <- ifelse((x$lci >0), "pos.sig", 
                    ifelse(x$uci < 0, "neg.sig", "neutral"))
    x$sig <- as.factor(x$sig)
  } else  if (case == 2) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("Gross et al. 2017", dim(x)[1])
    x$cite <- as.factor(x$cite)
    x$sig <- ifelse((x$lci >0), "pos.sig", 
                    ifelse(x$uci < 0, "neg.sig", "neutral"))
    x$sig <- as.factor(x$sig)
  } else  if (case == 3) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("Jing et al. 2015", dim(x)[1])
    x$cite <- as.factor(x$cite)
    x$sig <- ifelse((x$lci >0), "pos.sig", 
                    ifelse(x$uci < 0, "neg.sig", "neutral"))
    x$sig <- as.factor(x$sig)
  } else  if (case == 4) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("Quero et al. 2013", dim(x)[1])
    x$cite <- as.factor(x$cite)
    x$sig <- ifelse((x$lci >0), "pos.sig", 
                    ifelse(x$uci < 0, "neg.sig", "neutral"))
    x$sig <- as.factor(x$sig)
  }
  return(x)
}

df.berd <- addCiteSigSigns(df.berd, 1)
df.gross <- addCiteSigSigns(df.gross, 2)
df.jing <- addCiteSigSigns(df.jing, 3)
df.quero <- addCiteSigSigns(df.quero, 4)

# ggplot
df.all <- rbind(df.berd, df.gross, df.jing, df.quero)
df.all$sig <- factor(df.all$sig, levels = c("pos.sig", "neutral"))
ggplot(df.all, aes(x = nfuncs, y = coef)) +
  geom_jitter(width = 0.2, size = 0.8, 
              aes(color = sig, shape = sig)) +
  geom_smooth(aes(color = sig), lwd = 0.6, method = "lm", se = FALSE) +
  geom_smooth(method = "lm", lwd = 0.6, color = "black", se = FALSE) +
  scale_color_manual(values = c("#fb6a4a", "gray")) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, col = "grey") +
  labs(x = "Number of functions", y = "Slope of estimate") +
  facet_wrap(~ cite, ncol = 2) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
ggsave("./outputs/Observations.pdf", width = 5.5, height = 4.5)


