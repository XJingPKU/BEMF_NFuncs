###########################################################
# The relationship between BEMF relationship and the number of functions
# 
# Data from the BIODEPTH, FunDivEUROPE and Jena
# 
# First created: 10/25/2017
#
# Contacts: Xin Jing <Xin.Jing@uvm.edu>
###########################################################
options(useFancyQuotes = FALSE)
rm(list = ls())

# load library
library(dplyr)
library(plyr)
library(lme4)
library(ggplot2)
library(gridExtra)

###########################################################
# BioDepth
# import data
biodepth <- read.csv("./data/all_biodepth.csv")
site.inf <- read.csv("./data/all_biodepth_site_vars.csv")

# data cleaning
biodepth <- left_join(biodepth, site.inf, by = "location")
df <- biodepth %>% 
  mutate(fb = factor(NewBlock)) %>% 
  select(location, NewBlock, latitude, longitude, altitude, 
         MAP, GST, pH, Diversity, biomassY3,
         root3, cotton3, wood3, N.Soil, N.g.m2, light3)
df[, -c(1:5, 9)] <- data.frame(apply(df[, -c(1:5, 9)], 2, function(x) {
  # x <- log(x + 1)
  x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

# linear-mixed effects models
vars4func <- c("biomassY3", "root3", "cotton3", 
               "N.Soil", "cotton3",
               "N.g.m2", "light3")
res.df <- NULL
for (i in 1:7) {
  N <- dim(combn(vars4func, i))[2]
  for (j in 1:N) {
    dat <- df[, combn(vars4func, i)[, j]]
    dat <- data.frame(dat)
    emf <- apply(dat, 1, function(x) {
      ret <- sum(x, na.rm = TRUE) / i
    })
    df2 <- cbind(df, emf)
    res <- matrix(NA, nrow = 1, ncol = 5)
    mod <- lmer(emf ~ log(Diversity) + (Diversity|location), REML = TRUE, 
                control = lmerControl(optCtrl = list(maxfun = 2e5)),
                data = df2)
    res[, 1] <- coef(summary(mod))[2, 1]
    res[, 2] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
    res[, 3] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
    res[, 4] <- i
    res[, 5] <- j
    cat(i, j, sep = "\n")
    res.df <- rbind(res.df, res)
  }
}
# write.csv(res.df, "./outputs/res.df.BioDepth.csv")
res.df <- data.frame(res.df)
names(res.df) <- c("coef", "lci", "uci", "nfuncs", "ncombns")
res.df$sig <- ifelse((res.df$lci >0), "pos.sig", 
                     ifelse(res.df$uci < 0, "neg.sig", "neutral"))
res.df$sig <- factor(res.df$sig,
                     levels = c("pos.sig", "neutral", "neg.sig"))

p1 <- ggplot(res.df, aes(x = nfuncs, y = coef)) +
  geom_jitter(width = 0.2, size = 0.8, 
              aes(color = sig, shape = sig)) +
  geom_smooth(aes(color = sig), lwd = 0.6, method = "lm", se = FALSE) +
  geom_smooth(method = "lm", lwd = 0.6, color = "black", se = FALSE) +
  scale_color_manual(values = c("#fb6a4a", "gray", "#4393c3")) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, col = "grey") +
  xlim(c(-2.0, 35)) +
  ylim(c(-0.65, 0.75)) +
  labs(x = "Number of functions", y = "Slope estimate") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        # legend.title = element_blank(),
        # legend.background = element_blank())
        legend.position = 'none')

###########################################################
# FunDivEUROPE
# import data
# fundiv.EU <- read.csv("./data/Ratcliffe_et_al_ELE_12849_ecosystem_function_variables.csv")
# site.inf <- read.csv("./data/Ratcliffe_et_al_ELE_12849_regional_context_variables.csv")
# 
# # data cleaning
# fundiv.EU <- left_join(fundiv.EU, site.inf, by = "region")
# df <- fundiv.EU
# df[, 7:32] <- data.frame(apply(df[, 7:32], 2, function(x) {
#   # x <- log(x + 1)
#   x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
# }))
# 
# set.seed(20171130)
# 
# vars4func <- names(df)[7:32]
# emf.dat <- df[vars4func]
# names(emf.dat) <- paste("EF", 1:26, sep = "")
# vars4func <- names(emf.dat)
# 
# # linear-mixed effects models
# res.df <- NULL
# for (i in 1:26) {
#   N <- dim(combn(vars4func, i))[2]
#   if (N > 500) J <- 1:500 else J <- 1:N
#   k <- 0
#   res1 <- numeric(length(J))
#   res2 <- numeric(length(J))
#   res3 <- numeric(length(J))
#   res4 <- numeric(length(J))
#   res5 <- numeric(length(J))
#   
#   for (j in J) {
#     # make an empty data for each loop
#     dat <- NULL
#     emf <- NULL
#     df2 <- NULL
#     mod <- NULL
#     
#     # calculate EMF
#     if (i == 1) {
#       dat <- emf.dat[, j]
#       emf <- dat
#     } else {
#       dat <- sample(emf.dat[vars4func], i)
#       emf <- apply(dat, 1, function(x) {
#         ret <- sum(x, na.rm = TRUE) / i
#       })
#     }
#     df2 <- cbind(df, emf)
#     
#     # fit the model
#     mod <- lmer(emf ~ log(target_species_richness) + 
#                   (target_species_richness|region) +
#                   (1|composition), REML = TRUE, 
#                 data = df2)
#     k <- k + 1
#     res1[k] <- coef(summary(mod))[2, 1]
#     res2[k] <- coef(summary(mod))[2, 1] - 1.96 * coef(summary(mod))[2, 2]
#     res3[k] <- coef(summary(mod))[2, 1] + 1.96 * coef(summary(mod))[2, 2]
#     res4[k] <- i
#     res5[k] <- j
#     cat(i, j, sep = "\n")
#   }
#   res <- matrix(NA, nrow = 1, ncol = 5)
#   res <- cbind(res1, res2, res3, res4, res5)
#   # return
#   res.df <- rbind(res.df, res)
# }
# write.csv(res.df, "./outputs/res.df.FunDivEUROPE.csv")
# res.df <- data.frame(res.df)
# names(res.df) <- c("coef", "lci", "uci", "nfuncs", "ncombns")
# res.df$sig <- ifelse((res.df$lci >0), "pos.sig", 
#                      ifelse(res.df$uci < 0, "neg.sig", "neutral"))
res.df <- read.csv("./outputs/combn_functions_FunDiv_excluding_drivers.csv")
res.df$sig <- factor(res.df$sig,
                     levels = c("pos.sig", "neutral", "neg.sig"))

p2 <- ggplot(res.df, aes(x = nfuncs, y = coef)) +
  geom_jitter(width = 0.2, size = 0.8, 
              aes(color = sig, shape = sig)) +
  geom_smooth(aes(color = sig), lwd = 0.6, method = "lm", se = FALSE) +
  geom_smooth(method = "lm", lwd = 0.6, color = "black", se = FALSE) +
  scale_color_manual(values = c("#fb6a4a", "gray", "#4393c3")) +
  ylim(c(-0.65, 0.75)) +
  xlim(c(-2.0, 35)) +
  # scale_x_continuous(breaks = seq(0, 28, 2)) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, col = "grey") +
  labs(x = "Number of functions", y = "Slope estimate") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        # legend.title = element_blank(),
        # legend.background = element_blank())
        legend.position = 'none')

###########################################################
# # import data
# jena.exp <- read.csv("./data/Meyer_etal_NEE_Jena_exp.csv")
# jena.funcs <- read.csv("./data/Meyer_etal_NEE_Jena_exp_funcs.csv")
# 
# # data cleaning
# jena <- left_join(jena.exp, jena.funcs, by = "plotcode")
# 
# df <- jena
# df[, 15:96] <- data.frame(apply(df[, 15:96], 2, function(x) {
#   # x <- log(x + 1)
#   x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
# }))
# 
# set.seed(20171130)
# 
# vars4func <- names(df)[15:96]
# emf.dat <- df[vars4func]
# names(emf.dat) <- paste("EF", 1:82, sep = "")
# vars4func <- names(emf.dat)
# 
# res.df <- NULL
# for (i in 1:82) {
#   # N <- dim(combn(vars4func, i))[2]
#   # if (N > 500) J <- 1:500 else J <- 1:N
#   ifelse(i == 1, J <- 1:82, 
#          ifelse(i == 81, J <- 1:82, 
#                 ifelse(i == 82, J <- 1, J <- 1:500)))
#   k <- 0
#   res1 <- numeric(length(J))
#   res2 <- numeric(length(J))
#   res3 <- numeric(length(J))
#   res4 <- numeric(length(J))
#   res5 <- numeric(length(J))
#   
#   for (j in J) {
#     # make an empty data for each loop
#     dat <- NULL
#     emf <- NULL
#     df2 <- NULL
#     mod <- NULL
#     
#     # calculate EMF
#     if (i == 1) {
#       dat <- emf.dat[, j]
#       emf <- dat
#     } else {
#       if (i == 81) {
#         dat <- emf.dat[, -j]
#         emf <- apply(dat, 1, function(x) {
#           ret <- sum(x, na.rm = TRUE) / 81
#         })
#       } else {
#         dat <- sample(emf.dat[vars4func], i)
#         emf <- apply(dat, 1, function(x) {
#           ret <- sum(x, na.rm = TRUE) / i
#         })
#       }
#     }
#     
#     df2 <- cbind(df, emf)
#     
#     # fit the model
#     mod <- lm(emf ~ block + log(sowndiv), 
#               data = df2)
#     k <- k + 1
#     res1[k] <- coef(summary(mod))[5, 1]
#     res2[k] <- coef(summary(mod))[5, 1] - 1.96 * coef(summary(mod))[5, 2]
#     res3[k] <- coef(summary(mod))[5, 1] + 1.96 * coef(summary(mod))[5, 2]
#     res4[k] <- i
#     res5[k] <- j
#     cat(i, j, sep = "\n")
#   }
#   res <- matrix(NA, nrow = 1, ncol = 5)
#   res <- cbind(res1, res2, res3, res4, res5)
#   # return
#   res.df <- rbind(res.df, res)
# }
# res.df <- data.frame(res.df)
# names(res.df) <- c("coef", "lci", "uci", "nfuncs", "ncombns")
# res.df$sig <- ifelse((res.df$lci >0), "pos.sig", 
#                      ifelse(res.df$uci < 0, "neg.sig", "neutral"))
# res.df$sig <- factor(res.df$sig,
#                      levels = c("pos.sig", "neutral", "neg.sig"))
# write.csv(res.df, "./outputs/combn_functions_Jena.csv")
res.df <- read.csv("./outputs/combn_functions_Jena.csv")
res.df$sig <- factor(res.df$sig,
                     levels = c("pos.sig", "neutral", "neg.sig"))

p3 <- ggplot(res.df, aes(x = nfuncs, y = coef)) +
  geom_jitter(width = 0.2, size = 0.8, 
              aes(color = sig, shape = sig)) +
  geom_smooth(method = "lm", lwd = 0.6, color = "black", se = FALSE) +
  # geom_smooth(aes(color = sig), lwd = 0.6, method = "lm", se = FALSE) +
  scale_color_manual(values = c("#fb6a4a", "gray", "#4393c3")) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, col = "grey") +
  labs(x = "Number of functions", y = "Slope estimate") +
  ylim(c(-0.65, 0.75)) +
  xlim(c(0.80, 82.5)) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        # legend.title = element_blank(),
        # legend.background = element_blank())
        legend.position = 'none')

pdf("./outputs/Experiments.pdf", width = 9, height = 6)
grid.arrange(arrangeGrob(p1, p2, ncol = 2),
                        p3, nrow = 2)
dev.off()

###########################################################