###########################################################
# Simulations: tweak what functions get included 
# 
# Xin Jing
# 8/3/2018
###########################################################

# clean the working environment
rm(list = ls())
# load library
library(TeachingDemos)
library(MASS)
library(plyr)
library(dplyr)
library(ggplot2)
library(magrittr)
# set seed
char2seed("correlation")


###########################################################
# simulate ecosystem functions with average correlation values of
# 0, 0.25, 0.5, 0.75 and 1
mean.corrs <- c(0, 0.25, 0.5, 0.75, 1)  # average correlation coefficients of functions in monoculture
no.functions <- 9  # number of functions simulated
mean.functions <- 0  # mean of functions with normal distribution
sd.functions <- 1  # sd of functions
no.communities <- 100  # number of artifical communities
richness.values <- rep(c(1:10), no.communities/10)  # richness levels
regional.richness.value <- 20  # regional richness levels
CF <- c(0, 3, 6, 20)  # maximum complementarity factor (CF)
r <- 0.25 # rate of CF grows to its maximum value

m1 <- df.sum <- list()
for (a in 1:100) {  # replications
  for(i in 1:length(mean.corrs)) {
    res.df <- NULL
    ecosystem.function <- matrix(NA, ncol = no.functions, nrow = no.communities)
    Sigma2 <- diag(no.functions)
    Sigma2[which(diag(no.functions) == 0)] <- mean.corrs[i]
    function.pool <- mvrnorm(n = regional.richness.value, 
                             mu = rep(mean.functions, no.functions),
                             Sigma = Sigma2,
                             tol = 1e-6, 
                             empirical = FALSE,
                             EISPACK = FALSE)
    function.pool <- apply(function.pool, 2, function(x) {x + abs(min(x))})
    for (j in 1:no.communities) {
      function.values <- function.pool[sample(c(1:regional.richness.value), richness.values[j]), ]
      if(richness.values[j] == 1) {
        for (k in 1:no.functions) {
          ecosystem.function[j, k] <- function.values[k]
        }
      } else {
        for (k in 1:no.functions) {
          ecosystem.function[j, k] <- mean(function.values[, k])
        }
      }
    }
    for (c in 1:length(CF)) {
      if (CF[c] == 0) {
        ecosystem.function2 <- ecosystem.function
      } else {
        ecosystem.function2 <- ecosystem.function * (CF[c] * ( 1 - ( 1 - 1/CF[c]) * exp(1 - richness.values^r)))
      }
      ecosystem.function2 <- scale(ecosystem.function2)
      ecosystem.function2 <- data.frame(ecosystem.function2)
      names(ecosystem.function2) <- paste("F", 1:no.functions, sep = "")
      vars4funcs <- names(ecosystem.function2)
      
      for (l in 1:no.functions) {
        N <- dim(combn(vars4funcs, l))[2]
        for (n in 1:N) {
          temp.dat <- ecosystem.function2[, combn(vars4funcs, l)[, n]]
          temp.dat <- data.frame(temp.dat)
          emf <- rowMeans(temp.dat)
          emf <- (emf - min(emf)) / (max(emf) - min(emf))
          df <- cbind(emf, richness.values)
          df <- data.frame(df)
          names(df) <- c("emf", "richness.values")
          mod <- lm(emf ~ log(richness.values), data = df)
          res <- matrix(NA, nrow = 1, ncol = 4)
          res[, 1] <- coef(summary(mod))[2, 1]
          res[, 2] <- l
          res[, 3] <- n
          res.df <- rbind(res.df, res)
        }
      }
    }
    res.df <- data.frame(res.df)
    names(res.df) <- c("coef", "nfuncs", "ncombns", "CFs")
    res.df$CFs <- rep(CF, each = dim(res.df)[1] / length(CF))
    m1[[i]] <- ddply(res.df, .(CFs), function(x) {
      mod <- coef(summary(lm(coef ~ nfuncs, data = x)))[2, ]
    })
  }
  df.sum[[a]] <- do.call("rbind", m1)
}
df <- do.call("rbind", df.sum)
df$cor.vars <- rep(rep(c(0, 0.25, 0.5, 0.75, 1), each = length(CF)), dim(df)[1] / 20)
df$lci <- df$Estimate - 1.96*df$`Std. Error`
df$uci <- df$Estimate + 1.96*df$`Std. Error`
df.sum2 <- ddply(df, .(cor.vars, CFs), function(x) {
  data.frame(
    Estimate.mean = mean(x$Estimate),
    Lower.limit.mean = mean(x$lci),
    Upper.limit.mean = mean(x$uci))
})

# df.sum2$CFs <- factor(df.sum2$CFs, levels = c("0", "3"),
                      # labels = c("Complementarity factor: 0", "Complementarity factor: 3"))
df.sum2[df.sum2$CFs %in% c(0, 3), ] %>% 
  mutate(CFs = factor(CFs, levels = c("0", "3"),
                      labels = c("Without Complementarity", "With Complementarity"))) %>% 
  ggplot(aes(cor.vars, Estimate.mean)) +
  geom_point(size = 3.5, color = "grey") +
  geom_errorbar(aes(ymin = Lower.limit.mean,
                    ymax = Upper.limit.mean), width = 0.025) +
  geom_hline(yintercept = 0, color = "gray", lty = 2) +
  facet_wrap(~ CFs) +
  labs(x = "Average correlation coefficients of functions", 
       y = "Average biodiversity effect on EMF") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave("./outputs/simulations_corr_matrix_1.2.pdf", width = 6.46, height = 4.1)

###########################################################
#                   End of Script                         #
###########################################################