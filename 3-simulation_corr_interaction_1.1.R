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
library(ggplot2)
# set seed
char2seed("correlation")

###########################################################
# simulate ecosystem functions with average correlation values of
# 0, 0.25, 0.5, 0.75 and 1
mean.corrs <- c(0, 0.25, 0.5, 0.75, 1)
no.functions <- 9
mean.functions <- 0
sd.functions <- 1
interaction.strength <- 0.05
# interaction.strength <- 0
no.communities <- 100
richness.values <- rep(c(1:10), no.communities/10)
regional.richness.value <- 20

m1 <- df.sum <- list()

for (a in 1:100) {  # replications
  res.df <- NULL
  # effect.size <- lower.limit <- upper.limit <- c()
  for(i in 1:length(mean.corrs)) {
    ecosystem.function <- ecosystem.function.int <- matrix(NA, ncol = no.functions, nrow = no.communities)
    Sigma2 <- diag(no.functions)
    Sigma2[which(diag(no.functions) == 0)] <- mean.corrs[i]
    function.pool <- mvrnorm(n = regional.richness.value, 
                             mu = rep(mean.functions, no.functions),
                             Sigma = Sigma2,
                             tol = 1e-6, 
                             empirical = FALSE,
                             EISPACK = FALSE)
    for (j in 1:no.communities) {
      function.values <- function.pool[sample(c(1:regional.richness.value), richness.values[j]), ]
      interaction.values <- interaction.strength * sample(abs(function.values), (richness.values[j] * richness.values[j]), replace = TRUE) 
      # if(richness.values[j] == 1) {interaction.values <- 0}
      if(richness.values[j] == 1) {
        for (k in 1:no.functions) {
          ecosystem.function[j, k] <- ecosystem.function.int[j, k] <- function.values[k]
        }
      } else {
        for (k in 1:no.functions) {
          ecosystem.function[j, k] <- mean(function.values[, k])
          ecosystem.function.int[j, k] <- mean(function.values[, k]) + sum(interaction.values)
        }
      }
    }
    ecosystem.function <- data.frame(scale(ecosystem.function))
    ecosystem.function.int <- data.frame(scale(ecosystem.function.int))
    # print(mean(cor(ecosystem.function)))
    # print(mean(cor(ecosystem.function.int)))
    names(ecosystem.function) <- names(ecosystem.function.int) <- paste("F", 1:no.functions, sep = "")
    vars4funcs <- names(ecosystem.function)
    for (l in 1:no.functions) {
      N <- dim(combn(vars4funcs, l))[2]
      for (n in 1:N) {
        temp.dat <- ecosystem.function[, combn(vars4funcs, l)[, n]]
        temp.dat <- data.frame(temp.dat)
        emf <- rowMeans(temp.dat)
        emf <- (emf - min(emf)) / (max(emf) - min(emf))
        mod <- lm(emf ~ log(richness.values))
        res <- matrix(NA, nrow = 1, ncol = 3)
        res[, 1] <- coef(summary(mod))[2, 1]
        res[, 2] <- l
        res[, 3] <- n
        res.df <- rbind(res.df, res)
        # repeat again
        temp.dat <- ecosystem.function.int[, combn(vars4funcs, l)[, n]]
        temp.dat <- data.frame(temp.dat)
        emf <- rowMeans(temp.dat)
        emf <- (emf - min(emf)) / (max(emf) - min(emf))
        mod <- lm(emf ~ log(richness.values))
        res <- matrix(NA, nrow = 1, ncol = 3)
        res[, 1] <- coef(summary(mod))[2, 1]
        res[, 2] <- l
        res[, 3] <- n
        res.df <- rbind(res.df, res)
      }
    }
  }
  res.df <- data.frame(res.df)
  names(res.df) <- c("coef", "nfuncs", "ncombns")
  res.df$cor.vars <- rep(c("0", "0.25", "0.50", "0.75", "1.00"),
                         each = dim(res.df)[1] / length(mean.corrs))
  res.df$int.fac <- rep(c("no.int", "int"),
                        dim(res.df)[1] / 2)
  m1[[a]] <- ddply(res.df, .(int.fac, cor.vars), function(x) {
    mod <- coef(summary(lm(coef ~ nfuncs, data = x)))[2, ]
  })
}

df <- do.call("rbind", m1)
df$lci <- df$Estimate - 1.96*df$`Std. Error`
df$uci <- df$Estimate + 1.96*df$`Std. Error`
df.sum <- ddply(df, .(int.fac, cor.vars), function(x) {
  data.frame(
    Estimate.mean = mean(x$Estimate),
    Lower.limit.mean = mean(x$lci),
    Upper.limit.mean = mean(x$uci))
})

df.sum$int.fac <- factor(df.sum$int.fac,
                         levels = c("no.int", "int"),
                         labels = c("Without Interaction", "With Interaction"))

ggplot(df.sum, aes(cor.vars, Estimate.mean)) +
  geom_point(size = 3.5, color = "grey") +
  geom_errorbar(aes(ymin = Lower.limit.mean,
                    ymax = Upper.limit.mean), width = 0.06) +
  geom_hline(yintercept = 0, color = "gray", lty = 2) +
  facet_wrap(~ int.fac) +
  labs(x = "Average correlation coefficients of functions", 
       y = "Average biodiversity effect on EMF") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
ggsave("./outputs/Figure5_simulations_corr_matrix.pdf", width = 6.46, height = 4.1)

###########################################################
#                   End of Script                         #
###########################################################