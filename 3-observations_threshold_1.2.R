###########################################################
# Observations
# By threshold approach
# 
# First created: 6/20/2018
#
# Contacts: Xin Jing <Xin.Jing@uvm.edu>
###########################################################

rm(list = ls())

# load library
library(TeachingDemos)
library(tidyverse)
library(multifunc)
library(cowplot)

char2seed("multifunc")

###########################################################
# Jena grassland biodiversity experiment
# Data from Meyer et al. 2018
# Full cite: Meyer, Sebastian T., et al. "Biodiversity–multifunctionality relationships depend on identity and number of measured functions." Nature ecology & evolution 2.1 (2018): 44.
# import data
jena.dat <- read.csv("./data/Meyer_etal_NEE_Jena_exp.csv")
jena.funcs <- read.csv("./data/Meyer_etal_NEE_Jena_exp_funcs.csv")

# data cleaning
jena.dat <- left_join(jena.dat, jena.funcs, by = "plotcode")
vars4func.jena <- names(jena.dat)[15:96]

# Sensitive analysis
# number of functions considered in total using one identical function
# by soil microbial biomass
jena.dat %>% 
  mutate(BM_microbes.std = BM_microbes/mean(sort(BM_microbes, decreasing = TRUE)[1:6])) %>% 
  ggplot(aes(sowndiv, BM_microbes.std)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm", size = 0.6) +
  labs(x = "Biodiversity (species richness)",
       y = "Soil microbial biomass") +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank())
ggsave("./outputs/observations_Jena_SMB.pdf", 
       height = 5.5, width = 6.0)


df.sim <- jena.dat %>% 
  dplyr::select(sowndiv, BM_microbes)
for (i in 3:81) {
  df.sim[i] <- df.sim$BM_microbes
}
vars4func.sim <- names(df.sim)[2:81]
LinearSlopes <- LinearSlopes.std <- NULL
funcsMaxed <- NULL
for(i in seq(5, 80, 10)) {
  simLinearSlopesSub <- simLinearSlopesSub.std <- NULL
  for (j in 1) {
    vars4func.sim.temp <- sample(vars4func.sim, i)
    simThresh <- getFuncsMaxed(df.sim, vars4func.sim.temp, threshmin = 0.05,
                               threshmax = 0.99, prepend = c("sowndiv"), 
                               maxN = 6)
    funcsMaxed <- rbind(funcsMaxed, simThresh)
    simLinearSlopes <- getCoefTab(funcMaxed ~ sowndiv, data = simThresh,
                                  coefVar = "sowndiv")
    simLinearSlopes$nFuncs <- i
    simLinearSlopesSub[[j]] <- simLinearSlopes
    # simThresh <- simThresh %>% 
    #   group_by(thresholds) %>% 
    #   mutate(funcMaxed.std = (funcMaxed - mean(funcMaxed)) / sd(funcMaxed)) %>% 
    #   ungroup() %>% 
    #   na.omit()
    simThresh$funcMaxed.std <- with(simThresh, funcMaxed/nFunc)
    simLinearSlopes.std <- getCoefTab(funcMaxed.std ~ sowndiv, data = simThresh,
                                      coefVar = "sowndiv")
    simLinearSlopes.std$nFuncs <- i
    simLinearSlopesSub.std[[j]] <- simLinearSlopes.std
  }
  LinearSlopes[[i]] <- do.call(rbind, simLinearSlopesSub)
  LinearSlopes.std[[i]] <- do.call(rbind, simLinearSlopesSub.std)
}
LinearSlopesAll <- do.call(rbind, LinearSlopes)
LinearSlopesAll$std <- "No"
LinearSlopesAll.std <- do.call(rbind, LinearSlopes.std)
LinearSlopesAll.std$std <- "Yes"
df.all <- rbind(LinearSlopesAll, LinearSlopesAll.std)
df.all %>% 
  filter(nFuncs %in% c(5, 15, 25, 35, 45)) %>% 
  mutate(std = factor(std, levels = c("No", "Yes"),
                      labels = c("Without standardization",
                                 "With standardization"))) %>%  
  ggplot(aes(x = thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate),
             alpha = 0.6, shape = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_grid(std ~ nFuncs, scale = "free_y") +
  labs(x = "Thresholds (%)",
       y = "Change in number of functions\nper addition of one species") +
  # guides(color = guide_legend("Number of functions\nin total considered")) +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
ggsave("./outputs/observations_multithresh_funcinTotal_SMB.pdf", 
       height = 5.5, width = 10)

# plot slopes
funcsMaxed$nFuncs <- rep(seq(5, 80, 10), each = dim(funcsMaxed)[1] / length(seq(5, 80, 10)))
pSlope <- funcsMaxed %>% 
  filter(nFuncs %in% c(5, 15, 25, 35, 45)) %>% 
  mutate(percent = 100 * thresholds) %>% 
  ggplot(aes(x = sowndiv, y = funcMaxed, group = percent)) +
  stat_smooth(method = "glm", 
              # method.args = list(family = quasipoisson(link = "log")),
              fill = NA, aes(color = percent), lwd = 0.6) +
  scale_color_gradient(name = "Percent of \nMaximum", 
                       low = "blue", high = "red") +
  ylim(c(-2, 50)) +
  labs(x = "Biodiversity (species richness)",
       y = "Change in number of functions\nper addition of one species") +
  facet_grid(~ nFuncs) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

pSlope.std <- funcsMaxed %>% 
  filter(nFuncs %in% c(5, 15, 25, 35, 45)) %>% 
  mutate(percent = 100 * thresholds) %>% 
  ggplot(aes(x = sowndiv, y = funcMaxed / nFunc, group = percent)) +
  stat_smooth(method = "glm", 
              # method.args = list(family = quasipoisson(link = "log")),
              fill = NA, aes(color = percent), lwd = 0.6) +
  scale_color_gradient(name = "Percent of \nMaximum", 
                       low = "blue", high = "red") +
  ylim(c(-0.05, 1.1)) +
  facet_grid(~ nFuncs) +
  labs(x = "Biodiversity (species richness)",
       y = "Change in percentage of functions\nper addition of one species") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
plot_grid(pSlope, pSlope.std, ncol = 1)
ggsave("./outputs/observations_multithresh_funcinTotal_SMBSlopes.pdf", 
       height = 5.5, width = 10)

# Sensitive analysis
# change the number of functions in total considered, but use subset variables
# the number of functions in total considered are 75, 60, 45, 30, 25, 20, 15, 10, 5 
# repeat 500 times at each number of functions in total considered
# run this section will take time
# LinearSlopes <- LinearSlopes.std <- NULL
# for(k in 1:500) {
#   print(k)
#   jenaLinearSlopesSub <- jenaLinearSlopesSub.std <- NULL
#   for(i in c(75, 60, 45, 30, 25, 20, 15, 10, 5)) {
#       if (i == 75) {
#         vars4func.jena.temp <- sample(vars4func.jena, i)
#       } else vars4func.jena.temp <- sample(vars4func.jena.temp, i)
#       
#       jenaThresh <- getFuncsMaxed(jena.dat, vars4func.jena.temp, threshmin = 0.05,
#                                   threshmax = 0.99, prepend = c("sowndiv"), maxN = 6)
#       jenaLinearSlopes <- getCoefTab(funcMaxed ~ sowndiv, data = jenaThresh,
#                                      coefVar = "sowndiv")
#       jenaLinearSlopes$nFuncs <- i
#       jenaLinearSlopesSub <- rbind(jenaLinearSlopesSub, jenaLinearSlopes)
#       jenaThresh$funcMaxed.std <- with(jenaThresh, funcMaxed / nFunc)
#       jenaLinearSlopes.std <- getCoefTab(funcMaxed.std ~ sowndiv, data = jenaThresh,
#                                          coefVar = "sowndiv")
#       jenaLinearSlopes.std$nFuncs <- i
#       jenaLinearSlopesSub.std <- rbind(jenaLinearSlopesSub.std, jenaLinearSlopes.std)
#   }
#   LinearSlopes <- rbind(LinearSlopes, jenaLinearSlopesSub)
#   LinearSlopes.std <- rbind(LinearSlopes.std, jenaLinearSlopesSub.std)
# }
# 
# LinearSlopesAll <- LinearSlopes
# LinearSlopesAll$std <- "No"
# LinearSlopesAll.std <- LinearSlopes.std
# LinearSlopesAll.std$std <- "Yes"
# df.all <- rbind(LinearSlopesAll, LinearSlopesAll.std)
# df.all %>% 
#   filter(nFuncs %in% c(5, 10, 15, 20, 25, 30)) %>%
#   mutate(std = factor(std, levels = c("No", "Yes"),
#                       labels = c("Without standardization",
#                                  "With standardization"))) %>%  
#   ggplot(aes(x = thresholds)) +
#   geom_point(aes(x = thresholds * 100, y = Estimate),
#              size = 0.2, shape = 1) +
#   geom_hline(yintercept = 0, lty = 2) +
#   facet_grid(std ~ nFuncs, scale = "free_y") +
#   labs(x = "Thresholds (%)",
#        y = "Change in number of functions\nper addition of one species") +
#   guides(color = guide_legend("Number of functions\nin total considered\nout of 82 functions")) +
#   theme_bw(base_size = 12.5) +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         legend.position = "none")
# ggsave("./outputs/observations_multithresh_funcinTotal.pdf", 
#        height = 5.5, width = 10)
# 
# df.all %>% 
#   filter(nFuncs %in% c(5, 10, 15, 20, 25, 30, 45, 60, 75)) %>% 
#   mutate(std = factor(std, levels = c("No", "Yes"),
#                       labels = c("Without standardization",
#                                  "With standardization"))) %>%  
#   group_by(std, nFuncs, thresholds) %>% 
#   summarise(mu = mean(Estimate),
#             lci = quantile(Estimate, 0.025),
#             uci = quantile(Estimate, 0.975)) %>% 
#   ggplot(aes(thresholds, mu)) +
#   geom_ribbon(fill = "grey50", aes(thresholds * 100,
#                                    ymin = lci,
#                                    ymax = uci)) +
#   geom_point(aes(thresholds * 100, mu), shape = 1, size = 0.8) +
#   geom_hline(yintercept = 0, lty = 2, color = "gray") +
#   facet_grid(std ~ nFuncs, scales = "free_y") +
#   scale_x_continuous(breaks = c(10, 50, 90)) +
#   labs(x = "Thresholds (%)",
#        y = "Change in number of functions\nper addition of one species") +
#   theme_bw(base_size = 14.5) +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         legend.position = "none")
# ggsave("./outputs/observations_multithresh_funcinTotal_CI.pdf", 
#        height = 5.0, width = 10)

# multiple-threshold analysis
jenaThresh <- getFuncsMaxed(jena.dat, vars4func.jena, threshmin = 0.05,
                            threshmax = 0.99, prepend = c("sowndiv"), maxN = 6)
jenaLinearSlopes <- getCoefTab(funcMaxed ~ sowndiv, data = jenaThresh,
                               groupVar = "thresholds", coefVar = "sowndiv")
ggplot(jenaLinearSlopes, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  theme_bw()
jenaThresh$funcMaxed.std <- with(jenaThresh, funcMaxed / nFunc)
jenaLinearSlopes.std <- getCoefTab(funcMaxed.std ~ sowndiv, data = jenaThresh,
                                   groupVar = "thresholds", coefVar = "sowndiv")
ggplot(jenaLinearSlopes.std, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  theme_bw()

###########################################################
# Global drylands 
# Data from Maestre et al. 2012
# Full cite: Maestre, Fernando T., et al. "Plant species richness and ecosystem multifunctionality in global drylands." Science 335.6065 (2012): 214-218.

# load data
maestre.dat <- read.csv("./data/Maestre_Global_drylands_final_236_28_5_2018.csv")

# Z-scoring transformation for each function
vars4func.maestre <- c("TON", "BGL", "FOS", "P.HCL", "ORC", "AMO", "NIT",
                       "AMI", "PRO", "PHE", "ARO", "HEX", "PEN", "NTR")

maestreThresh <- getFuncsMaxed(maestre.dat, vars4func.maestre, threshmin = 0.05,
                               threshmax = 0.99, prepend = c("SR"), maxN = 60)
maestreLinearSlopes <- getCoefTab(funcMaxed ~ SR, data = maestreThresh,
                                  coefVar = "SR")
ggplot(maestreLinearSlopes, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  theme_bw()
maestreThresh$funcMaxed.std <- with(maestreThresh, funcMaxed / nFunc)
maestreLinearSlopes.std <- getCoefTab(funcMaxed.std ~ SR, data = maestreThresh,
                                      coefVar = "SR")
ggplot(maestreLinearSlopes.std, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  theme_bw()
###########################################################
# Tibetan grasslands
# Data from Jing et al. 2015
# Full cite: Jing, Xin, et al. "The links between ecosystem multifunctionality and above-and belowground biodiversity are mediated by climate." Nature communications 6 (2015).

# import data
tibet.dat <- read.csv("./data/jing_etal_Tibet.csv")

# data transformation
vars4func.tibet <- c("AGB", "PTN", "PTP", "BGB",
                     "SSOC", "SSTN", "TAN", "SSTP")
tibetThresh <- getFuncsMaxed(tibet.dat, vars4func.tibet, threshmin = 0.05,
                               threshmax = 0.99, prepend = c("SR"), maxN = 60)
tibetLinearSlopes <- getCoefTab(funcMaxed ~ SR, data = tibetThresh,
                                  coefVar = "SR")
ggplot(tibetLinearSlopes, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  theme_bw()
tibetThresh$funcMaxed.std <- with(tibetThresh, funcMaxed / nFunc)
tibetLinearSlopes.std <- getCoefTab(funcMaxed.std ~ SR, data = tibetThresh,
                                      coefVar = "SR")
ggplot(tibetLinearSlopes.std, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  theme_bw()
###########################################################
# BioDepth
# Data from Spehn et al. 2015
# Full cite: Spehn, E. M., et al. "Ecosystem effects of biodiversity manipulations in European grasslands." Ecological monographs 75.1 (2005): 37-63.
# import data
biodepth.dat <- read.csv("./data/all_biodepth.csv")
vars4func.biodepth <- c("biomassY3", "root3", 
                        "N.Soil", "cotton3",
                        "N.g.m2", "light3")

biodepthThresh <- getFuncsMaxed(biodepth.dat, vars4func.biodepth, threshmin = 0.05,
                               threshmax = 0.99, prepend = c("Diversity"), maxN = 60)
biodepthLinearSlopes <- getCoefTab(funcMaxed ~ Diversity, data = biodepthThresh,
                                  coefVar = "Diversity")
ggplot(biodepthLinearSlopes, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  theme_bw()
biodepthThresh$funcMaxed.std <- with(biodepthThresh, funcMaxed / nFunc)
biodepthLinearSlopes.std <- getCoefTab(funcMaxed.std ~ Diversity, 
                                       data = biodepthThresh,
                                       coefVar = "Diversity")
ggplot(biodepthLinearSlopes.std, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  theme_bw()

###########################################################
# FunDivEUROPE
# Data from Ratcliffe et al. 2017
# Full cite: Ratcliffe, Sophia, et al. "Biodiversity and ecosystem functioning relations in European forests depend on environmental context." Ecology letters 20.11 (2017): 1414-1426.
# import data
fundiv.dat <- read.csv("./data/Ratcliffe_et_al_ELE_12849_ecosystem_function_variables.csv")
vars4func.fundiv <- names(fundiv.dat)[7:32]

fundivThresh <- getFuncsMaxed(fundiv.dat, vars4func.fundiv, threshmin = 0.05,
                               threshmax = 0.99, prepend = c("target_species_richness"), maxN = 60)
fundivLinearSlopes <- getCoefTab(funcMaxed ~ target_species_richness, data = fundivThresh,
                                  coefVar = "target_species_richness")
ggplot(fundivLinearSlopes, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  theme_bw()
fundivThresh$funcMaxed.std <- with(fundivThresh, funcMaxed / nFunc)
fundivLinearSlopes.std <- getCoefTab(funcMaxed.std ~ target_species_richness, data = fundivThresh,
                                      coefVar = "target_species_richness")
ggplot(fundivLinearSlopes.std, aes(x = Thresholds)) +
  geom_point(aes(x = thresholds * 100, y = Estimate)) +
  theme_bw()

df <- rbind(jenaLinearSlopes, jenaLinearSlopes.std,
            biodepthLinearSlopes, biodepthLinearSlopes.std,
            fundivLinearSlopes, fundivLinearSlopes.std,
            tibetLinearSlopes, tibetLinearSlopes.std,
            maestreLinearSlopes, maestreLinearSlopes.std)
df$labs <- rep(c("Jena", "Jena",
                 "Biodepth", "Biodepth",
                 "Fundiv", "Fundiv",
                 "Tibet", "Tibet",
                 "Drylands", "Drylands"),
               each = 95)
df$std <- rep(c("No", "Yes"), each = 95)
p.obsStd <- df %>% 
  # filter(labs != "Fundiv") %>%
  mutate(labs = factor(labs,
                       levels = c("Jena", "Fundiv", "Biodepth", 
                                  "Tibet", "Drylands"),
                       labels = c("Jena grassland", "European forests",
                                  "European grasslands", "Tibetan grasslands",
                                  "Global drylands"))) %>%
  filter(std != "No") %>% 
  ggplot(aes(x = Thresholds)) +
  geom_ribbon(fill = "grey50",
              aes(thresholds * 100,
                  ymin = Estimate - 1.96*`Std. Error`,
                  ymax = Estimate + 1.96*`Std. Error`)) +
  geom_point(shape = 1, size = 0.8, aes(x = thresholds * 100, y = Estimate)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray") +
  facet_grid( ~labs, scales = "fixed") +
  scale_x_continuous(breaks = c(10, 50, 90)) +
  labs(x = "Thresholds (%)",
       y = "Change in percentage of functions\nper addition of one species") +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave("./outputs/observations_multithresh_std.pdf", 
       height = 4.5, width = 10)

p.obs <- df %>% 
  # filter(labs != "Fundiv") %>%
  mutate(labs = factor(labs,
                       levels = c("Jena", "Fundiv", "Biodepth", 
                                  "Tibet", "Drylands"),
                       labels = c("Jena grassland", "European forests",
                                  "European grasslands", "Tibetan grasslands",
                                  "Global drylands"))) %>% 
  filter(std != "Yes") %>% 
  ggplot(aes(x = Thresholds)) +
  geom_ribbon(fill = "grey50",
              aes(thresholds * 100,
                  ymin = Estimate - 1.96*`Std. Error`,
                  ymax = Estimate + 1.96*`Std. Error`)) +
  geom_point(shape = 1, size = 0.8, aes(x = thresholds * 100, y = Estimate)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray") +
  facet_grid( ~labs, scales = "fixed") +
  scale_x_continuous(breaks = c(10, 50, 90)) +
  labs(x = "Thresholds (%)",
       y = "Change in number of functions\nper addition of one species") +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave("./outputs/observations_multithresh_nostd.pdf", 
       height = 4.5, width = 10)
plot_grid(p.obs, p.obsStd, ncol = 1, align = TRUE)
ggsave("./outputs/observations_multithresh_obs.pdf", 
              height = 6.8, width = 9.5)

###########################################################
#                   End of Script                         #
###########################################################