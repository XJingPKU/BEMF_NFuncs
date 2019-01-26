###########################################################
# Simutations
# The following code are modified by Xin Jing based on Gamfeldt and Roger 2017
# Full cite: Gamfeldt, Lars, and Fabian Roger. 
# "Revisiting the biodiversity-ecosystem multifunctionality relationship." 
# Nature ecology & evolution 1.7 (2017): 168.
###########################################################

rm(list = ls())
# load library
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
# load functions
source("./R/Multifunc_simulations_functions.R")

###########################################################
# data simulation
specnum <- 12  # the number of species
funcnum <- 9  # the number of functions
distribution = "runif"
maxrep <- choose(specnum, floor(specnum/2))
#maxrep <- 500

FuncMat <- FunctionValue(specnum, funcnum, distribution, min = 0, max = 1)

func.names <- as.character( unique( FuncMat$Functions))

SpecMat <- SpeciesMatrix(specnum = specnum, maxrep = maxrep)

CF = 3
r = 0.25

# empty dataframe to store results
Slope_res_mean <- Slope_res_sum <- Slope_res_unit <- data.frame(Estimate = numeric(),
                                                                `Std. Error` = numeric(),
                                                                `t value` = numeric(),    
                                                                `Pr(>|t|)` = numeric(),
                                                                nfunc = numeric(),
                                                                ncomp = numeric())

# loop over all possible number of functions with complementarity
for (l in 0:funcnum) {
  
  set.seed(999)
  
  # choose method = average if no functions with complementarity and method = comp otherwise
  if(l == 0) {
    method = "av"
  }  else {
    method = "comp"
    compfunc = func.names[1:l]
  }
  
  # draw function values and calculate mean function for all richness levels
  AvFunc <- AverageFunction(SpecMat, FuncMat,
                            method = method, 
                            CF = CF, 
                            compfunc = compfunc,
                            r = r)
  
  # standardize functions
  AvFunc <- AvFunc %>% 
    select(Richness, one_of(func.names)) %>% 
    mutate_at(vars(one_of(func.names)), function(x) {x / max(x)})
  
  # loop over all subsets of function of size 1:funcnum
  for (i in seq_len(funcnum)) { 
    
    # all poosibel combination of i out of funcnum functions
    func_comb <- combn(func.names, i)
    
    # loop over all function combinations of size i
    for ( k  in seq_len(ncol(func_comb))) { 
      
      # calculate mean function
      AvFunc_mean <- AvFunc %>%
        select(Richness, one_of(func_comb[ ,k])) %>% 
        mutate(meanFunction = rowMeans(.[func_comb[ ,k]]))
      
      # calculate sum function
      AvFunc_sum <- AvFunc %>%
        select(Richness, one_of(func_comb[ ,k])) %>% 
        mutate(meanFunction = rowSums(.[func_comb[ ,k]]))
      
      # calculate scaled function
      AvFunc_unit <- AvFunc %>%
        select(Richness, one_of(func_comb[ ,k])) %>% 
        mutate(meanFunction = rowMeans(.[func_comb[ ,k]])) %>% 
        mutate(meanFunction = (meanFunction - min(meanFunction)) / (max(meanFunction) - min(meanFunction)))
      
      fitLinearModel <- function(dat, Slope_res) {
        # fit linear model
        mod <- lm(meanFunction ~ Richness, data = dat)
        
        # get slope estimate
        est <- summary(mod)$coefficients[2,]
        
        # store results
        Slope_res <- data.frame(t(est)) %>% 
          mutate(., nfunc = i) %>% 
          mutate(ncomp = l) %>% 
          rbind(Slope_res, .)
      }
      
      Slope_res_mean <- fitLinearModel(dat = AvFunc_mean, Slope_res = Slope_res_mean)
      Slope_res_sum <- fitLinearModel(dat = AvFunc_sum, Slope_res = Slope_res_sum)
      Slope_res_unit <- fitLinearModel(dat = AvFunc_unit, Slope_res = Slope_res_unit)
    }
  }
}

# relationship between the number of functions and the biodiversity effect on EMF
# showing three approaches including the averaging, summing and scaled by [0, 1]
df <- rbind(Slope_res_mean, Slope_res_sum, Slope_res_unit)
df$std <- gl(3, dim(df)[1] / 3, labels = c("Averaging approach",
                                           "Summing approach",
                                           "Scaling approach"))
df %>% 
  filter(ncomp %in% c(0,ceiling(funcnum/3),2*ceiling(funcnum/3),funcnum)) %>%
  ggplot(aes(x = nfunc, y = Estimate, colour = as.factor(ncomp)))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.75),
             alpha = 0.5, shape = 21)+
  geom_smooth(method = "lm", se = F, size = 0.5, 
              position = position_dodge(width = 0.5))+
  scale_color_brewer(guide = guide_legend(title = "Number of functions\nwith complementarity",
                                          nrow=2,byrow=TRUE),
                     palette = "Set1")+
  scale_x_continuous(breaks = seq(1, funcnum, 2))+
  facet_wrap(~ std, scales = "free", ncol = 3) +
  labs(y = "Biodiversity effect on EMF",
       x = "Number of functions considered")+
  theme_bw(base_size = 14.5)+
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave("./outputs/simulations_gamfeldt-Roger.pdf", height = 4.5, width = 8.0)

# Classify the biodiversity effect on EMF into three groups
# significantly positive, neutral and significantly negative
# using 95% confidence interval for significance test
Slope_res <- df
Slope_res$lci <- Slope_res$Estimate - 1.96 * Slope_res$Std..Error
Slope_res$uci <- Slope_res$Estimate + 1.96 * Slope_res$Std..Error
Slope_res$sig <- with(Slope_res, ifelse((lci >0), "pos.sig", 
                                        ifelse(uci < 0, "neg.sig", "neutral")))
Slope_res$sig <- factor(Slope_res$sig, levels = c("pos.sig", "neutral"))

# Plot
Slope_res %>% 
  filter(ncomp %in% c(0, 3, 6, 9)) %>%
  mutate(sig = factor(sig, levels = c("pos.sig", "neutral", "neg.sig"),
                      labels = c("Positive", "Neutral", "Negative"))) %>%
  ggplot(aes(x = nfunc, y = Estimate)) +
  geom_smooth(method = "lm", lwd = 0.6, color = "black", se = FALSE) +
  geom_jitter(width = 0.15, size = 1.2, 
              aes(color = sig, shape = sig)) +
  geom_smooth(aes(color = sig), lwd = 0.6, method = "lm", se = FALSE) +
  scale_color_manual(values = c("#fb6a4a", "gray")) +
  scale_x_continuous(breaks = seq(1, funcnum, 2))+
  # scale_y_continuous(limits = c(NA, 0.042))+
  facet_grid(std ~ ncomp, scales = "free") +
  labs(y = "Slope estimate",
       x = "Number of functions")+
  theme_bw(base_size = 16.5)+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank(),
        # legend.background = element_blank()
        legend.position = 'right')
# ggsave("./outputs/simulation_scaled_by_unit_and_Sig.pdf", height = 6.5, width = 9.0)

# the proportion/probability of functions with a significantly positive/neutral relationship with diversity
df.nfuncs.sig <- ddply(Slope_res, .(std, ncomp, nfunc, sig), nrow)
df.tot.obs <- ddply(df.nfuncs.sig, .(std, ncomp, nfunc), function(x) sum(x$V1))
df.nfuncs.sig <- merge(df.nfuncs.sig, df.tot.obs, by = c("std", "ncomp", "nfunc"))
df.nfuncs.sig$prop <- with(df.nfuncs.sig, V1.x / V1.y)

df.nfuncs.sig %>%
  filter(ncomp %in% c(0, 3, 6, 9)) %>%
  mutate(sig = factor(sig, levels = c("pos.sig", "neutral", "neg.sig"),
                      labels = c("Positive", "Neutral", "Negative"))) %>%
  ggplot(aes(x = nfunc, y = prop, color = sig)) +
  geom_point(aes(shape = sig), size = 1.5) +
  geom_smooth(method = "glm",
              # method.args = list(family = "binomial"),
              se = FALSE, size = 0.25) +
  scale_x_continuous(breaks = seq(1, funcnum, 2))+
  scale_y_continuous(breaks = seq(0, 1.15, 0.20)) +
  # ylim(c(-0.15, 1.140)) +
  scale_color_manual(values = c("#fb6a4a", "gray", "#4393c3")) +
  facet_grid(std ~ ncomp, scales = "free") +
  labs(x = "Number of functions", y = "Probability of positive or neutral \nbiodiversity effect on EMF") +
  theme_bw(base_size = 16.5) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank())
ggsave("./outputs/simulation_probability.pdf", height = 6.2, width = 8.8)

###########################################################
#                   End of Script                         #
###########################################################