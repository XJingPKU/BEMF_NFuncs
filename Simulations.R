###########################################################
# simutation
# The following code are modified by Xin Jing from Gamfeldt and Roger 2017
# Full cite: Gamfeldt, Lars, and Fabian Roger. "Revisiting the biodiversity-ecosystem multifunctionality relationship." Nature ecology & evolution 1.7 (2017): 168.

rm(list = ls())
options(useFancyQuotes = FALSE)
library(dplyr)
library(tidyr)
library(ggplot2)

source("./R/Multifunc_simulations_functions.R")

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
Slope_res <- data.frame(Estimate = numeric(),
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
  #mutate_at(vars(one_of(func.names)), function(x) {(x - min(x)) / (max(x) - min(x))})
  
  # loop over all subsets of function of size 1:funcnum
  for (i in seq_len(funcnum)) { 
    
    # all poosibel combination of i out of funcnum functions
    func_comb <- combn(func.names, i)
    
    # loop over all function combinations of size i
    for ( k  in seq_len(ncol(func_comb))) { 
      
      # calculate mean function
      AvFunc_temp <- AvFunc %>%
        select(Richness, one_of(func_comb[ ,k])) %>% 
        mutate(meanFunction = rowMeans(.[func_comb[ ,k]]))
      
      # fit linear model
      mod <- lm(meanFunction ~ Richness, data = AvFunc_temp)
      
      # get slope estimate
      est <- summary(mod)$coefficients[2,]
      
      # store results
      Slope_res <- data.frame(t(est)) %>% 
        mutate(., nfunc = i) %>% 
        mutate(ncomp = l) %>% 
        rbind(Slope_res, .)
    }
  }
}

# Plot: the same as Figure 2a
Slope_res %>% 
  filter(ncomp %in% c(0,ceiling(funcnum/3),2*ceiling(funcnum/3),funcnum)) %>%
  ggplot(aes(x = nfunc, y = Estimate, colour = as.factor(ncomp)))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.75),
             alpha = 0.5, shape = 21)+
  geom_smooth(method = "lm", se = F, size = 0.5, 
              position = position_dodge(width = 0.5))+
  scale_color_brewer(guide = guide_legend(title = "Number of functions\nwith complementarity",
                                          nrow=2,byrow=TRUE),
                     palette = "Set1")+
  scale_x_continuous(breaks = seq(1,funcnum,1))+
  scale_y_continuous(limits = c(NA, 0.038))+
  labs(y = "Slope estimate",
       x = "Number of functions considered")+
  theme_bw()+
  theme(legend.position = "bottom")

# using 95% confidence interval to test the significance of the slope
Slope_res$lci <- Slope_res$Estimate - 1.96 * Slope_res$Std..Error
Slope_res$uci <- Slope_res$Estimate + 1.96 * Slope_res$Std..Error
Slope_res$sig <- with(Slope_res, ifelse((lci >0), "pos.sig", 
                                        ifelse(uci < 0, "neg.sig", "neutral")))
Slope_res$sig <- factor(Slope_res$sig, levels = c("pos.sig", "neutral"))

# Plot
Slope_res %>% 
  filter(ncomp %in% c(0:5, 7:9)) %>%
  ggplot(aes(x = nfunc, y = Estimate)) +
  geom_smooth(method = "lm", lwd = 0.6, color = "black", se = FALSE) +
  geom_jitter(width = 0.15, size = 1.2, 
              aes(color = sig, shape = sig)) +
  geom_smooth(aes(color = sig), lwd = 0.6, method = "lm", se = FALSE) +
  scale_color_manual(values = c("#fb6a4a", "gray")) +
  scale_x_continuous(breaks = seq(1, funcnum, 2))+
  scale_y_continuous(limits = c(NA, 0.042))+
  facet_wrap(~ ncomp, ncol = 3) +
  labs(y = "Slope estimate",
       x = "Number of functions")+
  theme_bw(base_size = 16.5)+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        # legend.title = element_blank(),
        # legend.background = element_blank()
        legend.position = 'none')
ggsave("./outputs/Simulation.pdf", width = 7.5, height = 7.0)
# ggsave("./outputs/Simulation.pdf", width = 7.5, height = 5)
