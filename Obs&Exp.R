###########################################################
# Observations and expeiments
# to test whether BEMF relationship is related to the number of functions
rm(list = ls())
options(useFancyQuotes = FALSE)

library(ggplot2)
library(gridExtra)
library(cowplot)

# import data
df.berd <- read.csv("./outputs/res.df.berd.csv")
df.gross <- read.csv("./outputs/res.df.gross.csv")
df.jing <- read.csv("./outputs/res.df.jing.csv")
df.quero <- read.csv("./outputs/res.df.quero.csv")
df.biodepth <- read.csv("./outputs/res.df.BioDepth.csv")
df.fundiv <- read.csv("./outputs/res.df.FunDivEUROPE.csv")

# add short cites to each data set and significant signs
addCiteSigSigns <- function(x, case) {
  if (case == 1) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("(a) Berdugo et al. 2017", dim(x)[1])
    x$cite <- as.factor(x$cite)
    x$sig <- ifelse((x$lci >0), "pos.sig", 
                    ifelse(x$uci < 0, "neg.sig", "neutral"))
    x$sig <- as.factor(x$sig)
  } else  if (case == 2) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("(b) Gross et al. 2017", dim(x)[1])
    x$cite <- as.factor(x$cite)
    x$sig <- ifelse((x$lci >0), "pos.sig", 
                    ifelse(x$uci < 0, "neg.sig", "neutral"))
    x$sig <- as.factor(x$sig)
  } else  if (case == 3) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("(c) Jing et al. 2015", dim(x)[1])
    x$cite <- as.factor(x$cite)
    x$sig <- ifelse((x$lci >0), "pos.sig", 
                    ifelse(x$uci < 0, "neg.sig", "neutral"))
    x$sig <- as.factor(x$sig)
  } else  if (case == 4) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("(d) Quero et al. 2013", dim(x)[1])
    x$cite <- as.factor(x$cite)
    x$sig <- ifelse((x$lci >0), "pos.sig", 
                    ifelse(x$uci < 0, "neg.sig", "neutral"))
    x$sig <- as.factor(x$sig)
  } else if (case == 5) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("(e) Spehn et al. 2005", dim(x)[1])
    x$cite <- as.factor(x$cite)
    x$sig <- ifelse((x$lci >0), "pos.sig", 
                    ifelse(x$uci < 0, "neg.sig", "neutral"))
    x$sig <- as.factor(x$sig)
  } else if (case == 6) {
    names(x) <- c("X", "coef", "lci", "uci", "nfuncs", "ncombns")
    x$cite <- rep("(f) Ratcliffe et al. 2017", dim(x)[1])
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
df.biodepth <- addCiteSigSigns(df.biodepth, 5)
df.fundiv <- addCiteSigSigns(df.fundiv, 6)

# ggplot
df.all <- rbind(df.berd, df.gross, df.jing, df.quero, df.biodepth, df.fundiv)
df.all$sig <- factor(df.all$sig, levels = c("pos.sig", "neutral", "neg.sig"))
p1 <- ggplot(df.all, aes(x = nfuncs, y = coef)) +
  geom_jitter(width = 0.2, size = 0.8, 
              aes(color = sig, shape = sig)) +
  geom_smooth(aes(color = sig), lwd = 0.6, method = "lm", se = FALSE) +
  geom_smooth(method = "lm", lwd = 0.6, color = "black", se = FALSE) +
  scale_color_manual(values = c("#fb6a4a", "gray", "#4393c3")) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, col = "grey") +
  xlim(c(-1.5, 38)) +
  ylim(c(-0.65, 0.75)) +
  labs(x = "Number of functions", y = "Slope estimate") +
  facet_wrap(~ cite, ncol = 2) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')

# Jena
res.df <- read.csv("./outputs/combn_functions_Jena.csv")
res.df$sig <- factor(res.df$sig,
                     levels = c("pos.sig", "neutral", "neg.sig"))

p2 <- ggplot(res.df, aes(x = nfuncs, y = coef)) +
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

pdf("./outputs/Obs&Exp.pdf", width = 6.5, height = 8.5)
plot_grid(p1, p2, nrow = 2,
             rel_heights = c(3/4, 1/4))
dev.off()
