###########################################################
# Schematic diagram represents the averaging/summing approaches 
# Xin Jing
# 8/1/2018
###########################################################

# clean the working environment
rm(list = ls())
# load library
library(TeachingDemos)
library(reshape2)
library(ggplot2)
library(plyr)
library(magrittr)
# set seed
char2seed("averaing")

###########################################################
# generate a dataframe for averaging approach
div <- rep(1:10, each = 10)  # diversity (10 levels)
f1 <- scale(div + 1.5*rnorm(100, 0, 1))  # random norm values (ecosystem function 1)
f2 <- scale(div + 10.5*rnorm(100, 0, 1))  # random norm values (ecosystem function 2)
f12.av <- (f1 + f2) / 2  # average of the two functions
f12.sm <- (f1 + f2)  # sum of the two functions
f12.sc.av <- (f12.av - min(f12.av)) / (max(f12.av) - min(f12.av))  # scaled averaing multifunctionality metric
f12.sc.sm <- (f12.sm - min(f12.sm)) / (max(f12.sm) - min(f12.sm))  # scaled averaing multifunctionality metric

# combine the data 
df <- data.frame(div, f1, f2, f12.av, f12.sm, f12.sc.av, f12.sc.sm)
df.long <- melt(df, id.vars = "div")

# inspect the relationship between biodiversity and ecosystem functions
df.long %>% 
  ggplot(aes(div, value, color = variable)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ variable, ncol = 2, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank())

# extract regression coefficients
df.coefs <- ddply(df.long, .(variable), function(x) {
  mod <- lm(value ~ div, data = x)
  cof <- summary(mod)$coefficients[2]
})

# inspect regression coefficients
df.coefs$variable <- factor(df.coefs$variable, 
                            levels = c("f1", "f2", "f12.av", "f12.sm", "f12.sc.av", "f12.sc.sm"),
                            labels = c("F1", "F2", "(F1 + F2) / 2", "F1 + F2", "Scaled [(F1 + F2) / 2]", "Scaled [F1 + F2]"))
main <- df.coefs %>% 
  dplyr::filter(variable %in% c("F1", "F2", "(F1 + F2) / 2", "F1 + F2")) %>%
  ggplot(aes(variable, V1)) +
  geom_col(fill = "gray85") +
  geom_text(aes(x = variable, y = (V1 - 0.05), 
                label = round(V1, 1)),
            size = 8.5) +
  ylim(c(0, 0.6)) +
  labs(x = "Ecosystem functions", 
       y = "Biodiversity effect on EMF\n(Slope estimate)") +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank())

sub.f1 <- ggplotGrob(ggplot(df, aes(x = div, y = f1)) +
                       geom_point(shape = 1, color = "gray", size = 0.8) +
                       geom_smooth(method = "lm", color = "blue", size = 0.8, se = FALSE) +
                       labs(x = "Biodiversity", y = "F1") +
                       ylim(c(-2.5, 2.5)) +
                       theme_bw() +
                       theme(panel.grid = element_blank()))

sub.f2 <- ggplotGrob(ggplot(df, aes(x = div, y = f2)) +
                       geom_point(shape = 1, color = "gray", size = 0.8) +
                       geom_smooth(method = "lm", color = "blue", size = 0.8, se = FALSE) +
                       labs(x = "Biodiversity", y = "F2") +
                       ylim(c(-2.5, 2.5)) +
                       theme_bw() +
                       theme(panel.grid = element_blank()))

sub.f3 <- ggplotGrob(ggplot(df.long[df.long$variable %in% c("f1", "f2", "f12.av"), ], 
                            aes(x = div, y = value, color = variable)) +
                       geom_point(shape = 1, size = 0.8) +
                       geom_smooth(method = "lm", size = 0.8, se = FALSE) +
                       labs(x = "Biodiversity", y = "(F1 + F2) / 2") +
                       ylim(c(-2.5, 2.5)) +
                       scale_color_manual(values = c("grey", "grey", "blue")) +
                       theme_bw() +
                       theme(panel.grid = element_blank(),
                             legend.position = 'none'))

sub.f4 <- ggplotGrob(ggplot(df.long[df.long$variable %in% c("f1", "f2", "f12.sm"), ], 
                            aes(x = div, y = value, color = variable)) +
                       geom_point(shape = 1, size = 0.8) +
                       geom_smooth(method = "lm", size = 0.8, se = FALSE) +
                       labs(x = "Biodiversity", y = "F1 + F2") +
                       ylim(c(-5, 5)) +
                       scale_color_manual(values = c("grey", "grey", "blue")) +
                       theme_bw() +
                       theme(panel.grid = element_blank(),
                             legend.position = 'none'))

main + 
  annotation_custom(grob = sub.f1, xmin = 0.50, xmax = 1.450,
                    ymin = 0.41, ymax = 0.62) +
  annotation_custom(grob = sub.f2, xmin = 1.50, xmax = 2.450, 
                    ymin = 0.41, ymax = 0.62) +
  annotation_custom(grob = sub.f3, xmin = 2.50, xmax = 3.450,
                    ymin = 0.41, ymax = 0.62) +
  annotation_custom(grob = sub.f4, xmin = 3.50, xmax = 4.450,
                    ymin = 0.41, ymax = 0.62)

ggsave("./outputs/Schematic_diagram.pdf", width = 8.5, height = 5.49)

###########################################################
#                   End of Script                         #
###########################################################