###########################################################
# Global map of the empirical datasets
# Five datasets were used, including Jena grassland, European grasslands, 
# European forests, Tibetan grasslands and global drylands
# XJ
# November 16 2018
###########################################################

rm(list = ls())
# load library
library(ggmap)
library(ggplot2)
library(ggthemes)
library(cowplot)

###########################################################
# Plotting locations on to a world map
# load data
loc.inf <- read.csv("./data/site_inf.csv")
loc.inf$site <- factor(loc.inf$site, 
                       levels = c("Jena grassland", 
                                  "European forests",
                                  "European grasslands", 
                                  "Tibetan grasslands",
                                  "Global drylands"))

# load world map
world <- map_data(map = "world")

# set coordiates for subsets of the data
lon <- c(-120, -120, -90, -90, -120, 
         -74, -74, -62, -62, -74, 
         -8, -8, 30.5, 30.5, -8,
         89, 89, 102, 102, 89, 
         -10, -10, 45, 45, -10, 
         140.5, 140.5, 150, 150, 140.5)
lat <- c(20, 40, 40, 20, 20, 
         -50, -20, -20, -50, -50, 
         30, 68, 68, 30, 30,
         28, 40, 40, 28, 28, 
         -35, 20, 20, -35, -35, 
         -43, -27, -27, -43, -43)
df <- data.frame(lon = lon, lat = lat)
df$loc <- rep(c("mp.na", "mp.eu", "mp.ti", "mp.sa", "mp.af", "mp.au"), 
              each = 5)
df.labs <- data.frame(lon = c())


# world map
mp.world <- ggplot() +
  geom_map(data = world, map = world, 
           aes(group = group, map_id = region),
           fill = "white", color = "gray90", size = 0.3) +
  geom_point(data = loc.inf, aes(x = longitude, y = latitude, 
                                 color = site), 
             size = 0.9, alpha = 0.8) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", 
                               "#e7298a", "#66a61e")) +
  geom_path(data = df, aes(x = lon, y = lat, group = loc),
            color = "gray70", size = 0.5) +
  labs(x = "", y = "") +
  theme_map() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.72, 0.25),
        legend.key = element_blank(),
        legend.title = element_blank())


# regional map
# north american
mp.na <- ggplot() +
  geom_map(data = world, map = world, 
           aes(group = group, map_id = region),
           fill = "white", color = "gray78", size = 0.3) +
  geom_point(data = loc.inf, aes(x = longitude, y = latitude, 
                                 color = site), 
             size = 2.0, alpha = 0.66) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", 
                                 "#e7298a", "#66a61e")) +
  scale_x_continuous(limits = c(-120, -90), breaks = c(-115, -105, -95),
                     labels = c("115 W", "105 W", "95 W")) +
  scale_y_continuous(limits = c(20, 40), breaks = c(24, 30, 36),
                     labels = c("22 N", "30 N", "38 N")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        line = element_line(size = 0.25),
        panel.border = element_rect(size = 0.5))
# south american
mp.sa <- ggplot() +
  geom_map(data = world, map = world, 
           aes(group = group, map_id = region),
           fill = "white", color = "gray78", size = 0.3) +
  geom_point(data = loc.inf, aes(x = longitude, y = latitude, 
                                 color = site), 
             size = 2.0, alpha = 0.66) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", 
                               "#e7298a", "#66a61e")) +
  scale_x_continuous(limits = c(-74, -62), breaks = c(-72, -68, -64),
                     labels = c("72 W", "68 W", "64 W")) +
  scale_y_continuous(limits = c(-50, -20), breaks = c(-45, -35, -25),
                     labels = c("45 S", "35 S", "25 S")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        line = element_line(size = 0.25),
        panel.border = element_rect(size = 0.5))
# Europe
mp.eu <- ggplot() +
  geom_map(data = world, map = world, 
           aes(group = group, map_id = region),
           fill = "white", color = "gray78", size = 0.3) +
  geom_point(data = loc.inf, aes(x = longitude, y = latitude, 
                                 color = site), 
             size = 2.0, alpha = 0.66) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", 
                               "#e7298a", "#66a61e")) +
  scale_x_continuous(limits = c(-8, 30), breaks = c(-5, 10, 25),
                     labels = c("5 W", "10 E", "25 E")) +
  scale_y_continuous(limits = c(30, 68), breaks = c(38, 50, 62),
                     labels = c("38 N", "50 N", "62 N")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        line = element_line(size = 0.25),
        panel.border = element_rect(size = 0.5))
# tibet
mp.ti <- ggplot() +
  geom_map(data = world, map = world, 
           aes(group = group, map_id = region),
           fill = "white", color = "gray78", size = 0.3) +
  geom_point(data = loc.inf, aes(x = longitude, y = latitude, 
                                 color = site), 
             size = 2.0, alpha = 0.66) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", 
                               "#e7298a", "#66a61e")) +
  scale_x_continuous(limits = c(89, 102), breaks = c(90, 95, 100),
                     labels = c("90 E", "95 E", "100 E")) +
  scale_y_continuous(limits = c(28, 40), breaks = c(30, 34, 38),
                     labels = c("30 N", "34 N", "38 N")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        line = element_line(size = 0.25),
        panel.border = element_rect(size = 0.5))
# south africa
mp.af <- ggplot() +
  geom_map(data = world, map = world, 
           aes(group = group, map_id = region),
           fill = "white", color = "gray78", size = 0.3) +
  geom_point(data = loc.inf, aes(x = longitude, y = latitude, 
                                 color = site), 
             size = 2.0, alpha = 0.66) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", 
                               "#e7298a", "#66a61e")) +
  scale_x_continuous(limits = c(-10, 45), breaks = c(-5, 15, 35),
                     labels = c("5 W", "15 E", "35 E")) +
  scale_y_continuous(limits = c(-35, 20), breaks = c(-30, -10, 10),
                     labels = c("30 S", "10 S", "10 N")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        line = element_line(size = 0.25),
        panel.border = element_rect(size = 0.5))
# austrilia
mp.au <- ggplot() +
  geom_map(data = world, map = world, 
           aes(group = group, map_id = region),
           fill = "white", color = "gray78", size = 0.3) +
  geom_point(data = loc.inf, aes(x = longitude, y = latitude, 
                                 color = site), 
             size = 2.0, alpha = 0.66) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", 
                               "#e7298a", "#66a61e")) +
  scale_x_continuous(limits = c(140.5, 150), breaks = c(142.5, 145.0, 147.5),
                     labels = c("142.5 E", "145.0 E", "147.5 E")) +
  scale_y_continuous(limits = c(-43, -27), breaks = c(-40, -35, -30),
                     labels = c("40 S", "35 S", "30 S")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        line = element_line(),
        panel.border = element_rect(size = 0.5))


# combine global and regional map
first <- plot_grid(mp.na, mp.eu, mp.ti, ncol = 3,
                   rel_widths = c(1/3, 1/3, 1/3),
                   labels = c("b)", "c)", "d)"))
second <- plot_grid(mp.world, labels = "a)")
third <- plot_grid(mp.sa, mp.af, mp.au, ncol = 3,
                   rel_widths = c(1/3, 1/3, 1/3),
                   labels = c("e)", "f)", "g)"))

# save global map
pdf("./outputs/0-global_map.pdf", height = 8.8, width = 8.5)
plot_grid(first, second, third, nrow = 3,
          rel_heights = c(2/7, 3/7, 2/7))
dev.off()

###########################################################
#                   End of Script                         #
###########################################################