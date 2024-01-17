######################################################################################  
# PENETRATION CORRELATION ANALYSIS ----
######################################################################################  

# are the tissue samples with higher virus penetration correlated with those with more
# immune cells further down in the tissue? (are immune cell spatial distributions related to
# those cells hunting down viral particles?)

# install packages

library(tidyverse)
library(ggpubr)

# load data

load("data/processed/cell_and_virus_distance.rda")


## Initial Visualization ----

cell_and_virus_distance |>
  group_by(target_cell, donor, tissue, condition) |>
  filter(layer == 'Epithelium') |> # is there some subset of cells being affected?
  reframe(mean_baseline_distance_um = mean(distance_to_annotation_with_baseline_um),
            mean_depth = mean_depth) |> 
  distinct() |> 
  ggplot(aes(x = mean_depth, y = mean_baseline_distance_um, color = target_cell)) +
  geom_point() +
  scale_y_reverse() +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(title = 'i have no idea whats going on with innate immune cells',
       subtitle = 'cells in epithelium ONLY',
       x = 'mean viral penetration (um)') +
  facet_grid(condition*tissue ~ target_cell) +
  theme_pubr()+scale_color_npg()
  
