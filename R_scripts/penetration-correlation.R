######################################################################################  
# PENETRATION CORRELATION ANALYSIS ----
######################################################################################  

# are the tissue samples with higher virus penetration correlated with those with more
# immune cells further down in the tissue? (are immune cell spatial distributions related to
# those cells hunting down viral particles?)

# are any of these cells in the area of the virus? binary

# install packages

library(tidyverse)
library(ggpubr)
library(ggsci)

# load data

load("data/processed/cell_and_virus_distance.rda")
load("data/processed/distance.rda")


## Initial Visualization ----

correlation_plot_1 <- cell_and_virus_distance |>
  group_by(target_cell, donor, tissue, condition) |>
  #filter(layer == 'Epithelium') |> # is there some subset of cells being affected?
  reframe(mean_epithelium_distance_um = mean(distance_to_annotation_with_epithelium_um),
            mean_depth = mean_depth) |> 
  distinct() |> 
  ggplot(aes(x = mean_depth, y = mean_epithelium_distance_um, color = target_cell)) +
  geom_point() +
  #scale_y_reverse() +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(title = 'Correlation between Viral Penetration and Mean Distance to Epithelium',
       x = 'mean viral penetration (um)') +
  facet_grid(condition*tissue ~ target_cell) +
  theme_pubr()+scale_color_npg()

# basically this plot tells us that anything within 30um of the epithelium barrier is
# "in the area of the virus"
cell_and_virus_distance |>
  ggplot(aes(x = mean_depth)) +
  geom_freqpoly(bins = 7)

cell_and_virus_distance |>
  group_by(target_cell, donor, tissue, condition) |>
  mutate(in_virus_area = case_when(
    distance_to_annotation_with_epithelium_um < 40 ~ TRUE,
    distance_to_annotation_with_epithelium_um >= 40 ~ FALSE
  )) |>
  group_by(target_cell, donor, tissue, condition) |>
  summarize(n = n(),
            num_positive = sum(in_virus_area),
            pct_positive = num_positive/n) |>
  #View()
  ggplot(aes(y = target_cell, x = pct_positive, color = target_cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(size = num_positive)) +
  labs(x = "% of cells in Viral Range (< 40um from Epithelium) after 4h HIV",
       y = "Immunofluorescence Target") +
  facet_wrap(~condition*tissue, ncol = 4) +
  theme_pubr()+scale_color_npg()

# can we also do this with our uninfected data?

correlation_plot_2 <- distance |>
  #filter(group == "Non infected") |>
  group_by(target_cell, donor, tissue, condition) |>
  mutate(in_virus_area = case_when(
    distance_to_annotation_with_epithelium_um < 30 ~ TRUE,
    distance_to_annotation_with_epithelium_um >= 30 ~ FALSE
  )) |>
  group_by(target_cell, donor, tissue, condition, group) |>
  summarize(n = n(),
            num_positive = sum(in_virus_area),
            pct_positive = num_positive/n) |>
  ggplot(aes(y = group, x = pct_positive, color = target_cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(size = num_positive)) +
  labs(x = "% of cells in Viral Range (< 30um from Epithelium)",
       y = "Immunofluorescence Target") +
  facet_grid(condition*tissue ~ target_cell) +
  #facet_wrap(~target_cell, ncol = 7) +
  theme_pubr()+scale_color_npg()




  
save(correlation_plot_1, file = "figures/plots/correlation_plot_1.rda")
save(correlation_plot_2, file = "figures/plots/correlation_plot_2.rda")

