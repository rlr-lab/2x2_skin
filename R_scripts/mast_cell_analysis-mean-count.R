## MAST CELL ANALYSIS ------

# load packages

library(ggpubr)
library(ggsci)
library(gt)
library(emmeans)
library(tidyverse)

source("R_scripts/functions/pairwise_comparisons.R")

# load data

tryptase <- read_csv("data/R21_PathCore.csv") |>
  janitor::clean_names(replace=janitor:::mu_to_u) |> 
  # accidentally labelled some target cells as 118/119 instead of 117
  mutate(target_cell = case_when(
    target_cell == 'CD118' ~ 'CD117',
    target_cell == 'CD119' ~ 'CD117',
    TRUE ~ target_cell
  )) |>
  mutate(donor = factor(donor),
         replicate = factor(replicate),
         positive = as.numeric(as.character(positive)),
         percent_pos = (positive/total_cells)*100,
         percent_pos_by_area = percent_pos/log10(area_um_2))

load("data/processed/tryptase.rds")



################################################################################
# CELL COUNT ANALYSIS ----
################################################################################


## Initial Visualizations

# lets take a look at some variables--
  
tryptase %>%
  ggplot(aes(x = condition, y = positive, color = target_cell)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_pubr() + scale_color_lancet() + 
  facet_grid(group ~ tissue, scales = "free") +
  rotate_x_text(90)


tryptase %>%
  ggplot(aes(x = condition, y = percent_pos, color = target_cell)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_pubr() + scale_color_lancet() + 
  facet_grid(group ~ tissue, scales = "free") +
  rotate_x_text(90)

tryptase %>%
  ggplot(aes(x = condition, y = total_cells, color = target_cell)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_pubr() + scale_color_lancet() + 
  facet_grid(group ~ tissue, scales = "free") +
  rotate_x_text(90)

tryptase %>%
  ggplot(aes(x = condition, y = area_um_2, color = target_cell)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_pubr() + scale_color_lancet() + 
  facet_grid(group ~ tissue, scales = "free") +
  rotate_x_text(90)

tryptase %>%
  #filter(target_cell == "CD117") |>
  ggplot(aes(x = condition, y = percent_pos_by_area, color = target_cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_pubr() + scale_color_lancet() + 
  facet_grid(group ~ tissue)

tryptase %>%
  ggplot(aes(
    x = interaction(group, condition),
    y = percent_pos_by_area,
    color = target_cell
  )) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_pubr() + scale_color_lancet() + 
  facet_grid(. ~ tissue) + 
  rotate_x_text(45)



################################################################################
### uninfected ----
################################################################################

# here's a model with just the uninfected tissues-- what's going on in non-infected case?
  
non_infected <- tryptase |>
  filter(group == "Non Infected") |>
  group_by(condition, tissue, target_cell, donor) |>
  summarize(mean_percent_pos_by_area = mean(percent_pos_by_area))

non_infected_model <-
  lme4::lmer(
    log10(mean_percent_pos_by_area) ~ condition * tissue * target_cell + (1 | donor),
    data = non_infected
  )

# is the model any good?
qqnorm(resid(non_infected_model))
qqline(resid(non_infected_model))
hist(resid(non_infected_model))
plot(fitted(non_infected_model), resid(non_infected_model))
abline(h = 0)

summary(non_infected_model)

performance::check_model(non_infected_model)

# make pairwise comparisons

compare_all_pairs(non_infected_model, factors = c('tissue', 'condition', 'target_cell')) |>
  display_comparison_table()


# what ramon did:

comps_to_make <- c("tissue",
                   "condition",
                   "target_cell",
                   "condition | tissue | target_cell",
                   "target_cell | condition | tissue",
                   "tissue | condition | target_cell",
                   "condition | tissue")


#comps_to_make <- c("condition | tissue",
#                   "condition | tissue | target_cell",
#                   "condition | target_cell")


compare_pairs(non_infected_model, 
              comparisons = comps_to_make,
              p_adjustment = 'fdr') |>
  display_comparison_table(title = "Non Infected Samples Model")

# looks different between CD117 and Tryptase
# also maybe diff between glans/shaft for circumcised individuals for both 
# CD117 and tryptase

p_vals <- compare_pairs(non_infected_model, 
                        comparisons = comps_to_make,
                        p_adjustment = 'fdr') |>
  get_pvals_for_comparisons()

# graphs of interesting comparisons -- non infected samples

non_infected |>
  #filter(tissue == "Shaft") |>
  ggplot(aes(x = condition, y = percent_pos_by_area, color = condition)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_pvalue_manual(p_vals |> filter(group1 == 'Circumcised', tissue != '.',
                                      target_cell != '.'), 
                     label = 'p_adj = {round(p.value, 3)}', y.position = 0.8) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(tissue ~ target_cell) + 
  labs(title = "Non Infected Samples",
       y = expression ("% of cells expressing CD117 or Tryptase per"~mu*m^2)) +
  rotate_x_text(45)

non_infected |>
  #filter(tissue == "Shaft") |>
  ggplot(aes(x = target_cell, y = mean_percent_pos_by_area, color = target_cell)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_pvalue_manual(p_vals |> filter(group1 == 'CD117', tissue != '.',
                                      condition != '.'), 
                     label = 'p_adj = {round(p.value, 5)}', y.position = 0.6) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(tissue ~ condition) + 
  labs(title = "Non Infected Samples",
       y = expression ("mean % of cells expressing CD117 or Tryptase per"~mu*m^2)) +
  rotate_x_text(45)

non_infected |>
  #filter(tissue == "Shaft") |>
  ggplot(aes(x = tissue, y = mean_percent_pos_by_area, color = tissue)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_pvalue_manual(p_vals |> filter(group1 == 'Glans', target_cell != '.',
                                      condition != '.'), 
                     label = 'p_adj = {round(p.value, 5)}', y.position = 0.6) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(target_cell ~ condition) + 
  labs(title = "Non Infected Samples",
       y = expression ("mean % of cells expressing CD117 or Tryptase per"~mu*m^2)) +
  rotate_x_text(45)

################################################################################
### uninfected CD117 ----
################################################################################

# here's a model with just the uninfected (original) tissues and the CD117

non_infected_cd117 <- tryptase |>
  filter(group == "Non Infected") |>
  filter(target_cell == "CD117")

non_infected_model <-
  lme4::lmer(
    log10(percent_pos_by_area) ~ condition * tissue + (1 | donor),
    data = non_infected_cd117
  )

# is the model any good?
qqnorm(resid(non_infected_model))
qqline(resid(non_infected_model))
hist(resid(non_infected_model))
plot(fitted(non_infected_model), resid(non_infected_model))
abline(h = 0)

summary(non_infected_model)


# make pairwise comparisons

compare_all_pairs(non_infected_model, factors = c('tissue', 'condition')) |>
  display_comparison_table()


p_vals <- compare_pairs(non_infected_model, 
                        comparisons = comps_to_make,
                        p_adjustment = 'fdr') |>
  get_pvals_for_comparisons()

# graphs of interesting comparisons -- non infected samples

non_infected |>
  #filter(tissue == "Shaft") |>
  ggplot(aes(x = tissue, y = percent_pos_by_area, color = tissue)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_pvalue_manual(p_vals |> filter(group1 == 'Glans', condition != '.'), 
                     label = 'p.adj.signif', y.position = 0.8, hide.ns = T,
                     tip.length = 0) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(condition ~ target_cell) + 
  labs(title = "Non Infected Samples",
       y = expression ("% of cells expressing CD117 or Tryptase per"~mu*m^2)) +
  rotate_x_text(45)





################################################################################
### infected ----
################################################################################

# We can then compare that to the infected case:
  
infected <- tryptase |>
  filter(group == "4h HIV")

infected_model <-
  lme4::lmer(
    log10(percent_pos_by_area) ~ condition * tissue * target_cell + (1|donor),
    data = infected
  )

# is the model any good?
summary(infected_model)
qqnorm(resid(infected_model))
qqline(resid(infected_model))
hist(resid(infected_model))
plot(fitted(infected_model), resid(infected_model))
abline(h = 0)

performance::check_model(infected_model)

# yeah it's okay

# make pairwise comparisons
compare_all_pairs(infected_model, factors = c('tissue', 'condition', 'target_cell')) |>
  display_comparison_table(title = "Infected Samples Model")

compare_pairs(infected_model, 
              comparisons = comps_to_make,
              p_adjustment = 'fdr') |>
  display_comparison_table(title = "Infected Samples Model")

# same differences as the infected model

p_vals <- compare_pairs(infected_model, 
                        comparisons = comps_to_make,
                        p_adjustment = 'fdr') |>
  get_pvals_for_comparisons()

infected |>
  ggplot(aes(x = tissue, y = percent_pos_by_area, color = tissue)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_pvalue_manual(p_vals |> filter(group1 == 'Glans', condition != '.'), 
                     label = 'p.adj.signif', y.position = 0.6, hide.ns = T,
                     tip.length = 0) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(condition ~ target_cell) + 
  labs(title = "Infected Samples",
       y = expression ("% of cells expressing CD117 or Tryptase per"~mu*m^2)) +
  rotate_x_text(45)


################################################################################
## paired samples: to do glans vs shaft -----
################################################################################

# additionally, we can look just at samples that are paired between the two sections

# get only donors that have paired infected samples
Donors_w4h <- unique(filter(tryptase, group == "4h HIV")$donor)

tryptase_paired <-
  tryptase %>% 
  filter(donor %in% Donors_w4h)

tryptase_paired %>%
  ggplot(aes(x = condition, y = percent_pos_by_area, color = group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(target_cell ~ tissue) + 
  rotate_x_text(45)

################################################################################
### glans ----
################################################################################

paired_glans <- tryptase_paired |> 
  filter(tissue == "Glans")

paired_glans_model <-
  lme4::lmer(
    log10(percent_pos_by_area) ~ condition * group * target_cell + (1|donor),
    data = paired_glans
  )

summary(paired_glans_model)

comps_to_make <- c("group",
                   "condition",
                   "target_cell",
                   "condition | group | target_cell",
                   "target_cell | condition | group",
                   "group | condition | target_cell")

comparisons_made <- compare_pairs(paired_glans_model, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')

compare_all_pairs(paired_glans_model, factors = c('group', 'condition', 'target_cell')) |>
  display_comparison_table()

display_comparison_table(comparisons_made)


# within glans, looks like we're getting signif differences for CD117 and Tryptase
# and across circumcision status between infected and not, 
# indicating that there is a shift in these immune cell populations between those two 

# also just in general

glans_p_vals <- compare_pairs(paired_glans_model, 
                        comparisons = comps_to_make,
                        p_adjustment = 'fdr') |>
  get_pvals_for_comparisons()

paired_glans |>
  ggplot(aes(x = group, y = percent_pos_by_area, color = group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_pvalue_manual(glans_p_vals |> filter(group1 == '4h HIV', condition != '.'), 
                     label = 'p.adj.signif', y.position = 0.5, hide.ns = T,
                     tip.length = 0) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(condition ~ target_cell) + 
  labs(title = "Differences within Glans Samples",
       y = expression ("% of cells expressing CD117 or Tryptase per"~mu*m^2)) +
  rotate_x_text(45)

################################################################################
### shaft ----
################################################################################

paired_shaft <- tryptase_paired |> 
  filter(tissue == "Shaft")

paired_shaft_model <-
  lme4::lmer(
    log10(percent_pos_by_area) ~ condition * group * target_cell + (1|donor),
    data = paired_shaft
  )

# is the model any good?
summary(paired_shaft_model)
qqnorm(resid(paired_shaft_model))
qqline(resid(paired_shaft_model))
hist(resid(paired_shaft_model))
plot(fitted(paired_shaft_model), resid(paired_shaft_model))
abline(h = 0)

performance::check_model(paired_shaft_model)



comps_to_make <- c("group",
                  "condition",
                  "target_cell",
                  "condition | group | target_cell",
                  "target_cell | condition | group",
                  "group | condition | target_cell")

comparisons_made <- compare_pairs(paired_shaft_model, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')

compare_all_pairs(paired_shaft_model, factors = c('group', 'condition', 'target_cell')) |>
  display_comparison_table(title = 'shaft model')

display_comparison_table(comparisons_made)

shaft_p_vals <- compare_pairs(paired_shaft_model, 
                              comparisons = comps_to_make,
                              p_adjustment = 'fdr') |>
  get_pvals_for_comparisons()

paired_shaft |>
  ggplot(aes(x = group, y = percent_pos_by_area, color = group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_pvalue_manual(shaft_p_vals |> filter(group1 == '4h HIV', condition != '.'), 
                     label = 'p.adj.signif', y.position = 0.5, hide.ns = T,
                     tip.length = 0) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(condition ~ target_cell) + 
  labs(title = "Differences within Shaft Samples",
       y = expression ("% of cells expressing CD117 or Tryptase per"~mu*m^2)) +
  rotate_x_text(45)


# this is where we see the difference in circumcised that isn't in uncircumcised

tryptase_paired |>
  ggplot(aes(x = group, y = percent_pos_by_area, color = group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_pubr() + scale_color_npg() + 
  #facet_grid(condition ~ tissue) + 
  rotate_x_text(45)


## mean 

# we can also look at the mean values of cells per area in each of these groups:
  
# but I'm not sure that we should?