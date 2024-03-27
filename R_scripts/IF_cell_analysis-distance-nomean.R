######################################################################################  
# TARGET CELL ANALYSIS ------
######################################################################################  

# load packages

library(ggpubr)
library(ggsci)
library(gt)
library(emmeans)
library(tidyverse)

source("R_scripts/functions/pairwise_comparisons.R")

# load data

load('data/processed/distance.rda')


################################################################################
# DISTANCE ANALYSIS ----
################################################################################


## Initial Visualization

# We can look at this data a few different ways:


distance |>
  ggplot(aes(x= group,
             y= distance_to_annotation_with_epithelium_um,
             color= target_cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha = 0.2) +
  theme_pubr() +
  scale_color_lancet() + 
  facet_grid(tissue*layer~target_cell, scales = "free") +
  rotate_x_text(90)


distance |>
  ggplot(aes(x=interaction(group,target_cell),
             y=distance_to_annotation_with_epithelium_um,
             color=target_cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge())+
  theme_pubr()+
  scale_y_reverse() +
  scale_color_lancet()+
  facet_grid(~tissue*layer, scales = "free")+
  rotate_x_text(45)

distance |>
  filter(target_cell=="CD4") |>
  ggplot(aes(x=group,
             y=distance_to_annotation_with_epithelium_um,
             color= group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(), alpha = 0.3) +
  theme_pubr()+
  labs(title = "CD4+ Cells Only") +
  scale_color_npg()

distance |>
  ggplot(aes(x=group,
             y=(distance_to_annotation_with_epithelium_um+0.1),
             color=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+
  scale_color_npg()+
  facet_grid(tissue~target_cell)+
  scale_y_log10()+
  rotate_x_text(45)

distance |>
  ggplot(aes(x=group,
             y=distance_to_annotation_with_epithelium_um+0.1,
             color=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+
  scale_color_npg()+
  facet_grid(tissue~target_cell)+
  scale_y_log10()+
  rotate_x_text(45)


## modeling

################################################################################
### glans ----
################################################################################



glans <- distance |> 
  filter(tissue == 'Glans')

# check distribution of outcome variable
ggplot(data = glans, aes(x = log10(distance_to_annotation_with_epithelium_um))) +
  geom_histogram()

model_glans <- lme4::lmer(log10(distance_to_annotation_with_epithelium_um) ~ condition*group*layer*target_cell+(1|donor), 
                          data = glans)

# is the model any good? -- better if we use log10 transform
summary(model_glans)
qqnorm(resid(model_glans))
qqline(resid(model_glans))
hist(resid(model_glans))
plot(fitted(model_glans), resid(model_glans))
abline(h = 0)

performance::check_model(model_glans)


# glance at all
compare_all_pairs(model_glans, factors = c('group', 'condition', 'target_cell', 'layer')) |>
  display_comparison_table_zvalue(title = 'glans model')


# here, we don't really care if certain cell types are more abundant than others, we 
# care if their abundance changes between conditions

comps_to_make <- c("group",
                   "layer",
                   "condition",
                   #"target_cell",
                   "condition | group",
                   "condition | group | target_cell",
                   "group | condition | target_cell",
                   "condition | group | target_cell | layer",
                   #"target_cell | condition | group | layer",
                   "group | condition | target_cell | layer")

comparisons_made <- compare_pairs(model_glans, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')


display_comparison_table_zvalue(comparisons_made, title = 'glans table')

# here we can see significant differences between CD3 as well as CD4 single positive 
# cells in uncircumcised epithelium between infected and non infected-- a significant difference was not detected in corresponding circumsized tissue:


p_vals_glans <- compare_pairs(model_glans, 
                        comparisons = comps_to_make,
                        p_adjustment = 'fdr') |>
  get_pvals_for_comparisons()

distance |>
  ggplot(aes(x = distance_to_annotation_with_epithelium_um, fill = layer)) +
  geom_histogram() +
  facet_wrap(~layer, nrow = 2)

distance |> 
  filter(tissue == "Glans") |>
  ggplot(aes(x = group, y = distance_to_annotation_with_epithelium_um, color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_pvalue_manual(p_vals_glans |> filter(group1 == '4h', layer == '.', condition != '.'), 
                     label = 'p.adj.signif', y.position = 3, hide.ns = T,
                     tip.length = 0) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(condition ~ target_cell) + 
  scale_y_log10() +
  scale_y_reverse() +
  labs(title = "Glans Samples",
       y = expression ("Distance to annotation epithelium in "~mu*m)) +
  rotate_x_text(45)

################################################################################
### shaft ----
################################################################################

shaft <- distance |> 
  filter(tissue == 'Shaft')


model_shaft <- lme4::lmer(log10(distance_to_annotation_with_epithelium_um) ~ condition*group*layer*target_cell+(1|donor), data = shaft)


# is the model any good? -- better if we use log10 transform
summary(model_shaft)
qqnorm(resid(model_shaft))
qqline(resid(model_shaft))
hist(resid(model_shaft))
plot(fitted(model_shaft), resid(model_shaft))
abline(h = 0)

performance::check_model(model_shaft)


# glance at all
compare_all_pairs(model_shaft, factors = c('group', 'condition', 'target_cell', 'layer')) |>
  display_comparison_table(title = 'shaft model')


# here, we don't really care if certain cell types are more abundant than others, we 
# care if their abundance changes between conditions

comps_to_make <- c("group",
                   "layer",
                   "condition",
                   "condition | group",
                   "condition | group | target_cell",
                   "group | condition | target_cell",
                   "condition | group | target_cell | layer",
                   "group | condition | target_cell | layer")

comparisons_made <- compare_pairs(model_shaft, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')

display_comparison_table_zvalue(comparisons_made, title = "shaft model")


p_vals_shaft <- compare_pairs(model_shaft, 
                        comparisons = comps_to_make,
                        p_adjustment = 'fdr') |>
  get_pvals_for_comparisons()

# make same comparison we did with glans
distance |> 
  filter(tissue == "Shaft") |>
  ggplot(aes(x = group, y = distance_to_annotation_with_epithelium_um, color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_pvalue_manual(p_vals_shaft |> filter(group1 == '4h', layer == '.', condition != '.'), 
                     label = 'p.adj.signif', y.position = 3, hide.ns = T,
                     tip.length = 0) +
  theme_pubr() + scale_color_npg() + 
  facet_grid(condition ~ target_cell) + 
  scale_y_reverse() +
  labs(title = "Shaft Samples",
       y = expression ("Distance to annotation epithelium in "~mu*m)) +
  rotate_x_text(45)


