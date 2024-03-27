## MAST CELL ANALYSIS ------

# load packages

library(ggpubr)
library(ggsci)
library(gt)
library(emmeans)
library(tidyverse)

source("R_scripts/functions/pairwise_comparisons.R")

# load data

distance <- read_csv("data/Target cells R21 DISTANCE.csv") |>
  janitor::clean_names(replace=janitor:::mu_to_u) |> 
  mutate(donor = factor(donor),
         replicate = factor(replicate))


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



### looking at mean 

# mean across replicates 

# why bother taking a mean here? -- this is probably my biggest question for Ramon
# maybe makes it easier to model idk

mean_dist <- distance |>
  dplyr::group_by(condition, group, layer, donor, tissue, target_cell) |>
  dplyr::summarize(mean_distance_to_annotation_with_epithelium_um = mean(distance_to_annotation_with_epithelium_um))

mean_dist |> 
  ggplot(aes(x=layer,
             y=mean_distance_to_annotation_with_epithelium_um,
             color= interaction(group,condition)))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+
  scale_color_npg()+
  facet_wrap(~target_cell, scales = "free")

mean_dist |> 
  ggplot(aes(x=layer,
             y=mean_distance_to_annotation_with_epithelium_um+1,
             color=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+
  scale_color_npg()+
  facet_grid(tissue*condition~target_cell)+
  scale_y_log10()+
  rotate_x_text(45)

mean_dist |> 
  ggplot(aes(x=group,
             y=mean_distance_to_annotation_with_epithelium_um+1,
             color=layer)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_pubr()+
  scale_color_jco()+
  facet_grid(tissue*condition~target_cell)+
  scale_y_log10()


## modeling

################################################################################
### glans ----
################################################################################



glans <- mean_dist |> 
  filter(tissue == 'Glans')

# check distribution of outcome variable
ggplot(data = glans, aes(x = log10(mean_distance_to_annotation_with_epithelium_um))) +
  geom_histogram()

model_glans <- lme4::lmer(log10(mean_distance_to_annotation_with_epithelium_um) ~ condition*group*layer*target_cell+(1|donor), 
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
  display_comparison_table(title = 'glans model')


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


display_comparison_table(comparisons_made, title = 'glans table')

# here we can see significant differences between CD3 as well as CD4 single positive 
# cells in uncircumsized epithelium between infected and non infected-- a significant difference was not detected in corresponding circumsized tissue:


# can we... ?
distance |>
  filter(layer == 'Connective' & (target_cell == 'CD4' | target_cell == 'CD3') & tissue == 'Glans') |>
  ggplot(aes(x = group,
             y = distance_to_annotation_with_epithelium_um,
             color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.4) +
  scale_y_log10() +
  stat_compare_means(
    comparisons = list(c('4h', 'Non infected')),
    method = "wilcox.test", 
    label = "p.signif") +
  facet_grid(rows = vars(target_cell), cols = vars(condition), scales = 'free') +
  labs(title = "CD3/CD4 Single Positive cells in Connective Tissue (Glans)") +
  theme_pubr() + scale_color_lancet() 

distance |>
  filter(layer == 'Connective' & (target_cell == 'CD4' | target_cell == 'CCR10') & tissue == 'Glans') |>
  ggplot(aes(x = condition,
             y = distance_to_annotation_with_epithelium_um,
             color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.4) +
  scale_y_log10() +
  stat_compare_means(
    comparisons = list(c('Circumcised', 'Uncircumcised')),
    method = "wilcox.test", 
    label = "p.signif") +
  facet_grid(rows = vars(target_cell), cols = vars(group), scales = 'free') +
  labs(title = "[not mean] CD4 Single Positive cells in Connective Tissue (Glans)") +
  theme_pubr() + scale_color_lancet() 

mean_dist |>
  filter(layer == 'Connective' & (target_cell == 'CD4' | target_cell == 'CD3') & tissue == 'Glans') |>
  ggplot(aes(x = group,
             y = mean_distance_to_annotation_with_epithelium_um,
             color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.4) +
  scale_y_log10() +
  stat_compare_means(
    comparisons = list(c('4h', 'Non infected')),
    method = "wilcox.test", 
    label = "p.signif") +
  facet_grid(rows = vars(target_cell), cols = vars(condition), scales = 'free') +
  labs(title = "[mean] CD4 Single Positive cells in Connective Tissue (Glans)") +
  theme_pubr() + scale_color_lancet()  

mean_dist |>
  filter(layer == 'Epithelium' & target_cell == 'CD4' & group == '4h') |>
  ggplot(aes(x = condition,
             y = mean_distance_to_annotation_with_epithelium_um,
             color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_y_log10() +
  theme_pubr() + scale_color_lancet() 


################################################################################
### shaft ----
################################################################################

shaft <- mean_dist |> 
  filter(tissue == 'Shaft')


model_shaft <- lme4::lmer(log10(mean_distance_to_annotation_with_epithelium_um) ~ condition*group*layer*target_cell+(1|donor), data = shaft)


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

display_comparison_table(comparisons_made, title = "shaft model")


# Now we can check out those by making plots-- I'm not too concerned here with cell types that were just differently abundant overall, but if that does end up being relevant somehow we can go back to it.
# god i am on the same page as past sean <3 

mean_dist |>
  filter(layer == 'Connective' & (target_cell == 'CD3') & tissue == 'Shaft' & condition == 'Circumcised') |>
  ggplot(aes(x = group,
             y = mean_distance_to_annotation_with_epithelium_um,
             color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.4) +
  scale_y_log10() +
  stat_compare_means(
    comparisons = list(c('4h', 'Non infected')),
    method = "wilcox.test", 
    label = "p.signif") +
  facet_grid(rows = vars(target_cell), cols = vars(condition), scales = 'free') +
  labs(title = "[mean] CD3 Single Positive cells in Connective Tissue (Shaft)",
       caption = 'pvalue is signif in emmeans model, might not mean much tho bc few data points') +
  theme_pubr() + scale_color_lancet() 

mean_dist |>
  filter(layer == 'Epithelium' & (target_cell == 'CD4+CCR10') & tissue == 'Shaft') |>
  ggplot(aes(x = group,
             y = mean_distance_to_annotation_with_epithelium_um,
             color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.4) +
  scale_y_log10() +
  stat_compare_means(
    comparisons = list(c('4h', 'Non infected')),
    method = "wilcox.test", 
    label = "p.signif") +
  facet_grid(rows = vars(target_cell), cols = vars(condition), scales = 'free') +
  labs(title = "[mean] CD4+CCR10+ Positive cells in Epithelial Tissue (Shaft)",
       caption = 'pvalue for circumcised is signif in emmeans model, might not mean much tho bc few data points') +
  theme_pubr() + scale_color_lancet() 

