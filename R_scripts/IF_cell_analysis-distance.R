## IMMUNOFLUORESCENCE ANALYSIS ------

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
# seems to be what tom wants <3

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
  ggplot(aes(x=group,
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
  filter(tissue == 'Glans') |>
  dplyr::rename("treat_group" = "group")
  

# check distribution of outcome variable
ggplot(data = glans, aes(x = log10(mean_distance_to_annotation_with_epithelium_um))) +
  geom_histogram()

model_glans <- lme4::lmer(log10(mean_distance_to_annotation_with_epithelium_um) ~ condition*treat_group*layer*target_cell+(1|donor), 
                          data = glans)

# is the model any good? -- better if we use log10 transform
summary(model_glans)
qqnorm(resid(model_glans))
qqline(resid(model_glans))
hist(resid(model_glans))
plot(fitted(model_glans), resid(model_glans))
abline(h = 0)

performance::check_model(model_glans)


# here, we don't really care if certain cell types are more abundant than others, we 
# care if their abundance changes between conditions

comps_to_make <- c("treat_group",
                   "layer",
                   "condition",
                   #"target_cell",
                   "condition | treat_group",
                   "condition | treat_group | target_cell",
                   "treat_group | condition | target_cell",
                   "condition | treat_group | target_cell | layer",
                   #"target_cell | condition | treat_group | layer",
                   "treat_group | condition | target_cell | layer")

comparisons_made <- compare_pairs(model_glans, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')


display_comparison_table(comparisons_made, title = 'glans table')


p_vals <- comparisons_made |>
  get_pvals_for_comparisons()

# here we can see significant differences between CD3 as well as CD4 single positive 
# cells in uncircumcised epithelium between infected and non infected-- a significant difference was not detected in corresponding circumsized tissue:


glans_if_plot <- glans |> 
  ggplot(aes(x=treat_group,
             y=mean_distance_to_annotation_with_epithelium_um+1,
             color=treat_group))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_log10() +
  geom_point(position = position_jitterdodge())+
  facet_grid(layer*condition~target_cell)+
  ggprism::add_pvalue(p_vals |> filter(contrast == '4h - Non infected', layer != '.', p.value < 0.05), 
                      label = 'p.adj = {round(p.value, 4)}', y.position = 2000, color = "black", bracket.size = 0.1) +
  theme_pubr() +
  labs(title = "Glans Samples") +
  coord_cartesian(ylim = c(1,3000)) +
  rotate_x_text(45)

save(glans_if_plot, file = "figures/plots/glans_if_plot.rda")


################################################################################
### shaft ----
################################################################################

shaft <- mean_dist |> 
  filter(tissue == 'Shaft') |>
  dplyr::rename("treat_group" = "group")


model_shaft <- lme4::lmer(log10(mean_distance_to_annotation_with_epithelium_um) ~ condition*treat_group*layer*target_cell+(1|donor), data = shaft)


# is the model any good? -- better if we use log10 transform
summary(model_shaft)
qqnorm(resid(model_shaft))
qqline(resid(model_shaft))
hist(resid(model_shaft))
plot(fitted(model_shaft), resid(model_shaft))
abline(h = 0)

performance::check_model(model_shaft)


# here, we don't really care if certain cell types are more abundant than others, we 
# care if their abundance changes between conditions

comps_to_make <- c("treat_group",
                   "layer",
                   "condition",
                   "condition | treat_group",
                   "condition | treat_group | target_cell",
                   "treat_group | condition | target_cell",
                   "condition | treat_group | target_cell | layer",
                   "treat_group | condition | target_cell | layer")

comparisons_made <- compare_pairs(model_shaft, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')

display_comparison_table(comparisons_made, title = "shaft model")

################################################################################
### by layer ----
################################################################################

# idrc about layer?? like that is the outcome variable LOL

connective <- mean_dist |> 
  dplyr::filter(layer == "Connective") |>
  dplyr::rename("treat_group" = "group")


# check distribution of outcome variable
ggplot(data = connective, aes(x = log10(mean_distance_to_annotation_with_epithelium_um))) +
  geom_histogram()

model_connective <- lme4::lmer(log10(mean_distance_to_annotation_with_epithelium_um) ~ condition*treat_group*target_cell*tissue+(1|donor), 
                          data = connective)

# is the model any good? -- better if we use log10 transform
summary(model_connective)
qqnorm(resid(model_connective))
qqline(resid(model_connective))
hist(resid(model_connective))
plot(fitted(model_connective), resid(model_connective))
abline(h = 0)

performance::check_model(model_connective)


# glance at all
compare_all_pairs(model_connective, factors = c('treat_group', 'condition', 'target_cell', "tissue")) |>
  display_comparison_table(title = 'connective model')


# here, we don't really care if certain cell types are more abundant than others, we 
# care if their abundance changes between conditions

comps_to_make <- c("treat_group",
                   "tissue",
                   "condition",
                   #"target_cell",
                   "condition | treat_group",
                   "condition | treat_group | target_cell",
                   "treat_group | condition | target_cell",
                   "condition | treat_group | target_cell | tissue",
                   "tissue | treat_group | target_cell | condition",
                   #"target_cell | condition | treat_group | layer",
                   "treat_group | condition | target_cell | tissue")

comparisons_made <- compare_pairs(model_connective, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')


display_comparison_table(comparisons_made, title = 'connective table')


p_vals <- comparisons_made |>
  get_pvals_for_comparisons()

# here we can see significant differences between CD3 as well as CD4 single positive 
# cells in uncircumcised epithelium between infected and non infected-- a significant difference was not detected in corresponding circumsized tissue:


connective_layer_plot <- connective |> 
  ggplot(aes(x=treat_group,
             y=mean_distance_to_annotation_with_epithelium_um+1,
             color=treat_group))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_log10() +
  geom_point(position = position_jitterdodge())+
  facet_grid(tissue*condition~target_cell)+
  ggprism::add_pvalue(p_vals |> filter(contrast == '4h - Non infected', tissue != '.', p.value < 0.05), 
                      label = 'p.adj = {round(p.value, 4)}', y.position = 2000, color = "black", bracket.size = 0.1) +
  theme_pubr() +
  coord_cartesian(ylim = c(1,3000)) +
  labs(title = "Cells in the Connective Layer") +
  rotate_x_text(45)

save(connective_layer_plot, file = "figures/plots/connective_layer_plot.rda")

################################################################################
### by layer: epithelium ----
################################################################################

# idrc about layer?? like that is the outcome variable LOL

epithelium <- mean_dist |> 
  dplyr::filter(layer == "Epithelium") |>
  dplyr::rename("treat_group" = "group")


# check distribution of outcome variable
ggplot(data = epithelium, aes(x = log10(mean_distance_to_annotation_with_epithelium_um))) +
  geom_histogram()

model_epithelium <- lme4::lmer(log10(mean_distance_to_annotation_with_epithelium_um) ~ condition*treat_group*target_cell*tissue+(1|donor), 
                               data = epithelium)

# is the model any good? -- better if we use log10 transform
summary(model_epithelium)
qqnorm(resid(model_epithelium))
qqline(resid(model_epithelium))
hist(resid(model_epithelium))
plot(fitted(model_epithelium), resid(model_epithelium))
abline(h = 0)

performance::check_model(model_epithelium)


# here, we don't really care if certain cell types are more abundant than others, we 
# care if their abundance changes between conditions

comps_to_make <- c("treat_group",
                   "tissue",
                   "condition",
                   #"target_cell",
                   "condition | treat_group",
                   "condition | treat_group | target_cell",
                   "treat_group | condition | target_cell",
                   "condition | treat_group | target_cell | tissue",
                   "tissue | treat_group | target_cell | condition",
                   #"target_cell | condition | treat_group | layer",
                   "treat_group | condition | target_cell | tissue")

comparisons_made <- compare_pairs(model_epithelium, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')


display_comparison_table(comparisons_made, title = 'epithelium table')

# nothing. great!

################################################################################
### uncircumcised ----
################################################################################

# idrc about layer?? like that is the outcome variable LOL

uncut <- mean_dist |> 
  dplyr::filter(condition == "Uncircumcised") |>
  dplyr::rename("treat_group" = "group")


# check distribution of outcome variable
ggplot(data = uncut, aes(x = log10(mean_distance_to_annotation_with_epithelium_um))) +
  geom_histogram()

model_uncut <- lme4::lmer(log10(mean_distance_to_annotation_with_epithelium_um) ~ layer*treat_group*target_cell*tissue+(1|donor), 
                               data = uncut)

# is the model any good? -- better if we use log10 transform
summary(model_uncut)
qqnorm(resid(model_uncut))
qqline(resid(model_uncut))
hist(resid(model_uncut))
plot(fitted(model_uncut), resid(model_uncut))
abline(h = 0)

performance::check_model(model_uncut)


# glance at all
compare_all_pairs(model_uncut, factors = c('treat_group', 'layer', 'target_cell', "tissue")) |>
  display_comparison_table(title = 'uncircumcised model')


# here, we don't really care if certain cell types are more abundant than others, we 
# care if their abundance changes between conditions

comps_to_make <- c("treat_group",
                   "tissue",
                   #"layer",
                   #"target_cell",
                   #"layer | treat_group",
                   #"layer | treat_group | target_cell",
                   "treat_group | tissue | target_cell",
                   #"layer | treat_group | target_cell | tissue",
                   "tissue | treat_group | target_cell | layer",
                   #"target_cell | condition | treat_group | layer",
                   "treat_group | layer | target_cell | tissue"
                   )

comparisons_made <- compare_pairs(model_uncut, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')


display_comparison_table(comparisons_made, title = 'Uncircumcised table')


p_vals <- comparisons_made |>
  get_pvals_for_comparisons()


uncut_if_plot <- uncut |> 
  ggplot(aes(x=treat_group,
             y=mean_distance_to_annotation_with_epithelium_um+1,
             color=treat_group))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_log10() +
  geom_point(position = position_jitterdodge())+
  facet_grid(layer*tissue~target_cell)+
  ggprism::add_pvalue(p_vals |> filter(contrast == '4h - Non infected', tissue != '.', p.value < 0.05), 
                      label = 'p.adj = {round(p.value, 4)}', y.position = 2000, color = "black", bracket.size = 0.1) +
  theme_pubr() +
  coord_cartesian(ylim = c(1,3000)) +
  labs(title = "Uncircumcised Samples") +
  rotate_x_text(45)

save(uncut_if_plot, file = "figures/plots/uncut_if_plot.rda")


################################################################################
### circumcised ----
################################################################################


cut <- mean_dist |> 
  dplyr::filter(condition == "Circumcised") |>
  dplyr::rename("treat_group" = "group")


# check distribution of outcome variable
ggplot(data = cut, aes(x = log10(mean_distance_to_annotation_with_epithelium_um))) +
  geom_histogram()

model_cut <- lme4::lmer(log10(mean_distance_to_annotation_with_epithelium_um) ~ layer*treat_group*target_cell*tissue+(1|donor), 
                          data = cut)

# is the model any good? -- better if we use log10 transform
summary(model_cut)
qqnorm(resid(model_cut))
qqline(resid(model_cut))
hist(resid(model_cut))
plot(fitted(model_cut), resid(model_cut))
abline(h = 0)

performance::check_model(model_cut)


# glance at all
compare_all_pairs(model_cut, factors = c('treat_group', 'layer', 'target_cell', "tissue")) |>
  display_comparison_table(title = 'circumcised model')


# here, we don't really care if certain cell types are more abundant than others, we 
# care if their abundance changes between conditions

comps_to_make <- c("treat_group",
                   "tissue",
                   #"layer",
                   #"target_cell",
                   #"layer | treat_group",
                   #"layer | treat_group | target_cell",
                   "treat_group | tissue | target_cell",
                   #"layer | treat_group | target_cell | tissue",
                   "tissue | treat_group | target_cell | layer",
                   #"target_cell | condition | treat_group | layer",
                   "treat_group | layer | target_cell | tissue"
)

comparisons_made <- compare_pairs(model_cut, 
                                  comparisons = comps_to_make,
                                  p_adjustment = 'fdr')


display_comparison_table(comparisons_made, title = 'Circumcised table')


p_vals <- comparisons_made |>
  get_pvals_for_comparisons()

# here, only salient difference is between glans/shaft 
cut_if_plot <- cut |> 
  ggplot(aes(x=tissue,
             y=mean_distance_to_annotation_with_epithelium_um+1,
             color=tissue))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_log10() +
  geom_point(position = position_jitterdodge())+
  facet_grid(treat_group~target_cell)+
  ggprism::add_pvalue(p_vals |> filter(contrast == 'Glans - Shaft', treat_group != '.', p.value < 0.05), 
                      label = 'p.adj = {round(p.value, 4)}', y.position = 2000, color = "black", bracket.size = 0.1) +
  theme_pubr() +
  coord_cartesian(ylim = c(1,3000)) +
  labs(title = "Circumcised Samples") +
  rotate_x_text(45)

save(cut_if_plot, file = "figures/plots/uncut_if_plot.rda")