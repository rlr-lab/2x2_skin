######################################################################################  
# VIRUS COUNT ANALYSIS ----
######################################################################################  

# install packages

library(ggpubr)
library(pscl)
library(glmmTMB)
library(tidyverse)
library(emmeans)
library(fitdistrplus)
library(performance)


source("R_scripts/functions/pairwise_comparisons.R")

# load data

load("data/processed/virus_counts.rda")

######################################################################################  
## distribution check

# based on some distribution evaluation things, it looks like lognormal is gonna fit our data the best:
  
plotdist(virus_counts$virions, histo = TRUE, demp = TRUE)
descdist(virus_counts$virions, discrete=FALSE, boot=500)
# i guess this is a mostly lognormal distribution??
fnb<-fitdist(virus_counts$virions, "nbinom")
fp<-fitdist(virus_counts$virions, "pois")
plot.legend <- c( "nbinom", "pois")
denscomp(list(fnb,fp), legendtext = plot.legend)
qqcomp(list(fnb,fp), legendtext = plot.legend)
# based on qqplot-- poisson sucks
cdfcomp(list(fnb,fp), legendtext = plot.legend)
ppcomp(list(fnb,fp), legendtext = plot.legend)
# idk what a pp plot is either but im just gonna say poisson has not been killing it here

gofstat(list(fnb,fp))

###########################################  
# Virion Counts ----
###########################################  

#so let's make a lognormal model!


initial_model <- lme4::lmer(virions ~ condition*tissue+(1|donor), virus_counts)

plot(performance::check_distribution(initial_model))


# ZERO INflated binomial
zinbm0<-glmmTMB(virions ~ condition*tissue+(1|donor), 
                family="nbinom2",
                ziformula=~1,data = virus_counts)

# negative binomial
zinbm1<-glmmTMB(virions ~ tissue*condition+(1|donor), 
                family="nbinom2",
                data = virus_counts)

# zero inflation, controlling for tissue
zinbm2<-glmmTMB(virions ~ condition*tissue+(1|donor), 
                family="nbinom2",
                ziformula=~tissue,data = virus_counts)

anova(zinbm0,zinbm1,zinbm2)
performance::model_performance(zinbm0)
performance::check_model(zinbm0)
performance::compare_performance(zinbm0,zinbm1,zinbm2)

######################################################################################  
# Then we can compare EMMs of our best model, zinbm0 :

# *SEAN:* is binom best here? idk lol... i need to figure out how this works

t <- pairs(emmeans(zinbm0, ~ condition))
s <- pairs(emmeans(zinbm0, ~ tissue))
t.s <- pairs(emmeans(zinbm0, ~ tissue | condition))
s.t <- pairs(emmeans(zinbm0, ~ condition | tissue))

p_vals <- rbind(t.s, s.t,adjust="fdr") |>
  as_tibble() |>
  mutate(group = str_split(contrast, pattern = " - ")) |>
  unnest_wider(group, names_sep = "")

p_vals_tissue <- p_vals |>
  filter(condition == ".") |>
  rename("supp" = "tissue") |>
  mutate(p.value = round(p.value, digits= 2)) |>
  mutate(p.signif = case_when(
    p.value > 0.05 ~ 'ns',
    p.value < 0.01 ~ '*'
  ))

p_vals_condition <- p_vals |>
  filter(tissue == ".") |>
  rename("supp" = "condition") |>
  mutate(p.value = round(p.value, digits= 2)) |>
  mutate(p.signif = case_when(
    p.value > 0.05 ~ 'ns',
    p.value < 0.05 ~ '*'
  ))

check_model(zinbm0)


virus_counts |>
  ggplot(aes(x = tissue, y= virions, color = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  scale_y_log10() +
  theme_pubr()+scale_color_npg() +
  facet_wrap(~condition) +
  # need p values on here
  stat_pvalue_manual(p_vals_condition, label = 'p.signif', y.position = c(100))


ggboxplot(virus_counts, x = "condition", 
          y = "virions",
          combine = TRUE,
          color = "tissue", palette = "npg",
          ylab = "virions", 
          add = "jitter",                              # Add dotplot
          add.params = list(binwidth = 0.2))+
  scale_y_log10() +
  #stat_pvalue_manual(p_vals_tissue, label = 'p.value', y.position = c(2.1,1), hide.ns = TRUE) +
  stat_pvalue_manual(p_vals_condition, label = 'p.signif', x = 'supp', y.position = c(200, 200), hide.ns = FALSE) +
  labs(title = 'Comparison of virion count in different tissue samples',
    caption = 'p values computed using pairwise comparisons in emmeans and adjusted by FDR')

###########################################  
## looking at means
###########################################  

# look at some means of viral counts instead of just looking at all of them:

mean_virus_counts <-virus_counts %>%
  group_by(condition, tissue, donor) %>%
  summarise(mean_virions = mean(virions)) 

mean_model <- lme4::lmer(log10(mean_virions) ~ condition*tissue+(1|donor), 
                    data = mean_virus_counts)

plot(check_distribution(mean_model))
t <- pairs(emmeans(mean_model, ~ condition))
s <- pairs(emmeans(mean_model, ~ tissue))
t.s <- pairs(emmeans(mean_model, ~ tissue | condition))
s.t <- pairs(emmeans(mean_model, ~ condition | tissue))
rbind(t.s, s.t,adjust="fdr")

# no significance!

check_model(mean_model)
qqnorm(resid(mean_model))
qqline(resid(mean_model))
hist(resid(mean_model))
plot(fitted(mean_model),resid(mean_model))
abline(h=0)

mean_virus_counts <- mean_virus_counts %>%
  ggplot(aes(x=condition,y=mean_virions,color=tissue))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+scale_color_npg()+
  facet_grid(.~tissue, scales = "free")

save(mean_virus_counts, file = "figures/plots/mean_virus_counts.rda")
# not really seeing any significant differences when we just look at virion counts overall-- is there something we can do to analyze the depths that these virions penetrated?
# this maybe makes sense tho because I think % penetrating is much more interesting anyways?
# why were we modeling counts lmao

################################################
# Penetrators ----
################################################
Penetrators <- virus_counts |>
  dplyr::select(condition, tissue, percent_of_penetrating_virions)


## penetrator eda ----

#  lets go back to basic EDA-- what does our virion number data even look like

virus_counts |>
  ggplot(aes(x = percent_of_penetrating_virions, fill = condition)) +
  geom_histogram(bins = 100) +
  facet_grid(condition~tissue) +
  theme_pubr()+scale_color_npg()


# ok right this is why we fit with a negative binomial. I am ok with this i guess-- it's not really like we can just throw out the zeros...

#But what if we JUST look at the zeros?
# no this doesnt make sense. unless maybe it does

# Here, let's look at the percent of each category that have ZERO penetration (not sure if i can do this-- ask ramon)

totals <- virus_counts |>
  dplyr::group_by(tissue, condition) |>
  dplyr::summarize(totals = n())

rejected_virions_plot <- virus_counts |>
  filter(percent_of_penetrating_virions == 0 & virions != 0) |>
  dplyr::group_by(tissue, condition) |>
  dplyr::summarize(counts = n()) |>
  left_join(y = totals, by = c('tissue', 'condition')) |>
  mutate(percent_rejected = (counts*100)/totals) |>
  ggplot(aes(x = tissue, y = percent_rejected, fill = condition)) +
  geom_col(position = 'dodge', width = 0.6)+
  labs(y = "percent of images with virions \nbut with no penetrating virions") +
  theme_pubclean()+scale_color_npg()

save(rejected_virions_plot, file ="figures/plots/rejected_virions_plot.rda")

rejected_virions <- virus_counts |>
  mutate(rejected_virions = virions - penetrating_virions,
         rejected_ratio = rejected_virions/penetrating_virions)

# what if i plot counts of rejected virions? does this make sense? 

rejected_virions |>
  filter(virions != 0 & penetrating_virions != 0) |>
  ggplot(aes(x = tissue, y = rejected_ratio, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2))

# yeah ok i think what im doing here is garBAGE but the model should be good enough because that takes into account the fact that we're gonna have zeros or whatever and also its giving me p values so. fun investigation but nothing incredibly shocking here. we should get same results doing percent penetrating and percent rejected bc they're proxies for one another.

# maybe see higher of pct of penetrating virions in uncircumcised?

# NEED TO GET A P VALUE ON THIS!!
virus_counts%>%
  ggplot(aes(x=condition,y=percent_of_penetrating_virions,color=condition))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+scale_color_npg()+facet_grid(.~tissue,scales = "free")

################################
## Penetrator modeling ----
################################


lm_penetrators <- gamlss::gamlss((percent_of_penetrating_virions/100) ~ condition*tissue,
                                 # beta inflated something something
                                 family = gamlss.dist::BEINF(), 
                                 data = Penetrators)
qqnorm(resid(lm_penetrators))
qqline(resid(lm_penetrators))
hist(resid(lm_penetrators))
plot(fitted(lm_penetrators),resid(lm_penetrators))
abline(h=0)

plot(performance::check_distribution(lm_penetrators))

# Check on this:
# binned_residuals(lm_penetrators, residuals = "response")
# check_model(lm_penetrators)

comparisons <- compare_all_pairs(lm_penetrators, 
                  factors = c('tissue', 'condition'),
                  p_adjustment = 'fdr') 

comparisons |>
  display_comparison_table_zvalue()

p_vals <- comparisons |> get_pvals_for_comparisons()


# here we see some significant differences between circumcised and uncircumcised-- specifically in shaft, but maybe in both? Yes definitely-- although there might be something funky about how I'm manipulating p values here by choosing which comparisons I want to make. 

Penetrators |>
  ggplot(aes(x = condition, y= percent_of_penetrating_virions, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  theme_pubr()+scale_color_npg() +
  stat_pvalue_manual(p_vals[4,], label = 'p.adj.signif', y.position = c(110)) +
  labs(title = "Uncircumcised samples are more likely to have higher numbers of penetrating virions") 



penetrators_nomean_plot <- Penetrators |>
  ggplot(aes(x = condition, y= percent_of_penetrating_virions, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  theme_pubr()+scale_color_npg() +
  facet_wrap(~tissue) +
  stat_pvalue_manual(p_vals[5:6,], label = 'p = {round(p.value, 4)}', y.position = c(105,105)) +
  coord_cartesian(ylim = c(0,110))

save(penetrators_nomean_plot, file = "figures/plots/penetrators_nomean_plot.rda")

# Ok this is cool... but there are a lot of 100%s and 0%s in either case. This might be worth investigating-- who is at these extremes?


# are values the same when we do mean?

mean_penetrators <-virus_counts %>%
  group_by(condition, tissue, donor) %>%
  summarise(mean_percent_of_penetrating_virions = mean(percent_of_penetrating_virions)) 

lm_mean_penetrators <- gamlss::gamlss((mean_percent_of_penetrating_virions/100) ~ condition*tissue,
                                 # beta inflated something something
                                 family = gamlss.dist::BEINF(), 
                                 data = mean_penetrators)

qqnorm(resid(lm_mean_penetrators))
qqline(resid(lm_mean_penetrators))
hist(resid(lm_mean_penetrators))
plot(fitted(lm_mean_penetrators),resid(lm_mean_penetrators))
abline(h=0)

plot(performance::check_distribution(lm_mean_penetrators))

# Check on this:
# binned_residuals(lm_mean_penetrators, residuals = "response")
# check_model(lm_mean_penetrators)

comparisons <- compare_all_pairs(lm_mean_penetrators, 
                                 factors = c('tissue', 'condition'),
                                 p_adjustment = 'fdr') 


comparisons |>
  display_comparison_table_zvalue()

p_vals <- comparisons |> get_pvals_for_comparisons()


# here we see some significant differences between circumcised and uncircumcised-- specifically in shaft, but maybe in both? Yes definitely-- although there might be something funky about how I'm manipulating p values here by choosing which comparisons I want to make. 

mean_penetrators |>
  ggplot(aes(x = condition, y= mean_percent_of_penetrating_virions, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  theme_pubr()+scale_color_npg() +
  stat_pvalue_manual(p_vals[4,], label = 'p.adj.signif', y.position = c(110)) +
  labs(title = "Uncircumcised samples are more likely to have higher numbers of penetrating virions") 



mean_penetrators_plot <- mean_penetrators |>
  ggplot(aes(x = condition, y= mean_percent_of_penetrating_virions, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  theme_pubr()+scale_color_npg() +
  facet_wrap(~tissue) +
  stat_pvalue_manual(p_vals[5:6,], label = 'p = {round(p.value, 4)}', y.position = c(105,105)) +
  coord_cartesian(ylim = c(0,110))

save(mean_penetrators_plot, file = "figures/plots/mean_penetrators_plot.rda")

## VERSION WITH P HACKING

comparisons <- compare_all_pairs(lm_mean_penetrators, 
                                 factors = c('tissue', 'condition'),
                                 p_adjustment = 'none') 

p_vals <- comparisons |> get_pvals_for_comparisons()

# when we start p hacking, the data looks exactly like you want it to!


mean_penetrators_plot_phack <- mean_penetrators |>
  ggplot(aes(x = condition, y= mean_percent_of_penetrating_virions, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  theme_pubr()+scale_color_npg() +
  facet_wrap(~tissue) +
  stat_pvalue_manual(p_vals[5:6,], label = 'p = {round(p.value, 4)}', y.position = c(105,105)) +
  coord_cartesian(ylim = c(0,110))

save(mean_penetrators_plot_phack, file = "figures/plots/mean_penetrators_plot_phack.rda")

################################################
# Depth ----
################################################


# analysis of how far the virions penetrated -- doesn't really seem to be a difference here

# just take a peek!

depth_distribution <- Depth |>
  ggplot(aes(y = tissue, color = condition, x = penetration_depth)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 0.1)+
  geom_boxplot(outlier.shape = NA, fill = NA)+
  theme_pubr()+scale_color_npg() +
  labs(title = 'Penetration Depth Distributions are Similar for Each Group')

save(depth_distribution, file = "figures/plots/depth_distribution.rda")

mean_depth <- Depth %>%
  group_by(condition, tissue, donor) %>%
  summarise(mean_depth = mean(penetration_depth))


mean_depth%>%
  ggplot(aes(x=condition,y=mean_depth,color=tissue))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+scale_color_npg()+
  facet_grid(.~tissue,scales = "free")

fl<-fitdist(mean_depth$mean_depth, "lnorm")
fn<-fitdist(mean_depth$mean_depth, "norm")
fg<-fitdist(mean_depth$mean_depth, "gamma")
plot.legend <- c("normal", "gamma", "lognorm")
denscomp(list(fn,fg,fl), legendtext = plot.legend)
qqcomp(list(fn,fg,fl), legendtext = plot.legend)
cdfcomp(list(fn,fg,fl), legendtext = plot.legend)
ppcomp(list(fn,fg,fl), legendtext = plot.legend)
gofstat(list(fn,fg,fl))

# I don't know if you actually need the random effects term if you're doing mean depth here
lm_mean_depth <- lme4::lmer(log10(mean_depth) ~ condition*tissue + (1|donor), data = mean_depth)
plot(check_distribution(lm_mean_depth))

binned_residuals(lm_mean_depth)
check_model(lm_mean_depth)

comparisons <- compare_all_pairs(lm_mean_depth, 
                                 factors = c('tissue', 'condition'),
                                 p_adjustment = 'fdr') 


comparisons |>
  display_comparison_table()

p_vals <- comparisons |> get_pvals_for_comparisons()



mean_depth_plot <- mean_depth |>
  ggplot(aes(x=condition,y=mean_depth,color=tissue)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+scale_color_npg()+
  facet_grid(.~tissue,scales = "free") +
  labs(title = "Mean depth of penetrating virions does not vary by individual") +
  stat_pvalue_manual(p_vals[5:6,], label = 'p = {round(p.value, 4)}', y.position = c(20,20)) +
  coord_cartesian(ylim = c(0,22))
  
save(mean_depth_plot, file ="figures/plots/mean_depth_plot.rda")

# try again with random effects but without mean

fl<-fitdist(Depth$penetration_depth, "lnorm")
fn<-fitdist(Depth$penetration_depth, "norm")
fg<-fitdist(Depth$penetration_depth, "gamma")
plot.legend <- c("normal", "gamma", "lognorm")
denscomp(list(fn,fg,fl), legendtext = plot.legend)
qqcomp(list(fn,fg,fl), legendtext = plot.legend)
cdfcomp(list(fn,fg,fl), legendtext = plot.legend)
ppcomp(list(fn,fg,fl), legendtext = plot.legend)
gofstat(list(fn,fg,fl))

lm_depth <- lme4::lmer(log10(penetration_depth) ~ condition*tissue +(1|donor), data = Depth)
plot(check_distribution(lm_depth))

binned_residuals(lm_depth)
check_model(lm_depth)

t <- pairs(emmeans(lm_depth, ~ condition))
s <- pairs(emmeans(lm_depth, ~ tissue))
t.s <- pairs(emmeans(lm_depth, ~ tissue | condition))
s.t <- pairs(emmeans(lm_depth, ~ condition | tissue))
rbind(t.s,s.t,adjust="FDR")


qqnorm(resid(lm_depth))
qqline(resid(lm_depth))
hist(resid(lm_depth))
plot(fitted(lm_depth),resid(lm_depth))
abline(h=0)


Depth |>
  ggplot(aes(color = condition, x = penetration_depth)) +
  geom_freqpoly() +
  facet_wrap(~tissue) +
  theme_pubr()+scale_color_npg() 

#no significant p values here.. but also a pretty bad fit of the data with the model

# i just can't convince myself that there is a difference. 
# i got nothing <3 hey siri play not strong enough by boygenius

