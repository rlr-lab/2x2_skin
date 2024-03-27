

# install packages

library(ggpubr)
library(ggsci)
library(pscl)
library(glmmTMB)
library(tidyverse)
library(emmeans)
library(fitdistrplus)


source("R_scripts/functions/pairwise_comparisons.R")

# read in data
virus_counts <- read_csv("data/Virus count - R21.csv") |>
  janitor::clean_names(replace=janitor:::mu_to_u) |> 
  mutate(donor = factor(donor),
         image = factor(image),
         virions = as.numeric(as.character(virions)),
         penetrating_virions = as.numeric(as.character(penetrating_virions)),
         percent_of_penetrating_virions = as.numeric(as.character(percent_of_penetrating_virions)))

# making some tactical deletions-- no percent of penetrating greater than 100; no donor 1905; no NA virions

virus_counts <- virus_counts |> 
  filter(donor != "1905") |> 
  filter(percent_of_penetrating_virions < 101) |>
  filter(!is.na(virions))


################################################
# Depth ----
################################################


# analysis of how far the virions penetrated -- doesn't really seem to be a difference here


Depth <- virus_counts[,c(1:4,9:59)] %>%
  dplyr::rename("x9" = "penetration_depths_um") %>%
  pivot_longer(x9:x59,
               names_to = "particle",
               values_to = "penetration_depth") |>
  dplyr::select(!particle) |>
  drop_na(penetration_depth) |>
  mutate(penetration_depth = as.numeric(as.character(penetration_depth)))

# just take a peek!

Depth |>
  ggplot(aes(y = tissue, color = condition, x = penetration_depth)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 0.1)+
  geom_boxplot(outlier.shape = NA, fill = NA)+
  theme_pubr()+scale_color_npg() +
  labs(title = 'Penetration Depth Distributions are Similar for Each Group')

ggplot(Depth, aes(x = penetration_depth)) +
  geom_histogram()

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

plot(performance::check_distribution(lm_depth))

qqnorm(resid(lm_depth))
qqline(resid(lm_depth))
hist(resid(lm_depth))
plot(fitted(lm_depth),resid(lm_depth))
abline(h=0)


compare_all_pairs(lm_depth, factors = c('condition', 'tissue'),
                  p_adjustment = 'fdr') |>
  display_comparison_table_zvalue()



Depth |>
  ggplot(aes(color = condition, x = penetration_depth)) +
  geom_freqpoly() +
  facet_wrap(~tissue) +
  theme_pubr()+scale_color_npg() 

#no significant p values here.. maybe just no penetration difference after only 4 hours
# i just can't convince myself that there is a difference. 
# i got nothing <3 hey siri play not strong enough by boygenius

