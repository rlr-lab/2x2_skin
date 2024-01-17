######################################################################################  
# LOADING + CLEANING DATA ----
######################################################################################  


# load packages ----
library(tidyverse)

# does cell population change between circumsized and not? glans (<-) vs shaft (should be control?) ? 
# tissue infected in lab -> paired samples

######################################################################################  
# experiments x donor ----
######################################################################################  


# information about which experiments were performed with which donors
exp_donor <- readxl::read_excel("data/ExperimentsxDonor.xlsx")

######################################################################################  
# Information about tryptase/cd117 staining from pathology core
######################################################################################  

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

save(tryptase, file = 'data/processed/tryptase.rda')

######################################################################################  
# Target cells R21
######################################################################################  

distance <- read_csv("data/Target cells R21 DISTANCE.csv") |>
  janitor::clean_names(replace=janitor:::mu_to_u) |> 
  mutate(donor = factor(donor),
         replicate = factor(replicate))

save(distance, file = 'data/processed/distance.rda')

# counts and some summary stats about area of tissue section for above

# the thing is this data isn't super useful, and I don't think we found any differences
# previously-- if you're just looking at counts of immune cells, those should be
# randomly distributed based on the section of penis that was shaved off
target_total <- read_csv("data/Target cells R21 TOTAL COUNT.csv") |>
  janitor::clean_names(replace=janitor:::mu_to_u) 

######################################################################################  
# virus counts  ----
######################################################################################  

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

## depth

Depth<-virus_counts[,c(1:4,9:59)]

Depth<-Depth%>%
  dplyr::rename("x9" = "penetration_depths_um") %>%
  pivot_longer(x9:x59,
               names_to = "particle",
               values_to = "penetration_depth") |>
  dplyr::select(!particle) |>
  drop_na(penetration_depth) |>
  mutate(penetration_depth = as.numeric(as.character(penetration_depth)))

save(virus_counts, Depth, file = 'data/processed/virus_counts.rda')

######################################################################################  
# viral penetration ~ target cell distance ----
######################################################################################  

# try with just mean depth first-- this will be easiest

mean_depth <- Depth %>%
  group_by(condition, tissue, donor) %>%
  summarise(mean_depth = mean(penetration_depth))

cell_and_virus_distance <- distance |>
  filter(group == '4h') |> # only infected samples will have viral counts
  left_join(y = mean_depth, by = c('donor', 'tissue', 'condition')) |> # 14% missingness
  drop_na(mean_depth)

save(cell_and_virus_distance, file = 'data/processed/cell_and_virus_distance.rda')


