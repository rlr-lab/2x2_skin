### LOADING DATA ----

# load packages ----
library(tidyverse)

# does cell population change between circumsized and not? glans (<-) vs shaft (should be control?) ? 
# tissue infected in lab -> paired samples


# experiments x donor ----

exp_donor <- readxl::read_excel("data/ExperimentsxDonor.xlsx")


# R21_path core

r21 <- read_csv("data/R21_PathCore.csv")


# Target cells R21

target_dist <- read_csv("data/Target cells R21 DISTANCE.csv") |>
  janitor::clean_names(replace=janitor:::mu_to_u)

# counts and some summary stats about area of tissue section for above

target_total <- read_csv("data/Target cells R21 TOTAL COUNT.csv") |>
  janitor::clean_names(replace=janitor:::mu_to_u) 


# virus count R21 ----

virus <- read_csv("data/Virus count - R21.csv")
