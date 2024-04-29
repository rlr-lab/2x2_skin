## AUXILARY FUNCTION: PERFORM ALL PAIRWISE COMBINATIONS WITH EMMEANS ------

# load packages
library(tidyverse)
library(gt)
library(emmeans)


################################################################################

# FIRST: make a function that, given a level A and n other levels B, C, D, ...
# does all pairwise comparisons of levels of A alone, A by B, A by CD, A by BCD, etc
# then combines all into one emmGrid object with the specified p_adjustment

compare_all_to_a <- function(model, factors, compare_level, p_adjustment = "fdr"){
  
  n_comps <- (2^(length(factors) - 1))
  
  other_factors <- factors[factors != compare_level]

  # get a list of all possible combinations with this neat lapply/utils::combn usage
  list_of_combos <- lapply(seq_along(other_factors), combn, x = other_factors, simplify = FALSE) |>
    unlist(recursive = FALSE) 
  comps <- vector('list', length = n_comps)
  for (i in seq_along(comps)){
    if (i == 1){
      comps[i] <- pairs(emmeans::emmeans(model, specs = compare_level)) # for the first one, just compare A alone
    }
    else {
      comps[i] <- pairs(emmeans::emmeans(model, specs = compare_level, by = list_of_combos[i - 1][[1]]))
    }
  }
  # make a call of each of the strings
  # have to add p value adjustment to list of arguments for do.call
  comps <- c(comps, adjust = p_adjustment)
  # do call is a key player here-- lets us call the function emmmeans::rbind with all the arguments as strings
  do.call(rbind, comps)
}


# THEN do that for A, B, C, and D in turn

compare_all_pairs <- function(model, factors, p_adjustment = "fdr"){
  comps <- lapply(factors, compare_all_to_a, model = model, factors = factors, p_adjustment = p_adjustment)
  comps <- c(comps, adjust = p_adjustment)
  do.call(rbind, comps)
}

do_single_comparison <- function(model, comparison){
  if (stringr::str_detect(comparison, pattern = "\\|")){
    # do the comparison by the other variables
    input_vars <- stringr::str_split(comparison, pattern = '\\ \\|\\ ') |>
      as_vector() # here, input_vars[[1]] is the first variable that is comp level
    # get other variables from the A | B | C entry
    by_vars <- input_vars[input_vars != input_vars[[1]]] 
    
    pairs(emmeans::emmeans(model, specs = input_vars[[1]], by = by_vars))
  }
  else {
    # just do the comparison alone
    pairs(emmeans::emmeans(model, specs = comparison))
  }
}

compare_pairs <- function(model, comparisons, p_adjustment = "fdr"){
  # some function that takes in string of comparisons, does the comparisons
  
  comps <- lapply(comparisons, do_single_comparison, model = model)
  comps <- c(comps, adjust = p_adjustment)
  do.call(rbind, comps)
}

display_comparison_table <- function(pairwise_output, title = "summary table"){
  pairwise_output |>
    data.frame(stringsAsFactors = FALSE) |>
    mutate(p.value = round(p.value, 5)) |>
    mutate_at(c('estimate', 'SE', 'df', 't.ratio'), function(x){round(x, 2)}) |>
    rename(`p value` = p.value, `t-ratio` = t.ratio, `d.f.` = df) |>
    gt() |>
    cols_align('center') |>
    tab_source_note(
      source_note = "P value adjustment: FDR"
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        rows = `p value` < 0.05)
    ) |>
    tab_header(
      title = md(title)
    ) 
}

display_comparison_table_zvalue <- function(pairwise_output, title = "summary table"){
  pairwise_output |>
    data.frame(stringsAsFactors = FALSE) |>
    mutate(p.value = round(p.value, 5)) |>
    mutate_at(c('estimate', 'SE', 'df', 'z.ratio'), function(x){round(x, 2)}) |>
    rename(`p value` = p.value, `z-ratio` = z.ratio, `d.f.` = df) |>
    gt() |>
    cols_align('center') |>
    tab_source_note(
      source_note = "P value adjustment: FDR"
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        rows = `p value` < 0.05)
    ) |>
    tab_header(
      title = md(title)
    ) 
}

get_pvals_for_comparisons <- function(pairwise_output){
  pairwise_output |>
    as_tibble() |>
    mutate(group = str_split(contrast, pattern = " - ")) |>
    unnest_wider(group, names_sep = "") |>
    mutate(p.adj.signif = case_when(
      p.value < 0.0001 ~ '****',
      p.value < 0.001 ~ '***',
      p.value < 0.01 ~ '**',
      p.value < 0.05 ~ '*',
      p.value > 0.05 ~ 'ns'
    ))
}


################################################################################
## testing 


# fax <- c('tissue', 'condition', 'target_cell')


# compare_all_pairs(non_infected_model, fax)

#do_single_comparison(non_infected_model, "tissue")
#pairs(lsmeans(non_infected_model, ~ condition | tissue | target_cell))
#do_single_comparison(non_infected_model, "condition | tissue | target_cell")


#a <- pairs(lsmeans(non_infected_model, ~ tissue))
#b <- pairs(lsmeans(non_infected_model, ~ condition))
#c <- pairs(lsmeans(non_infected_model, ~ target_cell))
#b.a.c <- pairs(lsmeans(non_infected_model, ~ condition | tissue | target_cell))
#b.c.a <- pairs(lsmeans(non_infected_model, ~  target_cell | condition | tissue))
#rbind(a, b, c, b.c.a, b.a.c, adjust = "fdr")

#comparisons_to_make <- c("tissue", 
                         #"condition",
                         #"target_cell",
                         #"condition | tissue | target_cell",
                         #"target_cell | condition | tissue")

#compare_pairs(non_infected_model, comparisons = comparisons_to_make)

# it works! yay