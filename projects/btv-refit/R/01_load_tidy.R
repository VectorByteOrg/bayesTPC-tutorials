# 01_load_tidy.R
# Load and tidy BTV data from El Moustaid et al. (2021) Appendix A.6

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# Read raw data
raw <- read_excel("data/BTV_Joe_ph_ph.xlsx")

# Clean and standardize columns
dat <- raw %>%
  transmute(
    species = str_to_lower(interactor1species),
    T = as.numeric(interactor1temp),
    trait_name = originaltraitname,
    value = as.numeric(originaltraitvalue),
    replicates = as.numeric(replicates)  # for potential binomial likelihood
  ) %>% 
  filter(!is.na(T), !is.na(value))

# Map trait names to canonical codes
map_trait <- function(x) {
  x <- str_to_lower(x)
  case_when(
    str_detect(x, "biting") ~ "a",
    str_detect(x, "fecund") ~ "F", 
    str_detect(x, "vector competence") ~ "b",
    str_detect(x, "extrinsic incubation|eip") ~ "EIP",
    str_detect(x, "mortality") ~ "mu",
    str_detect(x, "development time|survival time") ~ "rho_time",
    str_detect(x, "survival rate|survival prob") ~ "p",
    TRUE ~ "other"
  )
}

# Distinguish survival stages (pE/pL/pP)
stage_from_name <- function(x) {
  x <- tolower(x)
  case_when(
    grepl("egg", x) ~ "pE",
    grepl("larv", x) ~ "pL", 
    grepl("pup", x) ~ "pP",
    TRUE ~ "p"  # fallback
  )
}

# Apply trait mapping
dat <- dat %>%
  mutate(
    code = map_trait(trait_name),
    p_stage = ifelse(code == "p", stage_from_name(trait_name), NA_character_),
    trait = coalesce(p_stage, case_match(code,
      "a" ~ "a",
      "F" ~ "F", 
      "b" ~ "b",
      "EIP" ~ "nu",
      "rho_time" ~ "rho",
      "mu" ~ "mu",
      "other" ~ "other"
    ))
  )

# Source data sanity functions
source("R/data_sanity.R")

# - ν = 1/EIP (1/day) - handle zero EIP values (no development)
nu_tbl <- dat %>% 
  filter(code == "EIP") %>%
  mutate(y = ifelse(value == 0, 0, 1/value)) %>%  # 0 EIP -> 0 development rate
  transmute(species, trait = "nu", T, y)

# - development times (ρE, ρL, ρP) stay as days (Quadratic on time)
rho_tbl <- dat %>% 
  filter(code == "rho_time") %>% 
  transmute(species, trait = "rho", T, y = value)

# - survival: split p into pE/pL/pP if stage info available
p_tbl <- dat %>% 
  filter(trait %in% c("p", "pE", "pL", "pP")) %>% 
  transmute(species, trait, T, y = value, replicates)

# - remaining: a, F, b, mu
rest_tbl <- dat %>%
  filter(trait %in% c("a", "F", "b", "mu")) %>%
  transmute(species, trait, T, y = value, replicates)

# Combine all tables
btv <- bind_rows(nu_tbl, rho_tbl, p_tbl, rest_tbl)

# Apply transformations to the combined data
btv <- transform_btv(btv)

# Validate the data
validate_rowwise(btv)

# Create output directory and save
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
write_csv(btv, "outputs/tidy_btv.csv")

# Print summary
cat("Data tidying complete!\n")
cat("Species found:", unique(btv$species), "\n")
cat("Traits found:", unique(btv$trait), "\n")
cat("Total observations:", nrow(btv), "\n")

# Show trait counts
cat("\nTrait counts:\n")
btv %>% count(trait) %>% print()
