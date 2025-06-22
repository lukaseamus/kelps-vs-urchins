# 1. Load data ####
require(tidyverse)
require(here)
grazing_prior_posterior <- 
  read_rds(here("Urchins", "RDS", "grazing_prior_posterior.rds"))
C_diff <- 
  read_rds(here("Biochemistry", "C_N", "RDS", "C_diff.rds"))
distance_diff <- 
  read_rds(here("Sinking", "RDS", "distance_diff.rds"))

# 2. Calculation ####
require(magrittr)
CSP <- grazing_prior_posterior %>%
  filter(Season == "Annual") %>%
  mutate(Defecation_prop = beta / 100) %>% # Convert to proportion
  select(starts_with("."), Defecation_prop) %>%
  full_join(
    C_diff %>%
      filter(parameter == "obs_new") %>%
      select(starts_with("."), prop) %>%
      rename(Carbon_prop = prop),
    by = c(".chain", ".iteration", ".draw")
  ) %>%
  full_join(
    distance_diff %>%
      filter(parameter == "mu") %>%
      select(starts_with("."), prop) %>%
      rename(Distance_prop = prop),
    by = c(".chain", ".iteration", ".draw")
  ) %>%
  mutate(CSP = Defecation_prop * Carbon_prop * Distance_prop) %T>%
  print()

CSP %>%
  summarise(mean = mean(CSP),
            sd = sd(CSP),
            n = n())

# 3. Visualisation ####