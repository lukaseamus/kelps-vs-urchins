# 1. Load data ####
require(tidyverse)
require(here)
require(magrittr)
# There is only enough data on Strongylocentrotus droebachiensis to calculate
# carbon sequestration potential (CSP). I am expressing CSP of faeces relative
# to the kelp baseline.

# Defecated fraction
defecation <- here("Urchins", "RDS", "meta_prior_posterior_s.rds") %>%
  read_rds() %>%
  filter(Species == "Strongylocentrotus droebachiensis") %>%
  droplevels() %T>%
  print()

# Carbon fraction
carbon <- here("Biochemistry", "RDS", "biochem_prior_posterior_s.rds") %>%
  read_rds() %>%
  filter(Compound == "Carbon" & Species == "Strongylocentrotus droebachiensis") %>%
  droplevels() %T>%
  print()

# Depth
depth <- here("Sinking", "RDS", "depth_diff.rds") %>%
  read_rds() %>%
  filter(parameter == "obs") %T>%
  print()

# 2. Calculation ####
# Combine data
CSP <- defecation %>%
  select(starts_with("."), Proportion) %>%
  rename(Mass = Proportion) %>%
  full_join(
    carbon %>%
      select(starts_with("."), Ratio) %>%
      rename(Carbon = Ratio)
  ) %>%
  full_join(
    depth %>%
      mutate(Ratio = 10^logratio) %>%
      select(starts_with("."), Ratio) %>%
      rename(Export = Ratio)
  ) %T>%
  print()

# Calculate carbon transfer through urchin
CSP %<>%
  mutate(Transfer = Mass * Carbon) %T>%
  print()
  
# Calculate CSP based on carbon transfer and export
CSP %<>%
  mutate(CSP = Transfer * Export) %T>%
  print()

# 3. Summary ####
require(glue)
CSP_summary <- CSP %>%
  pivot_longer(
    cols = !starts_with("."),
    names_to = "Measure",
    values_to = "Ratio" # All measures represent faeces / kelp
  ) %>%
  mutate(log_Ratio = log(Ratio)) %>%
  summarise(
    across(
      c(Ratio, log_Ratio),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Ratio < 1 ), # equivalent to log_Ratio < 0
    n = n(),
    .by = Measure
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~if_else(.x < 100, signif(.x, 2), signif(.x, 3))
    ),
    Ratio = glue( "{Ratio_median} ({log_Ratio_mean} ± {log_Ratio_sd})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

CSP_summary %>%
  write_csv(here("CSP", "CSP.csv"))

require(officer)
read_docx() %>%
  body_add_table(
    value = CSP_summary
  ) %>%
  print(target = here("CSP", "CSP.docx"))

# 4. Visualisation ####
CSP %>%
  ggplot() +
    geom_vline(xintercept = 1) +
    geom_density(aes(Transfer)) +
    scale_x_log10() +
    theme_minimal()

CSP %>%
  ggplot() +
    geom_vline(xintercept = 1) +
    geom_density(aes(Export)) +
    scale_x_log10() +
    theme_minimal()

CSP %>%
  ggplot() +
    geom_vline(xintercept = 1) +
    geom_density(aes(CSP)) +
    scale_x_log10() +
    theme_minimal()
# The putative benefit of sinking speed is mostly cancelled by minimal carbon transfer.
# If other factors such as microbial decomposition are also considered, it is likely
# that faecal carbon is completely mineralised before it can be exported.
