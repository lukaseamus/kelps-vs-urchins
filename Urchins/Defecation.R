# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

meta <- here("Urchins", "Defecation.csv") %>%
  read_csv() %>%
  mutate(
    Kelp = Plant %>% # Add kelp binary
      str_detect("Alaria|Laminaria|Saccharina|Macrocystis|Nereocystis|Costaria|Agaru|Neoagarum|Ecklonia") &
      !str_detect(coalesce(Notes, ""), "reconstituted"),
    Genus = Urchin %>% str_extract("^\\S+"), # Extract urchin genus
    Function = if_else(
      Genus %>% # Distinguish kelp-feeders as functional group
        str_detect("Strongylocentrotus|Mesocentrotus|Paracentrotus|Psammechinus|Parechinus"),
      "Kelpivore", "Other"
    )
  ) %T>%
  print(n = 105)

# 1.2 Simulate observations ####
# Many of the reported percentages have large standard deviations.
# The beta distribution cannot handle this and degenerates. A truncated
# normal is not the correct data distribution but can handle these cases.
# Rather than use the beta distribution where it works and the truncated
# normal elsewhere, I am consistently using the truncated normal. My
# downstream model is beta distributed, so I believe this is legitimate.
require(extraDistr) # R has no native truncated normal
set.seed(100)
meta %<>%
  rowwise() %>%
  mutate(
    Proportion = if( !is.na(Percentage_mean) & is.na(Percentage_sd) ) {
      list( Percentage_mean / 100 )
    } else if( !is.na(Percentage_mean) & !is.na(Percentage_sd) & 
               !is.na(Percentage_n) ) {
      list(
        rtnorm( 
          Percentage_n , Percentage_mean / 100 ,
          Percentage_sd / 100 , 0 , 1
        )
      )
    } else if( is.na(Percentage_mean) & !is.na(Defecation_mean) &
               is.na(Defecation_sd) ) {
      list( Defecation_mean / Consumption_mean )
    } else if( is.na(Percentage_mean) & !is.na(Absorption_mean) &
               is.na(Absorption_sd) ) {
      list( 1 - Absorption_mean / Consumption_mean )
    } else if( is.na(Percentage_mean) & !is.na(Defecation_mean) &
               !is.na(Defecation_sd) & !is.na(Defecation_n) ) {
      D <- rgamma(
        1e6 , Defecation_mean^2 / Defecation_sd^2 , 
        Defecation_mean / Defecation_sd^2
      )
      C <- rgamma(
        1e6 , Consumption_mean^2 / Consumption_sd^2 , 
        Consumption_mean / Consumption_sd^2
      )
      n <- min(Defecation_n, Consumption_n)
      list( rtnorm( n , mean(D/C) , sd(D/C) , 0 , 1 ) )
    } else if( is.na(Percentage_mean) & !is.na(Absorption_mean) &
               !is.na(Absorption_sd) & !is.na(Absorption_n) ) {
      A <- rgamma(
        1e6 , Absorption_mean^2 / Absorption_sd^2 , 
        Absorption_mean / Absorption_sd^2
      ) 
      C <- rgamma(
        1e6 , Consumption_mean^2 / Consumption_sd^2, 
        Consumption_mean / Consumption_sd^2
      )
      n <- min(Absorption_n, Consumption_n)
      list( 1 - rtnorm( n , mean(A/C) , sd(A/C) , 0 , 1 ) )
    } else {
      list(NULL)
    }
  ) %>%
  unnest(Proportion) %T>%
  print()

# 1.3 Add data from this study ####
experiment <- here("Urchins", "RDS", "grazing.rds") %>% 
  read_rds() %>%
  mutate(Proportion = Defecation / Consumption,
         Reference = "This study",
         Urchin = "Strongylocentrotus droebachiensis",
         Genus = Urchin %>% str_extract("^\\S+"),
         Plant = "Laminaria hyperborea",
         Kelp = TRUE, 
         Function = "Kelpivore") %>%
  select(Reference, Genus, Urchin, Plant, Kelp, Function, Proportion) %T>%
  print(n = 50)

defecation <- meta %>%
  select(Reference, Genus, Urchin, Plant, Kelp, Function, Proportion) %>%
  bind_rows(experiment) %>%
  mutate(Reference = Reference %>% fct(),
         Function = Function %>% fct(),
         Genus = Genus %>% fct(),
         Urchin = Urchin %>% fct()) %>%
  rename(Species = Urchin) %T>%
  print()

defecation %>% nrow() # 1368 observations
defecation %>% filter(Reference != "This study") %>% nrow() # 1324 observations
defecation %$% nlevels(Reference) # 32 references
defecation %$% nlevels(Genus) # 10 genera
defecation %$% nlevels(Species) # 17 species

defecation %$% range(Proportion)
# There are zeros and ones which is not allowed.
# Replace zeros with 1e-4 and ones with 1-1e-4:
defecation %<>%
  mutate(
    Proportion = case_when(
      Proportion == 0 ~ 1e-4,
      Proportion == 1 ~ 1 - 1e-4,
      TRUE ~ Proportion
    )
  ) %T>%
  print()

defecation %$% range(Proportion) # fixed

# 1.4 Data exploration ####
require(ggridges)
defecation %>%
  ggplot() +
  stat_density_ridges(aes(x = Proportion, fill = Kelp, y = Species), 
                      colour = NA, n = 2^10, rel_min_height = 0.001, bandwidth = 0.02, 
                      scale = 3, alpha = 0.6, from = 0, to = 1) +
  theme_minimal()

defecation %>%
  ggplot() +
  stat_density_ridges(aes(x = Proportion, fill = Kelp, y = Genus), 
                      colour = NA, n = 2^10, rel_min_height = 0.001, bandwidth = 0.02, 
                      scale = 3, alpha = 0.6, from = 0, to = 1) +
  theme_minimal()

defecation %>%
  ggplot() +
  stat_density_ridges(aes(x = Proportion, fill = Kelp, y = Function), 
                      colour = NA, n = 2^10, rel_min_height = 0.001, bandwidth = 0.02, 
                      scale = 3, alpha = 0.6, from = 0, to = 1) +
  theme_minimal()
 
defecation %>%
  count(Genus, Species, Kelp, Function) %>%
  print(n = 30)
# Including food source inn the model probably makes no sense, but I do want
# to distinguish functional groups in the prediction. Urchins that don't feed
# on kelp (mostly tropical species) can be summarised at their functional
# group level which kelp-feeding urchins will be reported at the species level.

defecation %>%
  group_by(Species) %>%
  summarise(References = n_distinct(Reference),
            Observations = n())
# 12 references and 612 observations on S. droebachiensis alone, 
# 4 references each on the other Strongylocentrotus species.
# Strongylocentrotus is definitely the most-studied genus.

# 2. Model ####
# 2.1 Prior simulation ####
tibble(n = 1:1e4,
       alpha = rnorm( 1e4 , qlogis(0.5) , 0.8 ), # global mean
       sigma_f = rtnorm( 1e4 , 0 , 0.5 , 0 ), # functional group variation
       sigma_g = rtnorm( 1e4 , 0 , 0.5 , 0 ), # genus variation
       sigma_s = rtnorm( 1e4 , 0 , 0.5 , 0 ), # species variation
       sigma_r = rtnorm( 1e4 , 0 , 0.5 , 0 ), # reference variation
       nu = rgamma( 1e4 , 30^2 / 20^2 , 30 / 20^2 )) %>% # likelihood precision
  mutate(
    f = rnorm( n() , alpha , sigma_f ),
    g = rnorm( n() , 0 , sigma_g ),
    s = rnorm( n() , 0 , sigma_s ),
    r = rnorm( n() , 0 , sigma_r ),
    mu = plogis( f + g + s + r ),
    p = rbeta( n() , mu * nu , (1 - mu) * nu )
  ) %>%
  pivot_longer(cols = c(mu, p),
               names_to = "parameter") %>%
  ggplot() +
    geom_density(aes(value), fill = "black") +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Essentially uniform

# 2.2 Stan model ####
require(cmdstanr)
meta_c_model <- here("Urchins", "Stan", "meta_c.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

meta_nc_model <- here("Urchins", "Stan", "meta_nc.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
meta_c_samples <- meta_c_model$sample(
          data = defecation %>%
            select(Proportion, Function, Genus, 
                   Species, Reference) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()

meta_nc_samples <- meta_nc_model$sample(
          data = defecation %>%
            select(Proportion, Function, Genus, 
                   Species, Reference) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()

# 2.3 Model diagnostics ####
# 2.3.1 R-hat and ESS ####
meta_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 35% of rhat above 1.001. rhat = 1.00 ± 0.000610.

meta_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 67% of rhat above 1.001. rhat = 1.00 ± 0.00256.

# 2.3.2 Chains ####


# 2.3.3 Pairs ####


# 2.4 Prior-posterior comparison ####

# 2.5 Prediction ####
# 2.5.1 Global parameters ####

# 2.5.2 Functional group parameters ####

# 2.5.3 Genus parameters ####

# 2.5.4 Species parameters ####


