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
# Including food source in the model probably makes no sense, but I do want
# to distinguish functional groups in the prediction. Urchins that don't feed
# on kelp (mostly tropical species) can be summarised at their functional
# group level which kelp-feeding urchins will be reported at the species level.

defecation %>%
  summarise(References = n_distinct(Reference),
            Observations = n(),
            .by = Species) %>%
  arrange(desc(Observations))
# 12 references and 612 observations on S. droebachiensis alone, 
# 4 references each on the other Strongylocentrotus species.
# Apart from Paracentrotus lividus and Psammechinus miliaris,
# other urchins only have one study. Strongylocentrotus is 
# definitely the most-studied genus.

# 2. Model ####
# 2.1 Prior simulation ####
tibble(n = 1:1e4,
       alpha = rnorm( 1e4 , qlogis(0.5) , 0.5 ), # global mean
       sigma_f = rtnorm( 1e4 , 0 , 0.5 , 0 ), # functional group variation
       sigma_g = rtnorm( 1e4 , 0 , 0.5 , 0 ), # genus variation
       sigma_s = rtnorm( 1e4 , 0 , 0.5 , 0 ), # species variation
       sigma_r = rtnorm( 1e4 , 0 , 0.5 , 0 ), # reference variation
       nu = rgamma( 1e4 , 20^2 / 10^2 , 20 / 10^2 )) %>% # likelihood precision
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
# Close to uniform

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

# Save draws
meta_c_samples$draws() %>%
  write_rds(here("Urchins", "RDS", "meta_c_samples.rds"))
meta_c_samples$draws(format = "df") %>%
  write_rds(here("Urchins", "RDS", "meta_c_samples_df.rds"))

meta_nc_samples$draws() %>%
  write_rds(here("Urchins", "RDS", "meta_nc_samples.rds"))
meta_nc_samples$draws(format = "df") %>%
  write_rds(here("Urchins", "RDS", "meta_nc_samples_df.rds"))

# 2.3 Model diagnostics ####
# 2.3.1 R-hat and ESS ####
meta_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 100% of rhat above 1.001. rhat = 1.02 ± 0.0160.

meta_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000133.
# Already the non-centred model is more promising.

# 2.3.2 Chains ####
require(bayesplot)
meta_c_chains <- meta_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off")

meta_nc_chains <- meta_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Non-centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off")

require(patchwork)
( ( meta_c_chains / meta_nc_chains ) +
  plot_layout(guides = "collect",
              heights = c(8/11, 1)) &
    theme(legend.position = "top",
          legend.justification = 0) ) %>%
  ggsave(filename = "meta_chains.pdf", path = here("Urchins", "Plots"),
         device = cairo_pdf, width = 80, height = 80, units = "cm")
# Non-centred chains are clearly better.

# 2.3.3 Pairs ####
meta_c_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("alpha", "f[1]", "f[2]", "g[1]", "g[2]", 
             "s[1]", "s[2]", "r[1]", "r[2]",
             "sigma_f", "sigma_g", "sigma_s", "sigma_r", 
             "nu"),
    grid_args = list(top = "Centred model")
  ) %>%
  ggsave(filename = "meta_c_pairs.png", path = here("Urchins", "Plots"),
         width = 50, height = 50, units = "cm", bg = "white")

meta_nc_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("alpha", "f[1]", "f[2]", "g[1]", "g[2]", 
             "s[1]", "s[2]", "r[1]", "r[2]",
             "sigma_f", "sigma_g", "sigma_s", "sigma_r", 
             "nu"),
    grid_args = list(top = "Non-centred model")
  ) %>%
  ggsave(filename = "meta_nc_pairs.png", path = here("Urchins", "Plots"),
         width = 50, height = 50, units = "cm", bg = "white")
# No posterior correlation in either case, but non-centred posteriors
# look more stable (less spiky).

# 2.4 Prior-posterior comparison ####
# 2.4.1 Sample prior ####
source("functions.R")
# prior_samples only works properly with non-centred models, but
# this is fine since priors are identical for both models
meta_prior <- prior_samples(
  model = meta_nc_model, 
  data = defecation %>%
    select(Proportion, Function, Genus, 
           Species, Reference) %>%
    compose_data()
)

# 2.4.2 Centred model ####
meta_c_prior_posterior_global <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_c_samples,
    parameters = c("alpha", "sigma_f", "sigma_g", 
                   "sigma_s", "sigma_r", "nu"),
    format = "long"
    ) %>%
  prior_posterior_plot() +
  labs(title = "Centred model")

meta_c_prior_posterior_f <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_c_samples,
    group = defecation %>% select(Function),
    parameters = c("f[Function]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Function", ridges = TRUE)

meta_c_prior_posterior_g <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_c_samples,
    group = defecation %>% select(Genus),
    parameters = c("g[Genus]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Genus", ridges = TRUE)

meta_c_prior_posterior_s <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_c_samples,
    group = defecation %>% select(Species),
    parameters = c("s[Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species", ridges = TRUE)

meta_c_prior_posterior_r <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_c_samples,
    group = defecation %>% select(Reference),
    parameters = c("r[Reference]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Reference", ridges = TRUE)

# 2.4.3 Non-centred model ####
meta_nc_prior_posterior_global <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    parameters = c("alpha", "sigma_f", "sigma_g", 
                   "sigma_s", "sigma_r", "nu"),
    format = "long"
    ) %>%
  prior_posterior_plot() +
  labs(title = "Non-centred model")

meta_nc_prior_posterior_f <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    group = defecation %>% select(Function),
    parameters = c("f[Function]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Function", ridges = TRUE)

meta_nc_prior_posterior_g <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    group = defecation %>% select(Genus),
    parameters = c("g[Genus]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Genus", ridges = TRUE)

meta_nc_prior_posterior_s <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    group = defecation %>% select(Species),
    parameters = c("s[Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species", ridges = TRUE)

meta_nc_prior_posterior_r <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    group = defecation %>% select(Reference),
    parameters = c("r[Reference]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Reference", ridges = TRUE)

# 2.4.4 Combined ####
meta_prior_posterior <- 
  ( ( meta_c_prior_posterior_global /
      meta_c_prior_posterior_f / 
      meta_c_prior_posterior_g / 
      meta_c_prior_posterior_s / 
      meta_c_prior_posterior_r ) +
      plot_layout(heights = c(4/32, 2/32, 10/32, 17/32, 1)) |
    ( meta_nc_prior_posterior_global /
        meta_nc_prior_posterior_f / 
        meta_nc_prior_posterior_g / 
        meta_nc_prior_posterior_s / 
        meta_nc_prior_posterior_r ) +
      plot_layout(heights = c(4/32, 2/32, 10/32, 17/32, 1)) ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", 
        legend.justification = 0)

meta_prior_posterior %>%
  ggsave(filename = "meta_prior_posterior.pdf", path = here("Urchins", "Plots"),
         device = cairo_pdf, width = 60, height = 60, units = "cm")
# Non-centred model posteriors are a lot smoother.
# Proceed with non-centred model.

# 2.5 Prediction ####
# 2.5.1 Global parameters ####
meta_prior_posterior_global <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    parameters = c("alpha", "sigma_f", "sigma_g",
                   "sigma_s", "sigma_r", "nu"),
    format = "short"
  ) %>%
  mutate( 
    # Calculate proportion for new functional groups, genera, species and studies
    mu = plogis(
      rnorm( n() , alpha , sigma_f ) +
        rnorm( n() , 0 , sigma_g ) +
        rnorm( n() , 0 , sigma_s ) +
        rnorm( n() , 0 , sigma_r )
    ),
    Proportion = rbeta( n() , mu * nu , (1 - mu) * nu )
  ) %T>%
  print()

meta_prior_posterior_global %>% # Save
  write_rds(here("Urchins", "RDS", "meta_prior_posterior_global.rds"))

# 2.5.2 Functional group parameters ####
meta_prior_posterior_f <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    group = defecation %>% select(Function),
    parameters = c("f[Function]", "sigma_g", "sigma_s",
                   "sigma_r", "nu"),
    format = "short"
  ) %>%
  mutate( 
    # Calculate proportion for existing functional groups, 
    # but new genera, species and studies
    mu = plogis(
      f + rnorm( n() , 0 , sigma_g ) +
        rnorm( n() , 0 , sigma_s ) +
        rnorm( n() , 0 , sigma_r )
    ),
    Proportion = rbeta( n() , mu * nu , (1 - mu) * nu )
  ) %>%
  filter(Function == "Kelpivore" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in Function
    Function = if_else(
      distribution == "prior", "Prior", Function
    ) %>% fct()
  ) %>%
  select(starts_with("."), Function, mu, Proportion) %T>%
  print()

meta_prior_posterior_f %>% # Save
  write_rds(here("Urchins", "RDS", "meta_prior_posterior_f.rds"))

# 2.5.3 Genus parameters ####
meta_prior_posterior_g <- defecation %>%
  distinct(Function, Genus) %>%
  # Genera are not matched to existing functional groups so this needs 
  # to be done manually
  left_join(
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = defecation %>% select(Function),
        parameters = c("f[Function]"),
        format = "short"
      ),
    relationship = "many-to-many" # Function is repeated across genera
  ) %>%
  left_join(
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = defecation %>% select(Genus),
        parameters = c("g[Genus]", "sigma_s", 
                       "sigma_r", "nu"),
        format = "short"
      )
  ) %>%
  mutate( 
    # Calculate proportion for existing genera within functional groups,
    # but new species and studies
    mu = plogis(
      f + g + rnorm( n() , 0 , sigma_s ) +
        rnorm( n() , 0 , sigma_r )
    ),
    Proportion = rbeta( n() , mu * nu , (1 - mu) * nu )
  ) %>%
  filter(Function == "Kelpivore" & Genus == "Strongylocentrotus" & 
           distribution == "prior" | distribution == "posterior") %>%
  mutate(
    Function = if_else(
      distribution == "prior", "Prior", Function
    ) %>% fct(),
    Genus = if_else(
      distribution == "prior", "Prior", Genus
    ) %>% fct()
  ) %>%
  select(starts_with("."), Function, Genus, mu, Proportion) %T>%
  print()

meta_prior_posterior_g %>% # Save
  write_rds(here("Urchins", "RDS", "meta_prior_posterior_g.rds"))

# 2.5.4 Species parameters ####
meta_prior_posterior_s <- defecation %>%
  distinct(Function, Genus, Species) %>%
  left_join(
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = defecation %>% select(Function),
        parameters = c("f[Function]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = defecation %>% select(Genus),
        parameters = c("g[Genus]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = defecation %>% select(Species),
        parameters = c("s[Species]", "sigma_r", "nu"),
        format = "short"
      )
  ) %>%
  mutate( 
    # Calculate proportion for new studies on existing species, 
    # within genera within functional groups
    mu = plogis(
      f + g + s + rnorm( n() , 0 , sigma_r )
    ),
    Proportion = rbeta( n() , mu * nu , (1 - mu) * nu )
  ) %>%
  filter(Function == "Kelpivore" & Genus == "Strongylocentrotus" & 
           Species == "Strongylocentrotus droebachiensis" &
           distribution == "prior" | distribution == "posterior") %>%
  mutate(
    Function = if_else(
      distribution == "prior", "Prior", Function
    ) %>% fct(),
    Genus = if_else(
      distribution == "prior", "Prior", Genus
    ) %>% fct(),
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(starts_with("."), Function, Genus, Species, 
         mu, Proportion) %T>%
  print()

meta_prior_posterior_s %>% # Save
  write_rds(here("Urchins", "RDS", "meta_prior_posterior_s.rds"))

# 2.5.5 Reference parameters ####
meta_prior_posterior_r <- defecation %>%
  distinct(Function, Genus, Species, Reference) %>%
  left_join(
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = defecation %>% select(Function),
        parameters = c("f[Function]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = defecation %>% select(Genus),
        parameters = c("g[Genus]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = defecation %>% select(Species),
        parameters = c("s[Species]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = defecation %>% select(Reference),
        parameters = c("r[Reference]", "nu"),
        format = "short"
      )
  ) %>%
  mutate( 
    # Calculate proportion for existing references within species
    mu = plogis( f + g + s + r ),
    Proportion = rbeta( n() , mu * nu , (1 - mu) * nu )
  ) %>%
  filter(Function == "Kelpivore" & Genus == "Strongylocentrotus" & 
           Species == "Strongylocentrotus droebachiensis" &
           Reference == "Mamelona & Pelletier 2005" &
           distribution == "prior" | distribution == "posterior") %>%
  mutate(
    Function = if_else(
      distribution == "prior", "Prior", Function
    ) %>% fct(),
    Genus = if_else(
      distribution == "prior", "Prior", Genus
    ) %>% fct(),
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct(),
    Reference = if_else(
      distribution == "prior", "Prior", Reference
    ) %>% fct()
  ) %>%
  select(starts_with("."), Function, Genus, Species, Reference,
         mu, Proportion) %T>%
  print()

meta_prior_posterior_r %>% # Save
  write_rds(here("Urchins", "RDS", "meta_prior_posterior_r.rds"))

# Clean up
rm(
  list = c(
    ls(pattern = "meta_c|meta_nc"),
    "meta_prior", "meta_prior_posterior"
  )
)
gc()

# 2.6 Parameter estimates ####
# 2.6.1 Global parameters ####
require(glue)
meta_summary_global <- meta_prior_posterior_global %>%
  mutate(Percentage = Proportion * 100) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Proportion < 0.5 ),
    n = n(),
    .by = distribution
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_mean} ± {mu_sd} ({mu_median})" ),
    Percentage = glue( "{Percentage_mean} ± {Percentage_sd} ({Percentage_median})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.2 Functional group parameters ####
meta_summary_f <- meta_prior_posterior_f %>%
  mutate(Percentage = Proportion * 100) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Proportion < 0.5 ),
    n = n(),
    .by = Function
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_mean} ± {mu_sd} ({mu_median})" ),
    Percentage = glue( "{Percentage_mean} ± {Percentage_sd} ({Percentage_median})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.3 Genus parameters ####
meta_summary_g <- meta_prior_posterior_g %>%
  mutate(Percentage = Proportion * 100) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Proportion < 0.5 ),
    n = n(),
    .by = c(Function, Genus)
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_mean} ± {mu_sd} ({mu_median})" ),
    Percentage = glue( "{Percentage_mean} ± {Percentage_sd} ({Percentage_median})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.4 Species parameters ####
meta_summary_s <- meta_prior_posterior_s %>%
  mutate(Percentage = Proportion * 100) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Proportion < 0.5 ),
    n = n(),
    .by = c(Function, Genus, Species)
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_mean} ± {mu_sd} ({mu_median})" ),
    Percentage = glue( "{Percentage_mean} ± {Percentage_sd} ({Percentage_median})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.5 Reference parameters ####
meta_summary_r <- meta_prior_posterior_r %>%
  mutate(Percentage = Proportion * 100) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Proportion < 0.5 ),
    n = n(),
    .by = c(Function, Genus, Species, Reference)
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_mean} ± {mu_sd} ({mu_median})" ),
    Percentage = glue( "{Percentage_mean} ± {Percentage_sd} ({Percentage_median})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print(n = 40)

# 2.6.6 Experiment parameters ####
experiment_summary <- here("Urchins", "RDS", "grazing_prior_posterior.rds") %>%
  read_rds() %>% # beta is already %
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( beta < 50 ), 
    n = n(),
    .by = Season
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    Percentage = glue( "{beta_mean} ± {beta_sd} ({beta_median})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.7 Table 1 ####
# 2.7.1 Combine summaries ####
meta_summary <- bind_rows(
  meta_summary_global %>%
    filter(distribution != "prior") %>%
    select(-distribution),
  meta_summary_f %>%
    filter(Function != "Prior"),
  meta_summary_g %>%
    filter(Genus %in% c("Strongylocentrotus", "Mesocentrotus")),
  meta_summary_s %>%
    filter(Function == "Kelpivore")
) %>%
  full_join(
    bind_rows(
      defecation %>%
        summarise(References = n_distinct(Reference),
                  Observations = n()),
      defecation %>%
        summarise(References = n_distinct(Reference),
                  Observations = n(),
                  .by = Function),
      defecation %>%
        summarise(References = n_distinct(Reference),
                  Observations = n(),
                  .by = c(Function, Genus)) %>%
        filter(Genus %in% c("Strongylocentrotus", "Mesocentrotus")),
      defecation %>%
        summarise(References = n_distinct(Reference),
                  Observations = n(),
                  .by = c(Function, Genus, Species)) %>%
        filter(Function == "Kelpivore")
    )
  ) %>%
  mutate(
    Urchin = case_when(
      is.na(Function) ~ "All urchins (17 species)",
      Function == "Kelpivore" & is.na(Genus) ~
        "Kelp-feeding urchins (8 species)",
      Function == "Other" ~ "Other urchins (9 species)",
      !is.na(Genus) & is.na(Species) ~ Genus,
      !is.na(Species) ~ Species
    )
  ) %>%
  select(Urchin, Percentage, P, References, Observations) %T>%
  print()

# 2.7.2 Add to experimental data ####
Table_1 <- experiment_summary %>%
  filter(Season != "Prior") %>%
  full_join(
    here("Urchins", "RDS", "grazing.rds") %>% 
      read_rds() %>%
      summarise(Observations = n(),
                .by = Season) %>%
      add_row(Season = "Annual" %>% fct(), 
              Observations = sum(.$Observations))
  ) %>%
  rename(Urchin = Season) %>%
  select(Urchin, Percentage, P, Observations) %>%
  bind_rows(meta_summary) %>%
  mutate(
    Urchin = Urchin %>% fct_relevel(
      "Spring", "Summer", "Autumn", "Annual",
      "Strongylocentrotus droebachiensis",
      "Strongylocentrotus purpuratus",
      "Strongylocentrotus intermedius",
      "Strongylocentrotus",
      "Mesocentrotus franciscanus",
      "Mesocentrotus nudus", "Mesocentrotus",
      "Paracentrotus lividus", "Psammechinus miliaris",
      "Parechinus angulosus", "Kelp-feeding urchins (8 species)",
      "Other urchins (9 species)"
    )
  ) %>%
  select(Urchin, Percentage, P, References, Observations) %>%
  arrange(Urchin) %T>%
  print()

# 2.7.3 Save table ####
Table_1 %>%
  write_csv(here("Tables", "Table_1.csv"))

require(officer)
read_docx() %>%
  body_add_table(
    value = Table_1 %>%
      mutate(
        P = P %>% # round down to zero
          replace_values(0.000062 ~ 0)
      )
  ) %>%
  print(target = here("Tables", "Table_1.docx"))

# 2.8 Figure 1 ####
# 2.8.1 Experiment ####
# Load data
grazing <- here("Urchins", "RDS", "grazing.rds") %>%
  read_rds()
grazing_prediction_summary <- here("Urchins", "RDS", "grazing_prediction_summary.rds") %>%
  read_rds()
grazing_prior_posterior <- here("Urchins", "RDS", "grazing_prior_posterior.rds") %>%
  read_rds()
consumption_prior_posterior <- here("Urchins", "RDS", "consumption_prior_posterior.rds") %>%
  read_rds()
defecation_prior_posterior <- here("Urchins", "RDS", "defecation_prior_posterior.rds") %>%
  read_rds()

# Define custom theme
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black", lineend = "square"),
                 legend.key = element_blank(),
                 legend.key.width = unit(.25, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "top",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing = unit(0.6, "cm"),
                 text = element_text(family = "Futura"))

require(geomtextpath)
Fig_1a_left <- ggplot() + 
    geom_textline(data = tibble(x = c(0, 5), y = c(0, 5)), aes(x, y),
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_point(data = grazing %>%
                 mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
               aes(Consumption, Defecation, shape = Season),
               size = 2.5, alpha = 0.4, colour = "#7030a5") +
    geom_line(data = grazing_prediction_summary %>%
                filter(!Season %in% c("Annual", "Prior")),
              aes(Consumption, mu, colour = Season)) +
    geom_ribbon(data = grazing_prediction_summary %>%
                  filter(!Season %in% c("Annual", "Prior")),
                aes(Consumption, ymin = mu.lower, ymax = mu.upper,
                    fill = Season, alpha = factor(.width))) +
    geom_text(data = grazing_prediction_summary %>%
                filter(!Season %in% c("Annual", "Prior")) %>% 
                summarise(x = max(Consumption),
                          y = max(mu.upper),
                          .by = Season),
              aes(x = x, y = y, label = Season),
              family = "Futura", nudge_y = c(0.3, 0.3, -1.5),
              size = 10, size.unit = "pt") +
    geom_density(data = consumption_prior_posterior %>%
                   filter(!Season %in% c("Annual", "Prior")),
                 aes(x = mu, y = after_stat(density) * 0.5, fill = Season),
                 n = 2^10, bw = 10*0.02, alpha = 0.7, colour = NA) +
    geom_density(data = defecation_prior_posterior %>%
                   filter(!Season %in% c("Annual", "Prior")),
                 aes(x = after_stat(density) * 0.35, y = mu, fill = Season),
                 n = 2^10, bw = 5*0.02, alpha = 0.7, colour = NA) +
    scale_shape_manual(values = c(16, 17, 15)) +
    scale_fill_manual(values = rep("#7030a5", 3), guide = "none") +
    scale_colour_manual(values = rep("#7030a5", 3), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(limits = c(0, 10), 
                       breaks = seq(0, 10, 2), 
                       oob = scales::oob_keep) +
    scale_y_continuous(limits = c(0, 5),
                       oob = scales::oob_keep) +
    labs(x = expression("Consumption (mg g"^-1*" d"^-1*")"),
         y = expression("Defecation (mg g"^-1*" d"^-1*")"),
         shape = "Strongylocentrotus droebachiensis") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 1, 0.2, 0.2, unit = "cm"),
          legend.title = element_text(size = 12, face = "italic", hjust = 0,
                                      margin = margin(r = .5, unit = "cm")))

Fig_1a_left

# Add defecation (%) to data.
grazing %<>% 
  mutate(Percentage = Defecation / Consumption * 100) %T>%
  print()

Fig_1a_right <- ggplot() +
  geom_jitter(data = grazing %>%
                mutate(Season = Season %>% fct_relevel("Autumn", "Summer")),
              aes(x = Percentage, y = Season %>% as.numeric() - 0.5,
                  shape = Season),
              colour = "#7030a5", alpha = 0.4, size = 2.5, height = 0.36) +
  stat_density_ridges(data = grazing_prior_posterior %>%
                        filter(!Season %in% c("Prior", "Annual")) %>%
                        mutate(Season = Season %>% 
                                 fct_drop() %>%
                                 fct_relevel("Autumn", "Summer")),
                      aes(x = beta, y = Season %>% as.numeric(), fill = Season), 
                      colour = NA, n = 2^10, rel_min_height = 0.001, bandwidth = 100*0.02, 
                      scale = 1.2, alpha = 0.7, from = 0, to = 100) +
  geom_text(data = grazing %>%
              mutate(Season = Season %>% fct_relevel("Autumn", "Summer"),
                     y = Season %>% as.numeric(),
                     x = if_else(Season == "Spring", 100, 0),
                     hjust = if_else(Season == "Spring", 1, 0)) %>%
              distinct(Season, x, y, hjust),
            aes(x = x, y = y, label = Season, hjust = hjust),
            vjust = -0.2, size = 10, size.unit = "pt", family = "Futura") +
  scale_shape_manual(values = c(15, 17, 16),
                     guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = rep("#7030a5", 3), guide = "none") +
  scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
  xlab("Defecation (%)") +
  coord_cartesian(ylim = c(0, 4.5), expand = FALSE, clip = "off") +
  mytheme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"))

Fig_1a_right

# 2.8.2 Meta-analysis ####
# I will plot predictions for kelp-feeding urchins only
# but will stratify by urchin family for clarity:
# Strongylocentrotidae contains Strongylocentrotus and Mesocentrotus,
# Parechinidae contains Paracentrotus, Psammechinus and Parechinus.

# Select relevant posteriors
meta_prior_posterior <- bind_rows(
  Urchin = meta_prior_posterior_g %>% filter(Genus == "Mesocentrotus"),
  Urchin = meta_prior_posterior_s %>% filter(Function == "Kelpivore" & 
                                      Genus != "Mesocentrotus"),
  Reference = meta_prior_posterior_r %>% filter(Function == "Kelpivore"),
  .id = "Level"
) %>%
  mutate(
    Urchin = if_else(
      Genus == "Mesocentrotus", "Mesocentrotus species", Species
    ) %>% fct(),
    Reference = if_else(
      is.na(Reference), "Global", Reference
    ) %>% fct()
  ) %>%
  select(-c(Function, Species)) %T>%
  print()

# Add Urchin variable to data
defecation %<>%
  mutate(
    Urchin = if_else(
      Genus == "Mesocentrotus", "Mesocentrotus species", Species
    )
  ) %T>%
  print()

# Add family
family <- defecation %>%
  filter(Function == "Kelpivore") %>%
  distinct(Genus) %>%
  mutate(
    Family = case_when(
      Genus %in% c("Strongylocentrotus", "Mesocentrotus") ~
        "Strongylocentrotidae",
      Genus %in% c("Paracentrotus", "Psammechinus", "Parechinus") ~
        "Parechinidae"
    )
  ) %T>%
  print()

defecation %<>%
  full_join(family) %T>%
  print()

meta_prior_posterior %<>%
  full_join(family) %T>%
  print()

# Build densities manually
meta_dens <- meta_prior_posterior %>%
  group_by(Level, Family, Urchin, Reference) %>%
  reframe(
    x = c(
      0, 
      density(
        mu, n = 2^10, 
        from = 0, to = 1, 
        bw = 1 * 0.02
      )$x,
      1
    ),
    y = c(
      0, 
      density(
        mu, n = 2^10, 
        from = 0, to = 1, 
        bw = 1 * 0.02
      )$y, 
      0
    )
  ) %>%
  group_by(Level, Family, Urchin, Reference) %>% 
  # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate(y = y * 0.1 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

# Assign y offset
ys <- defecation %>%
  filter(Function == "Kelpivore") %>%
  distinct(Family, Urchin) %>%
  mutate(
    y = c(
      "Mesocentrotus species" = 1,
      "Parechinus angulosus" = 1,
      "Strongylocentrotus intermedius" = 4.5,
      "Psammechinus miliaris" = 4.5,
      "Strongylocentrotus purpuratus" = 8,
      "Paracentrotus lividus" = 8,
      "Strongylocentrotus droebachiensis" = 11.5
    )[Urchin]
  ) %T>%
  print()

defecation %<>%
  full_join(ys) %T>%
  print()

meta_dens %<>%
  full_join(ys, by = c("Urchin", "Family")) %>%
  mutate(x = x * 100,
         y = y.x + y.y,
         Reference = Reference %>% fct_relevel("Global", after = Inf),
         Group = interaction(Urchin, Reference, sep = "_"),
         Highlight = Reference == "Mamelona & Pelletier 2005") %T>%
  print()

# Plot
Fig_1b <- ggplot() +
  geom_jitter(data = defecation %>% filter(Function == "Kelpivore"),
              aes(x = Proportion * 100, y = y - 0.5),
              colour = "#7030a5", alpha = 0.4, shape = 16, 
              size = 2.5, height = 0.36) +
  geom_line(data = meta_dens,
            aes(x, y, group = Group, colour = Level, alpha = Level)) +
  geom_text(data = ys, aes(2, y + 2, label = Urchin),
            family = "Futura", fontface = "italic", 
            size.unit = "pt", size = 12, hjust = 0) +
  geom_textsegment(data = meta_dens %>%  
                     filter(Highlight) %>% 
                     slice(which.max(y) + 22), 
                   aes(x, y, xend = x + 65, yend = y, label = Reference),
                   hjust = 1, family = "Futura", size = 3.5) +
  scale_colour_manual(values = c("#7030a5", "black"), guide = "none") +
  scale_alpha_manual(values = c(0.4, 1), guide = "none") +
  facet_wrap(~Family %>% fct_relevel("Strongylocentrotidae")) +
  xlab("Defecation (%)") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 14),
                  expand = FALSE, clip = "off") +
  mytheme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "cm"),
        plot.margin = margin(0, 0.5, 0.2, 0, unit = "cm"))

Fig_1b

# 2.8.3 Combine ####
Fig_1 <- ( Fig_1a_left | Fig_1a_right + guides(shape = "none") ) / Fig_1b +
  plot_annotation(tag_levels = list(c("a", "", "b"))) +
  plot_layout(heights = c(0.4, 1)) &
  theme(plot.tag = element_text(family = "Futura",
                                size = 15, face = "bold"),
        plot.tag.position = c(0, 1))

Fig_1 %>%
  ggsave(filename = "Fig_1.pdf", path = "Figures",
         device = cairo_pdf, width = 20, height = 22, units = "cm")

# Clean up
rm(list = ls())
gc()