# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

meta <- here("Biochemistry", "Biochemistry.csv") %>%
  read_csv() %>%
  mutate(
    Kelp = Food %>% # Add kelp binary
      str_detect("Egregia|Pterygophora|Alaria|Laminaria|Saccharina|Macrocystis|Nereocystis|Costaria|Agaru|Neoagarum|Ecklonia") &
      !str_detect(coalesce(Notes, ""), "reconstituted"),
    Genus = Urchin %>% str_extract("^\\S+"), # Extract urchin genus
    Function = if_else(
      Genus %>% # Distinguish kelp-feeders as functional group
        str_detect("Strongylocentrotus|Mesocentrotus|Paracentrotus|Psammechinus|Parechinus"),
      "Kelpivore", "Other"
    )
  ) %T>%
  print()

# 1.2 Simulate observations ####
# Unlike the defecation meta-analysis (Defecation.R), the amount of a
# given compound remaining in urchin faeces is not a true proportion
# because it can exceed 1 when urchin digestion enriches the faeces,
# e.g. as has been shown for calories and lipids (Dethier et al. 2019, 
# doi: 10.1016/j.jembe.2019.03.016). Instead it's a ratio and the data
# distribution is not a beta distribution, but a beta prime distribution.
# The beta or truncated normal distributions are therefore not useful.
# When plant and faecal content are given separately, I will simulate
# each as a gamma distribution and calculate the ratio. When the ratio
# is given with uncertainty I will draw from a beta prime distribution.
require(extraDistr) # R has no native beta prime
set.seed(100)
meta %<>%
  rowwise() %>%
  mutate(
    Ratio = if( !is.na(Percentage_mean) & is.na(Percentage_sd) ) {
      list( Percentage_mean / 100 ) # Convert % to ratio
    } else if( !is.na(Percentage_mean) & !is.na(Percentage_sd) & 
               !is.na(Percentage_n) ) {
      mean <- Percentage_mean / 100
      sd <- Percentage_sd / 100
      list( # This parameterisation is derived from mean * (1 + nu), 2 + nu.
        rbetapr(
          Percentage_n,
          mean * ( 1 + mean * ( 1 + mean ) / sd^2 ),
          2 + mean * ( 1 + mean ) / sd^2
        )
      )
    } else if( is.na(Percentage_mean) & !is.na(Faeces_mean) & 
               is.na(Faeces_sd) ) {
      list( Faeces_mean / Plant_mean ) # Simple ratio
    } else if( is.na(Percentage_mean) & !is.na(Faeces_mean) &
               !is.na(Faeces_sd) & !is.na(Faeces_n) ) {
      # Sometimes sample size is uneven
      n <- min(Faeces_n, Plant_n)
      # Simulate separately
      Faeces <- rgamma( 
        n, 
        Faeces_mean^2 / Faeces_sd^2 , 
        Faeces_mean / Faeces_sd^2
      )
      Plant <- rgamma(
        n, 
        Plant_mean^2 / Plant_sd^2, 
        Plant_mean / Plant_sd^2
      )
      # Calculate ratio
      list( Faeces / Plant )
    } else {
      list(NULL)
    }
  ) %>%
  unnest(Ratio) %T>%
  print()

# 1.3 Add data from this study ####
# The only data relevant to the meta-analysis are carbon and nitrogen.
experiment <- here("Biochemistry", "C_N", "RDS", "C_N.rds") %>% 
  read_rds() %>% select(-C_N) %>%
  pivot_wider(names_from = Treatment,
              values_from = c(C, N)) %>%
  mutate(Ratio_Nitrogen = N_Faeces / N_Kelp,
         Ratio_Carbon = C_Faeces / C_Kelp,
         Reference = "This study",
         Location = "Tromsø, Norway",
         Urchin = "Strongylocentrotus droebachiensis",
         Food = "Laminaria hyperborea",
         Unit = "%",
         Genus = Urchin %>% str_extract("^\\S+"),
         Kelp = TRUE, 
         Function = "Kelpivore") %>%
  # NAs need to be removed because in some cases we didn't match
  # kelp and faecal samples, so a ratio cannot be calculated.
  filter(!is.na(Ratio_Nitrogen) & !is.na(Ratio_Carbon)) %>%
  select(-c(Individual, starts_with("C_"), starts_with("N_"))) %>%
  pivot_longer(cols = c(Ratio_Nitrogen, Ratio_Carbon),
               names_to = "Compound",
               values_to = "Ratio",
               names_prefix = "Ratio_") %T>%
  print()

biochem <- meta %>%
  select(Reference, Location, Season, Genus, Urchin, Function, Food, Kelp,
         Compound, Unit, Ratio) %>%
  bind_rows(experiment) %>%
  mutate(Reference = Reference %>% fct(),
         Function = Function %>% fct(),
         Genus = Genus %>% fct(),
         Urchin = Urchin %>% fct(),
         Compound = Compound %>% fct()) %>%
  rename(Species = Urchin) %T>%
  print()

biochem %>% nrow() # 1461 observations
biochem %>% filter(Reference != "This study") %>% nrow() # 1393 observations
biochem %$% nlevels(Reference) # 14 references
biochem %$% nlevels(Genus) # 7 genera
biochem %$% levels(Genus) # Strongylocentrotus, Mesocentrotus, Psammechinus etc.
biochem %$% nlevels(Species) # 12 species
biochem %$% nlevels(Compound) # 6 biochemical measures
biochem %$% levels(Compound)

biochem %$% range(Ratio)
# There are zeros which is not allowed for true ratios.
# Replace zeros with 1e-4:
biochem %<>%
  mutate(
    Ratio = if_else(
      Ratio == 0, 1e-4, Ratio
    )
  ) %T>%
  print()

biochem %$% range(Ratio) # fixed

# 1.4 Data exploration ####
require(ggridges)
biochem %>%
  ggplot() +
    geom_vline(xintercept = 1) +
    stat_density_ridges(aes(x = Ratio, fill = Kelp, y = Species),
                        colour = NA, n = 2^10, rel_min_height = 0.001,
                        scale = 3, alpha = 0.6, from = 0, to = 5) +
    facet_grid(~Compound) +
    theme_minimal()

biochem %>%
  ggplot() +
    geom_vline(xintercept = 1) +
    stat_density_ridges(aes(x = Ratio, fill = Kelp, y = Genus), 
                        colour = NA, n = 2^10, rel_min_height = 0.001,
                        scale = 3, alpha = 0.6, from = 0, to = 5) +
    facet_grid(~Compound) +
    theme_minimal()

biochem %>%
  ggplot() +
    geom_vline(xintercept = 1) +
    stat_density_ridges(aes(x = Ratio, fill = Kelp, y = Function), 
                        colour = NA, n = 2^10, rel_min_height = 0.001,
                        scale = 3, alpha = 0.6, from = 0, to = 5) +
    facet_grid(~Compound) +
    theme_minimal()

biochem %>%
  ggplot() +
    geom_vline(xintercept = 1) +
    stat_density_ridges(aes(x = Ratio, y = Kelp), 
                        colour = NA, n = 2^10, rel_min_height = 0.001,
                        scale = 3, alpha = 0.6, from = 0, to = 5) +
    facet_grid(~Compound) +
    theme_minimal()

biochem %>%
  ggplot() +
    geom_vline(xintercept = 1) +
    geom_density(aes(Ratio)) +
    scale_x_continuous(limits = c(0, 5), oob = scales::oob_keep) +
    facet_grid(rows = vars(Compound)) +
    theme_minimal()

biochem %>% count(Genus, Species, Kelp, Function)
# Including food source in the model probably makes no sense, but I do want
# to distinguish functional groups in the prediction. Urchins that don't feed
# on kelp (mostly tropical species) can be summarised at their functional
# group level while kelp-feeding urchins will be reported at the species level.

biochem %>% count(Function, Compound)
# There is good representation of compound and of course this has to be
# included in the model.

biochem %>%
  summarise(Compounds = n_distinct(Compound),
            References = n_distinct(Reference),
            Observations = n(),
            .by = Species) %>%
  arrange(desc(Observations))
# 6 references and 456 observations on S. droebachiensis alone.

biochem %>%
  summarise(Species = n_distinct(Species),
            Food = n_distinct(Food),
            References = n_distinct(Reference),
            Observations = n(),
            .by = Compound) %>%
  arrange(desc(Observations))
# Good spread compounds.

biochem %>%
  summarise(Species = n_distinct(Species),
            Compounds = n_distinct(Compound),
            Observations = n(),
            .by = Reference) %>%
  arrange(desc(Observations))
# But most references only report on on species and a few compounds.

# 2. Model ####
# 2.1 Prior simulation ####
tibble(n = 1:1e4, # Model on the log scale
       alpha = rnorm( 1e4 , log(1) , 0.5 ), # global mean
       sigma_c = rtnorm( 1e4 , 0 , 0.5 , 0 ), # compound variation
       sigma_f = rtnorm( 1e4 , 0 , 0.5 , 0 ), # functional group variation
       sigma_g = rtnorm( 1e4 , 0 , 0.5 , 0 ), # genus variation
       sigma_s = rtnorm( 1e4 , 0 , 0.5 , 0 ), # species variation
       sigma_r = rtnorm( 1e4 , 0 , 0.5 , 0 ), # reference variation
       nu = rgamma( 1e4 , 20^2 / 10^2 , 20 / 10^2 )) %>% # likelihood precision
  mutate(
    c = rnorm( n() , alpha , sigma_c ),
    f = rnorm( n() , 0 , sigma_f ),
    g = rnorm( n() , 0 , sigma_g ),
    s = rnorm( n() , 0 , sigma_s ),
    r = rnorm( n() , 0 , sigma_r ),
    log_mu = c + f + g + s + r,
    mu = exp( log_mu ),
    r = rbetapr( n() , mu * (1 + nu) , 2 + nu )
  ) %>%
  pivot_longer(cols = c(log_mu, mu, r),
               names_to = "parameter") %>%
  ggplot() +
    geom_density(aes(value), fill = "black") +
    geom_vline(data = . %>% 
                 mutate(x = if_else(parameter == "log_mu", 0, 1)), 
               aes(xintercept = x), colour ="white") +
    scale_x_continuous(limits = c(NA, 10), oob = scales::oob_keep) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Even though the mu and r densities look like they favour values
# below 0, log_mu shows that prior probability is perfectly centered
# on a ratio of 1, which translates to an unbiased prior belief that 
# urchins do not biochemically transform with lots of uncertainty.

# 2.2 Stan model ####
require(cmdstanr)
biochem_c_model <- here("Biochemistry", "Stan", "biochem_c.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

biochem_nc_model <- here("Biochemistry", "Stan", "biochem_nc.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
biochem_c_samples <- biochem_c_model$sample(
          data = biochem %>%
            select(Ratio, Compound, Function, Genus, 
                   Species, Reference) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()

biochem_nc_samples <- biochem_nc_model$sample(
          data = biochem %>%
            select(Ratio, Compound, Function, Genus, 
                   Species, Reference) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()

# Save draws
biochem_c_samples$draws() %>%
  write_rds(here("Biochemistry", "RDS", "biochem_c_samples.rds"))
biochem_c_samples$draws(format = "df") %>%
  write_rds(here("Biochemistry", "RDS", "biochem_c_samples_df.rds"))

biochem_nc_samples$draws() %>%
  write_rds(here("Biochemistry", "RDS", "biochem_nc_samples.rds"))
biochem_nc_samples$draws(format = "df") %>%
  write_rds(here("Biochemistry", "RDS", "biochem_nc_samples_df.rds"))

# 2.3 Model diagnostics ####
# 2.3.1 R-hat and ESS ####
biochem_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 100% of rhat above 1.001. rhat = 1.03 ± 0.0121.

biochem_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000836.
# Already the non-centred model is more promising.

# 2.3.2 Chains ####
require(bayesplot) 
( biochem_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off") +
  theme(legend.position = "top", legend.justification = 0) ) %>%
  ggsave(filename = "biochem_c_chains.pdf", path = here("Biochemistry", "Plots"),
         device = cairo_pdf, width = 100, height = 100, units = "cm")

( biochem_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Non-centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off") +
  theme(legend.position = "top", legend.justification = 0) ) %>%
  ggsave(filename = "biochem_nc_chains.pdf", path = here("Biochemistry", "Plots"),
         device = cairo_pdf, width = 100, height = 100, units = "cm")
# Non-centred chains are clearly better.

# 2.3.3 Pairs ####
biochem_c_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("alpha", "c[1]", "c[4]", "f[1,2]", "f[4,1]", "g[1,1]", "g[4,5]", 
             "s[1,3]", "s[4,8]", "r[1,10]", "r[4,4]", "sigma_c", "sigma_f[1]", 
             "sigma_f[4]", "sigma_g[1]", "sigma_g[4]", "sigma_s[1]", "sigma_s[4]", 
             "sigma_r[1]", "sigma_r[4]", "log_nu_mu", "log_nu_sigma", "log_nu[1]", 
             "log_nu[4]"),
    grid_args = list(top = "Centred model")
  ) %>%
  ggsave(filename = "biochem_c_pairs.png", path = here("Biochemistry", "Plots"),
         width = 100, height = 100, units = "cm", bg = "white")

biochem_nc_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("alpha", "c[1]", "c[4]", "f[1,2]", "f[4,1]", "g[1,1]", "g[4,5]", 
             "s[1,3]", "s[4,8]", "r[1,10]", "r[4,4]", "sigma_c", "sigma_f[1]", 
             "sigma_f[4]", "sigma_g[1]", "sigma_g[4]", "sigma_s[1]", "sigma_s[4]", 
             "sigma_r[1]", "sigma_r[4]", "log_nu_mu", "log_nu_sigma", "log_nu[1]", 
             "log_nu[4]"),
    grid_args = list(top = "Non-centred model")
  ) %>%
  ggsave(filename = "biochem_nc_pairs.png", path = here("Biochemistry", "Plots"),
         width = 100, height = 100, units = "cm", bg = "white")
# No posterior correlation in either case. 
# Non-centred posteriors are much smoother.

# 2.4 Prior-posterior comparison ####
# 2.4.1 Sample prior ####
source("functions.R")
# prior_samples only works properly with non-centred models, but
# this is fine since priors are identical for both models
biochem_prior <- prior_samples(
  model = biochem_nc_model, 
  data = biochem %>%
    select(Ratio, Compound, Function, Genus, 
           Species, Reference) %>%
    compose_data()
)

# 2.4.2 Centred model ####
biochem_c_prior_posterior_global <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_c_samples,
    parameters = c("alpha", "sigma_c", "log_nu_mu", "log_nu_sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot() +
  labs(title = "Centred model")

biochem_c_prior_posterior_c <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_c_samples,
    group = biochem %>% select(Compound),
    parameters = c("c[Compound]", "sigma_f[Compound]", "sigma_g[Compound]", 
                   "sigma_s[Compound]", "sigma_r[Compound]", "log_nu[Compound]"),
    format = "long"
  ) %>%
  prior_posterior_plot(group_name = "Compound", ridges = TRUE)

biochem_c_prior_posterior_f <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_c_samples,
    group = biochem %>% select(Compound, Function),
    parameters = c("f[Compound, Function]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Compound", 
    second_group_name = "Function",
    ridges = TRUE
  )

biochem_c_prior_posterior_g <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_c_samples,
    group = biochem %>% select(Compound, Genus),
    parameters = c("g[Compound, Genus]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Compound", 
    second_group_name = "Genus",
    ridges = TRUE
  )

biochem_c_prior_posterior_s <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_c_samples,
    group = biochem %>% select(Compound, Species),
    parameters = c("s[Compound, Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Compound", 
    second_group_name = "Species",
    ridges = TRUE
  )

biochem_c_prior_posterior_r <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_c_samples,
    group = biochem %>% select(Compound, Reference),
    parameters = c("r[Compound, Reference]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Compound", 
    second_group_name = "Reference",
    ridges = TRUE
  )

# 2.4.3 Non-centred model ####
biochem_nc_prior_posterior_global <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_nc_samples,
    parameters = c("alpha", "sigma_c", "log_nu_mu", "log_nu_sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot() +
  labs(title = "Non-centred model")

biochem_nc_prior_posterior_c <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_nc_samples,
    group = biochem %>% select(Compound),
    parameters = c("c[Compound]", "sigma_f[Compound]", "sigma_g[Compound]", 
                   "sigma_s[Compound]", "sigma_r[Compound]", "log_nu[Compound]"),
    format = "long"
  ) %>%
  prior_posterior_plot(group_name = "Compound", ridges = TRUE)

biochem_nc_prior_posterior_f <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_nc_samples,
    group = biochem %>% select(Compound, Function),
    parameters = c("f[Compound, Function]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Compound", 
    second_group_name = "Function",
    ridges = TRUE
  )

biochem_nc_prior_posterior_g <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_nc_samples,
    group = biochem %>% select(Compound, Genus),
    parameters = c("g[Compound, Genus]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Compound", 
    second_group_name = "Genus",
    ridges = TRUE
  )

biochem_nc_prior_posterior_s <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_nc_samples,
    group = biochem %>% select(Compound, Species),
    parameters = c("s[Compound, Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Compound", 
    second_group_name = "Species",
    ridges = TRUE
  )

biochem_nc_prior_posterior_r <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_nc_samples,
    group = biochem %>% select(Compound, Reference),
    parameters = c("r[Compound, Reference]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Compound", 
    second_group_name = "Reference",
    ridges = TRUE
  )

# 2.4.4 Combined ####
biochem_prior_posterior <- 
  ( ( biochem_c_prior_posterior_global /
        biochem_c_prior_posterior_c / 
        biochem_c_prior_posterior_f / 
        biochem_c_prior_posterior_g / 
        biochem_c_prior_posterior_s / 
        biochem_c_prior_posterior_r ) +
      plot_layout(heights = c(6/84, 14/84, 12/84, 42/84, 72/84, 1)) |
    ( biochem_nc_prior_posterior_global /
        biochem_nc_prior_posterior_c / 
        biochem_nc_prior_posterior_f / 
        biochem_nc_prior_posterior_g / 
        biochem_nc_prior_posterior_s / 
        biochem_nc_prior_posterior_r ) +
      plot_layout(heights = c(6/84, 14/84, 12/84, 42/84, 72/84, 1)) ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", 
        legend.justification = 0)

biochem_prior_posterior %>%
  ggsave(filename = "biochem_prior_posterior.pdf", path = here("Biochemistry", "Plots"),
         device = cairo_pdf, width = 60, height = 100, units = "cm")
# Non-centred model posteriors are a lot smoother.
# Proceed with non-centred model.

# 2.5 Prediction ####
# 2.5.1 Global parameters ####
biochem_prior_posterior_global <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_nc_samples,
    parameters = c("alpha", "sigma_c", "log_nu_mu", "log_nu_sigma"),
    format = "short"
  ) %>%
  mutate( 
    # Calculate ratio for new compounds
    mu = exp( rnorm( n() , alpha , sigma_c ) ),
    nu = exp( rnorm( n() , log_nu_mu , log_nu_sigma ) ),
    Ratio = rbetapr( n() , mu * (1 + nu) , 2 + nu )
  ) %T>%
  print()

biochem_prior_posterior_global %>% # Save
  write_rds(here("Biochemistry", "RDS", "biochem_prior_posterior_global.rds"))

# 2.5.2 Compound parameters ####
biochem_prior_posterior_c <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_nc_samples,
    group = biochem %>% select(Compound),
    parameters = c("c[Compound]", "sigma_f[Compound]", "sigma_g[Compound]", 
                   "sigma_s[Compound]", "sigma_r[Compound]", "log_nu[Compound]"),
    format = "short"
  ) %>%
  mutate( 
    # Calculate ratio for new groups, genera, species and studies 
    # within existing compounds
    mu = exp(
      c + rnorm( n() , 0 , sigma_f ) +
        rnorm( n() , 0 , sigma_g ) +
        rnorm( n() , 0 , sigma_s ) +
        rnorm( n() , 0 , sigma_r )
    ),
    nu = exp( log_nu ),
    Ratio = rbetapr( n() , mu * (1 + nu) , 2 + nu )
  ) %>%
  filter(Compound == "Carbon" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in Compound
    Compound = if_else(
      distribution == "prior", "Prior", Compound
    ) %>% fct()
  ) %>%
  select(-distribution) %T>%
  print()

biochem_prior_posterior_c %>% # Save
  write_rds(here("Biochemistry", "RDS", "biochem_prior_posterior_c.rds"))

# 2.5.3 Functional group parameters ####
biochem_prior_posterior_f <- biochem_prior %>% 
  prior_posterior_draws(
    posterior_samples = biochem_nc_samples,
    group = biochem %>% select(Compound, Function),
    parameters = c("c[Compound]", "f[Compound, Function]", 
                   "sigma_g[Compound]", "sigma_s[Compound]", 
                   "sigma_r[Compound]", "log_nu[Compound]"),
    format = "short"
  ) %>%
  mutate( 
    # Calculate ratio for new genera, species and studies 
    # in existing functional groups per compound
    mu = exp(
      c + f + rnorm( n() , 0 , sigma_g ) +
        rnorm( n() , 0 , sigma_s ) +
        rnorm( n() , 0 , sigma_r )
    ),
    nu = exp( log_nu ),
    Ratio = rbetapr( n() , mu * (1 + nu) , 2 + nu )
  ) %>%
  filter(Compound == "Carbon" & Function == "Kelpivore" & 
           distribution == "prior" | distribution == "posterior") %>%
  mutate(
    Compound = if_else(
      distribution == "prior", "Prior", Compound
    ) %>% fct(),
    Function = if_else(
      distribution == "prior", "Prior", Function
    ) %>% fct()
  ) %>%
  select(starts_with("."), Compound, Function, mu, Ratio) %T>%
  print()

biochem_prior_posterior_f %>% # Save
  write_rds(here("Biochemistry", "RDS", "biochem_prior_posterior_f.rds"))

# 2.5.4 Genus parameters ####
biochem_prior_posterior_g <- biochem %>%
  distinct(Compound, Function, Genus) %>%
  # Genera are not matched to existing compounds or functional groups 
  # so this needs to be done manually
  left_join(
    biochem_prior %>% 
      prior_posterior_draws(
        posterior_samples = biochem_nc_samples,
        group = biochem %>% select(Compound, Function),
        parameters = c("c[Compound]", "f[Compound, Function]", "sigma_s[Compound]", 
                       "sigma_r[Compound]", "log_nu[Compound]"),
        format = "short"
      ),
    relationship = "many-to-many" # Compound and function are repeated across genera
  ) %>%
  left_join(
    biochem_prior %>% 
      prior_posterior_draws(
        posterior_samples = biochem_nc_samples,
        group = biochem %>% select(Compound, Genus),
        parameters = c("g[Compound, Genus]"),
        format = "short"
      )
  ) %>%
  mutate( 
    # Calculate ratio for new species and references in existing genera
    mu = exp(
      c + f + g + rnorm( n() , 0 , sigma_s ) +
        rnorm( n() , 0 , sigma_r )
    ),
    nu = exp( log_nu ),
    Ratio = rbetapr( n() , mu * (1 + nu) , 2 + nu )
  ) %>%
  filter(Compound == "Carbon" & Function == "Kelpivore" & 
           Genus == "Strongylocentrotus" & 
           distribution == "prior" | distribution == "posterior") %>%
  mutate(
    Compound = if_else(
      distribution == "prior", "Prior", Compound
    ) %>% fct(),
    Function = if_else(
      distribution == "prior", "Prior", Function
    ) %>% fct(),
    Genus = if_else(
      distribution == "prior", "Prior", Genus
    ) %>% fct()
  ) %>%
  select(starts_with("."), Compound, Function, Genus, mu, Ratio) %T>%
  print()

biochem_prior_posterior_g %>% # Save
  write_rds(here("Biochemistry", "RDS", "biochem_prior_posterior_g.rds"))

# 2.5.5 Species parameters ####
biochem_prior_posterior_s <- biochem %>%
  distinct(Compound, Function, Genus, Species) %>%
  left_join(
    biochem_prior %>% 
      prior_posterior_draws(
        posterior_samples = biochem_nc_samples,
        group = biochem %>% select(Compound, Function),
        parameters = c("c[Compound]", "f[Compound, Function]",
                       "sigma_r[Compound]", "log_nu[Compound]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    biochem_prior %>% 
      prior_posterior_draws(
        posterior_samples = biochem_nc_samples,
        group = biochem %>% select(Compound, Genus),
        parameters = c("g[Compound, Genus]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    biochem_prior %>% 
      prior_posterior_draws(
        posterior_samples = biochem_nc_samples,
        group = biochem %>% select(Compound, Species),
        parameters = c("s[Compound, Species]"),
        format = "short"
      )
  ) %>%
  mutate( 
    # Calculate ratio for new studies on existing species
    mu = exp(
      c + f + g + s + rnorm( n() , 0 , sigma_r )
    ),
    nu = exp( log_nu ),
    Ratio = rbetapr( n() , mu * (1 + nu) , 2 + nu )
  ) %>%
  filter(Compound == "Carbon" & Function == "Kelpivore" & 
           Genus == "Strongylocentrotus" & 
           Species == "Strongylocentrotus droebachiensis" &
           distribution == "prior" | distribution == "posterior") %>%
  mutate(
    Compound = if_else(
      distribution == "prior", "Prior", Compound
    ) %>% fct(),
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
  select(starts_with("."), Compound, Function, 
         Genus, Species, mu, Ratio) %T>%
  print()

biochem_prior_posterior_s %>% # Save
  write_rds(here("Biochemistry", "RDS", "biochem_prior_posterior_s.rds"))

# 2.5.6 Reference parameters ####
biochem_prior_posterior_r <- biochem %>%
  distinct(Compound, Function, Genus, Species, Reference) %>%
  left_join(
    biochem_prior %>% 
      prior_posterior_draws(
        posterior_samples = biochem_nc_samples,
        group = biochem %>% select(Compound, Function),
        parameters = c("c[Compound]", "f[Compound, Function]",
                       "log_nu[Compound]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    biochem_prior %>% 
      prior_posterior_draws(
        posterior_samples = biochem_nc_samples,
        group = biochem %>% select(Compound, Genus),
        parameters = c("g[Compound, Genus]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    biochem_prior %>% 
      prior_posterior_draws(
        posterior_samples = biochem_nc_samples,
        group = biochem %>% select(Compound, Species),
        parameters = c("s[Compound, Species]"),
        format = "short"
      ),
    relationship = "many-to-many"
  ) %>%
  left_join(
    biochem_prior %>% 
      prior_posterior_draws(
        posterior_samples = biochem_nc_samples,
        group = biochem %>% select(Compound, Reference),
        parameters = c("r[Compound, Reference]"),
        format = "short"
      )
  ) %>%
  mutate( 
    # Calculate ratio for existing references
    mu = exp( c + f + g + s + r ),
    nu = exp( log_nu ),
    Ratio = rbetapr( n() , mu * (1 + nu) , 2 + nu )
  ) %>%
  filter(Compound == "Carbon" & Function == "Kelpivore" & 
           Genus == "Strongylocentrotus" & 
           Species == "Strongylocentrotus droebachiensis" &
           Reference == "This study" &
           distribution == "prior" | distribution == "posterior") %>%
  mutate(
    Compound = if_else(
      distribution == "prior", "Prior", Compound
    ) %>% fct(),
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
  select(starts_with("."), Compound, Function, Genus, 
         Species, Reference, mu, Ratio) %T>%
  print()

biochem_prior_posterior_r %>% # Save
  write_rds(here("Biochemistry", "RDS", "biochem_prior_posterior_r.rds"))

# Clean up
rm(
  list = c(
    ls(pattern = "biochem_c|biochem_nc"),
    "biochem_prior", "biochem_prior_posterior"
  )
)
gc()

# 2.6 Parameter estimates ####
# 2.6.1 Global parameters ####
require(glue)
biochem_summary_global <- biochem_prior_posterior_global %>%
  mutate(
    log_mu = log(mu),
    log_Ratio = log(Ratio)
  ) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Ratio < 1 ),
    n = n(),
    .by = distribution
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_median} ({log_mu_mean} ± {log_mu_sd})" ),
    nu = glue( "{nu_mean} ± {nu_sd}" ),
    Ratio = glue( "{Ratio_median} ({log_Ratio_mean} ± {log_Ratio_sd})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.2 Compound parameters ####
biochem_summary_c <- biochem_prior_posterior_c %>%
  mutate(
    log_mu = log(mu),
    log_Ratio = log(Ratio)
  ) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Ratio < 1 ),
    n = n(),
    .by = Compound
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_median} ({log_mu_mean} ± {log_mu_sd})" ),
    nu = glue( "{nu_mean} ± {nu_sd}" ),
    Ratio = glue( "{Ratio_median} ({log_Ratio_mean} ± {log_Ratio_sd})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.3 Functional group parameters ####
biochem_summary_f <- biochem_prior_posterior_f %>%
  mutate(
    log_mu = log(mu),
    log_Ratio = log(Ratio)
  ) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Ratio < 1 ),
    n = n(),
    .by = c(Compound, Function)
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_median} ({log_mu_mean} ± {log_mu_sd})" ),
    Ratio = glue( "{Ratio_median} ({log_Ratio_mean} ± {log_Ratio_sd})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.4 Genus parameters ####
biochem_summary_g <- biochem_prior_posterior_g %>%
  mutate(
    log_mu = log(mu),
    log_Ratio = log(Ratio)
  ) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Ratio < 1 ),
    n = n(),
    .by = c(Compound, Function, Genus)
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_median} ({log_mu_mean} ± {log_mu_sd})" ),
    Ratio = glue( "{Ratio_median} ({log_Ratio_mean} ± {log_Ratio_sd})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.5 Species parameters ####
biochem_summary_s <- biochem_prior_posterior_s %>%
  mutate(
    log_mu = log(mu),
    log_Ratio = log(Ratio)
  ) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Ratio < 1 ),
    n = n(),
    .by = c(Compound, Function, Genus, Species)
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_median} ({log_mu_mean} ± {log_mu_sd})" ),
    Ratio = glue( "{Ratio_median} ({log_Ratio_mean} ± {log_Ratio_sd})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.6 Reference parameters ####
biochem_summary_r <- biochem_prior_posterior_r %>%
  mutate(
    log_mu = log(mu),
    log_Ratio = log(Ratio)
  ) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Ratio < 1 ),
    n = n(),
    .by = c(Compound, Function, Genus, Species, Reference)
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_median} ({log_mu_mean} ± {log_mu_sd})" ),
    Ratio = glue( "{Ratio_median} ({log_Ratio_mean} ± {log_Ratio_sd})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.6.7 Experiment parameters ####
experiment_summary <- bind_rows(
  Carbon = here("Biochemistry", "C_N", "RDS", "C_prior_posterior.rds") %>% 
    read_rds(),
  Nitrogen = here("Biochemistry", "C_N", "RDS", "N_prior_posterior.rds") %>% 
    read_rds(),
  Phenols = here("Biochemistry", "Phenol", "RDS", "phenol_prior_posterior.rds") %>% 
    read_rds(),
  Pigments = here("Biochemistry", "Pigments", "RDS", "total_prior_posterior.rds") %>% 
    read_rds(),
  .id = "Compound"
) %>%
  filter(Treatment != "Prior") %>%
  droplevels() %>%
  select(starts_with("."), Compound, Treatment, mu_new, obs_new) %>%
  pivot_wider(names_from = Treatment, values_from = c(mu_new, obs_new)) %>%
  mutate(
    mu = mu_new_Faeces / mu_new_Kelp,
    Ratio = obs_new_Faeces / obs_new_Kelp,
    log_mu = log(mu),
    log_Ratio = log(Ratio)
  ) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    P = mean( Ratio < 1 ), 
    n = n(),
    .by = Compound
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median"), "P"),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_median} ({log_mu_mean} ± {log_mu_sd})" ),
    Ratio = glue( "{Ratio_median} ({log_Ratio_mean} ± {log_Ratio_sd})" )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.7 Table 2 ####
# 2.7.1 Combine summaries ####
biochem_summary <- bind_rows(
  Experiment = experiment_summary,
  Metaanalysis = biochem_summary_s %>%
    filter(Function == "Kelpivore"),
  Metaanalysis = biochem_summary_f %>%
    filter(Function == "Other"),
  .id = "Study"
) %>%
  full_join(
    bind_rows(
      Experiment = experiment %>%
        summarise(Observations = n(),
                  .by = Compound) %>% # Same samples for all compounds
        add_row(Compound = c("Phenols", "Pigments"),
                Observations = .$Observations),
      Metaanalysis = biochem %>%
        summarise(References = n_distinct(Reference),
                  Observations = n(),
                  .by = c(Compound, Function, Genus, Species)) %>%
        filter(Function == "Kelpivore"),
      Metaanalysis = biochem %>%
        summarise(References = n_distinct(Reference),
                  Observations = n(),
                  .by = c(Compound, Function)) %>%
        filter(Function == "Other"),
      .id = "Study"
    )
  ) %>%
  mutate(
    Urchin = case_when(
      is.na(Function) ~ "Strongylocentrotus droebachiensis",
      Function == "Kelpivore" ~ Species,
      Function == "Other" ~ "Other urchins (7 species)"
    )
  ) %>%
  select(Study, Urchin, Compound, mu, Ratio, P, References, Observations) %T>%
  print(n = 30)

# 2.7.2 Arrange table ####
Table_2 <- biochem_summary %>%
  mutate(
    Urchin = Urchin %>% fct_relevel(
      "Strongylocentrotus droebachiensis",
      "Strongylocentrotus purpuratus",
      "Strongylocentrotus intermedius",
      "Mesocentrotus franciscanus",
      "Psammechinus miliaris"
    ),
    Compound = Compound %>% fct_relevel(
      "Carbon", "Nitrogen", "Phenols", "Pigments", 
      "Carbohydrates", "Proteins", "Lipids"
    )
  ) %>%
  arrange(Study, Urchin, Compound) %>%
  select(-c(Study, mu)) %T>%
  print(n = 30)

# 2.7.3 Save table ####
Table_2 %>%
  write_csv(here("Tables", "Table_2.csv"))

require(officer)
read_docx() %>%
  body_add_table(
    value = Table_2
  ) %>%
  print(target = here("Tables", "Table_2.docx"))

# 3. Figure 2 ####
# 3.1 Load data ####
# 3.1.1 Load carbon and nitrogen data ####
C_N <- 
  here("Biochemistry", "C_N", "RDS", "C_N.rds") %>%
  read_rds()
C_prior_posterior <- 
  here("Biochemistry", "C_N", "RDS", "C_prior_posterior.rds") %>%
  read_rds()
C_diff  <- 
  here("Biochemistry", "C_N", "RDS", "C_diff.rds") %>%
  read_rds()
N_prior_posterior <- 
  here("Biochemistry", "C_N", "RDS", "N_prior_posterior.rds") %>%
  read_rds()
N_diff <- 
  here("Biochemistry", "C_N", "RDS", "N_diff.rds") %>%
  read_rds()
C_N_prior_posterior <- 
  here("Biochemistry", "C_N", "RDS", "C_N_prior_posterior.rds") %>%
  read_rds()
C_N_diff <- 
  here("Biochemistry", "C_N", "RDS", "C_N_diff.rds") %>%
  read_rds()

# 3.1.2 Load phenol data ####
phenol <- 
  here("Biochemistry", "Phenol", "RDS", "phenol.rds") %>% 
  read_rds()
phenol_ID_dens <- 
  here("Biochemistry", "Phenol", "RDS", "ID_dens.rds") %>% 
  read_rds()
phenol_prior_posterior <- 
  here("Biochemistry", "Phenol", "RDS", "phenol_prior_posterior.rds") %>% 
  read_rds()
phenol_diff <- 
  here("Biochemistry", "Phenol", "RDS", "phenol_diff.rds") %>%
  read_rds()

# 3.1.3 Load pigments data ####
pigments <- 
  here("Biochemistry", "Pigments", "RDS", "pigments.rds") %>% 
  read_rds()
total_ID_dens <- 
  here("Biochemistry", "Pigments", "RDS", "total_ID_dens.rds") %>% 
  read_rds()
chlvspheo_ID_dens <- 
  here("Biochemistry", "Pigments", "RDS", "chlvspheo_ID_dens.rds") %>% 
  read_rds()
total_prior_posterior <- 
  here("Biochemistry", "Pigments", "RDS", "total_prior_posterior.rds") %>% 
  read_rds()
chlvspheo_prior_posterior <- 
  here("Biochemistry", "Pigments", "RDS", "chlvspheo_prior_posterior.rds") %>% 
  read_rds()
total_diff <- 
  here("Biochemistry", "Pigments", "RDS", "total_diff.rds") %>% 
  read_rds()
chlvspheo_diff <- 
  here("Biochemistry", "Pigments", "RDS", "chlvspheo_diff.rds") %>% 
  read_rds()

# 3.1.4 Rescale ID densities ####
phenol_ID_dens <- phenol %>%
  select(Treatment, Season, Individual, ID, Samples_Data) %>%
  unnest(cols = Samples_Data) %>%
  group_by(Treatment, Season, ID) %>% 
  # Starts from negative because some densities are not strictly positive.
  reframe(x = density(Concentration, n = 2^10, from = -0.1, to = 2)$x,
          y = density(Concentration, n = 2^10, from = -0.1, to = 2)$y) %>%
  group_by(ID) %>%
  mutate(y_area = y * 0.014 / ( sum(y) * ( x[2] - x[1] ) )) %>%
  ungroup() %>%
  filter(y > 0.1) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = c(x, x %>% rev()),
          y_area = c(y_area, -y_area %>% rev())) %>%
  ungroup()

total_ID_dens <- pigments %>%
  select(Treatment, Season, Individual, ID, Samples_Data) %>%
  unnest(cols = Samples_Data) %>%
  select(-c(Chlorophyll, Pheopigments)) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = density(Total, n = 2^10, from = 0, to = 3)$x,
          y = density(Total, n = 2^10, from = 0, to = 3)$y) %>%
  group_by(ID) %>%
  mutate(y_area = y * 0.02 / ( sum(y) * ( x[2] - x[1] ) )) %>%
  ungroup() %>%
  filter(y > 0.095) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = c(x, x %>% rev()),
          y_area = c(y_area, -y_area %>% rev())) %>%
  ungroup()

chlvspheo_ID_dens <- pigments %>%
  select(Treatment, Season, Individual, ID, Samples_Data) %>%
  unnest(cols = Samples_Data) %>%
  mutate(Proportion = Chlorophyll / ( Chlorophyll + Pheopigments ) * 100) %>%
  select(-c(Total, Chlorophyll, Pheopigments)) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = density(Proportion, n = 2^10, from = 0, to = 100)$x,
          y = density(Proportion, n = 2^10, from = 0, to = 100)$y) %>%
  group_by(ID) %>%
  mutate(y_area = y * 0.7 / ( sum(y) * ( x[2] - x[1] ) )) %>%
  ungroup() %>%
  filter(y > 0.003) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = c(x, x %>% rev()),
          y_area = c(y_area, -y_area %>% rev())) %>%
  ungroup()

# 3.2 Visualise ####
# 3.2.1 Carbon ####
# Define theme
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0, 1, 0, 0, unit = "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black", lineend = "square"),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.line.y = element_blank(),
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

# Prediction
require(ggridges)
Fig_2a_left_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = C, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.4, size = 2.5, height = 0.36) +
  stat_density_ridges(data = C_prior_posterior %>%
                        filter(Treatment != "Prior") %>% # comment out to show prior
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 60, rel_min_height = 0.001, 
                      bandwidth = 60*0.02, scale = 1.2, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 20),
                     oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE, order = 1)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15), # circle, triangle, square
                     # override grey fill legend shapes because they can be confused with prior
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Carbon content (%)") +
  coord_cartesian(ylim = c(0, 3), expand = FALSE, clip = "off") +
  mytheme

# Difference
require(geomtextpath)
# Fig_2a_left_bottom <- C_diff %>% 
#   filter(Parameter %in% c("mu_new", "obs_new")) %>%
#   mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
#   ggplot() +
#   stat_density_ridges(aes(x = Difference, y = Parameter,
#                           fill = if_else(after_stat(x) < 0,
#                                          "Faeces", "Kelp")), 
#                       geom = "density_ridges_gradient", n = 2^10,
#                       colour = NA, linewidth = 0, bandwidth = 1,
#                       from = -50, to = 50, rel_min_height = 0.001,
#                       scale = 1) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 17 + 1,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura", 
#                    size = 3.5, hjust = 0.8, vjust = 0,
#                    n = 2^10, bw = 1, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 17 + 2,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura", 
#                    size = 3.5, hjust = 0.75, vjust = 0,
#                    n = 2^10, bw = 1, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 17 + 1,
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura", 
#                    size = 3.5, hjust = 0.35, vjust = 0,
#                    n = 2^10, bw = 1, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 17 + 2,
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura", 
#                    size = 3.5, hjust = 0.32, vjust = 0,
#                    n = 2^10, bw = 1, text_only = TRUE) +
#   geom_vline(xintercept = 0) +
#   annotate("text", x = -50, y = c(1, 2), 
#            label = c("italic(tilde('y'))", "italic('µ')"),
#            hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
#            parse = TRUE) +
#   scale_x_continuous(limits = c(-50, 50), oob = scales::oob_keep,
#                      breaks = seq(-50, 50, 25),
#                      labels = scales::label_number(style_negative = "minus")) +
#   scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
#                     guide = "none") +
#   xlab("Δ carbon content (%)") +
#   coord_cartesian(expand = FALSE, clip = "off") +
#   mytheme

Fig_2a_left_bottom <- C_diff %>% 
  filter(parameter == "obs_new") %>%
  ggplot() +
  stat_density_ridges(aes(x = diff, y = 0,
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 100*0.02,
                      from = -50, to = 50, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 100*0.02, text_only = T) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.35, vjust = 0,
                   n = 2^10, bw = 100*0.02, text_only = T) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-50, 50), oob = scales::oob_keep,
                     breaks = seq(-50, 50, 25),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ carbon content (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 3.2.2 Nitrogen ####
# Prediction
Fig_2a_middle_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = N, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.4, size = 2.5, height = 0.36) +
  stat_density_ridges(data = N_prior_posterior %>%
                        filter(Treatment != "Prior") %>% # comment out to show prior
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 3, rel_min_height = 0.001, 
                      bandwidth = 3*0.02, scale = 1.2, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 3), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE, order = 1)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15),
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Nitrogen content (%)") +
  coord_cartesian(ylim = c(0, 3), expand = FALSE, clip = "off") +
  mytheme

# Difference
# Fig_2a_middle_bottom <- N_diff %>% 
#   filter(Parameter %in% c("mu_new", "obs_new")) %>%
#   mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
#   ggplot() +
#   stat_density_ridges(aes(x = Difference, y = Parameter, 
#                           fill = if_else(after_stat(x) < 0,
#                                          "Faeces", "Kelp")), 
#                       geom = "density_ridges_gradient", n = 2^10,
#                       colour = NA, linewidth = 0, bandwidth = 0.04,
#                       from = -2, to = 2, rel_min_height = 0.001,
#                       scale = 1) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 0.92 + 1,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura", 
#                    size = 3.5, hjust = 0.7, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 0.92 + 2,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura", 
#                    size = 3.5, hjust = 0.65, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 0.92 + 1,
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura", 
#                    size = 3.5, hjust = 0.3, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 0.92 + 2,
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura", 
#                    size = 3.5, hjust = 0.35, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_vline(xintercept = 0) +
#   annotate("text", x = -2, y = c(1, 2), 
#            label = c("italic(tilde('y'))", "italic('µ')"),
#            hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
#            parse = TRUE) +
#   scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
#                      labels = scales::label_number(style_negative = "minus")) +
#   scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
#                     guide = "none") +
#   xlab("Δ nitrogen content (%)") +
#   coord_cartesian(expand = FALSE, clip = "off") +
#   mytheme

Fig_2a_middle_bottom <- N_diff %>% 
  filter(Parameter == "obs_new") %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 4*0.02,
                      from = -2, to = 2, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.7, vjust = 0,
                   n = 2^10, bw = 4*0.02, text_only = TRUE) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.3, vjust = 0,
                   n = 2^10, bw = 4*0.02, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ nitrogen content (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 3.2.3 Carbon-nitrogen ####
# Prediction
Fig_2a_right_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = C_N, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.4, size = 2.5, height = 0.36) +
  stat_density_ridges(data = C_N_prior_posterior %>%
                        filter(Treatment != "Prior") %>% # comment out to show prior
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 120, rel_min_height = 0.001, 
                      bandwidth = 120*0.02, scale = 1.2, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 120), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE, order = 1)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15),
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Carbon-nitrogen ratio") +
  coord_cartesian(ylim = c(0, 3), expand = FALSE, clip = "off") +
  mytheme + # remove space for all right plots
  theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))

# Difference
# Fig_2a_right_bottom <- C_N_diff %>% 
#   filter(Parameter %in% c("mu_new", "obs_new")) %>%
#   mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
#   ggplot() +
#   stat_density_ridges(aes(x = Difference, y = Parameter, 
#                           fill = if_else(after_stat(x) < 0,
#                                          "Faeces", "Kelp")), 
#                       geom = "density_ridges_gradient", n = 2^10,
#                       colour = NA, linewidth = 0, bandwidth = 2.4,
#                       from = -120, to = 120, rel_min_height = 0.001,
#                       scale = 1) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 50 + 1,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura", 
#                    size = 3.5, hjust = 0.8, vjust = 0,
#                    n = 2^10, bw = 2.4, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 50 + 2,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura", 
#                    size = 3.5, hjust = 0.7, vjust = 0,
#                    n = 2^10, bw = 2.4, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 50 + 1,
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura", 
#                    size = 3.5, hjust = 0.41, vjust = 0,
#                    n = 2^10, bw = 2.4, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 50 + 2,
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura", 
#                    size = 3.5, hjust = 0.43, vjust = 0,
#                    n = 2^10, bw = 2.4, text_only = TRUE) +
#   geom_vline(xintercept = 0) +
#   annotate("text", x = -120, y = c(1, 2), 
#            label = c("italic(tilde('y'))", "italic('µ')"),
#            hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
#            parse = TRUE) +
#   scale_x_continuous(limits = c(-120, 120), breaks = seq(-120, 120, 60),
#                      oob = scales::oob_keep,
#                      labels = scales::label_number(style_negative = "minus")) +
#   scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
#                     guide = "none") +
#   xlab("Δ carbon-nitrogen ratio") +
#   coord_cartesian(expand = FALSE, clip = "off") +
#   mytheme

Fig_2a_right_bottom <- C_N_diff %>% 
  filter(Parameter == "obs_new") %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 240*0.02,
                      from = -120, to = 120, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 240*0.02, text_only = TRUE) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.41, vjust = 0,
                   n = 2^10, bw = 240*0.02, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-120, 120), breaks = seq(-120, 120, 60),
                     oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ carbon-nitrogen ratio") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme + # remove space for all right plots
  theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))

# 3.2.4 Phenol ####
# Prediction
Fig_2b_top <- ggplot() +
  geom_polygon(data = phenol_ID_dens %>% # Stratify by Treatment
                 mutate(y_area = y_area + if_else(Treatment == "Faeces", 0.5, 1.5)) %>%
                 group_by(ID) %>% # Jitter
                 mutate(y_area = y_area + runif( 1 , -0.36 , 0.36 )),
               aes(x = x, y = y_area, group = ID, 
                   fill = Treatment), 
               alpha = 0.4) +
  stat_density_ridges(data = phenol_prior_posterior %>%
                        filter(Treatment != "Prior") %>% # comment out to show prior
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), colour = NA, n = 2^10,
                      from = 0, to = 2, rel_min_height = 0.001, 
                      bandwidth = 2*0.02, scale = 2.6, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5),
                     labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1)),
                     oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = "none") +
  xlab("Phenolic content (%)") +
  coord_cartesian(ylim = c(0, 3), expand = FALSE, clip = "off") +
  mytheme

# Difference
# Fig_2b_bottom <- phenol_diff %>% 
#   filter(Parameter %in% c("mu_new", "obs_new")) %>%
#   mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
#   ggplot() +
#   stat_density_ridges(aes(x = Difference, y = Parameter, 
#                           fill = if_else(after_stat(x) < 0,
#                                          "Faeces", "Kelp")), 
#                       geom = "density_ridges_gradient", n = 2^10,
#                       colour = NA, linewidth = 0, bandwidth = 0.04,
#                       from = -2, to = 2, rel_min_height = 0.001,
#                       scale = 1) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 0.515 + 1,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura",
#                    size = 3.5, hjust = 0.84, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 0.515 + 2,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura",
#                    size = 3.5, hjust = 0.805, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 0.515 + 1, 
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura",
#                    size = 3.5, hjust = 0.32, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 0.515 + 2, 
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura",
#                    size = 3.5, hjust = 0.25, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_vline(xintercept = 0) +
#   annotate("text", x = -2, y = c(1, 2), 
#            label = c("italic(tilde('y'))", "italic('µ')"),
#            hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
#            parse = TRUE) +
#   scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
#                      breaks = seq(-2, 2, 1),
#                      labels = scales::label_number(style_negative = "minus")) +
#   scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
#                     guide = "none") +
#   xlab("Δ phenolic content (%)") +
#   coord_cartesian(expand = FALSE, clip = "off") +
#   mytheme

Fig_2b_bottom <- phenol_diff %>% 
  filter(parameter == "obs_new") %>%
  ggplot() +
  stat_density_ridges(aes(x = diff, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 4*0.02,
                      from = -2, to = 2, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura",
                   size = 3.5, hjust = 0.84, vjust = 0,
                   n = 2^10, bw = 4*0.02, text_only = TRUE) +
  geom_textdensity(aes(x = diff, y = after_stat(density), 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.3, vjust = 0,
                   n = 2^10, bw = 4*0.02, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                     breaks = seq(-2, 2, 1),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ phenolic content (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 3.2.5 Total pigment ####
# Prediction
require(magrittr)
Fig_2c_left_top <- ggplot() +
  geom_polygon(data = total_ID_dens %>% # Stratify by Treatment
                 mutate(y_area = y_area + if_else(Treatment == "Faeces", 0.5, 1.5)) %>%
                 group_by(ID) %>% # Jitter
                 mutate(y_area = y_area + runif( 1 , -0.36 , 0.36 )),
               aes(x = x, y = y_area, group = ID, 
                   fill = Treatment), 
               alpha = 0.4) +
  stat_density_ridges(data = total_prior_posterior %>%
                        filter(Treatment != "Prior") %>% # comment out to show prior
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), colour = NA, n = 2^10,
                      from = 0, to = 3, rel_min_height = 0.001, 
                      bandwidth = 3*0.02, scale = 1.2, alpha = 0.7) +
  geom_textpath(data = pigments %>%
                  select(ID, Treatment, Season, Samples_Data_Summary) %>%
                  unnest(cols = Samples_Data_Summary) %>%
                  filter(Season == "Spring" & Treatment == "Kelp") %>%
                  summarise(min = min(Total_mean - Total_sd),
                            max = max(Total_mean + Total_sd)) %$%
                  tibble(x = c(rep(min, 2), rep(max, 2)),
                         y = c(2.1, 2.2, 2.2, 2.1)),
                aes(x = x, y = y, label = "Spring"),
                colour = "#dabc23", linejoin = "mitre", lineend = "square",
                family = "Futura", size = 3.5, vjust = -0.1) +
  scale_x_continuous(limits = c(0, 3), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = "none") +
  xlab(expression("Total pigment (mg g"^-1*")")) +
  coord_cartesian(ylim = c(0, 3), expand = FALSE, clip = "off") +
  mytheme + # Adjust the title margin to counteract the superscript.
  theme(axis.title.x = element_text(margin = margin(t = 0)))

# Difference
# Fig_2c_left_bottom <- total_diff %>% 
#   filter(Parameter %in% c("mu_new", "obs_new")) %>%
#   mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
#   ggplot() +
#   stat_density_ridges(aes(x = Difference, y = Parameter, 
#                           fill = if_else(after_stat(x) < 0,
#                                          "Faeces", "Kelp")), 
#                       geom = "density_ridges_gradient", n = 2^10,
#                       colour = NA, linewidth = 0, bandwidth = 0.04,
#                       from = -2, to = 2, rel_min_height = 0.001,
#                       scale = 1) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 1.07 + 1,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura",
#                    size = 3.5, hjust = 0.71, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 1.07 + 2,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura",
#                    size = 3.5, hjust = 0.68, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 1.07 + 1, 
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura",
#                    size = 3.5, hjust = 0.35, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 1.07 + 2, 
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura",
#                    size = 3.5, hjust = 0.38, vjust = 0,
#                    n = 2^10, bw = 0.04, text_only = TRUE) +
#   geom_vline(xintercept = 0) +
#   annotate("text", x = -2, y = c(1, 2), 
#            label = c("italic(tilde('y'))", "italic('µ')"),
#            hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
#            parse = TRUE) +
#   scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
#                      breaks = seq(-2, 2, 1),
#                      labels = scales::label_number(style_negative = "minus")) +
#   scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
#                     guide = "none") +
#   xlab(expression("Δ total pigment (mg g"^-1*")")) +
#   coord_cartesian(expand = FALSE, clip = "off") +
#   mytheme + # The previous plot also affected this one in the same way.
#   theme(axis.title.x = element_text(margin = margin(t = 0)))

Fig_2c_left_bottom <- total_diff %>% 
  filter(Parameter == "obs_new") %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 4*0.02,
                      from = -2, to = 2, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura",
                   size = 3.5, hjust = 0.71, vjust = 0,
                   n = 2^10, bw = 4*0.02, text_only = TRUE) +
  geom_textdensity(aes(x = Difference, y = after_stat(density), 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.35, vjust = 0,
                   n = 2^10, bw = 4*0.02, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                     breaks = seq(-2, 2, 1),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab(expression("Δ total pigment (mg g"^-1*")")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme + # The previous plot also affected this one in the same way.
  theme(axis.title.x = element_text(margin = margin(t = 0)))

# 3.2.6 Chlorophyll vs. pheopigments ####
# Prediction
Fig_2c_right_top <- ggplot() +
  geom_polygon(data = chlvspheo_ID_dens %>% # Stratify by Treatment
                 mutate(y_area = y_area + if_else(Treatment == "Faeces", 0.5, 1.5)) %>%
                 group_by(ID) %>% # Jitter
                 mutate(y_area = y_area + runif( 1 , -0.36 , 0.36 )),
               aes(x = x, y = y_area, group = ID, 
                   fill = Treatment), 
               alpha = 0.4) +
  stat_density_ridges(data = chlvspheo_prior_posterior %>%
                        filter(Treatment != "Prior") %>% # comment out to show prior
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), colour = NA, n = 2^10,
                      from = 0, to = 100, rel_min_height = 0.001, 
                      bandwidth = 100*0.02, scale = 1.2, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = "none") +
  xlab(expression("Intact chlorophyll (%)")) +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme + # remove space for all right plots
  theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))

# Difference
# Fig_2c_right_bottom <- chlvspheo_diff %>% 
#   filter(Parameter %in% c("mu_new", "obs_new")) %>%
#   mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
#   ggplot() +
#   stat_density_ridges(aes(x = Difference, y = Parameter, 
#                           fill = if_else(after_stat(x) < 0,
#                                          "Faeces", "Kelp")), 
#                       geom = "density_ridges_gradient", n = 2^10,
#                       colour = NA, linewidth = 0, bandwidth = 1.6,
#                       from = -80, to = 80, rel_min_height = 0.001,
#                       scale = 1) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 42.5 + 1,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura",
#                    size = 3.5, hjust = 0.9, vjust = 0,
#                    n = 2^10, bw = 1.6, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 42.5 + 2,
#                        label = label_Kelp),
#                    colour = "#dabc23", family = "Futura",
#                    size = 3.5, hjust = 0.82, vjust = 0,
#                    n = 2^10, bw = 1.6, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
#                    aes(x = Difference, y = after_stat(density) * 42.5 + 1, 
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura",
#                    size = 3.5, hjust = 0.37, vjust = 0,
#                    n = 2^10, bw = 1.6, text_only = TRUE) +
#   geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
#                    aes(x = Difference, y = after_stat(density) * 42.5 + 2, 
#                        label = label_Faeces),
#                    colour = "#7030a5", family = "Futura",
#                    size = 3.5, hjust = 0.36, vjust = 0,
#                    n = 2^10, bw = 1.6, text_only = TRUE) +
#   geom_vline(xintercept = 0) +
#   annotate("text", x = -80, y = c(1, 2),
#            label = c("italic(tilde('y'))", "italic('µ')"),
#            hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
#            parse = TRUE) +
#   scale_x_continuous(limits = c(-80, 80), oob = scales::oob_keep,
#                      labels = scales::label_number(style_negative = "minus")) +
#   scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
#                     guide = "none") +
#   xlab("Δ intact chlorophyll (%)") +
#   coord_cartesian(expand = FALSE, clip = "off") +
#   mytheme

Fig_2c_right_bottom <- chlvspheo_diff %>% 
  filter(parameter == "obs_new") %>%
  ggplot() +
  stat_density_ridges(aes(x = diff, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 160*0.02,
                      from = -80, to = 80, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura",
                   size = 3.5, hjust = 0.9, vjust = 0,
                   n = 2^10, bw = 160*0.02, text_only = TRUE) +
  geom_textdensity(aes(x = diff, y = after_stat(density), 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.37, vjust = 0,
                   n = 2^10, bw = 160*0.02, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-80, 80), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ intact chlorophyll (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme + # remove space for all right plots
  theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))

# 3.2.7 Combine ####
require(patchwork)
Fig_2 <- ( ( Fig_2a_left_top | Fig_2a_middle_top | Fig_2a_right_top ) /
           ( ( Fig_2a_left_bottom | Fig_2a_middle_bottom | Fig_2a_right_bottom ) +
                # margin allows finer adjustment than plot_spacer
                theme(plot.margin = margin(t = 0.3, unit = "cm")) ) /
           plot_spacer() / # create space between the two grouped rows
           ( Fig_2b_top | Fig_2c_left_top | Fig_2c_right_top ) /
           ( ( Fig_2b_bottom | Fig_2c_left_bottom | Fig_2c_right_bottom ) +
                theme(plot.margin = margin(t = 0.3, unit = "cm")) ) ) +
  plot_annotation(tag_levels = list(c("a", rep("", 5), "b", "c", rep("", 4))),
                  theme = theme(plot.margin = margin(0.75, 0.45, 0.5, 0.5, unit = "cm"))) +
  plot_layout(heights = c(1, 0.3, 0.2, 1, 0.3),
              guides = "collect") &
  theme(plot.tag = element_text(family = "Futura",
                                size = 15, face = "bold"),
        plot.tag.position = c(-0.04, 1.06),
        legend.position = "top",
        legend.justification = 0)

Fig_2 %>%
  ggsave(filename = "Fig_2.pdf", device = cairo_pdf, path = "Figures", 
         height = 16, width = 20, units = "cm")

# 4. Figure 3 ####
# 4.1 Select relevant posteriors ####
biochem_prior_posterior <- bind_rows(
  Species = biochem_prior_posterior_s %>% 
    filter(Species == "Strongylocentrotus droebachiensis"),
  Species = biochem_prior_posterior_r %>% 
    filter(Species == "Strongylocentrotus droebachiensis"),
  Group = biochem_prior_posterior_f %>% filter(Function != "Prior"),
  Group = biochem_prior_posterior_r %>% filter(Function != "Prior"),
  .id = "Level"
) %>%
  mutate(
    Urchin = if_else(
      Level == "Species", Species, Function
    ) %>% fct_drop() %>%
      fct_relevel("Strongylocentrotus droebachiensis", "Kelpivore") %>%
      fct_recode("Other urchins (7 species)" = "Other",
                 "Kelp-feeding urchins (5 species)" = "Kelpivore"),
    Reference = if_else( # Ensure reference remains within species
      is.na(Reference), "Global", interaction(Species, Reference, sep = "_")
    ) %>% fct_relevel("Global", after = Inf),
    Compound = Compound %>% 
      fct_drop() %>% # Remove lingering Prior
      fct_relevel("Carbon", "Nitrogen", "Carbohydrates",
                  "Proteins", "Lipids") %>%
      fct_recode(Carbs = "Carbohydrates")
  ) %T>%
  print()

# 4.2 Prepare data for plotting ####
# Values for S. droebachiensis need to be duplicated.
biochem_dup <- bind_rows(
  Species = biochem %>% 
    filter(Species == "Strongylocentrotus droebachiensis"),
  Group = biochem,
  .id = "Level"
) %>%
  mutate(
    Urchin = if_else(
      Level == "Species", Species, Function
    ) %>% fct_drop() %>%
      fct_relevel("Strongylocentrotus droebachiensis", "Kelpivore") %>%
      fct_recode("Other urchins (7 species)" = "Other",
                 "Kelp-feeding urchins (5 species)" = "Kelpivore"),
    Compound = Compound %>%
      fct_relevel("Carbon", "Nitrogen", "Carbohydrates",
                  "Proteins", "Lipids") %>%
      fct_recode(Carbs = "Carbohydrates")
  ) %T>%
  print()

# 4.3 Build densities manually ####
biochem_dens <- biochem_prior_posterior %>%
  mutate(
    max = case_when(
      Compound %in% c("Lipids", "Calories") ~ 6,
      Compound == "Proteins" ~ 4,
      Compound %in% c("Carbs", "Nitrogen") ~ 3,
      Compound == "Carbon" ~ 2
    )
  ) %>%
  group_by(Compound, Level, Urchin, Reference) %>%
  reframe(
    x = c(
      0, 
      density(
        mu, n = 2^10, 
        from = 0, to = max[1], 
        bw = 2 * 0.02
      )$x,
      max[1]
    ),
    y = c(
      0, 
      density(
        mu, n = 2^10, 
        from = 0, to = max[1], 
        bw = 2 * 0.02
      )$y, 
      0
    )
  ) %>%
  group_by(Compound, Level, Urchin, Reference) %>% 
  # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate(y = y * 0.2 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

# 4.4 Visualise ####
require(ggh4x)
Fig_3 <- ggplot() +
  geom_jitter(data = biochem_dup %>% 
                filter(Compound == "Carbs" & Ratio <= 3 | # Remove outliers 
                         Compound %in% c("Lipids", "Calories") & Ratio <= 6 |
                         Compound %in% c("Carbon", "Nitrogen", "Proteins")),
              aes(x = Ratio, y = -0.5, colour = Ratio < 1),
              alpha = 0.4, shape = 16, size = 2.5, height = 0.36) +
  geom_line(data = biochem_dens %>% filter(Reference != "Global" & x < 1),
            aes(x, y, group = Reference), colour = "#dabc23", alpha = 0.4) +
  geom_line(data = biochem_dens %>% filter(Reference != "Global" & x >= 1),
            aes(x, y, group = Reference), colour = "#7030a5", alpha = 0.4) +
  geom_line(data = biochem_dens %>% filter(Reference == "Global"), 
            aes(x, y), lineend = "square") +
  geom_text(data = biochem_dup %>% 
              filter(Compound == "Carbon") %>% 
              distinct(Compound, Urchin) %>%
              mutate(face = if_else(Urchin == "Strongylocentrotus droebachiensis",
                                    "italic", "plain")),
            aes(0.115, 2, label = Urchin, fontface = face),
            family = "Futura", size.unit = "pt", size = 12, hjust = 0) +
  annotate("segment", x = 1, xend = 1, y = -1, yend = 1.6) +
  scale_colour_manual(values = c("#7030a5", "#dabc23"),
                      labels = c("FALSE" = "Faeces",
                                 "TRUE"  = "Food"),
                      guide = guide_legend(reverse = TRUE, nrow = 1)) +
  scale_alpha_manual(values = c(0.4, 1), guide = "none") +
  facet_grid(Urchin ~ Compound, scales = "free") +
  facetted_pos_scales(
    x = list(
      Compound == "Carbon" ~ scale_x_continuous(limits = c(0, 2), breaks = 0:2),
      Compound %in% c("Nitrogen", "Carbs") ~ 
        scale_x_continuous(limits = c(0, 3), breaks = 0:3),
      Compound %in% c("Proteins", "Lipids", "Calories") ~ 
        scale_x_continuous(breaks = seq(0, 6, 2))
    )
  ) +
  xlab("Ratio (Faeces / Food)") +
  coord_cartesian(ylim = c(-1, 2.2), expand = FALSE, clip = "off") +
  mytheme +
  theme(strip.text = element_text(face = "bold"),
        strip.text.y = element_blank(),
        plot.margin = margin(0, 0.5, 0.2, 0.5, unit = "cm"),
        legend.position = c(0.25, -0.1061))

Fig_3

Fig_3 %>%
  ggsave(filename = "Fig_3.pdf", device = cairo_pdf, path = "Figures", 
         height = 12, width = 20, units = "cm")

# Clean up
rm(list = ls())
gc()
