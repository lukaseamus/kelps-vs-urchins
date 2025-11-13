# 1. Prepare data ####
# 1.1 Load raw data ####
require(tidyverse)
require(magrittr)
require(here)

enzymes <- here("Microbes", "Enzymes", "Raw") %>%
  list.files(pattern = "\\.csv$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Name = Path %>% basename() %>% str_remove("\\..*$"),
    Date = Name %>% str_remove("^X") %>% str_split_i("_", 1) %>% ymd(),
    Fluorophore = Name %>% str_split_i("_", 2),
    Plate = Name %>% str_split_i("_", 3) %>% as.numeric(),
    Data = Path %>%
      map2(
        Date,
        ~ if(.y < dmy(190923)) {
          .x %>% # This is the data wrangling required for BMG Omega files.
            read_csv(skip = 9) %>%
            unite(Well, `Well\nRow`, `Well\nCol`, sep = "") %>%
            select(-Content) %>%
            pivot_longer(cols = -Well) %>%
            pivot_wider(names_from = Well,
                        values_from = value) %>%
            rename(Time = NANA) %>%
            mutate(
              name = case_when(
                name %>% str_detect("Raw Data") ~ "Fluorescence",
                name %>% str_detect("Temperature") ~ "Temperature"
              )
            ) %>%
            pivot_longer(cols = -c(Time, name),
                         names_to = "Well",
                         values_to = "value") %>%
            pivot_wider(names_from = name,
                        values_from = value)
        } else {
          .x %>% # And this is for BioTek Synergy H1 files.
            read_csv(skip = 30, col_select = -1) %>%
            rename(Temperature = starts_with("T°")) %>%
            mutate(
              Time = Time / dminutes(), # Convert to decimal minutes
              across(where(is.character), ~ .x %>% na_if("OVRFLW")), # Check overflow
              across(where(is.character), as.numeric)
            ) %>%
            pivot_longer(cols = -c(Time, Temperature),
                         names_to = "Well",
                         values_to = "Fluorescence")
        }
      )
    ) %>%
  select(-Path) %T>%
  print(n = 33)

enzymes$Data[[1]] 
enzymes$Data[[30]]

# 1.2 Load metadata ####
meta <- here("Microbes", "Enzymes", "Meta.csv") %>% 
  read_csv() %>%
  mutate(Date = Date %>% dmy()) %>%
  group_by(Date, Fluorophore, Plate) %>%
  nest(.key = "Metadata") %T>%
  print(n = 33)

# 1.3 Join metadata to data ####
# Higher-level join
enzymes %<>%
  full_join(meta, by = c("Date", "Fluorophore", "Plate")) %T>%
  print(n = 33)

enzymes$Metadata[[1]]
enzymes$Metadata[[30]]

# Lower-level join
enzymes %<>% 
  mutate(
    Data = Data %>%
      map2(
        Metadata, 
        ~ .x %>% 
          full_join(.y, by = "Well")
      )
  ) %>%
  select(-Metadata) %T>%
  print(n = 33)

enzymes$Data[[1]]
enzymes$Data[[30]]

# 1.4 Separate standard and samples ####
enzymes %<>% 
  mutate(
    Standard_Data = Data %>% 
      map(
        ~ .x %>%
          filter(Content == "Standard") %>%
          # Group distinct standard curves
          group_by(Treatment, Quenched) %>%
          mutate(Group = cur_group_id() %>% 
                   as.character() %>%
                   fct()) %>%
          ungroup()
      ), 
    Samples_Data = Data %>%
      map(
        ~ .x %>% 
          filter(Content %in% c("Sample", "Blank"))
      )
  ) %>%
  select(-Data) %T>%
  print()

enzymes$Standard_Data[[1]]
enzymes$Standard_Data[[30]]

enzymes$Samples_Data[[1]]
enzymes$Samples_Data[[30]]

# 1.5 Temperature ####
enzymes %<>%
  rowwise() %>% # rowwise is easier for nested plotting
  mutate(
    Temperature_Plot = 
      list(
        Samples_Data %>%
          distinct(Time, Temperature) %>%
            ggplot(aes(Time, Temperature)) +
              geom_point(shape = 16, size = 0.5, alpha = 0.5) +
              ggtitle(Name) +
              theme_minimal()
        )
    ) %>%
  ungroup() # undo rowwise

require(patchwork)
enzymes %$% 
  wrap_plots(Temperature_Plot) %>%
  ggsave(filename = "Temperature_Plot.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 50, height = 25)
# Temperature is measured at 0.1°C resolution on both readers.
# BMG Omega seems to be a lot more consistent with a range of 
# just 0.2°C, which BioTek Syngergy H1 heated up by several
# degrees in all cases. I was aware of this and tried to 
# counteract with a fan, but wasn't able to control temperature.

temp <- enzymes %>%
  mutate(
    Temperature = Samples_Data %>%
      map(
        ~ .x %>% distinct(Time, Temperature)
      )
  ) %>%
  unnest(Temperature) %>%
  select(Date, Fluorophore, Plate, Temperature) %T>%
  print()

temp %>%
  group_by(Date, Fluorophore, Plate) %>%
  summarise(
    across(
      Temperature, 
      list(
        mean = mean, 
        sd = sd, 
        min = min, 
        max = max
      )
    )
  ) %>%
  print(n = 33)

temp %>%
  group_by(Date < dmy(190923)) %>%
  summarise(
    across(
      Temperature, 
      list(
        mean = mean, 
        sd = sd, 
        min = min, 
        max = max
      )
    )
  ) %>%
  print()
# Samples run on BioTek Synergy H1 were on average 1.4°C warmer.
# All samples were run in the range 24.9 to 27.8°C.

# 2. Standard curve models ####
# 2.1 Visualise data ####
enzymes %<>%
  rowwise() %>%
  mutate(
    Standard_Plot_Time = 
      list(
        Standard_Data %>%
            ggplot(aes(Time, Fluorescence,
                       colour = Concentration)) +
              geom_point(shape = 16, size = 0.5, alpha = 0.5) +
              facet_grid(~ Treatment + Quenched) +
              ggtitle(Name) +
              theme_minimal()
        )
    ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Standard_Plot_Time) %>%
  ggsave(filename = "Standard_Time.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# Various different trajectories with time, so no clear model
# is apparent. But typically quenching leads to an exponential
# decline in fluorescence and in the absence of quenching 
# fluorescence is fairly linear. One solution would be to fit
# a standard curve for every timepoint and convert each sample
# fluorescence value with its time-specific standard curve.
# However, this assumes that discs were added to sample and
# standard wells at the same time, which was not the case. 
# Due to this temporal mismatch it is better to fit a single 
# standard curve, treating timepoints as replicates for each
# concentration. If this fails because of the amount of data,
# the nearest equivalent would be to calculate mean and 
# standard deviation for each concentration across time to 
# factor in temporal uncertainty. Let's have a look.

enzymes %<>%
  rowwise() %>%
  mutate(
    Standard_Plot_Concentration = 
      list(
        Standard_Data %>%
            ggplot(aes(Concentration, Fluorescence,
                       colour = Quenched, group = Group)) +
              geom_point(shape = 16, size = 0.5, alpha = 0.5) +
              geom_smooth(se = FALSE) +
              ggtitle(Name) +
              theme_minimal()
        )
    ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Standard_Plot_Concentration) %>%
  ggsave(filename = "Standard_Concentration.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# Fluorescence always saturates with increasing concentration. This is
# likely due to the inner filter effect. I will have to test which 
# saturating model fits best.

# 2.2 Prior simulation ####
# There are four saturating models that can be tested: rectangular hyperbola, 
# exponential saturation, hyperbolic tangent, and piece-wise linear.
# F0 (fluorescence intercept), Fmax (fluorescence maximum) and beta (linear 
# relationship between concentration and fluorescence) are necessarily positive. 
# F0 is expected to be very near zero but no other information is available, 
# so an exponential distribution is best. Fmax depends on gain adjustment, which
# was always automatically done, so is variable but is generally expected to be
# around 150 thousand. beta is expected to be near Fmax / 150 µM, so around
# 1000 fluorescence units µM^-1. Both need a gamma distribution. 
# A truncated normal likelihood ensures positive predictions but doesn't mess 
# with the nonlinear model like a gamma likelihood would.

require(extraDistr) # R doesn't have a built-in truncated normal distribution.
tibble(n = 1:1e3,
       F0 = rexp( 1e3 , 0.01 ),
       Fmax = rgamma( 1e3 , 15e4^2 / 8e4^2 , 15e4 / 8e4^2 ),
       beta = rgamma( 1e3 , 1e3^2 / 500^2 , 1e3 / 500^2 ),
       sigma = rexp( 1e3, 1e-4 )) %>%
  expand_grid(c = seq(0, 150, 2)) %>%
  mutate(mu_rh = F0 + Fmax * beta * c / ( Fmax + beta * c ),
         mu_es = F0 + Fmax * ( 1 - exp( -beta * c / Fmax ) ),
         mu_ht = F0 + Fmax * tanh( beta * c / Fmax ),
         mu_pl = if_else(c <= Fmax / beta, F0 + beta * c, F0 + Fmax),
         F_rh = rtnorm( n = n() , mean = mu_rh , sd = sigma , a = 0 ),
         F_es = rtnorm( n = n() , mean = mu_es , sd = sigma , a = 0 ),
         F_ht = rtnorm( n = n() , mean = mu_ht , sd = sigma , a = 0 ),
         F_pl = rtnorm( n = n() , mean = mu_pl , sd = sigma , a = 0 ),) %>%
  pivot_longer(cols = c(starts_with("mu"), starts_with("F_")),
               names_to = "par_mod",
               values_to = "value") %>%
  separate(par_mod, into = c("par", "mod")) %>%
  mutate(mod = mod %>% fct_relevel("rh", "es", "ht"),
         par = par %>% fct_relevel("mu")) %>%
  ggplot(aes(c, value, group = n)) +
    geom_line(alpha = 0.01) +
    facet_grid(par ~ mod) +
    coord_cartesian(expand = FALSE) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Covers all reasonable possibilities.

# 2.3 Model selection ####
# 2.3.1 Summarise data ####
# To effectively compare the models I need to calculate the log likelihood
# and LOO, but this is not possible for the enormous datasets I have, so in
# the first instance I will summarise data across time as mean.
enzymes %<>%
  mutate(
    Standard_Data_Summary = Standard_Data %>%
      map(
        ~ .x %>%
          group_by(Concentration, Group, Quenched, Treatment) %>% 
          summarise(Fluorescence = mean(Fluorescence),
                    n = n()) %>%
          ungroup()
      )
  )

enzymes$Standard_Data_Summary[[1]]

# enzymes %<>%
#   mutate(
#     Standard_Data_Subset = if_else(
#       row_number() %in% sample(n(), 8), # subsample 8 datasets
#       Standard_Data, list(NA)
#     )
#   )
# 
# enzymes$Standard_Data_Subset

# 2.3.1 Load models ####
require(cmdstanr)
standard_rh_model <- here("Microbes", "Enzymes", "Stan", "standard_rh.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

standard_es_model <- here("Microbes", "Enzymes", "Stan", "standard_es.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

standard_ht_model <- here("Microbes", "Enzymes", "Stan", "standard_ht.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

standard_pl_model <- here("Microbes", "Enzymes", "Stan", "standard_pl.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

# 2.3.3 Run models ####
require(tidybayes)
enzymes %<>%
  mutate(
    Standard_rh_Samples = Standard_Data_Summary %>%
      map(
        ~ standard_rh_model$sample(
          data = .x %>%
            select(Concentration, Fluorescence, Group) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_es_Samples = Standard_Data_Summary %>%
      map(
        ~ standard_es_model$sample(
          data = .x %>%
            select(Concentration, Fluorescence, Group) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_ht_Samples = Standard_Data_Summary %>%
      map(
        ~ standard_ht_model$sample(
          data = .x %>%
            select(Concentration, Fluorescence, Group) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_pl_Samples = Standard_Data_Summary %>%
      map(
        ~ standard_pl_model$sample(
          data = .x %>%
            select(Concentration, Fluorescence, Group) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      )
  )

# enzymes %<>%
#   mutate(
#     Standard_rh_loglik = Standard_Data_Subset %>%
#       map(
#         ~ if( is_tibble(.x) ) {
#           standard_rh_model$sample(
#             data = .x %>%
#               select(Concentration, Fluorescence, Group) %>%
#               compose_data(),
#             chains = 4, # Low number of chains
#             parallel_chains = parallel::detectCores(),
#             iter_warmup = 500, # Low number of iterations
#             iter_sampling = 500
#           )
#         } else {
#           NA
#         }
#       ),
#     Standard_es_loglik = Standard_Data_Subset %>%
#       map(
#         ~ if( is_tibble(.x) ) {
#           standard_es_model$sample(
#             data = .x %>%
#               select(Concentration, Fluorescence, Group) %>%
#               compose_data(),
#             chains = 4,
#             parallel_chains = parallel::detectCores(),
#             iter_warmup = 500,
#             iter_sampling = 500
#           )
#         } else {
#           NA
#         }
#       ),
#     Standard_ht_loglik = Standard_Data_Subset %>%
#       map(
#         ~ if( is_tibble(.x) ) {
#           standard_ht_model$sample(
#             data = .x %>%
#               select(Concentration, Fluorescence, Group) %>%
#               compose_data(),
#             chains = 4,
#             parallel_chains = parallel::detectCores(),
#             iter_warmup = 500,
#             iter_sampling = 500
#           )
#         } else {
#           NA
#         }
#       ),
#     Standard_pl_loglik = Standard_Data_Subset %>%
#       map(
#         ~ if( is_tibble(.x) ) {
#           standard_pl_model$sample(
#             data = .x %>%
#               select(Concentration, Fluorescence, Group) %>%
#               compose_data(),
#             chains = 4,
#             parallel_chains = parallel::detectCores(),
#             iter_warmup = 500,
#             iter_sampling = 500
#           )
#         } else {
#           NA
#         }
#       )
#   )

# 2.3.4 LOO ####
require(loo)
enzymes %<>%
  rowwise() %>%
  mutate(
    loo_comparison =
      list(
        loo_compare(
          list(
            rh = Standard_rh_Samples$loo(cores = parallel::detectCores()),
            es = Standard_es_Samples$loo(cores = parallel::detectCores()),
            ht = Standard_ht_Samples$loo(cores = parallel::detectCores()),
            pl = Standard_pl_Samples$loo(cores = parallel::detectCores())
          )
        )        
      )
  ) %>%
  ungroup()

enzymes$loo_comparison

loo <- enzymes$loo_comparison %>%
  imap(~ .x %>% 
         as.data.frame() %>%
         rownames_to_column(var = "model") %>%
         mutate(dataset = .y)) %>%
  bind_rows() %>%
  as_tibble() %T>%
  print()

loo %>%
  ggplot(aes(elpd_loo, model, colour = model, group = dataset)) +
    geom_pointrange(aes(xmin = elpd_loo - se_elpd_loo,
                        xmax = elpd_loo + se_elpd_loo),
                    position = position_dodge(width = 0.8),
                    alpha = 0.5, shape = 16) +
    scale_colour_discrete(guide = "none") +
    theme_minimal()
# Rectangular hyperbola and exponential saturation have among the highest
# ELPD but also have a greater spread than hyperbolic tangent for instance.

# Summarise ELPD across datasets
loo %>%
  mutate( 
    draws = map2(
      elpd_loo, se_elpd_loo, # Monte Carlo
      ~ rnorm( 1e4 , mean = .x , sd = .y )
    )
  ) %>%
  unnest(draws) %>%
  group_by(model) %>%
  summarise(elpd_loo_mean = mean(draws),
            elpd_loo_sd = sd(draws)) %>%
  arrange(desc(elpd_loo_mean))
# In reality ELPD is so variable that the rectangular hyperbola, exponential
# saturation and hyperbolic tangent are practically indistinguishable.

# enzymes %<>%
#   mutate(
#     Standard_rh_loo = Standard_rh_loglik %>%
#       map(
#         ~ if( .x %>% inherits("CmdStanMCMC") ) {
#           .x$loo(cores = 4)
#         } else {
#           NA
#         }
#       ),
#     Standard_es_loo = Standard_es_loglik %>%
#       map(
#         ~ if( .x %>% inherits("CmdStanMCMC") ) {
#           .x$loo(cores = 4)
#         } else {
#           NA
#         }
#       ),
#     Standard_ht_loo = Standard_ht_loglik %>%
#       map(
#         ~ if( .x %>% inherits("CmdStanMCMC") ) {
#           .x$loo(cores = 4)
#         } else {
#           NA
#         }
#       ),
#     Standard_pl_loo = Standard_pl_loglik %>%
#       map(
#         ~ if( .x %>% inherits("CmdStanMCMC") ) {
#           .x$loo(cores = 4)
#         } else {
#           NA
#         }
#       )
#   )

# 2.3.5 Rhat ####
enzymes %<>%
  mutate(
    summary_rh = Standard_rh_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_rh = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_rh = mean(rhat),
                    rhat_sd_rh = sd(rhat))
      ),
    summary_es = Standard_es_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_es = sum(rhat_check) / length(rhat),
                    rhat_mean_es = mean(rhat),
                    rhat_sd_es = sd(rhat))
      ),
    summary_ht = Standard_ht_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_ht = sum(rhat_check) / length(rhat),
                    rhat_mean_ht = mean(rhat),
                    rhat_sd_ht = sd(rhat))
      ),
    summary_pl = Standard_pl_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_pl = sum(rhat_check) / length(rhat),
                    rhat_mean_pl = mean(rhat),
                    rhat_sd_pl = sd(rhat))
      )
  ) %>%
  unnest(cols = c(summary_rh, summary_es, summary_ht, summary_pl))

enzymes %>% 
  select(Name, ends_with("rh")) %>%
  print(n = 33)
enzymes %>% 
  select(Name, ends_with("_es")) %>%
  print(n = 33)
enzymes %>% 
  select(Name, ends_with("ht")) %>%
  print(n = 33)
enzymes %>% 
  select(Name, ends_with("pl")) %>%
  print(n = 33)
# Only rectangular hyperbola has no rhat above 1.001, but only
# piecewise linear has more than one rhat above 1.001.

# 2.3.6 Chains ####
require(bayesplot)
enzymes %<>%
  rowwise() %>%
  mutate(
    Standard_rh_Chains =
      list(
        Standard_rh_Samples$draws(format = "df") %>%
          mcmc_rank_overlay(pars = c("F0[1]", "F0[2]", 
                                     "Fmax[1]", "Fmax[2]", 
                                     "beta[1]", "beta[2]")) +
          ggtitle(Name)
      ),
    Standard_es_Chains =
      list(
        Standard_es_Samples$draws(format = "df") %>%
          mcmc_rank_overlay(pars = c("F0[1]", "F0[2]", 
                                     "Fmax[1]", "Fmax[2]", 
                                     "beta[1]", "beta[2]")) +
          ggtitle(Name)
      ),
    Standard_ht_Chains =
      list(
        Standard_ht_Samples$draws(format = "df") %>%
          mcmc_rank_overlay(pars = c("F0[1]", "F0[2]", 
                                     "Fmax[1]", "Fmax[2]", 
                                     "beta[1]", "beta[2]")) +
          ggtitle(Name)
      ),
    Standard_pl_Chains =
      list(
        Standard_ht_Samples$draws(format = "df") %>%
          mcmc_rank_overlay(pars = c("F0[1]", "F0[2]", 
                                     "Fmax[1]", "Fmax[2]", 
                                     "beta[1]", "beta[2]")) +
          ggtitle(Name)
      )
    ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Standard_rh_Chains) %>%
  ggsave(filename = "Standard_rh_Chains.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
enzymes %$%
  wrap_plots(Standard_es_Chains) %>%
  ggsave(filename = "Standard_es_Chains.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
enzymes %$% 
  wrap_plots(Standard_ht_Chains) %>%
  ggsave(filename = "Standard_ht_Chains.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
enzymes %$% 
  wrap_plots(Standard_pl_Chains) %>%
  ggsave(filename = "Standard_pl_Chains.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# All chains look good.

# 2.3.7 Pairs ####
enzymes %<>%
  mutate(
    Standard_rh_Pairs = Standard_rh_Samples %>%
      map(
        ~ .x$draws(format = "df") %>% # mcmc_pairs does not allow adding titles
          mcmc_pairs(pars = c("F0[1]", "Fmax[1]", "beta[1]"))
      ),
    Standard_es_Pairs = Standard_es_Samples %>%
      map(
        ~ .x$draws(format = "df") %>%
          mcmc_pairs(pars = c("F0[1]", "Fmax[1]", "beta[1]"))
      ),
    Standard_ht_Pairs = Standard_ht_Samples %>%
      map(
        ~ .x$draws(format = "df") %>%
          mcmc_pairs(pars = c("F0[1]", "Fmax[1]", "beta[1]"))
      ),
    Standard_pl_Pairs = Standard_pl_Samples %>%
      map(
        ~ .x$draws(format = "df") %>%
          mcmc_pairs(pars = c("F0[1]", "Fmax[1]", "beta[1]"))
      )
  )

enzymes %$% 
  wrap_plots(Standard_rh_Pairs) %>%
  ggsave(filename = "Standard_rh_Pairs.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
enzymes %$% 
  wrap_plots(Standard_es_Pairs) %>%
  ggsave(filename = "Standard_es_Pairs.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
enzymes %$% 
  wrap_plots(Standard_ht_Pairs) %>%
  ggsave(filename = "Standard_ht_Pairs.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
enzymes %$% 
  wrap_plots(Standard_pl_Pairs) %>%
  ggsave(filename = "Standard_pl_Pairs.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# Some correlation between Fmax and beta for rectangular hyperbola,
# and exponential saturation, but slightly weaker for hyperbolic tangent.
# Some bimodal distributions and L-shaped correlations for piecewise 
# linear models, suggesting that no reliable prediction is possible.

# 2.3.8 Prior-posterior comparison ####
source("functions.R")
enzymes %<>%
  mutate(
    Standard_rh_Prior = Standard_Data_Summary %>%
      map(
        ~ prior_samples(model = standard_rh_model,
                        data = .x %>%
                          select(Concentration, Fluorescence, Group) %>%
                          compose_data())
      ),
    Standard_es_Prior = Standard_Data_Summary %>%
      map(
        ~ prior_samples(model = standard_es_model,
                        data = .x %>%
                          select(Concentration, Fluorescence, Group) %>%
                          compose_data())
      ),
    Standard_ht_Prior = Standard_Data_Summary %>%
      map(
        ~ prior_samples(model = standard_ht_model,
                        data = .x %>%
                          select(Concentration, Fluorescence, Group) %>%
                          compose_data())
      ),
    Standard_pl_Prior = Standard_Data_Summary %>%
      map(
        ~ prior_samples(model = standard_pl_model,
                        data = .x %>%
                          select(Concentration, Fluorescence, Group) %>%
                          compose_data())
      )
  )

enzymes %<>%
  mutate(
    Standard_rh_Prior_Posterior = pmap(
      list(Standard_rh_Prior, Standard_rh_Samples, Standard_Data_Summary),
      ~ prior_posterior_draws(prior_samples = ..1,
                              posterior_samples = ..2,
                              group = ..3 %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "short")
    ),
    Standard_es_Prior_Posterior = pmap(
      list(Standard_es_Prior, Standard_es_Samples, Standard_Data_Summary),
      ~ prior_posterior_draws(prior_samples = ..1,
                              posterior_samples = ..2,
                              group = ..3 %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "short")
    ),
    Standard_ht_Prior_Posterior = pmap(
      list(Standard_ht_Prior, Standard_ht_Samples, Standard_Data_Summary),
      ~ prior_posterior_draws(prior_samples = ..1,
                              posterior_samples = ..2,
                              group = ..3 %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "short")
    ),
    Standard_pl_Prior_Posterior = pmap(
      list(Standard_pl_Prior, Standard_pl_Samples, Standard_Data_Summary),
      ~ prior_posterior_draws(prior_samples = ..1,
                              posterior_samples = ..2,
                              group = ..3 %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "short")
    )
  )

enzymes %<>%
  rowwise() %>%
  mutate(
    Standard_rh_Prior_Posterior_Plot =
      list(
        prior_posterior_draws(prior_samples = Standard_rh_Prior,
                              posterior_samples = Standard_rh_Samples,
                              group = Standard_Data_Summary %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "long") %>%
          prior_posterior_plot(group_name = "Group", ridges = FALSE) +
          ggtitle(Name)
      ),
    Standard_es_Prior_Posterior_Plot =
      list(
        prior_posterior_draws(prior_samples = Standard_es_Prior,
                              posterior_samples = Standard_es_Samples,
                              group = Standard_Data_Summary %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "long") %>%
          prior_posterior_plot(group_name = "Group", ridges = FALSE) +
          ggtitle(Name)
      ),
    Standard_ht_Prior_Posterior_Plot =
      list(
        prior_posterior_draws(prior_samples = Standard_ht_Prior,
                              posterior_samples = Standard_ht_Samples,
                              group = Standard_Data_Summary %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "long") %>%
          prior_posterior_plot(group_name = "Group", ridges = FALSE) +
          ggtitle(Name)
      ),
    Standard_pl_Prior_Posterior_Plot =
      list(
        prior_posterior_draws(prior_samples = Standard_pl_Prior,
                              posterior_samples = Standard_pl_Samples,
                              group = Standard_Data_Summary %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "long") %>%
          prior_posterior_plot(group_name = "Group", ridges = FALSE) +
          ggtitle(Name)
      )
  ) %>%
  ungroup()

enzymes %$% 
  ( wrap_plots(Standard_rh_Prior_Posterior_Plot) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") ) %>%
  ggsave(filename = "Standard_rh_Prior_Posterior.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 100)

enzymes %$% 
  ( wrap_plots(Standard_es_Prior_Posterior_Plot) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom") ) %>%
  ggsave(filename = "Standard_es_Prior_Posterior.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 100)

enzymes %$% 
  ( wrap_plots(Standard_ht_Prior_Posterior_Plot) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom") ) %>%
  ggsave(filename = "Standard_ht_Prior_Posterior.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 100)

enzymes %$% 
  ( wrap_plots(Standard_pl_Prior_Posterior_Plot) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom") ) %>%
  ggsave(filename = "Standard_pl_Prior_Posterior.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 100)
# The posterior for Fmax is mostly shifted higher than the prior
# for rectangular hyperbola, suggesting the fluorescence limit may not
# be as well parameterised as for exponential saturation and hyperbolic
# tangent.

# 2.3.9 Prediction ####
enzymes %<>%
  mutate(
    Standard_rh_Prediction = Standard_rh_Prior_Posterior %>%
      map2(
        Standard_Data_Summary,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration",
                            group_name = "Group",
                            length = 50) %>%
          mutate(mu = F0 + Fmax * beta * Concentration / 
                      ( Fmax + beta * Concentration ),
                 obs = rtnorm( n = n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration, Group) %>%
          mean_qi(mu, obs, .width = c(.5, .8, .9))
      ),
    Standard_es_Prediction = Standard_es_Prior_Posterior %>%
      map2(
        Standard_Data_Summary,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration",
                            group_name = "Group",
                            length = 50) %>%
          mutate(mu = F0 + Fmax * ( 1 - exp( -beta * Concentration / Fmax ) ),
                 obs = rtnorm( n = n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration, Group) %>%
          mean_qi(mu, obs, .width = c(.5, .8, .9))
      ),
    Standard_ht_Prediction = Standard_ht_Prior_Posterior %>%
      map2(
        Standard_Data_Summary,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration",
                            group_name = "Group",
                            length = 50) %>%
          mutate(mu = F0 + Fmax * tanh( beta * Concentration / Fmax ),
                 obs = rtnorm( n = n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration, Group) %>%
          mean_qi(mu, obs, .width = c(.5, .8, .9))
      ),
    Standard_pl_Prediction = Standard_pl_Prior_Posterior %>%
      map2(
        Standard_Data_Summary,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration",
                            group_name = "Group",
                            length = 50) %>%
          mutate(mu = if_else(Concentration <= Fmax / beta, 
                              F0 + beta * Concentration,
                              F0 + Fmax),
                 obs = rtnorm( n = n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration, Group) %>%
          mean_qi(mu, obs, .width = c(.5, .8, .9))
      )
  )

enzymes %<>%
  mutate(
    Standard_rh_Prior_Posterior = Standard_rh_Prior_Posterior %>%
      map2(
        Standard_Data_Summary,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      ),
    Standard_es_Prior_Posterior = Standard_es_Prior_Posterior %>%
      map2(
        Standard_Data_Summary,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      ),
    Standard_ht_Prior_Posterior = Standard_ht_Prior_Posterior %>%
      map2(
        Standard_Data_Summary,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      ),
    Standard_pl_Prior_Posterior = Standard_pl_Prior_Posterior %>%
      map2(
        Standard_Data_Summary,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      ),
    Standard_rh_Prediction = Standard_rh_Prediction %>%
      map2(
        Standard_Data_Summary,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      ),
    Standard_es_Prediction = Standard_es_Prediction %>%
      map2(
        Standard_Data_Summary,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      ),
    Standard_ht_Prediction = Standard_ht_Prediction %>%
      map2(
        Standard_Data_Summary,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      ),
    Standard_pl_Prediction = Standard_pl_Prediction %>%
      map2(
        Standard_Data_Summary,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      )
  )

enzymes %<>%
  rowwise() %>% 
  mutate(
    Standard_rh_Prediction_Plot =
      list(
        Standard_rh_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data_Summary, aes(Concentration, Fluorescence),
                       shape = 16) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper,
                            alpha = factor(.width))) +
            geom_ribbon(data = . %>% filter(distribution == "prior", .width == 0.9),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            facet_grid(~ Quenched + Treatment) +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      ),
    Standard_es_Prediction_Plot =
      list(
        Standard_es_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data_Summary, aes(Concentration, Fluorescence),
                       shape = 16) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper,
                            alpha = factor(.width))) +
            geom_ribbon(data = . %>% filter(distribution == "prior", .width == 0.9),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            facet_grid(~ Quenched + Treatment) +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      ),
    Standard_ht_Prediction_Plot =
      list(
        Standard_ht_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data_Summary, aes(Concentration, Fluorescence),
                       shape = 16) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper,
                            alpha = factor(.width))) +
            geom_ribbon(data = . %>% filter(distribution == "prior", .width == 0.9),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            facet_grid(~ Quenched + Treatment) +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      ),
    Standard_pl_Prediction_Plot =
      list(
        Standard_pl_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data_Summary, aes(Concentration, Fluorescence),
                       shape = 16) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper,
                            alpha = factor(.width))) +
            geom_ribbon(data = . %>% filter(distribution == "prior", .width == 0.9),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            facet_grid(~ Quenched + Treatment) +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Standard_rh_Prediction_Plot) %>%
  ggsave(filename = "Standard_rh_Prediction.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
enzymes %$% 
  wrap_plots(Standard_es_Prediction_Plot) %>%
  ggsave(filename = "Standard_es_Prediction.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
enzymes %$% 
  wrap_plots(Standard_ht_Prediction_Plot) %>%
  ggsave(filename = "Standard_ht_Prediction.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
enzymes %$% 
  wrap_plots(Standard_pl_Prediction_Plot) %>%
  ggsave(filename = "Standard_pl_Prediction.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# Except for piecewise linear these are all good fits.

# Plot only the mean of mu for all models for better comparison.
enzymes %<>%
  rowwise() %>% 
  mutate(
    Standard_Comp_Prediction = 
      list(
        bind_rows(Standard_rh_Prediction,
                  Standard_es_Prediction,
                  Standard_ht_Prediction,
                  Standard_pl_Prediction) %>%
          mutate(model = rep(c("rh", "es", "ht", "pl"), each = n() / 4))
        ),
    Standard_Comp_Prediction_Plot =
      list(
        Standard_Comp_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data_Summary, 
                       aes(Concentration, Fluorescence),
                       shape = 16) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu, colour = model)) +
            facet_grid(~ Quenched + Treatment) +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Standard_Comp_Prediction_Plot) %>%
  ggsave(filename = "Standard_Comp_Prediction.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# Saturating models clearly fit much better than the piecwise linear model 
# but it is still hard to tell which of them fits best. Sometimes rectangular 
# hyperbola seems to fit best but other times the steeper hyperbolic tangent
# fits best. That's why loo was also very similar. Exponential saturation had
# a marginally better loo and is intermediate in steepness between rectangular
# hyperbola and hyperbolic tangent, causing it to fit better where the others
# over- or undershoot the data, so I'll go with that.

# Clean up
enzymes %<>%
  select(-c(contains("_rh"), contains("_es"), contains("_ht"), contains("_pl"),
            # Note that contains("_pl") also gets rid of plots
            loo_comparison, Standard_Comp_Prediction, Standard_Data_Summary)) %T>%
  print()

rm(loo)

# 2.4 Run optimal model ####
# Re-load optimal model (exponential saturation) without generated quantities block
standard_model <- here("Microbes", "Enzymes", "Stan", "standard.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

# Run model on full data
enzymes %<>%
  mutate(
    Standard_Samples = Standard_Data %>%
      map(
        ~ standard_model$sample(
          data = .x %>%
            select(Concentration, Fluorescence, Group) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      )
  ) %T>%
  print(n = 33)

enzymes$Standard_Samples[[1]]

# 2.5 Model checks ####
# 2.5.1 Rhat ####
enzymes %<>%
  mutate(
    summary = Standard_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean = mean(rhat),
                    rhat_sd = sd(rhat))
      )
  ) %>% 
  unnest(cols = summary) %T>%
  print(n = 33)
# No rhat > 1.001, all mean rhats = 1.00

# 2.5.2 Chains ####
enzymes %<>%
  rowwise() %>%
  mutate(
    Standard_Chains =
      list(
        Standard_Samples$draws(format = "df") %>%
          mcmc_rank_overlay(pars = c("F0[1]", "F0[2]", 
                                     "Fmax[1]", "Fmax[2]", 
                                     "beta[1]", "beta[2]")) +
          ggtitle(Name)
      )
    ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Standard_Chains) %>%
  ggsave(filename = "Standard_Chains.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# Chains look good

# 2.5.3 Pairs ####
enzymes %<>%
  mutate(
    Standard_Pairs = Standard_Samples %>%
      map(
        ~ .x$draws(format = "df") %>%
          mcmc_pairs(pars = c("F0[1]", "Fmax[1]", "beta[1]"))
      )
  )

enzymes %$% 
  wrap_plots(Standard_Pairs) %>%
  ggsave(filename = "Standard_Pairs.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# Some correlation between Fmax and beta but no concerning.

# 2.6 Prior-posterior comparison ####
# 2.6.1 Sample prior ####
enzymes %<>%
  mutate(
    Standard_Prior = Standard_Data %>%
      map(
        ~ prior_samples(model = standard_model,
                        data = .x %>%
                          select(Concentration, Fluorescence, Group) %>%
                          compose_data())
      )
  )

# 2.6.2 Combine prior and posterior ####
enzymes %<>%
  mutate(
    Standard_Prior_Posterior = pmap(
      list(Standard_Prior, Standard_Samples, Standard_Data),
      ~ prior_posterior_draws(prior_samples = ..1,
                              posterior_samples = ..2,
                              group = ..3 %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "short")
    )
  )

# 2.6.3 Plot comparison ####
enzymes %<>%
  rowwise() %>% 
  # rowwise does the same as pmap albeit less efficiently 
  # but makes adding plot elements such as titles easier
  mutate(
    Standard_Prior_Posterior_Plot =
      list(
        prior_posterior_draws(prior_samples = Standard_Prior,
                              posterior_samples = Standard_Samples,
                              group = Standard_Data %>% select(Group),
                              parameters = c("F0[Group]", "Fmax[Group]", 
                                             "beta[Group]", "sigma"),
                              format = "long") %>%
          prior_posterior_plot(group_name = "Group", ridges = FALSE) +
          ggtitle(Name)
      )
  ) %>%
  ungroup()

enzymes %$% 
  ( wrap_plots(Standard_Prior_Posterior_Plot) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") ) %>%
  ggsave(filename = "Standard_Prior_Posterior.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 100)

# 2.7 Prediction ####
# 2.7.1 Calculate prediction ####
enzymes %<>%
  mutate(
    Standard_Prediction = Standard_Prior_Posterior %>%
      map2(
        Standard_Data,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration",
                            group_name = "Group",
                            length = 50) %>%
          mutate(mu = F0 + Fmax * ( 1 - exp( -beta * Concentration / Fmax ) ),
                 obs = rtnorm( n = n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration, Group) %>%
          mean_qi(mu, obs, .width = c(.5, .8, .9))
      )
  )

# 2.7.2 Recover group names ####
enzymes %<>%
  mutate(
    Standard_Prior_Posterior = Standard_Prior_Posterior %>%
      map2(
        Standard_Data,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      ),
    Standard_Prediction = Standard_Prediction %>%
      map2(
        Standard_Data,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      )
  )

# 2.7.3 Plot prediction ####
enzymes %<>%
  rowwise() %>% 
  mutate(
    Standard_Prediction_Plot =
      list(
        Standard_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, 
                       aes(Concentration, Fluorescence),
                       shape = 16, alpha = 0.01) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper,
                            alpha = factor(.width))) +
            geom_ribbon(data = . %>% filter(distribution == "prior", .width == 0.9),
                        aes(Concentration, ymin = obs.lower, ymax = obs.upper),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            facet_grid(~ Quenched + Treatment) +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Standard_Prediction_Plot) %>%
  ggsave(filename = "Standard_Prediction.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)

# Clean up
enzymes %<>%
  select(-c(contains("rhat"), contains("Plot"), Standard_Samples, 
            Standard_Chains, Standard_Pairs, Standard_Prior)) %T>%
  print()

# Good time to save progress
enzymes %>%
  write_rds(here("Microbes", "Enzymes", "RDS", "enzymes.rds"))

# Load progress
enzymes <- here("Microbes", "Enzymes", "RDS", "enzymes.rds") %>%
  read_rds() %T>%
  print(n = 33)

# 3. Sample data conversion ####
# 3.1 Remove priors ####
enzymes %<>%
  mutate(
    Standard_Posterior = Standard_Prior_Posterior %>%
      map(~ .x %>% filter(distribution == "posterior"))
  ) %>%
  select(-Standard_Prior_Posterior) %T>%
  print()

# 3.2 Inverse prediction ####
# 3.2.1 Rearrange parameters ####
# The optimal standard curve model takes the form F = F0 + Fmax * 
# ( 1 - exp( -beta * C / Fmax ) ), so it predicts fluorescence with
# concentration. I want to do the inverse, that is predict concentration
# with fluorescence. Solving for concentration, I get C = -Fmax / beta *
# log( 1 - ( F - F0 ) / Fmax ).

# 3.2.2 Calculate inverse prediction ####
enzymes %<>%
  mutate(
    Standard_Inverse_Prediction = Standard_Posterior %>%
      map2(
        Standard_Data,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Fluorescence",
                            group_name = "Group",
                            length = 50) %>%
          # Note that obs cannot be predicted here because sigma maps 
          # onto Fluorescence not Concentration.
          mutate(mu = -Fmax / beta * log( 1 - ( Fluorescence - F0 ) / Fmax )) %>%
          # The log expression is prone to undefined values so we need to remove
          # them before summarising, otherwise the summary becomes NaN.
          filter(is.finite(mu)) %>%
          group_by(distribution, Fluorescence, Group) %>%
          mean_qi(mu, .width = c(.5, .8, .9))
      )
  )

# What about those NaNs?
# In 9 plates mu is undefined for some combinations of parameter values
# because the fluorescence I am trying to predict concentration with 
# exceeds the asymptote in those cases. This is generally not a problem 
# since it just restrains probability space for the predicted concentration.

# 3.2.3 Recover group names ####
enzymes %<>%
  mutate(
    Standard_Inverse_Prediction = Standard_Inverse_Prediction %>%
      map2(
        Standard_Data,
        ~ .x %>%
          full_join(
            .y %>% distinct(Treatment, Quenched, Group),
            by = "Group"
          )
      )
  )

# 3.2.4 Plot inverse prediction ####
enzymes %<>%
  rowwise() %>% 
  mutate(
    Standard_Inverse_Prediction_Plot =
      list(
        Standard_Inverse_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, 
                       aes(Fluorescence, Concentration),
                       shape = 16, alpha = 0.01) +
            geom_line(aes(Fluorescence, mu)) +
            geom_ribbon(aes(Fluorescence, ymin = .lower, ymax = .upper,
                            alpha = factor(.width))) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            facet_grid(~ Quenched + Treatment) +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Standard_Inverse_Prediction_Plot) %>%
  ggsave(filename = "Standard_Inverse_Prediction.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# Looks fine, except for panel 2 in plot 230921_4-MU_1 which cannot
# predict above fluorescence of
enzymes$Standard_Inverse_Prediction[[12]] %>%
  filter(Quenched == TRUE) %$%
  max(Fluorescence) # 11152.61
# which is not high enough to capture all sample fluorescence (see below), 
# effectively reducing the amount of data. This may be because the
# fluorophore standard and released substrate fluorophore are somehow 
# quenched differently. But its not a problem since we have a time series
# for each well and can still estimate a linear rate from the data within
# the standard curve range.

# 3.3 Visualise sample data ####
enzymes %<>%
  rowwise() %>%
  mutate(
    Samples_Fluorescence_Plot = 
      list(
        Samples_Data %>%
            ggplot(aes(Time, Fluorescence,
                       colour = Concentration)) +
              geom_point(shape = 16, size = 0.5, alpha = 0.5) +
              facet_grid(~ Treatment + Quenched) +
              ggtitle(Name) +
              theme_minimal()
        )
    ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Samples_Fluorescence_Plot) %>%
  ggsave(filename = "Samples_Fluorescence.pdf",
         path = here("Microbes", "Enzymes", "Plots"),
         device = cairo_pdf, units = "cm",
         width = 100, height = 50)
# The inner filter effect can also clearly be seen for sample data
# which plateau over time. There are two alternative ways to proceed:

# 1. Assume linearity of fluorescence over time and estimate the rate
# using a linear model. Then convert the distribution for the slope
# (F t^-1) using the standard curve parameters to get C t^-1.

# 2. Predict the distribution for C for each value of F, which should
# linearise the time series. Then summarise C as mean and s.d. for each
# timepoint and fit a linear model with a measurement error term to account 
# for uncertainty in C. The slope is C t^-1.

# Approach 1 may work for shorter time series for which linearity can 
# be assumed, but this is clearly not the case here. Also it makes a
# standard curve over the full spectrum of concentrations redundant
# because the only fluorescence that needs to be converted is the slope
# which can be chosen to be arbitrarily small by changing the unit of
# time. Finally, the unit of time and hence magnitude of the linear rate
# affects the result because quenching is nonlinear, so a quenched and
# unquenched rate may not be different if the chosen unit is s^-1 but they
# are certainly different if it is h^-1. Hence approach 2 is better here, 
# albeit more computationally intensive.

# 3.4 Re-nest data ####
enzymes %<>%
  select(Name, Date, Fluorophore, Plate, Samples_Data, Standard_Posterior) %>% 
  mutate(
    Samples_Data_Higher = Samples_Data %>%
      map(~ .x %>% distinct(Content, Treatment, Quenched)),
    Samples_Data_Lower = Samples_Data %>%
      map(~ .x %>% select(-c(Content, Treatment, Quenched)))
  ) %>%
  select(-Samples_Data) %>%
  unnest(cols = Samples_Data_Higher) %>%
  rename(Samples_Data = Samples_Data_Lower) %T>%
  print(n = 95)

###########

# 3.5 Filter relevant standard parameters ####
enzymes %>%
  mutate(
    Standard_Posterior_Filtered =
      pmap(
        list(Standard_Posterior, Name, Quenched, Treatment),
        ~ if(..2 == "230929_4-MU_1") {
          ..1 %>% filter(Quenched == ..3 & Treatment == ..4)
        } else {
          ..1 %>% filter(Quenched == ..3)
        }
      )
  ) %>%
  print(n = 95)


enzymes_unnested <- enzymes %>%
  select(-c(Standard_Data, Standard_Prediction, Standard_Inverse_Prediction,
            Standard_Inverse_Prediction_Plot, Samples_Fluorescence_Plot)) %>%
  unnest(cols = Samples_Data) %T>%
  print()
# It worked, but patience is required.

# 3.6 Calculate concentration ####
# I need to calculate mean and standard deviation of concentration
# given the standard parameters and their uncertainty. I would use 
# left_join() typically but that would result in the number of rows
# in enzymes_unnested multiplied by 8e4 posterior samples which is
# more than is allowed by tibble. I need a rowwise oepration.

enzymes_unnested %>%
  slice(750000:750100) %>%
  rowwise() %>%
  mutate(
    Concentration = list(
      ( # 230929_4-MU_1 only has quenched standards for two kelps
        if(Name != "230929_4-MU_1" & Quenched == TRUE) {
          Standard_Posterior %>% filter(Quenched == TRUE)
        } else if(Name != "230929_4-MU_1" & Quenched == FALSE) {
          Standard_Posterior %>% filter(Quenched == FALSE)
        } else if(Name == "230929_4-MU_1" & !is.na(Treatment) & 
                  Treatment == "Alaria esculenta") {
          Standard_Posterior %>% filter(Treatment == "Alaria esculenta")
        } else if(Name == "230929_4-MU_1" & !is.na(Treatment) & 
                  Treatment == "Laminaria digitata") {
          Standard_Posterior %>% filter(Treatment == "Laminaria digitata")
        }
      ) %>% (
        if(!is.null(.)) {
          mutate(Concentration = -Fmax / beta * log( 1 - ( Fluorescence - F0 ) / Fmax )) %>%
            filter(!is.nan(Concentration)) %>%
            summarise(Fmax = mean(Fmax),
                      beta = mean(beta),
                      F0 = mean(F0),
                      Concentration_mean = mean(Concentration), 
                      Concentration_sd = sd(Concentration),
                      Concentration_n = n())
        } else {
          tibble(Fmax = NA,
                 beta = NA,
                 F0 = NA,
                 Concentration_mean = NA, 
                 Concentration_sd = NA,
                 Concentration_n = NA)
        }
      )
    )
  ) %>%
  ungroup() %>% 
  unnest(cols = Concentration) %>%
  select(Name, Fluorescence, Fmax, beta, F0, starts_with("Concentration"))




phenol %<>%
  mutate(
    Samples_Data = Technical_Prediction %>%
      map2(
        Standard_ht_Prior_Posterior,
        ~ .x %>%
          select(-c(mu, sigma)) %>%
          full_join(
            .y %>% 
              filter(distribution == "posterior") %>%
              select(-c(distribution, sigma)),
            by = c(".chain", ".iteration", ".draw")
            ) %>% 
          mutate( # Here's where the magic happens!
            Concentration = Amax / beta * atanh( ( Absorbance - A0 ) / Amax )
            )
      )
  )

phenol$Samples_Data %>% 
  map(~ .x %>% filter(is.nan(Concentration)))
# Very few NaNs produced (at mots a few per sample). That's good enough.
# Still enough samples and it adds further regularisation.

# Remove NaNs.
phenol %<>%
  mutate(
    Samples_Data = Samples_Data %>%
      map(
        ~ .x %>% filter(!is.nan(Concentration))
      )
  )

phenol$Samples_Data %>% 
  map_lgl(~ .x %$% any(is.nan(Concentration))) %>%
  any()
# No more NaNs.

# Multiply 50% diluted samples by 2.
phenol %<>%
  mutate(
    Samples_Data = if_else(ID %>% str_detect("50%"),
                           Samples_Data %>%
                             map(
                               ~ .x %>%
                                 mutate(Concentration = Concentration * 2)
                             ),
                           Samples_Data),
    ID = if_else(ID %>% str_detect("50%"),
                 ID %>% str_remove("_50%"),
                 ID) %>% fct()
    )

# Calculate summary.
phenol %<>%
  mutate(
    Samples_Data_Summary = Samples_Data %>%
      map(
        ~ .x %>%
          summarise(Concentration_mean = mean(Concentration),
                    Concentration_sd = sd(Concentration))
      )
  )

# 4. Enzyme activity rate models ####

# 5. Experimental models ####
# 5.1 Urchin grazing ####

# 5.2 Enzyme kinetics ####


