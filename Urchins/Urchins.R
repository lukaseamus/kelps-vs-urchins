# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

biomass <- here("Urchins", "Biomass.csv") %>%
  read_csv(col_types = list("f", "f")) %>%
  mutate(Sporophyte = case_when(
           Season == "Autumn" ~ Sporophyte,
           Season == "Spring" ~ Sporophyte + 15,
           Season == "Summer" ~ Sporophyte + 30,
         ),
         Deployment = Deployment %>% dmy_hm(),
         Retrieval = Retrieval %>% dmy_hm(),
         Days = Deployment %--% Retrieval / ddays(),
         Treatment = case_when(
           Tank %>% str_detect("U") ~ "Grazed",
           Tank %>% str_detect("C") ~ "Control",
           Tank %>% str_detect("M") ~ "Mechanical"
         ) %>% fct(), 
         Loss = Initial - Final, # Biomass loss (g)
         Lamina_ratio = Lamina_dry / Lamina_wet, # Dry-wet mass ratios
         Faeces_ratio = Faeces_dry / Faeces_wet) %T>%
  print(n = 105)

urchins <- here("Urchins", "Urchins.csv") %>%
  read_csv(col_types = list("f", "f")) %>%
  mutate(Sporophyte = case_when(
           Season == "Autumn" ~ Tank %>% 
             str_extract("\\d+") %>% 
             as.numeric(),
           Season == "Spring" ~ Tank %>% 
             str_extract("\\d+") %>% 
             as.numeric() + 15,
           Season == "Summer" ~ Tank %>% 
             str_extract("\\d+") %>% 
             as.numeric() + 30,
         ) %>% str_c() %>% fct(),
         Deployment = Deployment %>% dmy_hm(),
         Retrieval = Retrieval %>% dmy_hm(),
         Days = Deployment %--% Retrieval / ddays(),
         Growth = (Final - Initial) / Initial * 100 / Days) %T>% # RGR (% d^-1)
  print(n = 226)

# 1.2 Extract grazing data ####
grazing <- biomass %>%
  select(Season, Sporophyte, Treatment, Loss, Faeces, Days, 
         Lamina_ratio, Faeces_ratio) %>%
  filter(Treatment != "Mechanical") %>%
  pivot_wider(names_from = Treatment,
              values_from = -c(Season, Sporophyte, Treatment)) %>%
  select(where(~ !all(is.na(.)))) %>%
  rename(Faeces = Faeces_Grazed, Faeces_ratio = Faeces_ratio_Grazed) %T>%
  print(n = 45)
  
# 1.3 Join data ####
grazing %<>%
  left_join(urchins %>%
              group_by(Season, Sporophyte) %>%
              summarise(Urchins_Initial = sum(Initial),
                        Urchins_Final = sum(Final),
                        Urchins_Diameter = mean(Diameter),
                        Urchins_Days = unique(Days)) %>%
              mutate(Sporophyte = Sporophyte %>% as.numeric()),
            by = c("Season", "Sporophyte")) %T>%
  print(n = 45)

# 1.4 Calculate consumption, defecation and absorption ####
grazing %<>%
  mutate(Consumption = if_else(
            Season == "Autumn", # For autumn I have dry-wet mass ratios for all treatments.
            ( Loss_Grazed / Days_Grazed * Lamina_ratio_Grazed - 
                Loss_Control / Days_Control * Lamina_ratio_Control ) * 1e3 /
              Urchins_Initial,
            ( Loss_Grazed / Days_Grazed - Loss_Control / Days_Control ) *
              Lamina_ratio_Grazed * 1e3 / Urchins_Initial
            ), # Consumption (mg_dry g_urchin^-1 d^-1)
         Defecation = Faeces / Days_Grazed * Faeces_ratio * 1e3 / Urchins_Initial,
         # Defecation (mg_dry g_urchin^-1 d^-1
         Absorption = Consumption - Defecation) # Absorption (mg_dry g_urchin^-1 d^-1)

# 1.5 Calculate urchin biomass incorporation ####
# I define incorporation as the proportion of absorbed biomass that is added as 
# growth rather than lost to respiration. Urchin mass was measuredas wet weight 
# because I did not want to kill urchins, so urchin growth cannot be determined 
# in terms of dry mass. Conversely, absorption cannot be determined in terms of 
# wet mass because faeces hold far more water than kelp. It is therefore not 
# possible to directly estimate incorporation. One would either have to convert 
# urchin wet mass to dry mass using published ratios, which introduces uncertainty,
# or assume that proportional absorption based on dry mass can be used to convert
# wet mass consumption.

# Here are urchin mass ratio data I downloaded from GitHub at 
# https://github.com/jmschuster/Urchin-physiology-habitat/blob/main/TPC_data.csv and
# https://github.com/jmschuster/Urchin-physiology-habitat/blob/main/Acute_data.csv.
urchin_mass <- bind_rows(
  here("Urchins", "External", "TPC_data.csv") %>%
    read_csv(col_types = list("f")) %>%
    select(Experiment, Wet_weight_g, `Empty weigh boat`, `Dry weight w boat`),
  here("Urchins", "External", "Acute_data.csv") %>%
    read_csv(col_types = list("f")) %>%
    select(Experiment, Wet_weight_g, `Empty weigh boat`, `Dry weight w boat`)
  ) %>%
  mutate(Dry = `Dry weight w boat` - `Empty weigh boat`) %>%
  rename(Wet = Wet_weight_g) %>%
  drop_na() %>%
  group_by(Experiment) %>%
  distinct(Wet, Dry)

urchin_ratio <- urchin_mass %>%
  mutate(Ratio = Dry / Wet) %$%
  mean(Ratio) %T>%
  print()
# About 34% of Strongylocentrotus droebachiensis wet mass is dry mass.

grazing %>% # In the first instance I convert mg to g, urchin wet to dry mass and proportion to percentage.
  mutate(Absorption_percent_dry = Absorption / 1e3 / urchin_ratio * 100,
         # In the second instance I multiply consumed wet mass per day by the absorbed dry proportion. 
         Absorption_percent_wet = ( Loss_Grazed / Days_Grazed - Loss_Control / Days_Control ) * 
           (Absorption / Consumption) / Urchins_Initial * 100) %>%
  ggplot(aes(Absorption_percent_dry, Absorption_percent_wet)) +
    geom_point() +
    geom_abline(slope = 1) + 
    theme_minimal()
# Far from identical. Best not make any of these assumptions and just look at 
# urchin growth in relation to consumption, both of which can be derived from
# wet mass as relative consumption rate (RCR) and relative growth rate (RGR).

grazing %<>%
  mutate(RCR = ( Loss_Grazed / Days_Grazed - Loss_Control / Days_Control ) / 
           Urchins_Initial * 100, # Consumption (% d^-1)
         RGR = ( Urchins_Final - Urchins_Initial ) / Urchins_Initial * 100 / 
           Urchins_Days) # Growth (% d^-1)

# 2. Defecation model ####
# 2.1 Prior simulation ####
# The relationship between consumption and defecation is the defecated proportion.
# This is available in the literature from Sauchyn & Scheibling 2009, doi
# 10.3354/meps08296, who report 64.8 ± 3% absorption for the S. droebachiensis.
# So my prior for the slope beta is centred on 0.352. beta has to be between 0
# and 1 because there cannot be a negative relationship between consumption and
# defecation and defecation cannot exceed consumption. The y intercept could
# be positive (defecation of previous consumption without new consumption) or
# negative (consumption without defecation due to digestive lag), so it's more
# informative to parameterise alpha as the defecation at mean consumption by
# centring the predictor variable consumption. Sauchyn & Scheibling 2009 report
# lab consumption and defecation rates for clean kelp as 14.87 and 2.34 mg g^-1 d^-1.
# alpha can then only take positive values and should be centred on 2.34. I will 
# use partial pooling across seasons for both parameters.

grazing %<>% # Centre predictor variable consumption
  mutate(Consumption_c = Consumption - mean(Consumption))

tibble(n = 1:1e3, # simulate hierachical prior
       alpha_mu = rgamma( 1e3 , 2.34^2 / 0.8^2 , 2.34 / 0.8^2 ),
       alpha_theta = rexp( 1e3 , 1 ),
       alpha = rgamma( 1e3 , alpha_mu / alpha_theta , 1 / alpha_theta ),
       beta_mu = rbeta( 1e3 , 0.352 * 8 , (1 - 0.352) * 8 ),
       beta_nu = rgamma( 1e3 , 30^2 / 20^2 , 30 / 20^2 ),
       beta = rbeta( 1e3 , beta_mu * beta_nu , (1 - beta_mu) * beta_nu )) %>% # %$% hist(beta)
  expand_grid(Consumption_c = grazing %$% 
                seq(min(Consumption_c), max(Consumption_c), length.out = 50)) %>%
  mutate(Defecation = alpha + beta * Consumption_c) %>%
  ggplot(aes(Consumption_c, Defecation, group = n)) +
    geom_hline(yintercept = grazing %>% 
                 drop_na(Defecation) %$% 
                 range(Defecation)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(ylim = c(0, 4), expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Covers all reasonable possibilities.

# 2.2 Stan model ####
require(cmdstanr)
grazing_model <- here("Urchins", "Stan", "grazing.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
grazing_samples <- grazing_model$sample(
          data = grazing %>%
            select(Consumption_c, Defecation, Season) %>%
            drop_na() %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

# 2.3 Model checks ####
# Rhat
grazing_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000683.

# Chains
require(bayesplot)
grazing_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are great.

# Pairs
grazing_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "beta_mu", "sigma"))

grazing_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1]", "beta[1]", "sigma"))
# No correlation.

# 2.4 Prior-posterior comparison ####
source("functions.R")
grazing_prior <- prior_samples(
  model = grazing_model,
  data = grazing %>%
    select(Consumption_c, Defecation, Season) %>%
    drop_na() %>%
    compose_data()
  )

grazing_prior %>% 
  prior_posterior_draws(
    posterior_samples = grazing_samples,
    group = grazing %>%
      select(Consumption_c, Defecation, Season) %>%
      drop_na() %>%
      select(Season),
    parameters = c("alpha_mu", "alpha_theta", "beta_mu", "beta_nu",
                   "alpha[Season]", "beta[Season]", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Season", ridges = FALSE)
# Posteriors are not overly constrained by priors.

# 2.5 Prediction ####
# 2.5.1 Priors and posteriors for hyperparameters ####
grazing_prior_posterior_hyper <- grazing_prior %>% 
  prior_posterior_draws(
    posterior_samples = grazing_samples,
    parameters = c("alpha_mu", "alpha_theta", 
                   "beta_mu", "beta_nu", "sigma"),
    format = "short"
  ) %>% # Calculate predictions for new seasons
  mutate(alpha_new = rgamma( n() , alpha_mu / alpha_theta , 1 / alpha_theta ),
         beta_new = rbeta( n() , beta_mu * beta_nu , (1 - beta_mu) * beta_nu ))

# 2.5.2 Priors and posteriors for seasonal parameters ####
grazing_prior_posterior_season <- grazing_prior %>% 
  prior_posterior_draws(
    posterior_samples = grazing_samples,
    parameters = c("alpha[Season]", "beta[Season]", "sigma"),
    format = "short"
  )

########################################################

# 2.5.3 Spread across predictor ####

# Spread across predictor
GP_predictions <- GP_prior_posterior %>%
  left_join(GP_data %>%
              group_by(Species) %>%
              summarise(min = min(nP_mean),
                        max = max(nP_mean)),
            by = "Species") %>%
  mutate(min = if_else(is.na(min),
                       GP_data %$% min(nP_mean),
                       min),
         max = if_else(is.na(max),
                       GP_data %$% max(nP_mean),
                       max)) %>%
  rowwise() %>%
  mutate(nP_mean = list( seq(min, max, length.out = 100) )) %>%
  select(-c(min, max)) %>%
  unnest(nP_mean) %>%
  mutate(mu = G0 + GP * nP_mean,
         obs = rnorm( n() , mu , sigma ))

GP_predictions_summary <- GP_predictions %>%
  group_by(Species, nP_mean) %>%
  reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
          obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
  unnest(c(mu, obs), names_sep = "_")


# 2.8.3 Remove redundant prior ####
C_prior_posterior %<>% # priors are identical for both treatments ->
  filter(!(Treatment == "Faeces" & distribution == "prior")) %>% # remove one
  mutate(Treatment = if_else(distribution == "prior", # add Prior to treatment
                             "Prior", Treatment) %>% fct()) %>%
  select(-distribution)

# 2.8.3 Plot predictions ####
require(ggridges) # ggridges is better than ggdist for limiting distribution ranges
C_prior_posterior %>%
  pivot_longer(cols = c(mu, obs, mu_new, obs_new), 
               values_to = "C", names_to = "Level") %>%
  filter(Level %in% c("mu_new", "obs_new")) %>%
  ggplot(aes(C, Treatment, alpha = Level)) +
    geom_density_ridges(from = 0, to = 1) +
    scale_alpha_manual(values = c(0.8, 0.2)) +
    scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.8.4 Convert to percentage ####
C_prior_posterior %<>%
  mutate(mu = mu * 100,
         obs = obs * 100,
         mu_new = mu_new * 100,
         obs_new = obs_new * 100)

C_prior_posterior %>%
  pivot_longer(cols = c(mu, obs, mu_new, obs_new), 
               values_to = "C", names_to = "Level") %>%
  filter(Level %in% c("mu_new", "obs_new")) %>%
  ggplot(aes(C, Treatment, alpha = Level)) +
    geom_density_ridges(from = 0, to = 100) +
    scale_alpha_manual(values = c(0.8, 0.2)) +
    scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.8.5 Calculate difference ####
C_diff <- C_prior_posterior %>%
  filter(Treatment != "Prior") %>%
  droplevels() %>%
  select(-c(alpha_t, sigma_s, sigma_i, nu)) %>%
  pivot_wider(names_from = Treatment, values_from = c(mu, obs, mu_new, obs_new)) %>%
  mutate(mu = mu_Kelp - mu_Faeces, # calculate differences
         obs = obs_Kelp - obs_Faeces,
         mu_new = mu_new_Kelp - mu_new_Faeces,
         obs_new = obs_new_Kelp - obs_new_Faeces) %>%
  select(.chain, .iteration, .draw, mu, obs, mu_new, obs_new) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = "Parameter",
               values_to = "Difference")

# 2.8.6 Plot difference ####
C_diff %>%
  ggplot(aes(Difference, Parameter)) +
    geom_density_ridges(from = -30, to = 30) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-30, 30), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.8.7 Summarise difference ####
C_diff_summary <- C_diff %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P = mean(Difference > 0),
            n = length(Difference)) %T>%
  print()

# 2.8.8 Add labels to C_diff ####
C_diff %<>%
  left_join(C_diff_summary %>%
              select(Parameter, P), 
            by = "Parameter") %>%
  mutate(label_Kelp = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_Faeces = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))






