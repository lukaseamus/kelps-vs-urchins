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
# Urchins in tank with sporophyte 44 died, so there are no defecation rates and
# grazing rates are unreliable.

grazing %<>%
  drop_na(Faeces) %T>%
  print(n = 44)

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

# 1.5 Visualise ####
grazing %>%
  ggplot(aes(Consumption, Defecation)) +
    geom_point() +
    theme_minimal()
# There's one clear outlier. Let's investigate.

grazing %>%
  ggplot(aes(Consumption, Defecation, label = Sporophyte)) +
  geom_label() +
  theme_minimal()
# Sporophyte 39 is the issue. Check data.

grazing %>%
  filter(Sporophyte %in% 35:45)
# The grazing rate is not out of the ordinary. It's the autogenic control 
# growth rate: 88.4 g in 9.92 days. And the lamina mass ratio also seems
# off (0.41 instead of around 0.25). Probably measurement error, so I'll 
# take means of the other summer values to calculate consumption for no. 39.

Loss_Day_Control_Summer_mean <- grazing %>%
  filter(Season == "Summer") %>%
  mutate(Loss_Day_Control = Loss_Control / Days_Control) %$%
  mean(Loss_Day_Control)
  
Lamina_ratio_Grazed_Summer_mean <- grazing %>%
  filter(Season == "Summer") %$%
  mean(Lamina_ratio_Grazed)

# 1.6 Re-calculate ####
grazing %<>%
  mutate(Consumption = case_when(
            Season == "Autumn" ~
              ( Loss_Grazed / Days_Grazed * Lamina_ratio_Grazed - 
                  Loss_Control / Days_Control * Lamina_ratio_Control ) * 1e3 /
              Urchins_Initial,
            Season != "Autumn" & Sporophyte != 39 ~
              ( Loss_Grazed / Days_Grazed - Loss_Control / Days_Control ) *
              Lamina_ratio_Grazed * 1e3 / Urchins_Initial,
            Sporophyte == 39 ~
              ( Loss_Grazed / Days_Grazed - Loss_Day_Control_Summer_mean ) *
              Lamina_ratio_Grazed_Summer_mean * 1e3 / Urchins_Initial
            ),
         Defecation = Faeces / Days_Grazed * Faeces_ratio * 1e3 / Urchins_Initial,
         Absorption = Consumption - Defecation)

grazing %>%
  ggplot(aes(Consumption, Defecation, shape = Season)) +
  geom_point() +
  theme_minimal()
# Fixed!

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
         Absorption_percent_wet = if_else(Sporophyte == 39,
                                          ( Loss_Grazed / Days_Grazed - Loss_Day_Control_Summer_mean ) * 
                                            (Absorption / Consumption) / Urchins_Initial * 100,
                                          ( Loss_Grazed / Days_Grazed - Loss_Control / Days_Control ) * 
                                            (Absorption / Consumption) / Urchins_Initial * 100)) %>%
  ggplot(aes(Absorption_percent_dry, Absorption_percent_wet)) +
    geom_point() +
    geom_abline(slope = 1) + 
    theme_minimal()
# Far from identical. Best not make any of these assumptions and just look at 
# urchin growth in relation to consumption, both of which can be derived from
# wet mass as relative consumption rate (RCR) and relative growth rate (RGR).

grazing %<>%
  mutate(RCR = if_else(Sporophyte == 39,
                       ( Loss_Grazed / Days_Grazed - Loss_Day_Control_Summer_mean ) / 
                         Urchins_Initial * 100,
                       ( Loss_Grazed / Days_Grazed - Loss_Control / Days_Control ) / 
                         Urchins_Initial * 100), # Consumption (% d^-1)
         RGR = ( Urchins_Final - Urchins_Initial ) / Urchins_Days /
           Urchins_Initial * 100 ) # Growth (% d^-1)

grazing %>%
  ggplot(aes(RCR, RGR, label = Sporophyte)) +
  geom_label() +
  theme_minimal()
# No clear trend, likely because food conversion cannot be measured over such
# a short experimental timespan.

# 2. Grazing model ####
# 2.1 Prior simulation ####
# The relationship between consumption and defecation is the defecated proportion.
# This is the complement to absorption and there are several reports of one or
# the other in the literature on Strongylocentrotus droebachiensis fed kelp. 
# Miller & Mann 1973 (doi: 10.1007/BF00348685) report average defecation of 34%, 
# Vadas 1977 (doi: 10.2307/1942173) reports 35%, Larson et al. 1980 (doi: 
# 10.1007/BF00396982) report 41%, and Sauchyn & Scheibling 2009 (doi: 
# 10.3354/meps08296) report 35%, but Mamelona & Pelletier 2005 (doi: 
# 10.1016/j.jembe.2004.08.026), report much higher values of around 74%. So my 
# prior for beta has to be between 0 and 1 because there cannot be a negative 
# relationship between consumption and defecation and defecation cannot exceed 
# consumption. In addition it should be centred on 44% (mean of above estimates),
# i.e. 0.44, allowing plenty of variation and giving even values as high as 0.8
# substantial probability. There is no intercept in this model since steady
# state consumption and defecation are assumed because urchins were fed 
# continuously and defecation showed no temporal trend. Faeces that are produced 
# from food eaten prior to the experiment are balanced by faeces that are not 
# yet produced from food eaten during the experiment. I will use partial pooling 
# across seasons for beta.

# Simulate hierachical prior
tibble(n = 1:1e3, 
       # The mean is pretty well-defined around 0.44 for S. droebachiensis.
       beta_mu = rbeta( 1e3 , 0.44 * 15 , (1 - 0.44) * 15 ), 
       # But seasonal variation is large (Fuji 1962, doi: 10.18960/seitai.12.5_181).
       beta_nu = rgamma( 1e3 , 10^2 / 7^2 , 10 / 7^2 ), # low nu = large variation
       beta = rbeta( 1e3 , beta_mu * beta_nu , (1 - beta_mu) * beta_nu )) %>% # %$% hist(beta)
  expand_grid(Consumption = grazing %$% 
                seq(min(Consumption), max(Consumption), length.out = 50)) %>%
  mutate(Defecation = beta * Consumption) %>%
  ggplot(aes(Consumption, Defecation, group = n)) +
    geom_hline(yintercept = grazing %$% 
                 range(Defecation)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(ylim = c(0, 5), expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Covers all possibilities.

# 2.2 Stan model ####
require(cmdstanr)
grazing_model <- here("Urchins", "Stan", "grazing.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
grazing_samples <- grazing_model$sample(
          data = grazing %>% 
            select(Season, Consumption, Defecation) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )

# 2.3 Model checks ####
# Rhat
grazing_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000689.

# Chains
require(bayesplot)
grazing_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are great.

# Pairs
grazing_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("beta_mu", "beta_nu"))

grazing_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("beta_mu", "sigma"))

grazing_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("beta[1]", "sigma"))
# No correlation.

# 2.4 Prior-posterior comparison ####
source("functions.R")
grazing_prior <- prior_samples(
  model = grazing_model,
  data = grazing %>% 
    select(Season, Consumption, Defecation) %>%
    compose_data()
  )

grazing_prior %>% 
  prior_posterior_draws(
    posterior_samples = grazing_samples,
    group = grazing %>% select(Season),
    parameters = c("beta_mu", "beta_nu", "beta[Season]", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Season", ridges = FALSE)
# Posteriors are not overly constrained by priors.

# 2.5 Prediction ####
# 2.5.1 Priors and posteriors for hyperparameters ####
grazing_prior_posterior_hyper <- grazing_prior %>% 
  prior_posterior_draws(
    posterior_samples = grazing_samples,
    parameters = c("beta_mu", "beta_nu", "sigma"),
    format = "short"
  ) %>% # Calculate predictions for new seasons, i.e. annual.
  mutate(beta = rbeta( n() , beta_mu * beta_nu , (1 - beta_mu) * beta_nu ))

# 2.5.2 Priors and posteriors for seasonal parameters ####
grazing_prior_posterior_season <- grazing_prior %>% 
  prior_posterior_draws(
    posterior_samples = grazing_samples,
    parameters = c("beta[Season]", "sigma"),
    format = "short"
  ) %>% 
  # Since I want only one grouping variable, there is redundancy in distribution.
  filter(!(Season %in% c("Summer", "Spring") & 
             distribution == "prior")) %>% # Remove two redundant priors.
  mutate(Season = if_else(distribution == "prior", # Add Prior to Season.
                          "Prior", Season) %>% fct()) %>%
  select(-distribution)

# 2.5.3 Combine seasonal and hyperparameters ####
grazing_prior_posterior <- grazing_prior_posterior_hyper %>%
  rename(Season = distribution) %>%
  mutate(Season = if_else(Season == "posterior",
                          "Annual", "Annual prior") %>% fct()) %>%
  select(-c(beta_mu, beta_nu)) %>%
  bind_rows(grazing_prior_posterior_season) 

require(ggridges)
grazing_prior_posterior %>%
  filter(Season %in% c("Prior", "Annual prior")) %>%
  ggplot(aes(beta, Season)) +
    geom_density_ridges()
# Priors are identical because the same beast is simulated in 
# grazing_prior_posterior_hyper and sampled in grazing_prior_posterior_season.
# Prior looks smoother than Annual prior though, so I'll remove Annual prior.

grazing_prior_posterior %<>%
  filter(Season != "Annual prior") %>% 
  droplevels()

# 2.5.4 Predict across predictor range ####
require(truncnorm)
grazing_prediction <- grazing_prior_posterior %>%
  spread_continuous(data = grazing, predictor_name = "Consumption",
                    group_name = "Season") %>%
  mutate(mu = beta * Consumption,
         Defecation = rtruncnorm( n() , mean = mu , sd = sigma , a = 0 ))

# 2.5.5 Summarise predictions ####
grazing_prediction_summary <- grazing_prediction %>%
  group_by(Season, Consumption) %>%
  mean_qi(mu, Defecation, .width = c(.5, .8, .9))

# 2.5.6 Convert beta to percentage defecation ####
grazing_prior_posterior %<>%
  mutate(beta = beta * 100) %T>%
  print()

# Summarise predictions
grazing_prior_posterior %>%
  group_by(Season) %>%
  summarise(beta_mean = mean(beta),
            beta_sd = sd(beta),
            n = n())

# 2.5.7 Calculate seasonal differences ####
grazing_diff <- grazing_prior_posterior %>%
  filter(!Season %in% c("Prior", "Annual")) %>%
  droplevels() %>%
  select(-sigma) %>%
  pivot_wider(names_from = Season, values_from = beta) %>%
  mutate(Spring_Summer = Spring - Summer, # calculate differences
         Spring_Autumn = Spring - Autumn,
         Autumn_Summer = Autumn - Summer) %>%
  select(-c(Autumn, Spring, Summer)) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = "Contrast",
               values_to = "Difference")

# 2.5.8 Summarise difference ####
grazing_diff_summary <- grazing_diff %>%
  group_by(Contrast) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P = mean(Difference > 0),
            n = length(Difference)) %T>%
  print()

# 2.5.9 Add labels to grazing_diff ####
grazing_diff %<>%
  left_join(grazing_diff_summary %>%
              select(Contrast, P), 
            by = "Contrast") %>%
  mutate(label_right = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_left = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))

# 3. Consumption model ####
# 3.1 Prior simulation ####
# Sauchyn & Scheibling 2009 (doi: 10.3354/meps08296) estimated consumption the
# same way as me (dry mass per urchin wet mass per day) and report around 15
# mg g^-1 d^-1 in the lab. Let's see how this compares to other studies which
# report consumption as dry mass per urchin per day.
urchins %$% mean(Initial) # My average urchin weighed 49 g and there were 5 = 
urchins %$% mean(Initial) * 5 # 244 g
urchins %$% mean(Diameter) # The average test diameter was 49 mm.
# Mamelona & Pelletier 2005 (doi: 10.1016/j.jembe.2004.08.026) report 
# consumption between 40 and 250 mg urchin^-1 d^-1 for similarly sized urchins.
# This transaltes to
40 / urchins %$% mean(Initial) * 5
250 / urchins %$% mean(Initial) * 5
# 4 to 26 mg g^-1 d^-1
# As it were the midpoint between these limits is also 15 mg g^-1 d^-1.
# So my prior should be centred on 15 but must span the entire range.

tibble(n = 1:1e5,
       alpha_mu = rnorm( 1e5 , log(15) , 1 ),
       alpha_sigma = rexp( 1e5 , 2 ),
       alpha = rnorm( 1e5 , alpha_mu , alpha_sigma ),
       theta = rexp( 1e5 , 1 ),
       mu = exp(alpha),
       y_new = rgamma( 1e5 , mu / theta , 1 / theta )) %>%
  pivot_longer(cols = c(mu, y_new), 
               names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, parameter)) +
    geom_vline(xintercept = grazing %$% range(Consumption)) +
    geom_density_ridges(alpha = 0.5, from = 0, to = 26,
                        colour = NA) +
    scale_x_continuous(limits = c(0, 26), oob = scales::oob_keep) +
    coord_cartesian(expand = FALSE) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks good.

# 3.2 Stan model ####
# Because multilevel modelling here works with normal hierarchical priors,
# there are two possible parameterisations: centred and non-centred.

consumption_c_model <- here("Urchins", "Stan", "consumption_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

consumption_nc_model <- here("Urchins", "Stan", "consumption_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

consumption_c_samples <- consumption_c_model$sample(
          data = grazing %>% 
            select(Season, Consumption) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
          adapt_delta = 0.99 # force sampler to slow down
        )

consumption_nc_samples <- consumption_nc_model$sample(
          data = grazing %>% 
            select(Season, Consumption) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
          adapt_delta = 0.99
        )

# 3.3 Model checks ####
# Rhat
consumption_c_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# About 20-30% of rhat above 1.001. rhat = 1.00 ± 0.00114.

consumption_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000128.

consumption_c_samples$summary() %>%
  left_join(consumption_nc_samples$summary(),
            by = "variable") %>%
  rename(rhat_c = rhat.x, rhat_nc = rhat.y) %>%
  ggplot(aes(rhat_c, rhat_nc)) +
  geom_abline(slope = 1) +
  geom_point() +
  theme_minimal() +
  theme(panel.grid = element_blank())
# The non-centred model has better rhat.

# Chains
consumption_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are fine.

consumption_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are great.

# Pairs
consumption_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_sigma", "theta"))

consumption_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_sigma", "theta"))
# No correlations. Models look similar.

# 3.4 Prior-posterior comparison ####
consumption_c_prior <- prior_samples(
  model = consumption_c_model,
  data = grazing %>% 
    select(Season, Consumption) %>%
    compose_data(),
  adapt_delta = 0.99
  )

consumption_nc_prior <- prior_samples(
  model = consumption_nc_model,
  data = grazing %>% 
    select(Season, Consumption) %>%
    compose_data(),
  adapt_delta = 0.99
)

consumption_c_prior %>% 
  prior_posterior_draws(
    posterior_samples = consumption_c_samples,
    group = grazing %>% select(Season),
    parameters = c("alpha_mu", "alpha_sigma", "alpha[Season]", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Season", ridges = FALSE)

consumption_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = consumption_nc_samples,
    group = grazing %>% select(Season),
    parameters = c("alpha_mu", "alpha_sigma", "alpha_z[Season]",
                   "alpha[Season]", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Season", ridges = FALSE)
# Posteriors are not overly constrained by priors. Models
# are similar. I'll pick the non-centred model based on
# better rhat and chains.

# 3.5 Prediction ####
# 3.5.1 Priors and posteriors for hyperparameters ####
consumption_prior_posterior_hyper <- consumption_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = consumption_nc_samples,
    parameters = c("alpha_mu", "alpha_sigma", "theta"),
    format = "short"
  ) %>% # Calculate predictions for new seasons, i.e. annual.
  mutate(alpha = rnorm( n() , alpha_mu , alpha_sigma ))

# 3.5.2 Priors and posteriors for seasonal parameters ####
consumption_prior_posterior_season <- consumption_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = consumption_nc_samples,
    parameters = c("alpha[Season]", "theta"),
    format = "short"
  ) %>% 
  # Since I want only one grouping variable, there is redundancy in distribution.
  filter(!(Season %in% c("Summer", "Spring") & 
             distribution == "prior")) %>% # Remove two redundant priors.
  mutate(Season = if_else(distribution == "prior", # Add Prior to Season.
                          "Prior", Season) %>% fct()) %>%
  select(-distribution)

# 3.5.3 Combine seasonal and hyperparameters ####
consumption_prior_posterior <- consumption_prior_posterior_hyper %>%
  rename(Season = distribution) %>%
  mutate(Season = if_else(Season == "posterior",
                          "Annual", "Annual prior") %>% fct()) %>%
  select(-c(alpha_mu, alpha_sigma)) %>%
  bind_rows(consumption_prior_posterior_season) 

consumption_prior_posterior %>%
  filter(Season %in% c("Prior", "Annual prior")) %>%
  ggplot(aes(alpha, Season)) +
    geom_density_ridges()
# Priors are identical because the same beast is simulated in 
# grazing_prior_posterior_hyper and sampled in grazing_prior_posterior_season.
# Annual prior looks smoother than Prior, so I'll remove the latter.

consumption_prior_posterior %<>%
  filter(Season != "Prior") %>% 
  mutate(Season = Season %>% 
           fct_drop() %>% 
           fct_recode(Prior = "Annual prior"))

# 3.5.4 Predict new observations ####
consumption_prior_posterior %<>%
  mutate(mu = exp(alpha),
         obs = rgamma( n() , mu / theta , 1 / theta ))

# 3.5.5 Calculate seasonal differences ####
consumption_diff <- consumption_prior_posterior %>%
  filter(!Season %in% c("Prior", "Annual")) %>%
  droplevels() %>%
  select(-c(alpha, theta)) %>%
  pivot_wider(names_from = Season, values_from = c(mu, obs)) %>%
  mutate(mu_Spring_Summer = mu_Spring - mu_Summer, # calculate differences
         mu_Autumn_Spring = mu_Autumn - mu_Spring,
         mu_Autumn_Summer = mu_Autumn - mu_Summer,
         obs_Spring_Summer = obs_Spring - obs_Summer,
         obs_Autumn_Spring = obs_Autumn - obs_Spring,
         obs_Autumn_Summer = obs_Autumn - obs_Summer) %>%
  select(-c(mu_Autumn, mu_Spring, mu_Summer,
            obs_Autumn, obs_Spring, obs_Summer)) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = c("Parameter", "Contrast"),
               values_to = "Difference",
               names_pattern = "^([^_]+)_(.+)$")

# 3.5.6 Summarise difference ####
consumption_diff_summary <- consumption_diff %>%
  group_by(Parameter, Contrast) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P = mean(Difference > 0),
            n = length(Difference)) %T>%
  print()

# 3.5.7 Add labels to grazing_diff ####
consumption_diff %<>%
  left_join(consumption_diff_summary %>%
              select(Parameter, Contrast, P), 
            by = c("Parameter", "Contrast")) %>%
  mutate(label_right = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_left = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))

# 4. Defecation model ####
# 4.1 Prior simulation ####
# Sauchyn & Scheibling 2009 (doi: 10.3354/meps08296) estimated defecation the
# same way as me (dry mass per urchin wet mass per day) and report around 2.3
# mg g^-1 d^-1 in the lab. Since my proportional defecation compares well with
# this study, this seems a reasonable central tendency to go with.

tibble(n = 1:1e5,
       alpha_mu = rnorm( 1e5 , log(2.3) , 0.5 ),
       alpha_sigma = rexp( 1e5 , 2 ),
       alpha = rnorm( 1e5 , alpha_mu , alpha_sigma ),
       theta = rexp( 1e5 , 2 ),
       mu = exp(alpha),
       y_new = rgamma( 1e5 , mu / theta , 1 / theta )) %>%
  pivot_longer(cols = c(mu, y_new), 
               names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, parameter)) +
    geom_vline(xintercept = grazing %$% range(Defecation)) +
    geom_density_ridges(alpha = 0.5, from = 0, to = 5,
                        colour = NA) +
    scale_x_continuous(limits = c(0, 5), oob = scales::oob_keep) +
    coord_cartesian(expand = FALSE) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks good.

# 4.2 Stan model ####
# Because multilevel modelling here works with normal hierarchical priors,
# there are two possible parameterisations: centred and non-centred.

defecation_c_model <- here("Urchins", "Stan", "defecation_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

defecation_nc_model <- here("Urchins", "Stan", "defecation_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

defecation_c_samples <- defecation_c_model$sample(
          data = grazing %>% 
            select(Season, Defecation) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )

defecation_nc_samples <- defecation_nc_model$sample(
          data = grazing %>% 
            select(Season, Defecation) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )

# 4.3 Model checks ####
# Rhat
defecation_c_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000517.

defecation_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000236.

defecation_c_samples$summary() %>%
  left_join(defecation_nc_samples$summary(),
            by = "variable") %>%
  rename(rhat_c = rhat.x, rhat_nc = rhat.y) %>%
  ggplot(aes(rhat_c, rhat_nc)) +
  geom_abline(slope = 1) +
  geom_point() +
  theme_minimal() +
  theme(panel.grid = element_blank())
# Models are pretty much identical, but I lean towards
# the centred version.

# Chains
defecation_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are great.

defecation_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are good.

# Pairs
defecation_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_sigma", "theta"))

defecation_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_sigma", "theta"))
# No correlations. Models look similar.

# 4.4 Prior-posterior comparison ####
defecation_c_prior <- prior_samples(
  model = defecation_c_model,
  data = grazing %>% 
    select(Season, Defecation) %>%
    compose_data(),
  adapt_delta = 0.99 # force sampler to go slower for smooth priors
  )

defecation_nc_prior <- prior_samples(
  model = defecation_nc_model,
  data = grazing %>% 
    select(Season, Defecation) %>%
    compose_data()
)

defecation_c_prior %>% 
  prior_posterior_draws(
    posterior_samples = defecation_c_samples,
    group = grazing %>% select(Season),
    parameters = c("alpha_mu", "alpha_sigma", "alpha[Season]", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Season", ridges = FALSE)

defecation_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = defecation_nc_samples,
    group = grazing %>% select(Season),
    parameters = c("alpha_mu", "alpha_sigma", "alpha_z[Season]",
                   "alpha[Season]", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Season", ridges = FALSE)
# Posteriors are not overly constrained by priors. Models
# are similar. I'll pick the centred model based on better 
# rhat and chains.

# 4.5 Prediction ####
# 4.5.1 Priors and posteriors for hyperparameters ####
defecation_prior_posterior_hyper <- defecation_c_prior %>% 
  prior_posterior_draws(
    posterior_samples = defecation_c_samples,
    parameters = c("alpha_mu", "alpha_sigma", "theta"),
    format = "short"
  ) %>% # Calculate predictions for new seasons, i.e. annual.
  mutate(alpha = rnorm( n() , alpha_mu , alpha_sigma ))

# 4.5.2 Priors and posteriors for seasonal parameters ####
defecation_prior_posterior_season <- defecation_c_prior %>% 
  prior_posterior_draws(
    posterior_samples = defecation_c_samples,
    parameters = c("alpha[Season]", "theta"),
    format = "short"
  ) %>% 
  # Since I want only one grouping variable, there is redundancy in distribution.
  filter(!(Season %in% c("Summer", "Spring") & 
             distribution == "prior")) %>% # Remove two redundant priors.
  mutate(Season = if_else(distribution == "prior", # Add Prior to Season.
                          "Prior", Season) %>% fct()) %>%
  select(-distribution)

# 4.5.3 Combine seasonal and hyperparameters ####
defecation_prior_posterior <- defecation_prior_posterior_hyper %>%
  rename(Season = distribution) %>%
  mutate(Season = if_else(Season == "posterior",
                          "Annual", "Annual prior") %>% fct()) %>%
  select(-c(alpha_mu, alpha_sigma)) %>%
  bind_rows(defecation_prior_posterior_season) 

defecation_prior_posterior %>%
  filter(Season %in% c("Prior", "Annual prior")) %>%
  ggplot(aes(alpha, Season)) +
    geom_density_ridges()
# Priors are identical because the same beast is simulated in 
# defecation_prior_posterior_hyper and sampled in defecation_prior_posterior_season.
# Annual prior looks smoother than Prior, so I'll remove the latter.

defecation_prior_posterior %<>%
  filter(Season != "Prior") %>% 
  mutate(Season = Season %>% 
           fct_drop() %>% 
           fct_recode(Prior = "Annual prior"))

# 4.5.4 Predict new observations ####
defecation_prior_posterior %<>%
  mutate(mu = exp(alpha),
         obs = rgamma( n() , mu / theta , 1 / theta ))

# 4.5.5 Calculate seasonal differences ####
defecation_diff <- defecation_prior_posterior %>%
  filter(!Season %in% c("Prior", "Annual")) %>%
  droplevels() %>%
  select(-c(alpha, theta)) %>%
  pivot_wider(names_from = Season, values_from = c(mu, obs)) %>%
  mutate(mu_Spring_Summer = mu_Spring - mu_Summer, # calculate differences
         mu_Spring_Autumn = mu_Spring - mu_Autumn,
         mu_Autumn_Summer = mu_Autumn - mu_Summer,
         obs_Spring_Summer = obs_Spring - obs_Summer,
         obs_Spring_Autumn = obs_Spring - obs_Autumn,
         obs_Autumn_Summer = obs_Autumn - obs_Summer) %>%
  select(-c(mu_Autumn, mu_Spring, mu_Summer,
            obs_Autumn, obs_Spring, obs_Summer)) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = c("Parameter", "Contrast"),
               values_to = "Difference",
               names_pattern = "^([^_]+)_(.+)$")

# 4.5.6 Summarise difference ####
defecation_diff_summary <- defecation_diff %>%
  group_by(Parameter, Contrast) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P = mean(Difference > 0),
            n = length(Difference)) %T>%
  print()

# 4.5.7 Add labels to grazing_diff ####
defecation_diff %<>%
  left_join(defecation_diff_summary %>%
              select(Parameter, Contrast, P), 
            by = c("Parameter", "Contrast")) %>%
  mutate(label_right = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_left = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))


# 5. Food conversion model ####
# 6. Relative consumption model ####
# 7. Relative growth model ####

# 8. Visualisation ####
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
    # manually limit 1:1 line to range (geom_textabline isn't clipped)
    geom_textline(data = tibble(x = c(0, 5), y = c(0, 5)), aes(x, y),
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_point(data = grazing %>%
                 mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
               aes(Consumption, Defecation, shape = Season),
               size = 2, alpha = 0.5, colour = "#7030a5") +
    geom_line(data = grazing_prediction_summary %>%
                filter(!Season %in% c("Annual", "Prior")),
              aes(Consumption, mu, colour = Season)) +
    geom_ribbon(data = grazing_prediction_summary %>%
                  filter(!Season %in% c("Annual", "Prior")),
                aes(Consumption, ymin = mu.lower, ymax = mu.upper,
                    fill = Season, alpha = factor(.width))) +
    geom_text(data = grazing_prediction_summary %>%
                filter(!Season %in% c("Annual", "Prior")) %>% 
                group_by(Season) %>% 
                summarise(x = max(Consumption), y = max(mu)),
              aes(x = x, y = y, label = Season),
              hjust = 0, nudge_x = 0.1, family = "Futura", 
              size = 4.2, colour = "#7030a5") +
    geom_density(data = consumption_prior_posterior %>%
                   filter(!Season %in% c("Annual", "Prior")),
                 aes(x = mu, y = after_stat(density) * 0.5, fill = Season),
                 alpha = 0.6, colour = NA) +
    geom_density(data = defecation_prior_posterior %>%
                   filter(!Season %in% c("Annual", "Prior")),
                 aes(x = after_stat(density) * 0.35, y = mu, fill = Season),
                 alpha = 0.6, colour = NA) +
    scale_shape_manual(values = c(16, 17, 15)) +
    # geom_ribbon doesn't allow grouping by Season and .width, so I had to hack
    scale_fill_manual(values = rep("#7030a5", 3), guide = "none") +
    scale_colour_manual(values = rep("#7030a5", 3), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(breaks = seq(0, 10, 2)) +
    labs(x = expression("Consumption (mg g"^-1*" d"^-1*")"),
         y = expression("Defecation (mg g"^-1*" d"^-1*")")) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 5),
                    expand = FALSE, clip = "off") +
    mytheme
Fig_1a_left

# Add defecation (%) to data.
grazing %<>% mutate(Defecation_percent = Defecation / Consumption * 100)

Fig_1a_right <- ggplot() +
  geom_jitter(data = grazing %>%
                mutate(Season = Season %>% fct_relevel("Autumn", "Summer")),
              aes(x = Defecation_percent, y = Season %>% as.numeric() + 0.5,
                  shape = Season),
              colour = "#7030a5", alpha = 0.5, size = 2, height = 0.35) +
  stat_density_ridges(data = grazing_prior_posterior %>%
                        mutate(Season = Season %>% 
                                 fct_relevel("Annual", "Autumn", "Summer", "Spring")),
                      aes(x = beta, y = Season %>% as.numeric(), fill = Season), 
                      colour = NA, n = 2^10, rel_min_height = 0.001, bandwidth = 0.8, 
                      scale = 3, alpha = 0.6, from = 0, to = 100) +
  annotate("text", x = 0, y = 1:5, hjust = 0, vjust = -0.8, size = 4.2,
           label = c("Annual", "Autumn", "Summer", "Spring", "Prior"),
           colour = c("#363538", rep("#7030a5", 3), "#b5b8ba"), 
           family = "Futura") +
  scale_shape_manual(values = c(15, 17, 16),
                     guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c("#363538", rep("#7030a5", 3), "#b5b8ba"),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
  xlab("Defecation (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())
Fig_1a_right

# 9. Save relevant data ####
grazing %>%
  write_rds(here("Urchins", "RDS", "grazing.rds"))
grazing_prior_posterior %>% 
  write_rds(here("Urchins", "RDS", "grazing_prior_posterior.rds"))
grazing_prediction_summary %>%
  write_rds(here("Urchins", "RDS", "grazing_prediction_summary.rds"))
consumption_prior_posterior %>% 
  write_rds(here("Urchins", "RDS", "consumption_prior_posterior.rds"))
defecation_prior_posterior %>% 
  write_rds(here("Urchins", "RDS", "defecation_prior_posterior.rds"))
grazing_diff %>% 
  write_rds(here("Urchins", "RDS", "grazing_diff.rds"))
consumption_diff %>% 
  write_rds(here("Urchins", "RDS", "consumption_diff.rds"))
defecation_diff %>% 
  write_rds(here("Urchins", "RDS", "defecation_diff.rds"))