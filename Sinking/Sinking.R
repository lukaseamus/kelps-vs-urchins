# 1. Sinking speed ####
# 1.1 Load data ####
require(tidyverse)
require(here)
require(magrittr)
speed <- here("Sinking", "Speed.csv") %>% read_csv() %>%
  mutate(Tissue = case_when(
    Type2 %>% str_detect("Blade") ~ "Kelp",
    Type2 %>% str_detect("Poo") ~ "Faeces",
    TRUE ~ NA
  ) %>% fct(),
  # Calculate sinking speed in cm s^-1
  Speed = `drop length (m)` / `Sinking time (s)` %>% as.numeric() * 100
  ) %>%
  rename(Mass = "WW (g)", Area = "Area (mm2)") %>%
  select(Tissue, Speed, Mass, Area) %>%
  drop_na(Tissue, Speed) %T>%
  print(n = 346)
# Warning because some NAs were specified as characters (e.g. "no good").
# Some observations are missing area, some mass, and some both.

speed %>%
  count(Tissue)
# Way fewer faecal than kelp measurments

# Add meta-analysis data
set.seed(100) # Set seed for simulation of observations from mean ± s.e.m.
meta_speed <- here("Sinking", "Meta_speed.csv") %>% read_csv() %>%
  rowwise() %>% # Simulate from a gamma distribution to avoid negative values
  mutate(Speed = list( rgamma( N , 
                               Mean^2 / ( SEM * sqrt(N) )^2 , 
                               Mean / ( SEM * sqrt(N) )^2 ) )) %>%
  unnest(Speed) # Speed is already give in cm s^-1

speed %<>%
  bind_rows(
    meta_speed %>% 
      select(Speed) %>%
      mutate(Tissue = "Faeces" %>% fct(), 
             Mass = NA, Area = NA)
  ) %T>%
  print(n = 520)

speed %>%
  count(Tissue)
# Much more equal

# 1.2 Prior simulation ####
speed %$% range(Speed)
# Since this is a right-skewed distribution I am picking a 
# mean of around 10 cm s^-1.
require(ggridges)
tibble(n = 1:1e5,
       alpha = rnorm( 1e5 , log(10) , 1 ),
       theta = rexp( 1e5 , 1 ),
       mu = exp( alpha ),
       y_new = rgamma( 1e5 , mu / theta , 1 / theta )) %>%
  pivot_longer(cols = c(mu, y_new), 
               names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, parameter)) +
    geom_density_ridges(alpha = 0.5, colour = NA,
                        from = 0, to = 25) +
    scale_x_continuous(limits = c(0, 25), oob = scales::oob_keep) +
    coord_cartesian(expand = FALSE) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 1.3 Stan model ####
require(cmdstanr)
speed_model <- here("Sinking", "Stan", "speed.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
speed_samples <- speed_model$sample(
  data = speed %>%
    select(Tissue, Speed) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
)

# 1.4 Model checks ####
# Rhat
speed_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000448.

# Chains
require(bayesplot)
speed_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are well-mixed.

# Pairs
speed_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1]", "theta[1]"))

speed_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[2]", "theta[2]"))
# No correlation.

# 1.5 Prior-posterior comparison ####
source("functions.R")
speed_prior <- prior_samples(
  model = speed_model,
  data = speed %>%
    select(Tissue, Speed) %>%
    compose_data()
)

speed_prior %>% 
  prior_posterior_draws(
    posterior_samples = speed_samples,
    group = speed %>% select(Tissue),
    parameters = c("alpha[Tissue]", "theta[Tissue]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Tissue", ridges = FALSE)

# 1.6 Prediction ####
# Extract priors and posteriors
speed_prior_posterior <- speed_prior %>% 
  prior_posterior_draws(
    posterior_samples = speed_samples,
    group = speed %>% select(Tissue),
    parameters = c("alpha[Tissue]", "theta[Tissue]"),
    format = "short"
  )

# Predict mu and new observations
speed_prior_posterior %<>%
  mutate(
    mu = exp( alpha ),
    obs = rgamma( n() , mu / theta , 1 / theta )
  )

# Wrangle Tissue and distribution into one variable
speed_prior_posterior %<>% # priors are identical for both treatments ->
  filter(!(Tissue == "Faeces" & distribution == "prior")) %>% # remove one
  mutate(Tissue = if_else(distribution == "prior", # add Prior to Tissue
                          "Prior", Tissue) %>% fct()) %>%
  select(-distribution)

# Calculate differences and proportions
speed_diff <- speed_prior_posterior %>%
  filter(Tissue != "Prior") %>%
  droplevels() %>%
  select(-c(alpha, theta)) %>%
  pivot_wider(names_from = Tissue, values_from = c(mu, obs)) %>%
  mutate(mu_diff = mu_Kelp - mu_Faeces, # Calculate absolute differences
         obs_diff = obs_Kelp - obs_Faeces,
         mu_reldiff = (mu_Kelp - mu_Faeces) / mu_Kelp, # Calculate relative differences
         obs_reldiff = (obs_Kelp - obs_Faeces) / obs_Kelp,
         mu_logratio = log10(mu_Faeces / mu_Kelp), # Calculate log10 ratio (order of magnitude)
         obs_logratio = log10(obs_Faeces / obs_Kelp)) %>%
  select(starts_with("."), ends_with("diff"), 
         ends_with("reldiff"), ends_with("logratio")) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = c("parameter", "difference"),
               values_to = "value",
               names_sep = "_") %>%
  pivot_wider(values_from = value, names_from = difference)

# Summarise differences
speed_diff_summary <- speed_diff %>%
  group_by(parameter) %>%
  summarise(mean = mean(diff),
            sd = sd(diff),
            P = mean(diff > 0),
            n = n()) %T>%
  print()

speed_diff %>%
  group_by(parameter) %>%
  summarise(mean = mean(reldiff),
            sd = sd(reldiff),
            P = mean(reldiff > 0),
            n = n())

speed_diff %>%
  group_by(parameter) %>%
  summarise(mean = mean(logratio),
            sd = sd(logratio),
            P = mean(logratio < 0),
            n = n())

# Add labels to speed_diff
speed_diff %<>%
  left_join(speed_diff_summary %>%
              select(parameter, P), 
            by = "parameter") %>%
  mutate(label_Kelp = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_Faeces = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))

# 2. Distance ####
# This can be skipped since I am not using distance in my estimation of 
# carbon sequestration potential.
# 2.1 Load data ####
distance <- here("Sinking", "Distance.csv") %>% read_csv() %>%
  mutate(Tissue = if_else(particle == "feces", 
                          "Faeces", "Kelp") %>% fct()) %>%
  rename(Distance = cumulativedistancekm) %>%
  select(Tissue, Distance)

# 2.2 Prior simulation ####
distance %$% range(Distance)

tibble(n = 1:1e5,
       alpha = rnorm( 1e5 , log(100) , 1 ),
       theta = rexp( 1e5 , 1 ),
       mu = exp( alpha ),
       y_new = rgamma( 1e5 , mu / theta , 1 / theta )) %>%
  pivot_longer(cols = c(mu, y_new), 
               names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, parameter)) +
    geom_density_ridges(alpha = 0.5, colour = NA,
                        from = 0, to = 356) +
    scale_x_continuous(limits = c(0, 356), oob = scales::oob_keep) +
    coord_cartesian(expand = FALSE) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.3 Stan model ####
distance_model <- here("Sinking", "Stan", "distance.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

distance_samples <- distance_model$sample(
  data = distance %>%
    select(Tissue, Distance) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
)

# 2.4 Model checks ####
# Rhat
distance_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000103.

# Chains
distance_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are well-mixed. Structure in log probability likely due 
# to strong data structure.

# Pairs
distance_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1]", "theta[1]"))

distance_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[2]", "theta[2]"))
# Weak correlation. Not a problem.

# 2.5 Prior-posterior comparison ####
distance_prior <- prior_samples(
  model = distance_model,
  data = distance %>%
    select(Tissue, Distance) %>%
    compose_data()
)

distance_prior %>% 
  prior_posterior_draws(
    posterior_samples = distance_samples,
    group = distance %>% select(Tissue),
    parameters = c("alpha[Tissue]", "theta[Tissue]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Tissue", ridges = FALSE)
# Standard exponential(1) prior on theta is very tight but data are 
# strong enough to increase likelihood uncertainty, so not a problem.

# 2.6 Prediction ####
# Extract priors and posteriors
distance_prior_posterior <- distance_prior %>% 
  prior_posterior_draws(
    posterior_samples = distance_samples,
    group = distance %>% select(Tissue),
    parameters = c("alpha[Tissue]", "theta[Tissue]"),
    format = "short"
  )

# Predict mu and new observations
distance_prior_posterior %<>%
  mutate(
    mu = exp( alpha ),
    obs = rgamma( n() , mu / theta , 1 / theta )
  )

# Wrangle Tissue and distribution into one variable
distance_prior_posterior %<>% # priors are identical for both treatments ->
  filter(!(Tissue == "Faeces" & distribution == "prior")) %>% # remove one
  mutate(Tissue = if_else(distribution == "prior", # add Prior to Tissue
                          "Prior", Tissue) %>% fct()) %>%
  select(-distribution)
  
# Calculate differences and proportions
distance_diff <- distance_prior_posterior %>%
  filter(Tissue != "Prior") %>%
  droplevels() %>%
  select(-c(alpha, theta)) %>%
  pivot_wider(names_from = Tissue, values_from = c(mu, obs)) %>%
  mutate(mu_diff = mu_Kelp - mu_Faeces, # Calculate absolute differences
         obs_diff = obs_Kelp - obs_Faeces,
         mu_reldiff = (mu_Kelp - mu_Faeces) / mu_Kelp, # Calculate relative differences
         obs_reldiff = (obs_Kelp - obs_Faeces) / obs_Kelp,
         mu_prop = mu_Faeces / mu_Kelp, # Calculate proportions
         obs_prop = obs_Faeces / obs_Kelp) %>%
  select(starts_with("."), ends_with("diff"), ends_with("reldiff"), ends_with("prop")) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = c("parameter", "difference"),
               values_to = "value",
               names_sep = "_") %>%
  pivot_wider(values_from = value, names_from = difference)

# Summarise difference
distance_diff_summary <- distance_diff %>%
  group_by(parameter) %>%
  summarise(mean = mean(diff),
            sd = sd(diff),
            P = mean(diff > 0),
            n = n()) %T>%
  print()

distance_diff %>%
  group_by(parameter) %>%
  summarise(mean = mean(reldiff),
            sd = sd(reldiff),
            P = mean(reldiff > 0),
            n = n())

distance_diff %>%
  group_by(parameter) %>%
  summarise(mean = mean(prop),
            sd = sd(prop),
            n = n())

# Add labels to speed_diff
distance_diff %<>%
  left_join(distance_diff_summary %>%
              select(parameter, P), 
            by = "parameter") %>%
  mutate(label_Kelp = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_Faeces = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%")) %T>%
  print()

# 3. Depth ####
# 3.1 Load data ####
depth <- here("Sinking", "Depth.csv") %>% read_csv() %>%
  mutate(Speed = Speed * 100, # Convert m s^-1 to cm s^-1
         Speed_c = Speed - mean(Speed))

depth %>%
  ggplot(aes(Depth, Speed %>% factor())) +
    geom_density_ridges(from = 0, to = 500) +
    scale_x_continuous(limits = c(0, 500),
                       oob = scales::oob_keep) +
    theme_minimal()

# 3.2 Prior simulation ####
# Here I want a predictive model of final depth given sinking
# speed. Since depth cannot be negative and it is extremely
# right-skewed, a gamma likelihood is best.
depth %$% median(Depth)

tibble(n = 1:1e3,
       alpha = rnorm( 1e3 , log(50) , 1 ), # Average depth around 50 m
       beta = rnorm( 1e3 , -1 , 0.5 ), # Slight indication of expected negative slope
       theta = rexp( 1e3 , 1 )) %>%
  expand_grid(Speed_c = depth %$% 
                seq(min(Speed_c), max(Speed_c), length.out = 50)) %>%
  mutate(mu = exp( alpha + beta * Speed_c ), # Exponential link function
         Depth = rgamma( n() , mu / theta , 1 / theta )) %>%
  ggplot(aes(Speed_c, Depth, group = n)) +
    geom_hline(yintercept = depth %$% 
                 range(Depth)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(ylim = c(0, 3e3), expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks fine

# 3.3 Stan model ####
speed2depth_model <- here("Sinking", "Stan", "speed2depth.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

# Do not attempt to run this if you have a feeble computer!
speed2depth_samples <- speed2depth_model$sample(
  data = depth %>%
    select(Depth, Speed_c) %>%
    group_by(Speed_c) %>%
    slice_sample(n = 1e4) %>% # sample 1e4 observations from each speed
    compose_data(), # so there are 4e4 observations in total
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
)
# Takes ages, even with only slightly more than 1% of the data
speed2depth_samples$summary()

# n is huge (3.5 million), so here is a rare case where maximum likelihood 
# estimation can work better. I'll simply use the built-in glm function:
speed2depth_mle <- glm(Depth ~ Speed, family = Gamma(link = "log"),
                       data = depth) # Easy breezy!
summary(speed2depth_mle)

# 3.4 Model checks ####
depth %>%
  mutate(resid = resid(speed2depth_mle)) %>%
  ggplot(aes(resid, Speed %>% factor())) +
    geom_density_ridges() +
    theme_minimal()
# Residuals are fairly homogenous but not normal. Can't really
# do any better with maximum likelihood estimation in this case.

# 3.5 Prediction ####
# 3.5.1 Simulate posteriors ####
speed2depth_posterior <- summary(speed2depth_mle)$coefficients %>%
  as_tibble() %>%
  reframe(alpha = rnorm( 8e4 , Estimate[1] , `Std. Error`[1] ),
          beta = rnorm( 8e4 , Estimate[2] , `Std. Error`[2] ),
          phi = summary(speed2depth_mle)$dispersion) %T>%
  print()
  
speed2depth_prediction <- speed2depth_posterior %>%
  expand_grid(Speed = c(seq(0, 5, length.out = 1e3), 6)) %>%
  mutate(mu = exp( alpha + beta * Speed ),
         Depth = rgamma( n() , 1 / phi , 1 / (mu * phi) )) %T>%
  print()

speed2depth_prediction_summary <- speed2depth_prediction %>%
  group_by(Speed) %>%
  mean_qi(mu, Depth, .width = c(.5, .8, .9)) %T>%
  print()

# 3.5.1 Predict depth from speed ####
depth_posterior <- speed_prior_posterior %>%
  filter(Tissue != "Prior") %>%
  droplevels() %>%
  rename(Speed = obs) %>%
  select(starts_with("."), Tissue, Speed) %>%
  full_join(speed2depth_posterior %>%
              mutate(.draw = 1:8e4),
            by = ".draw") %>%
  mutate(mu = exp( alpha + beta * Speed ),
         Depth = rgamma( n() , 1 / phi , 1 / (mu * phi) )) %T>%
  print()

# Calculate differences and proportions
depth_diff <- depth_posterior %>%
  rename(obs = Depth) %>%
  select(starts_with("."), Tissue, mu, obs) %>%
  pivot_wider(names_from = Tissue, values_from = c(mu, obs)) %>%
  mutate(mu_diff = mu_Kelp - mu_Faeces, # Calculate absolute differences
         obs_diff = obs_Kelp - obs_Faeces,
         mu_reldiff = (mu_Kelp - mu_Faeces) / mu_Kelp, # Calculate relative differences
         obs_reldiff = (obs_Kelp - obs_Faeces) / obs_Kelp,
         mu_logratio = log10(mu_Faeces / mu_Kelp), # Calculate log10 ratio (order of magnitude)
         obs_logratio = log10(obs_Faeces / obs_Kelp)) %>%
  select(starts_with("."), ends_with("diff"), 
         ends_with("reldiff"), ends_with("logratio")) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = c("parameter", "difference"),
               values_to = "value",
               names_sep = "_") %>%
  pivot_wider(values_from = value, names_from = difference) %T>%
  print()

# Summarise differences
depth_diff_summary <- depth_diff %>%
  group_by(parameter) %>%
  summarise(mean = mean(diff),
            sd = sd(diff),
            P = mean(diff > 0),
            n = n()) %T>%
  print()

depth_diff %>%
  group_by(parameter) %>%
  summarise(mean = mean(reldiff),
            sd = sd(reldiff),
            P = mean(reldiff > 0),
            n = n())
# Summary in terms of mean and sd is inadequate due to heavy tails.

depth_diff %>%
  group_by(parameter) %>%
  summarise(mean = mean(logratio),
            sd = sd(logratio),
            P = mean(logratio < 0),
            n = n())
# Summary in terms of mean and sd is inadequate due to heavy tails.

# Add labels to speed_diff
depth_diff %<>%
  left_join(depth_diff_summary %>%
              select(parameter, P), 
            by = "parameter") %>%
  mutate(label_Kelp = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_Faeces = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))

# 4. Visualisation ####
# 4.1 Define custom theme ####
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

# 4.2 Speed ####
Fig_S5a_top <- ggplot() +
  geom_jitter(data = speed %>%
                mutate(Tissue = Tissue %>% fct_relevel("Faeces", "Kelp")),
              aes(x = Speed, y = Tissue %>% as.numeric() - 0.5, 
                  colour = Tissue), 
              alpha = 0.3, size = 1, height = 0.4, shape = 16) +
  stat_density_ridges(data = speed_prior_posterior %>%
                        mutate(Tissue = Tissue %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs, y = Tissue %>% as.numeric(), fill = Tissue), 
                      colour = NA, n = 2^10,
                      from = 0, to = 25, rel_min_height = 0.001, 
                      bandwidth = 0.25, scale = 2, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 6),
                     oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                      guide = "none") +
  xlab(expression("Sinking speed (cm s"^-1*")")) +
  coord_cartesian(ylim = c(0, 3), expand = FALSE, clip = "off") +
  mytheme + # Adjust the title margin to counteract the superscript.
  theme(axis.title.x = element_text(margin = margin(b = -20)))
Fig_S5a_top

require(geomtextpath)
Fig_S5a_middle <- speed_diff %>% 
  filter(parameter == "obs") %>%
  ggplot() +
  stat_density_ridges(aes(x = diff, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.32,
                      from = -16, to = 16, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.75, vjust = 0,
                   n = 2^10, bw = 0.32, text_only = TRUE) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 0.32, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -16, y = 0, 
           label = "italic(tilde('y'))['Kelp']*' − '*italic(tilde('y'))['Faeces']",
           hjust = 0, vjust = -1.6, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-16, 16), oob = scales::oob_keep,
                     breaks = seq(-16, 16, 8),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.8), alpha("#dabc23", 0.8)),
                    guide = "none") +
  xlab(expression("Δ sinking speed (cm s"^-1*")")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme + # Adjust the title margin to counteract the superscript.
  theme(axis.title.x = element_text(margin = margin(b = -20)))
Fig_S5a_middle

Fig_S5a_bottom <- speed_diff %>% 
  filter(parameter == "obs") %>%
  ggplot() +
  stat_density_ridges(aes(x = logratio, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.06,
                      from = -3, to = 3, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = logratio, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.24, vjust = 0,
                   n = 2^10, bw = 0.06, text_only = TRUE) +
  geom_textdensity(aes(x = logratio, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.7, vjust = 0,
                   n = 2^10, bw = 0.06, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = 3, y = 0, 
           label = "log[10]*' '*frac(italic(tilde('y'))['Faeces'],italic(tilde('y'))['Kelp'])",
           hjust = 1, vjust = -0.1, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1.5),
                     oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus",
                                                   accuracy = c(1, 0.1, 1, 0.1, 1))) +
  scale_fill_manual(values = c(alpha("#dabc23", 0.8), alpha("#7030a5", 0.8)),
                    guide = "none") +
  xlab("Relative sinking speed") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme
Fig_S5a_bottom

# 4.3 Distance ####
# Distance is not included in the estimation so will not be included in the supplement.
ggplot() +
  geom_jitter(data = distance %>%
                mutate(Tissue = Tissue %>% fct_relevel("Faeces", "Kelp")),
              aes(x = Distance, y = Tissue %>% as.numeric() - 0.5, 
                  colour = Tissue), 
              alpha = 0.2, size = 2, height = 0.4, shape = 16) +
  stat_density_ridges(data = distance_prior_posterior %>%
                        mutate(Tissue = Tissue %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs, y = Tissue %>% as.numeric(), fill = Tissue), 
                      colour = NA, n = 2^10,
                      from = 0, to = 400, rel_min_height = 0.001, 
                      bandwidth = 2, scale = 1.2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 400), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                      guide = "none") +
  xlab("Export distance (km)") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

distance_diff %>% 
  filter(Parameter == "obs") %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 4,
                      from = -400, to = 400, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.75, vjust = 0,
                   n = 2^10, bw = 4, text_only = TRUE) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.2, vjust = 0,
                   n = 2^10, bw = 4, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-400, 400), oob = scales::oob_keep,
                     breaks = seq(-400, 400, 200),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#dabc23", 0.6)),
                    guide = "none") +
  xlab("Δ export distance (km)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 4.4 Depth ####
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

Fig_S5b <- ggplot() + 
    # geom_jitter(data = depth, # Plotting 3.5 million points makes the pdf nearly unreadable
    #            aes(Speed, Depth),
    #            size = 0.5, alpha = 0.002, 
    #            shape = 16, width = 0.08) +
    geom_violin(data = depth,
                aes(Speed, Depth, group = Speed),
                width = 0.6, fill = NA) +
    geom_line(data = speed2depth_prediction_summary,
              aes(Speed, mu)) +
    geom_ribbon(data = speed2depth_prediction_summary,
                aes(Speed, ymin = Depth.lower, ymax = Depth.upper,
                    alpha = factor(.width))) +
    geom_density(data = speed_prior_posterior %>%
                   filter(Tissue != "Prior"),
                 aes(x = obs, y = after_stat(density) * -400, fill = Tissue),
                 alpha = 0.8, colour = NA, bw = 0.08) +
    geom_density(data = depth_posterior,
                 aes(x = after_stat(density) * -250, y = Depth, fill = Tissue),
                 alpha = 0.8, colour = NA, bw = 30) +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_fill_manual(values = c("#dabc23", "#7030a5"), guide = "none") +
    scale_x_continuous(position = "top", limits = c(NA, 6), breaks = seq(0, 6, by = 2), 
                       oob = scales::oob_keep) +
    scale_y_reverse(limits = c(3000, NA), breaks = seq(3000, 0, by = -500),
                    oob = scales::oob_keep) +
    labs(x = expression("Sinking speed (cm s"^-1*")"),
         y = "Export depth (m)") +
    coord_cartesian(xlim = c(0, 6), ylim = c(3000, 0),
                    expand = FALSE, clip = "off") +
    mytheme +
    # The central plot controls the top position of all plots so needs to be
    # tweaked with text margins.
    theme(axis.title.x.top = element_text(hjust = 1, margin = margin(t = -20)),
          legend.justification = 1)
Fig_S5b

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

Fig_S5c_top <- ggplot() +
  stat_density_ridges(data = depth_posterior %>%
                        mutate(Tissue = Tissue %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = Depth, y = Tissue %>% as.numeric() - 1, fill = Tissue), 
                      colour = NA, n = 2^10,
                      from = 0, to = 1000, rel_min_height = 0.001, 
                      bandwidth = 10, scale = 3, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 1000), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23"),
                    guide = "none") +
  scale_colour_manual(values = c("#7030a5", "#dabc23"),
                      guide = "none") +
  xlab("Export depth (m)") +
  coord_cartesian(ylim = c(0, 3), expand = FALSE, clip = "off") +
  mytheme +
  theme(axis.title.x = element_text(margin = margin(t = 3, b = -20)))
Fig_S5c_top

Fig_S5c_middle <- depth_diff %>% 
  filter(parameter == "obs") %>%
  ggplot() +
  stat_density_ridges(aes(x = diff, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 20,
                      from = -1000, to = 1000, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.76, vjust = 0,
                   n = 2^10, bw = 20, text_only = TRUE) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.22, vjust = 0,
                   n = 2^10, bw = 20, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  # annotate("text", x = -1000, y = 0.002, 
  #          label = "italic(tilde('y'))['Kelp']*' − '*italic(tilde('y'))['Faeces']",
  #          hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
  #          parse = TRUE) +
  scale_x_continuous(limits = c(-1000, 1000), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus", 
                                                   big.mark = "")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.8), alpha("#dabc23", 0.8)),
                    guide = "none") +
  xlab("Δ export depth (m)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme +
  theme(axis.title.x = element_text(margin = margin(t = 3, b = -20)))
Fig_S5c_middle

Fig_S5c_bottom <- depth_diff %>% 
  filter(parameter == "obs") %>%
  ggplot() +
  stat_density_ridges(aes(x = logratio, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.2,
                      from = -10, to = 10, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = logratio, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.38, vjust = 0,
                   n = 2^10, bw = 0.2, text_only = TRUE) +
  geom_textdensity(aes(x = logratio, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 0.2, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  # annotate("text", x = -10, y = 0.01, 
  #          label = "log[10]*' '*frac(italic(tilde('y'))['Faeces'],italic(tilde('y'))['Kelp'])",
  #          hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
  #          parse = TRUE) +
  scale_x_continuous(limits = c(-10, 10), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#dabc23", 0.8), alpha("#7030a5", 0.8)),
                    guide = "none") +
  xlab("Relative export depth") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme
Fig_S5c_bottom

# 4.5 Figure S5 ####
require(patchwork)
layout <- c(
  area(1, 1), area(2, 1), area(3, 1),
  area(1, 2, 3, 2),
  area(1, 3), area(2, 3), area(3, 3)
)

Fig_S5 <- wrap_plots(
  Fig_S5a_top, Fig_S5a_middle, Fig_S5a_bottom,
  Fig_S5b,
  Fig_S5c_top, Fig_S5c_middle, Fig_S5c_bottom,
  design = layout
) +
  plot_layout(heights = c(1, 0.4, 0.4)) +
  plot_annotation(tag_levels = list(c("a", "", "", "b", "c", "", ""))) &
  theme(plot.tag = element_text(family = "Futura",
                                size = 15, face = "bold"))

Fig_S5 %>%
  ggsave(filename = "Fig_S5.pdf", path = "Figures",
         device = cairo_pdf, height = 14, width = 21, 
         units = "cm")

# 5. Save relevant data ####
speed %>%
  write_rds(here("Sinking", "RDS", "speed.rds"))
speed_prior_posterior %>% 
  write_rds(here("Sinking", "RDS", "speed_prior_posterior.rds"))
speed_diff %>% 
  write_rds(here("Sinking", "RDS", "speed_diff.rds"))
depth_posterior %>% 
  write_rds(here("Sinking", "RDS", "depth_posterior.rds"))
depth_diff %>% 
  write_rds(here("Sinking", "RDS", "depth_diff.rds"))