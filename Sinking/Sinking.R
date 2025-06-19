# 1. Sinking speed ####
# 1.1 Load data ####
require(tidyverse)
require(here)

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
  drop_na(Tissue, Speed)
# Warning because some NAs were specified as characters (e.g. "no good").
speed %>%
  print(n = 346)
# Some observations are missing area, some mass, and some both.

# 1.2 Prior simulation ####
require(magrittr)
speed %$% range(Speed)

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
# No rhat above 1.001. rhat = 1.00 ± 0.0000162.

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
  
# Calculate difference
speed_diff <- speed_prior_posterior %>%
  filter(Tissue != "Prior") %>%
  droplevels() %>%
  select(-c(alpha, theta)) %>%
  pivot_wider(names_from = Tissue, values_from = c(mu, obs)) %>%
  mutate(mu = mu_Kelp - mu_Faeces, # calculate differences
         obs = obs_Kelp - obs_Faeces) %>%
  select(.chain, .iteration, .draw, mu, obs) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = "Parameter",
               values_to = "Difference")

# Summarise difference
speed_diff_summary <- speed_diff %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P = mean(Difference > 0),
            n = length(Difference)) %T>%
  print()

# Add labels to speed_diff
speed_diff %<>%
  left_join(speed_diff_summary %>%
              select(Parameter, P), 
            by = "Parameter") %>%
  mutate(label_Kelp = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_Faeces = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))

# 2. Distance ####
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
  
# Calculate difference
distance_diff <- distance_prior_posterior %>%
  filter(Tissue != "Prior") %>%
  droplevels() %>%
  select(-c(alpha, theta)) %>%
  pivot_wider(names_from = Tissue, values_from = c(mu, obs)) %>%
  mutate(mu = mu_Kelp - mu_Faeces, # calculate differences
         obs = obs_Kelp - obs_Faeces) %>%
  select(.chain, .iteration, .draw, mu, obs) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = "Parameter",
               values_to = "Difference")

# Summarise difference
distance_diff_summary <- distance_diff %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P = mean(Difference > 0),
            n = length(Difference)) %T>%
  print()

# Add labels to speed_diff
distance_diff %<>%
  left_join(distance_diff_summary %>%
              select(Parameter, P), 
            by = "Parameter") %>%
  mutate(label_Kelp = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_Faeces = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))



# 3. Depth ####
# 3.2 Prior simulation ####
# 3.3 Stan model ####
# 3.4 Model checks ####
# 3.5 Prior-posterior comparison ####
# 3.6 Prediction ####

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
              alpha = 0.5, size = 2, height = 0.4, shape = 16) +
  stat_density_ridges(data = speed_prior_posterior %>%
                        mutate(Tissue = Tissue %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs, y = Tissue %>% as.numeric(), fill = Tissue), 
                      colour = NA, n = 2^10,
                      from = 0, to = 25, rel_min_height = 0.001, 
                      bandwidth = 0.1, scale = 1.2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 25), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                      guide = "none") +
  xlab(expression("Sinking speed (cm s"^-1*")")) +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme
Fig_S5a_top

require(geomtextpath)
Fig_S5a_bottom <- speed_diff %>% 
  filter(Parameter == "obs") %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.2,
                      from = -16, to = 16, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.7, vjust = 0,
                   n = 2^10, bw = 0.2, text_only = TRUE) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.3, vjust = 0,
                   n = 2^10, bw = 0.2, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  # annotate("text", x = -16, y = 0, 
  #          label = "italic(tilde('y'))",
  #          hjust = 0, vjust = 0, 
  #          family = "Futura", size = 3.5,
  #          parse = TRUE) +
  scale_x_continuous(limits = c(-16, 16), oob = scales::oob_keep,
                     breaks = seq(-16, 16, 8),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#dabc23", 0.6)),
                    guide = "none") +
  xlab(expression("Δ sinking speed (cm s"^-1*")")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme
Fig_S5a_bottom

# 4.3 Distance ####
Fig_S5b_top <- ggplot() +
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
Fig_S5b_top

Fig_S5b_bottom <- distance_diff %>% 
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
Fig_S5b_bottom

# 4.4 Depth ####

# 5. Save relevant data ####
speed %>%
  write_rds(here("Sinking", "RDS", "speed.rds"))
speed_prior_posterior %>% 
  write_rds(here("Sinking", "RDS", "speed_prior_posterior.rds"))
speed_diff %>% 
  write_rds(here("Sinking", "RDS", "speed_diff.rds"))
distance %>%
  write_rds(here("Sinking", "RDS", "distance.rds"))
distance_prior_posterior %>% 
  write_rds(here("Sinking", "RDS", "distance_prior_posterior.rds"))
distance_diff %>% 
  write_rds(here("Sinking", "RDS", "distance_diff.rds"))