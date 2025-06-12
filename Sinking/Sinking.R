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
  Speed = `drop length (m)` / `Sinking time (s)` %>% as.numeric()
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

require(ggdist)
tibble(n = 1:1e5,
       mu_log = rnorm( 1e5 , log(0.1) , 1.2 ),
       theta = rexp( 1e5 , 100 ),
       mu = exp( mu_log ),
       P = rgamma( 1e5 , mu / theta , 1 / theta )) %>%
  pivot_longer(cols = c(mu, P), 
               names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, parameter)) +
    geom_vline(xintercept = c(0, 0.25)) +
    stat_slab(alpha = 0.5, height = 2, n = 3e3) +
    coord_cartesian(expand = F,
                    xlim = c(0, 0.25)) +
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
# No rhat above 1.001. rhat = 1.00 ± 0.0000407.

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

# Calculate mu and y_new
speed_prior_posterior %<>%
  mutate(
    mu = exp( alpha ),
    y_new = rgamma( n() , mu / theta , 1 / theta )
  )

# Wrangle Tissue and distribution into one variable
speed_prior_posterior %<>% # priors are identical for both treatments ->
  filter(!(Tissue == "Faeces" & distribution == "prior")) %>% # remove one
  mutate(Tissue = if_else(distribution == "prior", # add Prior to Tissue
                          "Prior", Tissue) %>% fct()) %>%
  select(-distribution)

# Visualise prediction
require(ggridges) # ggridges is better than ggdist for limiting distribution ranges
speed_prior_posterior %>%
  pivot_longer(cols = c(mu, y_new), 
               values_to = "Speed", names_to = "Level") %>%
  filter(Level %in% c("mu", "y_new")) %>%
  ggplot(aes(Speed, Tissue, alpha = Level)) +
  geom_density_ridges(from = 0, to = 0.25) +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  scale_x_continuous(limits = c(0, 0.25), oob = scales::oob_keep) +
  theme_minimal() +
  theme(panel.grid = element_blank())
  
# Calculate difference
speed_diff <- speed_prior_posterior %>%
  filter(Tissue != "Prior") %>%
  droplevels() %>%
  select(-c(alpha, theta)) %>%
  pivot_wider(names_from = Tissue, values_from = c(mu, y_new)) %>%
  mutate(mu = mu_Kelp - mu_Faeces, # calculate differences
         y_new = y_new_Kelp - y_new_Faeces) %>%
  select(.chain, .iteration, .draw, mu, y_new) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = "Parameter",
               values_to = "Difference")

# Visualise difference
speed_diff %>%
  ggplot(aes(Difference, Parameter)) +
    geom_density_ridges(from = -0.1, to = 0.1) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-0.1, 0.1), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

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
  rename(Distance = Cumdistance) %>%
  select(Tissue, Distance)

distance %>%
  ggplot(aes(Distance/2e3, y = Tissue)) +
    geom_density_ridges(bandwidth = 10, from = 0)

# 2.2 Prior simulation ####
# 2.3 Stan model ####
# 2.4 Model checks ####
# 2.5 Prior-posterior comparison ####
# 2.6 Prediction ####



# 3. Depth ####
# 3.2 Prior simulation ####
# 3.3 Stan model ####
# 3.4 Model checks ####
# 3.5 Prior-posterior comparison ####
# 3.6 Prediction ####

# 4. Supplementary figure ####
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

Fig_S5a_top <- ggplot() +
  geom_jitter(data = speed %>%
                mutate(Tissue = Tissue %>% fct_relevel("Faeces", "Kelp")),
              aes(x = Speed, y = Tissue %>% as.numeric() - 0.5, 
                  colour = Tissue), 
              alpha = 0.5, size = 2, height = 0.4) +
  stat_density_ridges(data = speed_prior_posterior %>%
                        mutate(Tissue = Tissue %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = y_new, y = Tissue %>% as.numeric(), fill = Tissue), 
                      colour = NA, n = 2^10,
                      from = 0, to = 0.25, rel_min_height = 0.001, 
                      bandwidth = 0.001, scale = 1.2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 0.25), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                      guide = "none") +
  xlab(expression("Sinking speed (m s"^-1*")")) +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme
Fig_S5a_top

require(geomtextpath)
Fig_S5a_bottom <- speed_diff %>% 
  filter(Parameter == "y_new") %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.002,
                      from = -0.16, to = 0.16, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.7, vjust = 0,
                   n = 2^10, bw = 0.002, text_only = TRUE) +
  geom_textdensity(aes(x = Difference, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.3, vjust = 0,
                   n = 2^10, bw = 0.002, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -0.16, y = 0, 
           label = "italic(tilde('y'))",
           hjust = 0, vjust = 0, 
           family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-0.16, 0.16), oob = scales::oob_keep,
                     breaks = seq(-0.16, 0.16, 0.08),
                     labels = scales::label_number(style_negative = "minus",
                                                   accuracy = c(rep(0.01, 2), 1, rep(0.01, 2)))) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#dabc23", 0.6)),
                    guide = "none") +
  xlab(expression("Δ sinking speed (m s"^-1*")")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme
Fig_S5a_bottom


# 5. Save relevant data ####
speed %>%
  write_rds(here("Sinking", "RDS", "speed.rds"))
speed_prior_posterior %>% 
  write_rds(here("Sinking", "RDS", "speed_prior_posterior.rds"))
speed_diff %>% 
  write_rds(here("Sinking", "RDS", "speed_diff.rds"))
