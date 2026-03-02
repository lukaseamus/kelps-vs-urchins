# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

amphi <- here("Amphipods", "Amphipods.csv") %>%
  read_csv() %>%
  mutate(
    # There's a mix of datetimes and dates
    Deployment = Deployment %>% 
      parse_date_time(orders = c("dmy HM", "dmy")),
    Retrieval = Retrieval %>% 
      parse_date_time(orders = c("dmy HM", "dmy")),
    # Dates are automatically assigned time 00:00, but 12:00 is more reasonable
    Deployment = Deployment %>% 
      update(hour = if_else( hour(.) == 0 , 12 , hour(.) )),
    Retrieval = Retrieval %>% 
      update(hour = if_else( hour(.) == 0 , 12 , hour(.) )),
    # Calculate experimental duration in days
    Day = Deployment %--% Retrieval / ddays(),
    # Food really has 3 rather than 2 levels because kelp was either not
    # grazed or pre-grazed by urchins
    Food = case_when(
      Food == "Faeces" ~ "Faeces",
      Food == "Kelp" & Origin == "Control" ~ "Kelp",
      Food == "Kelp" & Origin == "Urchin" ~ "Grazed kelp"
    )
  ) %T>%
  print()

# 1.2 Check treatments ####
amphi %>%
  mutate(
    Treatment = if_else(
      Tank %>% str_detect("U"), "Urchin", "Control"
    ),
    Check = Origin == Treatment
  ) %>%
  summarise(Identical = all(Check))
  
amphi %>%
  mutate(
    Amphipod = Arena == "Amphipod",
    Grazer = !is.na(Grazer_initial),
    Check = Amphipod == Grazer
  ) %>%
  summarise(Identical = all(Check))

# 1.3 Amphipod size and length ####
# Separate data
amphis <- amphi %>%
  select(Season, Assay, Sporophyte, Day, Food, 
         starts_with("Grazer"), Notes) %>%
  drop_na(Grazer_initial) %>%
  rename(Initial = Grazer_initial,
         Final = Grazer_final,
         Length = Grazer_length) %>% # RGR (% d^-1)
  mutate(Growth = (Final - Initial) / Initial * 100 / Day) %T>%
  print()

# Visualise
amphis %>%
  ggplot(aes(Initial, Final)) +
    geom_point() +
    coord_cartesian(ylim = c(0, 0.2))

amphis %>%
  ggplot(aes(Initial, Growth)) +
    geom_point() +
    coord_cartesian(ylim = c(-10, 10))

amphis %>%
  drop_na(Length) %>%
  ggplot(aes(Length, Initial)) +
    geom_point()

# 1.4 No-choice experiments ####
# 1.4.1 Separate data ####
nochoice <- amphi %>%
  filter(Assay == "No-choice") %>%
  select(where(~ !all(is.na(.)))) %T>%
  print()

# 1.4.2 Check final dry mass ####
nochoice %>%
  filter(Food == "Faeces") %>%
  mutate(
    Remainder = Final_tube - Tube,
    Check = Remainder %>% near(Final_dry)
  ) %>%
  summarise(Identical = all(Check))

# 1.4.3 Calculate initial kelp dry mass ####
# For kelp, this is simply calculated by assuming
# that the dry-fresh mass ratio remains constant.
# Final dry and fresh mass were measured for each
# sample, yielding a sample-specific mass ratio.

nochoice %<>%
  mutate(Initial_dry = Initial * Final_dry / Final) %T>%
  print()

# For faeces this is not possible because final
# fresh mass cannot be determined (the small amount of
# remaining faeces needs to be carefully siphoned, 
# centrifuged and dried before weighing is possible).
# Initial faecal dry mass therefore needs to be 
# derived with a global mass ratio.

# 1.4.4 Load faecal mass data ####
# Load mass data collected specifically for this purpose as
# well as those form the urchin experiment (see Urchins.R).
faeces <- bind_rows(
  Amphipods = here("Amphipods", "Faeces.csv") %>%
    read_csv() %>% select(-Bag),
  Urchins = here("Urchins", "Biomass.csv") %>%
    read_csv() %>% select(Faeces_wet, Faeces_dry) %>%
    rename(Buoyant = Faeces_wet, Dry = Faeces_dry) %>%
    drop_na(),
  .id = "Experiment"
) %>%
  mutate(Experiment = Experiment %>% fct()) %T>%
  print()

faeces %>%
  ggplot(aes(Buoyant, Dry, colour = Experiment)) +
  geom_point()
# Looks like the relationship holds across datasets. 
# Let's look across orders of magnitude.

faeces %>%
  ggplot(aes(Buoyant, Dry, colour = Experiment)) +
    geom_point() +
    geom_abline(slope = 1, intercept = -1) +
    scale_x_log10() +
    scale_y_log10()
# Yes, looks like I can use both sets of data. 
# But the relationship breaks down for very low
# buoyant masses, likely due to measurement error.
faeces %>% arrange(Buoyant)
# The first sensible value is dry mass = 0.03 g
# at buoyant mass = 0.432 g. I cannot simply pick
# and choose data, so the most defensible approach
# is filtering out everything below this threshold.
faeces %<>% filter(!Buoyant < 0.432)
faeces %>% arrange(Buoyant)

faeces %>%
  ggplot(aes(Buoyant, Dry, colour = Experiment)) +
    geom_point() +
    geom_abline(slope = 1, intercept = -1) +
    scale_x_log10() +
    scale_y_log10()
# Still a few outliers but these are due to random
# rather than systematic error. They won't affect 
# inference of a dry-fresh mass ratio.

# 1.4.5 Faecal mass ratio model ####
# Prior simulation
require(extraDistr) # No truncated normal in base R
tibble(n = 1:1e3,
       # Hyperparameters
       logit_beta_mu = rnorm( 1e3 , qlogis(0.1) , 0.5 ),
       logit_beta_sigma = rtnorm( 1e3 , 0 , 0.5 , 0 ), # half-normal prior
       log_sigma_mu = rnorm( 1e3 , log(1) , 0.5 ),
       log_sigma_sigma = rtnorm( 1e3 , 0 , 0.5 , 0 ),
       logit_beta = rnorm( 1e3 , logit_beta_mu , logit_beta_sigma ),
       log_sigma = rnorm( 1e3 , log_sigma_mu , log_sigma_sigma )) %>%
  expand_grid(Buoyant = faeces %$% 
                seq(min(Buoyant), max(Buoyant), length.out = 50)) %>%
  mutate(
    mu = plogis(logit_beta) * Buoyant,
    sigma = exp(log_sigma),
    Dry = rtnorm( n() , mu , sigma , 0 )
  ) %>%
  pivot_longer(cols = c(mu, Dry),
               names_to = "parameter") %>%
  ggplot(aes(Buoyant, value, group = n)) +
    geom_abline() +
    geom_hline(yintercept = faeces %$% range(Dry)) +
    geom_line(alpha = 0.05) +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Based on simulation a 0.1 dry-fresh mass ratio (beta) seems sensible.

# Stan model
require(cmdstanr)
ratio_model <- here("Amphipods", "Stan", "ratio.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
ratio_samples <- ratio_model$sample(
          data = faeces %>% compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()
# Effective sample sizes look good.

# Model checks
ratio_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000239.

require(bayesplot)
ratio_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
ratio_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("logit_beta_mu", "logit_beta_sigma", 
                      "logit_beta[1]", "logit_beta[2]"))
ratio_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_sigma_mu", "log_sigma_sigma", 
                      "log_sigma[1]", "log_sigma[2]"))
# No correlation. Looks good overall.

# Prior-posterior comparison
source("functions.R")
ratio_prior <- prior_samples(
  model = ratio_model, 
  data = faeces %>% compose_data()
)

ratio_prior %>% 
  prior_posterior_draws(
    posterior_samples = ratio_samples,
    group = faeces %>% select(Experiment),
    parameters = c("logit_beta_mu", "logit_beta_sigma", 
                   "logit_beta[Experiment]", 
                   "log_sigma_mu", "log_sigma_sigma", 
                   "log_sigma[Experiment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Experiment")
# All looks good.

# 1.4.6 Figure S4 ####
# Get posteriors
ratio_prior_posterior <- ratio_prior %>% 
  prior_posterior_draws( # Hyperparameters
    posterior_samples = ratio_samples,
    parameters = c("logit_beta_mu", "logit_beta_sigma",
                   "log_sigma_mu", "log_sigma_sigma"),
    format = "short"
  ) %>%
  mutate(
    logit_beta = rnorm( n() , logit_beta_mu , logit_beta_sigma ),
    log_sigma = rnorm( n() , log_sigma_mu , log_sigma_sigma ),
    Experiment = "Global" %>% fct()
  ) %>% 
  select(-c(logit_beta_mu, logit_beta_sigma,
            log_sigma_mu, log_sigma_sigma)) %>%
  bind_rows(
    ratio_prior %>% 
      prior_posterior_draws( # Parameters
        posterior_samples = ratio_samples,
        group = faeces %>% select(Experiment),
        parameters = c("logit_beta[Experiment]", 
                       "log_sigma[Experiment]"),
        format = "short"
      )
  ) %>% 
  filter(Experiment == "Amphipods" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in Experiment
    Experiment = if_else(
      distribution == "prior", "Prior", Experiment
    ) %>% fct(),
    beta = plogis(logit_beta), # Calculate slope (mass ratio)
    sigma = exp(log_sigma) # Calculate likelihood sd
  ) %>%
  select(-distribution) %T>%
  print()

# Summarise
ratio_summary <- ratio_prior_posterior %>%
  summarise(
    across(
      c(beta, sigma),
      list(
        mean = mean,
        sd = sd,
        median = median
      )
    ),
    .by = Experiment
  ) %>%
  full_join(
    faeces %>% 
      count(Experiment) %>%
      add_row(Experiment = "Global" %>% fct(),
              n = sum(.$n))
  ) %T>%
  print()

# Predict across new buoyant masses
ratio_prediction <- ratio_prior_posterior %>%
  spread_continuous(data = faeces,
                    predictor_name = "Buoyant",
                    group_name = "Experiment",
                    length = 150) %>%
  mutate(
    mu = beta * Buoyant,
    Dry = rtnorm( n() , mu , sigma , 0 )
  ) %>%
  group_by(Experiment, Buoyant) %>%
  median_qi(mu, Dry, .width = c(.5, .8, .9)) %T>%
  print()

# Define theme
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
                 strip.text = element_text(size = 12, hjust = 0, face = "bold"),
                 panel.spacing = unit(1, "cm"),
                 text = element_text(family = "Futura"))

# Annotation
require(glue)
ratio_ann <- ratio_summary %>%
  filter(!Experiment %in% c("Prior", "Global")) %>%
  mutate( 
    across( where(is.numeric), ~signif(.x, 2) ),
    label = glue("Ratio = {beta_mean} ± {beta_sd} (n = {n})"),
    Experiment = Experiment %>%
      fct_drop() %>%
      fct_relevel("Amphipods")
  ) %>%
  select(Experiment, label) %T>%
  print()

ratio_ann_global <- ratio_summary %>%
  filter(Experiment == "Global") %>%
  mutate( 
    across( where(is.numeric), ~signif(.x, 2) ),
    label = glue("Ratio = {beta_median} ({beta_mean} ± {beta_sd}, n = {n})")
  ) %>%
  select(label) %T>%
  print()

# Plot
require(ggh4x)
Fig_S4a <- ratio_prediction %>% 
  filter(!Experiment %in% c("Global", "Prior")) %>%
  mutate(Experiment = Experiment %>%
           fct_drop() %>%
           fct_relevel("Amphipods")) %>%
  ggplot() +
    geom_point(data = faeces, aes(Buoyant, Dry),
               shape = 16, size = 2.5, alpha = 0.4, colour = "#7030a5") +
    geom_ribbon(aes(Buoyant, ymin = mu.lower, ymax = mu.upper,
                    alpha = factor(.width)), fill = "#7030a5") +
    geom_line(aes(Buoyant, mu), colour = "#7030a5") +
    geom_text(data = ratio_ann, aes(0, c(0.3, 10), label = label),
              size.unit = "pt", size = 12, family = "Futura",
              hjust = -0.0242, vjust = 1.2) +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    facet_wrap(
      ~ Experiment, scales = "free",
      labeller = labeller(
        Experiment = c(
        Amphipods = "Amphipod experiments",
        Urchins = "Urchin experiments"
        )
      )
    ) +
    facetted_pos_scales(
      x = list(
        Experiment == "Amphipods" ~
          scale_x_continuous(limits = c(0, 5)),
        Experiment == "Urchins" ~
          scale_x_continuous(limits = c(0, 120))
      ),
      y = list(
        Experiment == "Amphipods" ~
          scale_y_continuous(
            limits = c(0, 0.3),
            breaks = seq(0, 0.3, 0.1),
            labels = scales::label_number(accuracy = c(1, rep(0.1, 3)))
          ),
        Experiment == "Urchins" ~
          scale_y_continuous(
            limits = c(0, 10),
            breaks = seq(0, 10, 2.5),
            labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1))
          )
      )
    ) +
    labs(x = "Faecal buoyant mass (g)", y = "Faecal dry mass (g)") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme

Fig_S4a

Fig_S4b <- ratio_prediction %>% 
  filter(Experiment != "Prior") %>%
  mutate(Experiment = Experiment %>%
           fct_drop() %>%
           fct_relevel("Amphipods", "Urchins")) %>%
  ggplot() +
    geom_point(data = faeces, aes(Buoyant, Dry),
               shape = 16, size = 2.5, alpha = 0.4, colour = "#7030a5") +
    geom_ribbon(aes(Buoyant, ymin = mu.lower, ymax = mu.upper,
                    alpha = factor(.width), 
                    group = interaction(Experiment, .width)), 
                fill = "#7030a5") +
    geom_line(aes(Buoyant, mu, group = Experiment), colour = "#7030a5") +
    geom_text(data = ratio_ann_global, aes(0.3, 10, label = label),
              size.unit = "pt", size = 12, family = "Futura",
              hjust = -0.02, vjust = 1.2) +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_log10(breaks = c(0.3, 1, 3, 10, 30, 100),
                  labels = scales::label_number(accuracy = c(0.1, rep(1, 5)))) +
    scale_y_log10(breaks = 10^(-2:1), 
                  labels = scales::label_number(accuracy = c(0.01, 0.1, 1, 1))) +
    labs(x = "Faecal buoyant mass (g)", y = "Faecal dry mass (g)") +
    coord_cartesian(xlim = c(0.3, 100), ylim = 10^c(-2, 1),
                    expand = FALSE, clip = "off") +
    mytheme

Fig_S4b

require(patchwork)
Fig_S4 <- ( Fig_S4a + 
              theme(axis.title = element_blank(),
                    plot.margin = margin(0, 0.5, 1, 0.2, unit = "cm")) ) / 
          ( Fig_S4b +
              theme(plot.margin = margin(0, 0.5, 0.2, 0.2, unit = "cm")) )
Fig_S4

Fig_S4 %>%
  ggsave(filename = "Fig_S4.pdf", path = "Figures",
         device = cairo_pdf, width = 20, height = 15, units = "cm")

# 1.4.7 Calculate initial faecal dry mass ####
# I would usually propagate uncertainty by using the full posterior
# of the ratio, but since initial kelp dry mass is calculated without
# error, this would create a mix of measurement error and no error 
# which is difficult to model. Given that the fit is very tight, I
# think it's defensible to proceed with just the median ratio.
nochoice %>% 
  filter(Food == "Faeces") %$%
  range(Initial)
# Based on this range it makes sense to use the ratio derived
# from the amphipod experiments data which covers it:
faeces %>%
  filter(Experiment == "Amphipods") %$%
  range(Buoyant)

ratio_median <- ratio_summary %>%
  filter(Experiment == "Amphipods") %$% 
  beta_median %T>%
  print()

ratio_median %>% print(digits = 16)
# 0.06100929790540193

nochoice %<>%
  mutate(
    Initial_dry = if_else(
      Food == "Faeces",
      Initial * ratio_median,
      Initial_dry
    )
  ) %T>%
  print()

# 1.4.8 Calculate consumption rates ####
# Correct for autogenic control
nochoice %<>%
  # Biomass loss (g dry mass d^-1)
  mutate(Loss = (Initial_dry - Final_dry) / Day) %>%
  # Select relevant variables
  select(Season, Sporophyte, Arena, Food, 
         Initial_dry, Loss) %>%
  rename(Initial = Initial_dry) %>%
  pivot_wider(names_from = c(Food, Arena),
              values_from = c(Initial, Loss)) %>%
  mutate(
    # Consumption rates (g dry mass amphipod^-1 d^-1)
    Faeces_simple = Loss_Faeces_Amphipod - Loss_Faeces_Control,
    Kelp_simple = Loss_Kelp_Amphipod - Loss_Kelp_Control,
    # For grazed kelp there are cases with no autogenic control;
    # the next best autogenic control is ungrazed kelp
    `Grazed kelp_simple` = if_else(
      !is.na(`Loss_Grazed kelp_Control`),
      `Loss_Grazed kelp_Amphipod` - `Loss_Grazed kelp_Control`,
      `Loss_Grazed kelp_Amphipod` - Loss_Kelp_Control
    ),
    # Alternative with proportional autogenic loss
    Faeces_prop = Loss_Faeces_Amphipod - Loss_Faeces_Control *
      Initial_Faeces_Amphipod / Initial_Faeces_Control,
    Kelp_prop = Loss_Kelp_Amphipod - Loss_Kelp_Control *
      Initial_Kelp_Amphipod / Initial_Kelp_Control,
    `Grazed kelp_prop` = if_else(
      !is.na(`Loss_Grazed kelp_Control`),
      `Loss_Grazed kelp_Amphipod` - `Loss_Grazed kelp_Control` *
        `Initial_Grazed kelp_Amphipod` / `Initial_Grazed kelp_Control`,
      `Loss_Grazed kelp_Amphipod` - Loss_Kelp_Control *
        `Initial_Grazed kelp_Amphipod` / Initial_Kelp_Control
    )
  ) %>%
  select(Season, Sporophyte, 
         Faeces_simple, Kelp_simple, `Grazed kelp_simple`,
         Faeces_prop, Kelp_prop, `Grazed kelp_prop`) %>%
  pivot_longer(cols = -c(Season, Sporophyte)) %>%
  separate(name, into = c("Food", "Rate"), sep = "_") %>%
  mutate(Rate = if_else(Rate == "prop", "Consumption_prop", "Consumption")) %>%
  pivot_wider(names_from = Rate, values_from = value) %T>%
  print()

# Standardise by initial amphipod mass
nochoice %<>%
  left_join(
    amphis %>%
      filter(Assay == "No-choice") %>%
      select(Season, Sporophyte, Food, Initial, Growth) %>%
      rename(Amphipod = Initial)
  ) %>%
  mutate(
    # Mass-standardised rate (g dry mass g^-1 amphipod d^-1)
    Consumption_g = Consumption / Amphipod,
    Consumption_prop_g = Consumption_prop / Amphipod
  ) %T>%
  print()

# 1.4.9 Make sure sporophyte ID is unique ####
nochoice %>% 
  distinct(Season, Sporophyte) %>%
  print(n = 30)
# Currently sporophyte ID repeats across Season, but
# there were different individuals.

nochoice %<>%
  mutate(
    Individual = if_else(
      Season == "Summer", Sporophyte + 15, Sporophyte
    ) %>% factor(),
    Season = Season %>% fct(),
    Food = Food %>% fct()
  ) %>%
  select(-Sporophyte) %T>%
  print()

# Check for NAs
nochoice %>% filter(if_any(everything(), is.na))
# One NA in faecal consumption, another in amphipod growth

# Remove NA in consumption
nochoice %<>% drop_na(Consumption) %T>% print()

# Save data
nochoice %>%
  write_rds(here("Amphipods", "RDS", "nochoice.rds"))

# 1.5 Choice experiments ####
# 1.5.1 Separate data ####
choice <- amphi %>%
  filter(Assay == "Choice") %>%
  select(where(~ !all(is.na(.)))) %T>%
  print()

# 1.5.2 Calculate consumption rates ####
# Correct for autogenic control
choice %<>%
  # Biomass loss (g blotted fresh mass d^-1)
  mutate(Loss = (Initial - Final) / Day) %>%
  # Select relevant variables
  select(Season, Sporophyte, Arena, Food, Initial, Loss) %>%
  pivot_wider(names_from = c(Food, Arena),
              values_from = c(Initial, Loss)) %>%
  mutate(
    # Consumption rates (g blotted fresh mass amphipod^-1 d^-1)
    Faeces_simple = Loss_Faeces_Amphipod - Loss_Faeces_Control,
    `Grazed kelp_simple` = `Loss_Grazed kelp_Amphipod` - `Loss_Grazed kelp_Control`,
    # Alternative with proportional autogenic loss
    Faeces_prop = Loss_Faeces_Amphipod - Loss_Faeces_Control * 
      Initial_Faeces_Amphipod / Initial_Faeces_Control,
    `Grazed kelp_prop` = `Loss_Grazed kelp_Amphipod` - `Loss_Grazed kelp_Control` * 
      `Initial_Grazed kelp_Amphipod` / `Initial_Grazed kelp_Control`,
  ) %>%
  select(Season, Sporophyte, Faeces_simple, `Grazed kelp_simple`,
         Faeces_prop, `Grazed kelp_prop`) %>%
  pivot_longer(cols = -c(Season, Sporophyte)) %>%
  separate(name, into = c("Food", "Rate"), sep = "_") %>%
  mutate(Rate = if_else(Rate == "prop", "Consumption_prop", "Consumption")) %>%
  pivot_wider(names_from = Rate, values_from = value) %T>%
  print()

# Standardise by initial amphipod mass
choice %<>%
  left_join(
    amphis %>%
      filter(Assay == "Choice") %>%
      select(Season, Sporophyte, Food, Initial, Growth) %>%
      rename(Amphipod = Initial)
  ) %>%
  mutate(
    # Mass-standardised rate (g blotted fresh mass g^-1 amphipod d^-1)
    Consumption_g = Consumption / Amphipod,
    Consumption_prop_g = Consumption_prop / Amphipod
  ) %T>%
  print()

# 1.5.3 Make sure sporophyte ID is unique ####
choice %>%
  distinct(Season, Sporophyte)
# Sporophyte ID is already unique.

choice %<>%
  mutate(
    Individual = Sporophyte %>% factor(),
    Season = Season %>% fct(),
    Food = Food %>% fct()
  ) %>%
  select(-Sporophyte) %T>%
  print()

# Check for NAs
choice %>% filter(if_any(everything(), is.na))
# No NAs

# Save data
choice %>%
  write_rds(here("Amphipods", "RDS", "choice.rds"))

# 1.6 Data exploration ####
# 1.6.1 No-choice experiments ####
nochoice %>%
  ggplot(aes(Food, Consumption)) +
    geom_boxplot()

nochoice %>%
  ggplot(aes(Food, Consumption_g)) +
    geom_boxplot()

nochoice %>%
  ggplot(aes(Food, Consumption_prop)) +
    geom_boxplot()

nochoice %>%
  ggplot(aes(Food, Consumption_prop_g)) +
    geom_boxplot()

nochoice %>%
  ggplot(aes(Consumption_g, Consumption_prop_g, colour = Food)) +
    geom_point() +
    geom_abline()

nochoice %>%
  ggplot(aes(Consumption_prop, Consumption_prop_g, colour = Food)) +
    geom_point()
# Very similar. Mass-standardised rates are likely more accurate. 
# Proportional autogenic control is preferable when initial masses
# differ between grazed and control treatments as seems to be the
# case for some replicates in the kelp treatment. It doesn't make
# much of a difference because the initial mass ratio is ~1:
amphi %>% 
  filter(Assay == "No-choice") %>%
  select(Season, Sporophyte, Origin, Arena, Food, Initial) %>% 
  pivot_wider(names_from = Arena, values_from = Initial) %>%
  drop_na(Control) %>%
  ggplot(aes(Amphipod/Control)) + 
    geom_boxplot() +
    geom_vline(xintercept = 1)
# But it still may be more accurate, so I'll go with Consumption_prop_g.

# Negative values are an issue in all cases because they're impossible, 
# so must be due to measurement error. The smallest positive value is
nochoice %>% 
  filter(Consumption_prop_g > 0) %$% 
  min(Consumption_prop_g) # ~5e-4 g
# and balance resolution is 0.001 g. The experimental duration was about a week,
# so the daily rate would have an accuracy of about 1e-4 g. I will replace each 
# zero or negative value with 1e-4.

nochoice %<>%
  # Remove other measures of consumption
  select(-c(Consumption, Consumption_prop, Consumption_g)) %>%
  # Make Consumption_prop_g the main measure of consumption
  rename(Consumption = Consumption_prop_g) %>%
  # Replace negatives and zeros
  mutate(Consumption = if_else(Consumption <= 0, 1e-4, Consumption)) %T>%
  print()

nochoice %>%
  ggplot(aes(Food, Consumption)) +
    geom_boxplot()

# 1.6.2 Choice experiments ####
choice %>%
  ggplot(aes(Food, Consumption)) +
    geom_boxplot()

choice %>%
  ggplot(aes(Food, Consumption_g)) +
    geom_boxplot()

choice %>%
  ggplot(aes(Food, Consumption_prop)) +
    geom_boxplot()

choice %>%
  ggplot(aes(Food, Consumption_prop_g)) +
    geom_boxplot()

choice %>%
  ggplot(aes(Consumption_g, Consumption_prop_g, colour = Food)) +
    geom_point() +
    geom_abline()
# More similar than for no-choice experiment.

choice %>%
  ggplot(aes(Consumption_prop, Consumption_prop_g, colour = Food)) +
    geom_point()
# A bit more influence of amphipod mass.

amphi %>% 
  filter(Assay == "Choice") %>%
  select(Season, Sporophyte, Origin, Arena, Food, Initial) %>% 
  pivot_wider(names_from = Arena, values_from = Initial) %>% 
  ggplot(aes(Amphipod/Control)) + 
    geom_boxplot() +
    geom_vline(xintercept = 1)

# Again, I'll proceed with Consumption_prop_g.
choice %>% 
  filter(Consumption_prop_g > 0) %$% 
  min(Consumption_prop_g) # ~0.11 g
# This is blotted fresh mass so obviously much higher than
# for the no-choice experiment, but the balance accuracy
# and rationale remain the same.
choice %<>%
  select(-c(Consumption, Consumption_prop, Consumption_g)) %>%
  rename(Consumption = Consumption_prop_g) %>%
  mutate(Consumption = if_else(Consumption <= 0, 1e-4, Consumption)) %T>%
  print()

choice %>%
  ggplot(aes(Food, Consumption)) +
    geom_boxplot()

# 2. No-choice model ####
# 2.1 Prior simulation ####
# Gamma likelihood is a sensible choice.
tibble(n = 1:1e4,
       # Highly uncertain global intercept around 0.5 g
       alpha_mu = rnorm( 1e4 , log(0.5) , 2 ),
       # More variation between foods is expected
       alpha_sigma_f = rtnorm( 1e4 , 0 , 1 , 0 ), # Food
       alpha_sigma_s = rtnorm( 1e4 , 0 , 1 , 0 ), # Season
       alpha_sigma_fs = rtnorm( 1e4 , 0 , 1 , 0 ), # Food × Season
       alpha_sigma_i = rtnorm( 1e4 , 0 , 1 , 0 ), # Individual
       log_theta_mu = rnorm( 1e4 , log(0.5) , 1 ),
       log_theta_sigma = rtnorm( 1e4 , 0 , 1 , 0 )) %>%
  mutate(
    mu = exp(
      rnorm( n() , alpha_mu , alpha_sigma_f ) +
        rnorm( n() , 0 , alpha_sigma_s ) +
        rnorm( n() , 0 , alpha_sigma_fs ) +
        rnorm( n() , 0 , alpha_sigma_i )
    ),
    theta = exp( rnorm( n() , log_theta_mu , log_theta_sigma ) ),
    Consumption = rgamma( n() , mu / theta , 1 / theta )
  ) %>%
  pivot_longer(cols = c(mu, Consumption), 
               names_to = "parameter", values_to = "value") %>%
  ggplot() +
    geom_density(aes(value), fill = "black") +
    geom_vline(xintercept = nochoice %$% range(Consumption, na.rm = TRUE),
               colour ="white") +
    scale_x_continuous(limits = c(0, 3), oob = scales::oob_keep) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable

# 2.2 Stan model ####
nochoice_c_model <- here("Amphipods", "Stan", "nochoice_c.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

nochoice_nc_model <- here("Amphipods", "Stan", "nochoice_nc.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

nochoice_c_samples <- nochoice_c_model$sample(
          data = nochoice %>%
            select(Food, Season, Individual, Consumption) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()

nochoice_nc_samples <- nochoice_nc_model$sample(
          data = nochoice %>%
            select(Food, Season, Individual, Consumption) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()

# Save draws
nochoice_c_samples$draws() %>%
  write_rds(here("Amphipods", "RDS", "nochoice_c_samples.rds"))
nochoice_c_samples$draws(format = "df") %>%
  write_rds(here("Amphipods", "RDS", "nochoice_c_samples_df.rds"))

nochoice_nc_samples$draws() %>%
  write_rds(here("Amphipods", "RDS", "nochoice_nc_samples.rds"))
nochoice_nc_samples$draws(format = "df") %>%
  write_rds(here("Amphipods", "RDS", "nochoice_nc_samples_df.rds"))

# 2.3 Model diagnostics ####
# 2.3.1 R-hat and ESS ####
nochoice_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 98% of rhat above 1.001. rhat = 1.01 ± 0.00539.

nochoice_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000992.
# Already the non-centred model is more promising.

# 2.3.2 Chains ####
nochoice_c_chains <- nochoice_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off")

nochoice_nc_chains <- nochoice_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Non-centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off")

( ( nochoice_c_chains / nochoice_nc_chains ) +
  plot_layout(guides = "collect",
              heights = c(7/10, 1)) &
    theme(legend.position = "top",
          legend.justification = 0) ) %>%
  ggsave(filename = "nochoice_chains.pdf", path = here("Amphipods", "Plots"),
         device = cairo_pdf, width = 80, height = 80, units = "cm")
# Non-centred chains are clearly better.

# 2.3.3 Pairs ####
nochoice_c_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("alpha_mu", "alpha_f[1]", "alpha_f[2]", "alpha_f[3]", 
             "alpha_s[1]", "alpha_s[2]", "alpha_fs[1,1]", "alpha_fs[1,2]",
             "alpha_fs[2,1]", "alpha_fs[2,2]", "alpha_fs[3,1]", "alpha_fs[3,2]",
             "alpha_i[1]", "alpha_i[10]", "alpha_i[20]", "log_theta_mu", 
             "log_theta[1]", "log_theta[2]", "log_theta[3]"),
    grid_args = list(top = "Centred model")
  ) %>%
  ggsave(filename = "nochoice_c_pairs.png", path = here("Amphipods", "Plots"),
         width = 95, height = 95, units = "cm", bg = "white")

nochoice_nc_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("alpha_mu", "alpha_f[1]", "alpha_f[2]", "alpha_f[3]", 
             "alpha_s[1]", "alpha_s[2]", "alpha_fs[1,1]", "alpha_fs[1,2]",
             "alpha_fs[2,1]", "alpha_fs[2,2]", "alpha_fs[3,1]", "alpha_fs[3,2]",
             "alpha_i[1]", "alpha_i[10]", "alpha_i[20]", "log_theta_mu", 
             "log_theta[1]", "log_theta[2]", "log_theta[3]"),
    grid_args = list(top = "Non-centred model")
  ) %>%
  ggsave(filename = "nochoice_nc_pairs.png", path = here("Amphipods", "Plots"),
         width = 95, height = 95, units = "cm", bg = "white")
# Some positive posterior correlation within predictors and
# negative correlation between predictors. Nothing bad.
# Non-centred posteriors are smoother.

# 2.4 Prior-posterior comparison ####
# 2.4.1 Sample prior ####
# prior_samples only works properly with non-centred models, but
# this is fine since priors are identical for both models
nochoice_prior <- prior_samples(
  model = nochoice_nc_model, 
  data = nochoice %>%
    select(Food, Season, Individual, Consumption) %>%
    compose_data()
)

# 2.4.2 Centred model ####
nochoice_c_prior_posterior_global <- nochoice_prior %>% 
  prior_posterior_draws(
    posterior_samples = nochoice_c_samples,
    parameters = c("alpha_mu", "alpha_sigma_f", "alpha_sigma_s", 
                   "alpha_sigma_fs", "alpha_sigma_i", 
                   "log_theta_mu", "log_theta_sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot() +
  labs(title = "Centred model")

nochoice_c_prior_posterior_fs <- nochoice_prior %>% 
  prior_posterior_draws(
    posterior_samples = nochoice_c_samples,
    group = nochoice %>% select(Food, Season),
    parameters = c("alpha_f[Food]", "alpha_s[Season]",
                   "alpha_fs[Food, Season]"),
    format = "long"
  ) %>%
  prior_posterior_plot(group_name = "Food",
                       second_group_name = "Season")

nochoice_c_prior_posterior_i <- nochoice_prior %>% 
  prior_posterior_draws(
    posterior_samples = nochoice_c_samples,
    group = nochoice %>% select(Individual),
    parameters = c("alpha_i[Individual]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Individual", ridges = TRUE)

# 2.4.3 Non-centred model ####
nochoice_nc_prior_posterior_global <- nochoice_prior %>% 
  prior_posterior_draws(
    posterior_samples = nochoice_nc_samples,
    parameters = c("alpha_mu", "alpha_sigma_f", "alpha_sigma_s", 
                   "alpha_sigma_fs", "alpha_sigma_i", 
                   "log_theta_mu", "log_theta_sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot() +
  labs(title = "Non-centred model")

nochoice_nc_prior_posterior_fs <- nochoice_prior %>% 
  prior_posterior_draws(
    posterior_samples = nochoice_nc_samples,
    group = nochoice %>% select(Food, Season),
    parameters = c("alpha_f[Food]", "alpha_s[Season]",
                   "alpha_fs[Food, Season]"),
    format = "long"
  ) %>%
  prior_posterior_plot(group_name = "Food",
                       second_group_name = "Season")

nochoice_nc_prior_posterior_i <- nochoice_prior %>% 
  prior_posterior_draws(
    posterior_samples = nochoice_nc_samples,
    group = nochoice %>% select(Individual),
    parameters = c("alpha_i[Individual]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Individual", ridges = TRUE)

# 2.4.4 Combined ####
nochoice_prior_posterior <- 
  ( ( nochoice_c_prior_posterior_global /
        nochoice_c_prior_posterior_fs / 
        nochoice_c_prior_posterior_i ) +
      plot_layout(heights = c(0.4, 0.6, 1)) |
    ( nochoice_nc_prior_posterior_global /
        nochoice_nc_prior_posterior_fs / 
        nochoice_nc_prior_posterior_i ) +
      plot_layout(heights = c(0.4, 0.6, 1)) ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", 
        legend.justification = 0)

nochoice_prior_posterior %>%
  ggsave(filename = "nochoice_prior_posterior.pdf", 
         path = here("Amphipods", "Plots"),
         device = cairo_pdf, width = 40, height = 40, units = "cm")
# Non-centred model posteriors are smoother.
# Proceed with non-centred model.

# 2.5 Prediction ####
# 2.5.1 Parameter distributions ####
nochoice_prior_posterior <- bind_rows(
  # Global
  nochoice_prior %>% 
    prior_posterior_draws(
      posterior_samples = nochoice_nc_samples,
      parameters = c("alpha_mu", "alpha_sigma_f", "alpha_sigma_s",
                     "alpha_sigma_fs", "alpha_sigma_i", 
                     "log_theta_mu", "log_theta_sigma"),
      format = "short"
    ) %>%
    mutate( 
      # Calculate consumption for new foods
      mu = exp(
        rnorm( n() , alpha_mu , alpha_sigma_f ) +
          rnorm( n() , 0 , alpha_sigma_s ) +
          rnorm( n() , 0 , alpha_sigma_fs ) +
          rnorm( n() , 0 , alpha_sigma_i )
      ),
      theta = exp(
        rnorm( n() , log_theta_mu , log_theta_sigma )
      ),
      Consumption = rgamma( n() , mu / theta , 1 / theta ),
      Food = "New" %>% fct(),
      Season = "Annual" %>% fct()
    ),
  # Food
  nochoice_prior %>% 
    prior_posterior_draws(
      posterior_samples = nochoice_nc_samples,
      group = nochoice %>% select(Food),
      parameters = c("alpha_f[Food]", "alpha_sigma_s", 
                     "alpha_sigma_fs", "alpha_sigma_i", 
                     "log_theta[Food]"),
      format = "short"
    ) %>%
    mutate( 
      # Calculate consumption of existing foods 
      # for new seasons and sporophytes
      mu = exp(
        rnorm( n() , alpha_f , alpha_sigma_s ) +
          rnorm( n() , 0 , alpha_sigma_fs ) +
          rnorm( n() , 0 , alpha_sigma_i )
      ),
      theta = exp( log_theta ),
      Consumption = rgamma( n() , mu / theta , 1 / theta ),
      Season = "Annual" %>% fct()
    ),
  # Season
  nochoice_prior %>% 
    prior_posterior_draws(
      posterior_samples = nochoice_nc_samples,
      group = nochoice %>% select(Season),
      parameters = c("alpha_s[Season]", "alpha_sigma_f", 
                     "alpha_sigma_fs", "alpha_sigma_i",
                     "log_theta_mu", "log_theta_sigma"),
      format = "short"
    ) %>%
    mutate( 
      # Calculate consumption of existing seasons 
      # for new foods and sporophytes
      mu = exp(
        rnorm( n() , alpha_s , alpha_sigma_f ) +
          rnorm( n() , 0 , alpha_sigma_fs ) +
          rnorm( n() , 0 , alpha_sigma_i )
      ),
      theta = exp(
        rnorm( n() , log_theta_mu , log_theta_sigma )
      ),
      Consumption = rgamma( n() , mu / theta , 1 / theta ),
      Food = "New" %>% fct()
    )
) %>%
  filter(Food == "Kelp" & Season == "Annual" & distribution == "prior" |
           distribution == "posterior") %>%
  mutate(
    Food = if_else(
      distribution == "prior", "Prior", Food
    ) %>% fct(),
    Season = if_else(
      distribution == "prior", "Prior", Season
    ) %>% fct()
  ) %>%
  select(starts_with("."), Food, Season, mu, Consumption) %T>%
  print()

nochoice_prior_posterior %>% # Save
  write_rds(here("Amphipods", "RDS", "nochoice_prior_posterior.rds"))

nochoice_prior_posterior_interaction <- nochoice_prior %>% 
  prior_posterior_draws(
    posterior_samples = nochoice_nc_samples,
    group = nochoice %>% select(Food, Season),
    parameters = c("alpha_f[Food]", "alpha_s[Season]",
                   "alpha_fs[Food, Season]", "alpha_sigma_i", 
                   "log_theta[Food]"),
    format = "short"
  ) %>%
  mutate( 
    # Calculate consumption for each food and season combination 
    # across sporophytes
    mu = exp(
      alpha_f + alpha_s + alpha_fs + rnorm( n() , 0 , alpha_sigma_i )
    ),
    theta = exp( log_theta ),
    Consumption = rgamma( n() , mu / theta , 1 / theta )
  ) %>%
  filter(Food == "Kelp" & Season == "Spring" & distribution == "prior" |
           distribution == "posterior") %>%
  mutate(
    Food = if_else(
      distribution == "prior", "Prior", Food
    ) %>% fct(),
    Season = if_else(
      distribution == "prior", "Prior", Season
    ) %>% fct()
  ) %>%
  select(starts_with("."), Food, Season, mu, Consumption) %T>%
  print()

nochoice_prior_posterior_interaction %>% # Save
  write_rds(here("Amphipods", "RDS", "nochoice_prior_posterior_interaction.rds"))

# 2.5.2 Contrast distributions ####
# Food
nochoice_contrast_food <- nochoice_prior_posterior %>%
  filter(!Food %in% c("Prior", "New")) %>%
  droplevels() %>%
  pivot_wider(names_from = Food, values_from = c(mu, Consumption)) %>%
  mutate(
    # Means
    Faeces_Kelp_Diff_mu = mu_Faeces - mu_Kelp, # Calculate differences
    `Faeces_Grazed kelp_Diff_mu` = mu_Faeces - `mu_Grazed kelp`,
    `Kelp_Grazed kelp_Diff_mu` = mu_Kelp - `mu_Grazed kelp`,
    Faeces_Kelp_Ratio_mu = mu_Faeces / mu_Kelp, # Calculate ratios
    `Faeces_Grazed kelp_Ratio_mu` = mu_Faeces / `mu_Grazed kelp`,
    `Kelp_Grazed kelp_Ratio_mu` = mu_Kelp / `mu_Grazed kelp`,
    # Observations
    Faeces_Kelp_Diff_obs = Consumption_Faeces - Consumption_Kelp,
    `Faeces_Grazed kelp_Diff_obs` = Consumption_Faeces - `Consumption_Grazed kelp`,
    `Kelp_Grazed kelp_Diff_obs` = Consumption_Kelp - `Consumption_Grazed kelp`,
    Faeces_Kelp_Ratio_obs = Consumption_Faeces / Consumption_Kelp,
    `Faeces_Grazed kelp_Ratio_obs` = Consumption_Faeces / `Consumption_Grazed kelp`,
    `Kelp_Grazed kelp_Ratio_obs` = Consumption_Kelp / `Consumption_Grazed kelp`
  ) %>%
  select(starts_with("."), ends_with("mu"), ends_with("obs")) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = c("a", "b", "type", "Parameter"),
               values_to = "value",
               names_pattern = "^(.*)_(.*)_(.*)_(.*)$") %>%
  unite("Contrast", a, b, sep = " vs. ") %>%
  pivot_wider(values_from = value, names_from = type) %T>%
  print()

# Season
nochoice_contrast_season <- nochoice_prior_posterior %>%
  filter(!Season %in% c("Prior", "Annual")) %>%
  droplevels() %>%
  pivot_wider(names_from = Season, values_from = c(mu, Consumption)) %>%
  mutate(
    Summer_Spring_Diff_mu = mu_Summer - mu_Spring, # Calculate difference
    Summer_Spring_Ratio_mu = mu_Summer / mu_Spring, # Calculate ratio
    Summer_Spring_Diff_obs = Consumption_Summer - Consumption_Spring,
    Summer_Spring_Ratio_obs = Consumption_Summer / Consumption_Spring
  ) %>%
  select(starts_with("."), ends_with("mu"), ends_with("obs")) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = c("a", "b", "type", "Parameter"),
               values_to = "value",
               names_pattern = "^(.*)_(.*)_(.*)_(.*)$") %>%
  unite("Contrast", a, b, sep = " vs. ") %>%
  pivot_wider(values_from = value, names_from = type) %T>%
  print()

# 2.5.3 Estimates ####
# Something to be aware of: gamma introduced zeros
nochoice_prior_posterior %>%
  count(Consumption == 0)
nochoice_prior_posterior_interaction %>%
  count(Consumption == 0)
nochoice_contrast_food %>%
  count(Ratio == 0)
nochoice_contrast_season %>%
  count(Ratio == 0)
# A few Infs or NaNs were introduced because of this
nochoice_contrast_food %>%
  count(!is.finite(Ratio))
nochoice_contrast_season %>%
  count(!is.finite(Ratio))
# This makes summarising logs a little tricky

nochoice_summary <- bind_rows(
  Global = nochoice_prior_posterior,
  Interaction = nochoice_prior_posterior_interaction,
  .id = "Effect"
) %>%
  mutate(
    log_mu = log(mu), # Here I ensure there are no -Infs
    log_Consumption = if_else(Consumption == 0, NA, log(Consumption))
  ) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = T), # Ensure removal on NAs
        sd = ~sd(.x, na.rm = T),
        median = ~median(.x, na.rm = T)
      )
    ),
    n = n(),
    .by = c(Effect, Food, Season)
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median")),
      ~signif(.x, 2)
    ),
    mu = glue( "{mu_median} ({log_mu_mean} ± {log_mu_sd})" ),
    Consumption = glue( 
      "{Consumption_median} ({log_Consumption_mean} ± {log_Consumption_sd})" 
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

nochoice_contrast_summary <- bind_rows(
  Food = nochoice_contrast_food,
  Season = nochoice_contrast_season,
  .id = "Predictor"
) %>%
  mutate(
    log_Ratio = if_else(
      Ratio == 0 | !is.finite(Ratio), NA, log(Ratio)
    )
  ) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = T),
        sd = ~sd(.x, na.rm = T),
        median = ~median(.x, na.rm = T)
      )
    ),
    n = n(),
    P = mean( Diff > 0 ),
    .by = c(Predictor, Contrast, Parameter)
  ) %>%
  mutate(
    across(
      c(contains("mean"), contains("sd"), contains("median")),
      ~signif(.x, 2)
    ),
    Diff = glue( "{Diff_mean} ± {Diff_sd} ({Diff_median})" ),
    Ratio = glue( 
      "{Ratio_median} ({log_Ratio_mean} ± {log_Ratio_sd})" 
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 3. Absolute choice model ####

# 4. Relative choice model ####


# 5. Figure 4 ####
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


# Prediction
require(ggridges)
Fig_4a_top <- nochoice %>%
  mutate(Food = Food %>% fct_relevel("Faeces", "Grazed kelp"),
         Kelp = Food %>% str_detect("Kelp|kelp")) %>%
  ggplot() +
  geom_jitter(aes(x = Consumption, y = Food %>% as.numeric() - 0.5,
                  colour = Kelp, shape = Season),
              alpha = 0.4, size = 2.5, height = 0.36) +
  stat_density_ridges(data = nochoice_prior_posterior_interaction %>%
                        filter(Food != "Prior") %>%
                        droplevels() %>%
                        mutate(Food = Food %>% fct_relevel("Faeces", "Grazed kelp"),
                               Kelp = Food %>% str_detect("Kelp|kelp")),
                      aes(x = mu, y = Food %>% as.numeric(), 
                          group = interaction(Food, Season), fill = Kelp),
                      colour = NA, n = 2^10,
                      from = 0, to = 0.6, rel_min_height = 0.001,
                      bandwidth = 0.6*0.02, scale = 1.2, alpha = 0.7) +
  stat_density_ridges(data = nochoice_prior_posterior %>%
                        filter(!Food %in% c("Prior", "New")) %>%
                        droplevels() %>%
                        mutate(Food = Food %>% fct_relevel("Faeces", "Grazed kelp"),
                               Kelp = Food %>% str_detect("Kelp|kelp")),
                      aes(x = mu, y = Food %>% as.numeric(), group = Food), 
                      fill = NA, n = 2^10, from = 0, to = 0.6, 
                      rel_min_height = 0.001,
                      bandwidth = 0.6*0.02, scale = 1.2) +
  geom_text(aes(0.6, Food %>% as.numeric(), label = Food),
            family = "Futura", size.unit = "pt", size = 12, 
            hjust = 1, vjust = -0.5, check_overlap = TRUE) +
  scale_x_continuous(limits = c(0, 0.6), #breaks = seq(0, 60, 20),
                     oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23"),
                    labels = c("TRUE" = "Kelp",
                               "FALSE" = "Faeces"),
                    guide = guide_legend(reverse = TRUE, order = 1)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17), # circle, triangle
                     # override grey fill legend shapes
                     guide = guide_legend(override.aes = list(shape = c(1, 2)))) +
  xlab("Consumption") +
  coord_cartesian(ylim = c(0, 4.5), expand = FALSE, clip = "off") +
  mytheme

Fig_4a_top

# Difference
require(geomtextpath)
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



