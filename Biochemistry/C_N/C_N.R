# 1. Load data ####
# 1.1 Load raw data ####
require(tidyverse)
require(here)

C_N <- here("Biochemistry", "C_N", "C_N.csv") %>%
  read_csv() %>%
  mutate(Date = Date %>% dmy(),
         ID = ID %>% fct(),
         Method = Method %>% fct())
C_N

# 1.2 Visualise ####
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

Fig_S2a_top <- C_N %>%
  filter(year(Date) != 2025) %>%
  mutate(Oven = if_else(ID %>% str_detect("_oven"),
                        "Yes", "No") %>% fct(),
         ID = ID %>% str_remove("_oven") %>% fct()) %>%
  ggplot(aes(Mass, C, colour = Method, 
             shape = Oven, group = ID)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_path(alpha = 0.6, linewidth = 0.7) +
    scale_shape_manual(values = c(16, 17), guide = "none") +
    scale_colour_manual(values = c("#f5a54a", "#6a98b4"),
                        labels = c("Thermal conductivity detector",
                                   "Isotope ratio mass spectrometer")) +
    scale_x_continuous(labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1, 0.1))) +
    scale_y_continuous(breaks = seq(18, 42, 6)) +
    labs(x = "Sample mass (mg)", y = "Carbon content (%)") +
    coord_cartesian(xlim = c(1, 3.5), ylim = c(18, 42),
                    expand = FALSE, clip = "off") +
    mytheme
Fig_S2a_top

Fig_S2a_bottom <- C_N %>%
  filter(year(Date) != 2025) %>%
  mutate(Oven = if_else(ID %>% str_detect("_oven"),
                        "Yes", "No") %>% fct(),
         ID = ID %>% str_remove("_oven") %>% fct()) %>%
  ggplot(aes(Mass, N, colour = Method, 
             shape = Oven, group = ID)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_path(alpha = 0.6, linewidth = 0.7) +
    scale_shape_manual(values = c(16, 17), guide = "none") +
    scale_colour_manual(values = c("#f5a54a", "#6a98b4"),
                        labels = c("Thermal conductivity detector",
                                   "Isotope ratio mass spectrometer")) +
    scale_x_continuous(labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1, 0.1))) +
    scale_y_continuous(labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))) +
    labs(x = "Sample mass (mg)", y = "Nitrogen content (%)") +
    coord_cartesian(xlim = c(1, 3.5), ylim = c(0, 1.5),
                    expand = FALSE, clip = "off") +
    mytheme
Fig_S2a_bottom

# Differences in method are more obvious for nitrogen than for carbon. Generally,
# sample masses below 2 mg seem to be too low and the thermal conductivity
# detector (TCD) underestimates relative to the isotope ratio mass spectrometer 
# (IRMS). Oven re-drying generally increased the estimated elemental content.

# I trust IRMS more but there seems to be no difference between ~2.2 and ~3.2 mg
# sample mass, so I'll average across these higher sample masses. I only have oven 
# re-dried samples for a few cases so cannot include those. Unfortunately, I only
# have ~1.2 mg sample mass TCD estimates for sample U3f, so those values will need 
# to be adjusted using the mean difference between the ~1.2 mg TCD estimate and 
# the mean of the IRMS estimates.

Fig_S2b <- C_N %>%
  drop_na(d15N, d13C) %>%
  mutate(Treatment = if_else(ID %>% str_detect("f"),
                             "Faeces", "Kelp") %>% fct_relevel("Faeces"),
         d15N_confidence = if_else(is.na(Notes), 
                                   "High", "Low") %>% fct()) %>%
  ggplot() +
    geom_point(aes(d13C, d15N, colour = Treatment, shape = d15N_confidence),
               size = 3, alpha = 0.6) +
    geom_polygon(data = . %>% group_by(Treatment) %>%
                   slice(chull(d13C, d15N)),
                 aes(d13C, d15N, fill = Treatment, colour = Treatment),
                 alpha = 0.2, linewidth = 0.7) +
    scale_shape_manual(values = c(17, 16), guide = "none") +
    scale_fill_manual(values = c("#7030a5", "#c3b300"),
                      guide = guide_legend(reverse = TRUE)) +
    scale_colour_manual(values = c("#7030a5", "#c3b300"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_y_continuous(labels = scales::label_number(accuracy = c(0.1, 1, 0.1, 1, 0.1, 1))) +
    scale_x_continuous(breaks = seq(-22, -15, 1),
                       labels = scales::label_number(style_negative = "minus")) +
    labs(x = expression(delta^13*"C"["VPDB"]*" (‰)"), y = expression(delta^15*"N"["air"]*" (‰)")) +
    coord_cartesian(ylim = c(1.5, 4), xlim = c(-22, -15),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.position = "inside",
          legend.position.inside = c(0.85, 0.9))
Fig_S2b

# Faeces and Kelp would be hard to tease apart with stable isotope analysis.

require(patchwork)
Fig_S2 <- ( ( Fig_S2a_top + 
                theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank()) ) / 
            ( Fig_S2a_bottom + 
                theme(legend.position = "none") ) /
              Fig_S2b ) +
          plot_annotation(tag_levels = list(c("a", "", "b"))) &
          theme(plot.tag = element_text(family = "Futura", size = 15, face = "bold"))

Fig_S2 %>%
  ggsave(filename = "Fig_S2.pdf", device = cairo_pdf, path = "Figures", 
         height = 21, width = 16, units = "cm")

# 1.3 Prepare data ####
# Average IRMS estimates for samples with more than one repeat.
require(magrittr)
C_N %<>%
  filter(year(Date) != 2024) %>% # filter out third run
  left_join( # then join it by column with the other runs
    C_N %>% 
      filter(Date %>% year() == 2024) %>%
      select(ID, Method, N, C),
    by = c("ID", "Method") # join by ID so re-dried samples are filtered out
  ) %>% # re-dried estimates are removed by left_join
  mutate(N = (N.x + N.y) / 2, # calculate the mean of IRMS samples
         C = (C.x + C.y) / 2,
         N = coalesce(N, N.x), # merge IRMS means with remainder
         C = coalesce(C, C.x))
C_N

# Get difference between initial TDC and average IRMS.
diff <- C_N %>%
  filter(Date == "2023-11-23" | 
           Date == "2023-12-21" & Method == "IRMS") %>%
  select(ID, Method, N, C) %>%
  pivot_wider(names_from = Method, values_from = c(N, C)) %>%
  mutate(N = N_IRMS - N_TCD,
         C = C_IRMS - C_TCD) %>%
  drop_na() %>%
  summarise(N = mean(N),
            C = mean(C))
diff 
# In the first run (~1.2 mg, TCD) nitrogen was underestimated by 0.36% 
# and carbon by 1.6% relative to the second and third repeats (IRMS).

# Adjust U3f.
C_N %<>%
  mutate(N = if_else(ID == "U3f",
                     N + diff$N,
                     N),
         C = if_else(ID == "U3f",
                     C + diff$C,
                     C),
         Method = if_else(ID == "U3f",
                          "IRMS",
                          Method) %>% fct())

# Select relevant data.
C_N %<>%
  filter(Method == "IRMS") %>%
  select(Date, ID, N, C) 

# Add explanatory variables from ID.
C_N %<>%
  mutate(Treatment = if_else(ID %>% str_detect("f") | ID %>% str_detect("F"),
                             "Faeces", "Kelp") %>% fct(),
         Season = case_when(
                    Date %>% year() == 2023 ~ "Autumn",
                    Date %>% year() == 2025 &
                      ID %>% str_split_i(pattern = "_",
                                         i = 2) %>%
                      str_detect("1") ~ "Spring",
                    Date %>% year() == 2025 &
                      ID %>% str_split_i(pattern = "_",
                                         i = 2) %>%
                      str_detect("2") ~ "Summer"
                    ) %>% fct(),
         Individual = if_else(Date %>% year() == 2023,
                              ID %>% str_extract("\\d+"),
                              ID %>% str_split_i(pattern = "_",
                                                 i = 1) %>%
                                str_extract("\\d+")),
         Individual = case_when( # Make sporophyte number unique
                        Season == "Autumn" ~ Individual %>% as.numeric(),
                        Season == "Spring" ~ Individual %>% as.numeric() + 15,
                        Season == "Summer" ~ Individual %>% as.numeric() + 30
                        ) %>% str_c() %>% fct()
         )

# These are the data that will be analysed.
C_N %<>% 
  select(-c(Date, ID)) %T>%
  print(n = 79)

# 2. Carbon ####
# 2.1 Prior simulation ####
# For Laminaria hyperborea, carbon content is expected around 30%
# (28.67 ± 0.41% and 30.78 ± 0.92%, Wright et al., 10.1111/gcb.16299)
# but ranges from 25% to 40% across the North Atlantic (Filbee-Dexter et al., 
# 10.1371/journal.pbio.3001702). It is safer to use a beta likelihood
# than a gamma likelihood here because we're getting close to the 100%
# ceiling.

require(ggdist)
tibble(n = 1:1e5,
       mu_logit = rnorm( 1e5 , log( 0.3 / (1 - 0.3) ) , 0.3 ),
       # sigma = rexp( 1e5, 20 ),
       mu = 1 / ( 1 + exp(-mu_logit) ),
       # nu = ( mu * (1 - mu) / sigma^2 - 1 ),
       nu = rgamma( 1e5 , 30^2 / 20^2 , 30 / 20^2 ),
       sigma = sqrt( mu * ( 1 - mu ) / ( 1 + nu ) ),
       C = rbeta( 1e5 , mu * nu , (1 - mu) * nu )) %>% # %$% hist(sigma)
  pivot_longer(cols = c(mu, C), 
               names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, parameter)) +
    stat_slab(alpha = 0.5, height = 2, n = 3e3) +
    coord_cartesian(expand = F,
                    xlim = c(0, 1)) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 2.2 Stan models ####
require(cmdstanr)
C_c_model <- here("Biochemistry", "C_N", "Stan", "C_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

C_nc_model <- here("Biochemistry", "C_N", "Stan", "C_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
C_c_samples <- C_c_model$sample(
          data = C_N %>%
            select(-N) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

C_nc_samples <- C_nc_model$sample(
          data = C_N %>%
            select(-N) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

# 2.3 Model checks ####
# 2.3.1 Rhat ####
C_c_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Most rhat above 1.001. rhat = 1.00 ± 0.00169.

C_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000592.

# Plot comparison between centred and non-centred parameterisation.
C_nc_samples$summary() %>%
  left_join(C_c_samples$summary(),
            by = "variable") %>%
  rename(rhat_nc = rhat.x, rhat_c = rhat.y) %>%
  ggplot(aes(rhat_c, rhat_nc)) +
    geom_abline(slope = 1) +
    geom_point() +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Warning because z-scores are dropped as they have no equivalent in
# the centred parameteristion. The non-centred model is better.

# 2.3.2 Chains ####
require(bayesplot)
C_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "C_c_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "C_N", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains aren't great.

C_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "C_nc_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "C_N", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains are very good.

# 2.3.3 Pairs ####
C_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]", 
                      "nu[1]"))

C_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "nu[2]"))
# Correlation between alpha_t and alpha_s for both treatments.

C_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]",
                      "nu[1]"))

C_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "nu[2]"))
# Same as above.

# 2.4 Prior-posterior comparison ####
# 2.4.1 Sample priors ####
source("functions.R")
C_c_prior <- prior_samples(
  model = C_c_model,
  data = C_N %>%
    select(-N) %>%
    compose_data()
  )

C_nc_prior <- prior_samples(
  model = C_nc_model,
  data = C_N %>%
    select(-N) %>%
    compose_data()
)

# 2.4.2 Plot prior-posterior comparison ####
C_c_prior %>% 
  prior_posterior_draws(
    posterior_samples = C_c_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "nu[Treatment]"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# sigma_s seems a bit large, potentially because alpha_s is absorbing 
# some of alpha_t. Otherwise generally looks fine.

C_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = C_nc_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "z_s[Treatment, Season]", 
                   "z_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "nu[Treatment]"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# Same as above.

# 2.5 Stan models ####
# I can do better by using sum_to_zero_vector for alpha_s and alpha_i because
# I know a priori that they shouldn't contribute to the Treatment effect, so
# should be identifiable.

C_c_s2z_model <- here("Biochemistry", "C_N", "Stan", "C_c_s2z.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

C_nc_s2z_model <- here("Biochemistry", "C_N", "Stan", "C_nc_s2z.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

C_c_s2z_samples <- C_c_s2z_model$sample(
          data = C_N %>%
            select(-N) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

C_nc_s2z_samples <- C_nc_s2z_model$sample(
          data = C_N %>%
            select(-N) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

# 2.6 Model checks ####
# 2.6.1 Rhat ####
C_c_s2z_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Most rhat above 1.001. rhat = 1.01 ± 0.00333.

C_nc_s2z_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000813.

# Plot comparison between centred and non-centred parameterisation.
C_nc_s2z_samples$summary() %>%
  left_join(C_c_s2z_samples$summary(),
            by = "variable") %>%
  rename(rhat_nc = rhat.x, rhat_c = rhat.y) %>%
  ggplot(aes(rhat_c, rhat_nc)) +
    geom_abline(slope = 1) +
    geom_point() +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Warning because z-scores are dropped as they have no equivalent in
# the centred parameteristion. The non-centred model is better.

# 2.6.2 Chains ####
C_c_s2z_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "C_c_s2z_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "C_N", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains are slightly better than before.

C_nc_s2z_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "C_nc_s2z_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "C_N", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains are great.

# 2.6.3 Pairs ####
C_c_s2z_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]", 
                      "nu[1]"))

C_c_s2z_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "nu[2]"))
# Correlation gone. Looks good.

C_nc_s2z_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]",
                      "nu[1]"))

C_nc_s2z_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "nu[2]"))
# Looks optimal.

# 2.7 Prior-posterior comparison ####
# 2.7.1 Sample priors ####
C_c_s2z_prior <- prior_samples(
  model = C_c_s2z_model,
  data = C_N %>%
    select(-N) %>%
    compose_data()
  )

C_nc_s2z_prior <- prior_samples(
  model = C_nc_s2z_model,
  data = C_N %>%
    select(-N) %>%
    compose_data()
)

# 2.7.2 Plot prior-posterior comparison ####
C_c_s2z_prior %>% 
  prior_posterior_draws(
    posterior_samples = C_c_s2z_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "nu[Treatment]"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# Looks better than before. Much stronger alpha_t predictions.

C_nc_s2z_prior %>% 
  prior_posterior_draws(
    posterior_samples = C_nc_s2z_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "z_s[Treatment, Season]", 
                   "z_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "nu[Treatment]"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# Looks optimal.

# 2.8 Prediction ####
# 2.8.1 Combine relevant priors and posteriors ####
C_prior_posterior <- C_nc_s2z_prior %>% 
  prior_posterior_draws(
    posterior_samples = C_nc_s2z_samples,
    group = C_N %>%
      select(-N) %>%
      select(Treatment),
    parameters = c("alpha_t[Treatment]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "nu[Treatment]"),
    format = "short"
  )

# 2.8.2 Calculate predictions for new seasons and sporophytes ####
C_prior_posterior %<>%
  mutate(mu = 1 / ( 1 + exp( -alpha_t ) ), # inverse logit (logit = log(p / (1 - p)))
         obs = rbeta( n() , mu * nu , (1 - mu) * nu ),
         mu_new = 1 / ( 1 + exp( -( alpha_t + rnorm( n() , 0 , sigma_s ) + rnorm( n() , 0 , sigma_i ) ) ) ),
         obs_new = rbeta( n() , mu_new * nu , (1 - mu_new) * nu ))

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

# 3. Nitrogen ####
# 3.1 Prior simulation ####
# For Laminaria hyperborea, nitrogen content is expected around 1-2%
# (1.74%, Wright et al., 10.1111/gcb.16299) but ranges from 0.5% to 
# 3% across the North Atlantic (Filbee-Dexter et al., 
# 10.1371/journal.pbio.3001702). These are similarly low concentrations
# to phenolic content, so nitrogen can safely be estimated with a
# gamma likelihood, which is easier to parameterise.

tibble(n = 1:1e5,
       mu_log = rnorm( 1e5 , log(1.74) , 0.4 ),
       theta = rexp( 1e5, 5 ),
       mu = exp(mu_log),
       N = rgamma( 1e5 , mu / theta , 1 / theta )) %>%
  pivot_longer(cols = c(mu, N), 
               names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, parameter)) +
    stat_slab(alpha = 0.5, height = 2, n = 3e3) +
    coord_cartesian(expand = F,
                    xlim = c(0, 5)) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 3.2 Stan models ####
N_c_model <- here("Biochemistry", "C_N", "Stan", "N_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

N_nc_model <- here("Biochemistry", "C_N", "Stan", "N_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

N_c_samples <- N_c_model$sample(
          data = C_N %>%
            select(-C) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

N_nc_samples <- N_nc_model$sample(
          data = C_N %>%
            select(-C) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

# 3.3 Model checks ####
# 3.3.1 Rhat ####
N_c_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Most rhat above 1.001. rhat = 1.01 ± 0.00671.

N_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Less than 5% of rhat above 1.001. rhat = 1.00 ± 0.000217.

# Plot comparison between centred and non-centred parameterisation.
N_nc_samples$summary() %>%
  left_join(N_c_samples$summary(),
            by = "variable") %>%
  rename(rhat_nc = rhat.x, rhat_c = rhat.y) %>%
  ggplot(aes(rhat_c, rhat_nc)) +
    geom_abline(slope = 1) +
    geom_point() +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Warning because z-scores are dropped as they have no equivalent in
# the centred parameteristion. The non-centred model is better.

# 3.3.2 Chains ####
N_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "N_c_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "C_N", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains aren't great.

N_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "N_nc_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "C_N", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains are pretty good.

# 3.3.3 Pairs ####
N_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]", 
                      "theta[1]"))

N_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "theta[2]"))
# Correlation between alpha_t and alpha_s for both treatments.

N_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]",
                      "theta[1]"))

N_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "theta[2]"))
# Same as above.

# 3.4 Prior-posterior comparison ####
# 3.4.1 Sample priors ####
N_c_prior <- prior_samples(
  model = N_c_model,
  data = C_N %>%
    select(-C) %>%
    compose_data()
  )

N_nc_prior <- prior_samples(
  model = N_nc_model,
  data = C_N %>%
    select(-C) %>%
    compose_data()
)

# 3.4.2 Plot prior-posterior comparison ####
N_c_prior %>% 
  prior_posterior_draws(
    posterior_samples = N_c_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "theta[Treatment]"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# alpha_t is very uncertain because of non-identifiability with alpha_s.

N_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = N_nc_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "z_s[Treatment, Season]", 
                   "z_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "theta[Treatment]"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# Same as above.

# 3.5 Stan models ####
# I can do better by using sum_to_zero_vector for alpha_s and alpha_i because
# I know a priori that they shouldn't contribute to the Treatment effect, so
# should be identifiable.

N_c_s2z_model <- here("Biochemistry", "C_N", "Stan", "N_c_s2z.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

N_nc_s2z_model <- here("Biochemistry", "C_N", "Stan", "N_nc_s2z.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

N_c_s2z_samples <- N_c_s2z_model$sample(
          data = C_N %>%
            select(-C) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

N_nc_s2z_samples <- N_nc_s2z_model$sample(
          data = C_N %>%
            select(-C) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

# 3.6 Model checks ####
# 3.6.1 Rhat ####
N_c_s2z_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Most rhat above 1.001. rhat = 1.00 ± 0.00327.

N_nc_s2z_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Less than 1% of rhat above 1.001. rhat = 1.00 ± 0.000152.

# Plot comparison between centred and non-centred parameterisation.
N_nc_s2z_samples$summary() %>%
  left_join(N_c_s2z_samples$summary(),
            by = "variable") %>%
  rename(rhat_nc = rhat.x, rhat_c = rhat.y) %>%
  ggplot(aes(rhat_c, rhat_nc)) +
    geom_abline(slope = 1) +
    geom_point() +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Warning because z-scores are dropped as they have no equivalent in
# the centred parameteristion. The non-centred model is better.

# 3.6.2 Chains ####
N_c_s2z_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "N_c_s2z_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "C_N", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains are slightly better than before.

N_nc_s2z_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "N_nc_s2z_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "C_N", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains are great.

# 3.6.3 Pairs ####
N_c_s2z_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]", 
                      "theta[1]"))

N_c_s2z_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "theta[2]"))
# Correlation mostly gone.

N_nc_s2z_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]",
                      "theta[1]"))

N_nc_s2z_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "theta[2]"))
# Same as above.

# 3.7 Prior-posterior comparison ####
# 3.7.1 Sample priors ####
N_c_s2z_prior <- prior_samples(
  model = N_c_s2z_model,
  data = C_N %>%
    select(-C) %>%
    compose_data()
  )

N_nc_s2z_prior <- prior_samples(
  model = N_nc_s2z_model,
  data = C_N %>%
    select(-C) %>%
    compose_data()
)

# 3.7.2 Plot prior-posterior comparison ####
N_c_s2z_prior %>% 
  prior_posterior_draws(
    posterior_samples = N_c_s2z_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "theta[Treatment]"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# Looks fine but posteriors are a bit noisy.

N_nc_s2z_prior %>% 
  prior_posterior_draws(
    posterior_samples = N_nc_s2z_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "z_s[Treatment, Season]", 
                   "z_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "theta[Treatment]"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# Looks optimal.

# 3.8 Prediction ####
# 3.8.1 Combine relevant priors and posteriors ####
N_prior_posterior <- N_nc_s2z_prior %>% 
  prior_posterior_draws(
    posterior_samples = N_nc_s2z_samples,
    group = C_N %>%
      select(-C) %>%
      select(Treatment),
    parameters = c("alpha_t[Treatment]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "theta[Treatment]"),
    format = "short"
  )

# 3.8.2 Calculate predictions for new seasons and sporophytes ####
N_prior_posterior %<>%
  mutate(mu = exp( alpha_t ),
         obs = rgamma( n() , mu / theta , 1 / theta ),
         mu_new = exp( alpha_t + rnorm( n() , 0 , sigma_s ) + rnorm( n() , 0 , sigma_i ) ),
         obs_new = rgamma( n() , mu_new / theta , 1 / theta ))

# 3.8.3 Remove redundant prior ####
N_prior_posterior %<>% # priors are identical for both treatments ->
  filter(!(Treatment == "Faeces" & distribution == "prior")) %>% # remove one
  mutate(Treatment = if_else(distribution == "prior", # add Prior to treatment
                             "Prior", Treatment) %>% fct()) %>%
  select(-distribution)

# 3.8.3 Plot predictions ####
N_prior_posterior %>%
  pivot_longer(cols = c(mu, obs, mu_new, obs_new), 
               values_to = "N", names_to = "Level") %>%
  filter(Level %in% c("mu_new", "obs_new")) %>%
  ggplot(aes(N, Treatment, alpha = Level)) +
    geom_density_ridges(from = 0, to = 5) +
    scale_alpha_manual(values = c(0.8, 0.2)) +
    scale_x_continuous(limits = c(0, 5), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 3.8.4 Calculate difference ####
N_diff <- N_prior_posterior %>%
  filter(Treatment != "Prior") %>%
  droplevels() %>%
  select(-c(alpha_t, sigma_s, sigma_i, theta)) %>%
  pivot_wider(names_from = Treatment, values_from = c(mu, obs, mu_new, obs_new)) %>%
  mutate(mu = mu_Kelp - mu_Faeces, # calculate differences
         obs = obs_Kelp - obs_Faeces,
         mu_new = mu_new_Kelp - mu_new_Faeces,
         obs_new = obs_new_Kelp - obs_new_Faeces) %>%
  select(.chain, .iteration, .draw, mu, obs, mu_new, obs_new) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = "Parameter",
               values_to = "Difference")

# 3.8.5 Plot difference ####
N_diff %>%
  ggplot(aes(Difference, Parameter)) +
    geom_density_ridges(from = -1, to = 1) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-1, 1), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 3.8.6 Summarise difference ####
N_diff_summary <- N_diff %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P = mean(Difference > 0),
            n = length(Difference)) %T>%
  print()

# 3.8.7 Add labels to N_diff ####
N_diff %<>%
  left_join(N_diff_summary %>%
              select(Parameter, P), 
            by = "Parameter") %>%
  mutate(label_Kelp = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_Faeces = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))

# 4. Carbon-nitrogen ####
# 4.1 Descision not to model ####
# For Laminaria hyperborea, the carbon-nitrogen ratio is expected to be around 17
# (16.74 ± 0.46, Wright et al., 10.1111/gcb.16299) but ranges from 15 to 29
# across the North Atlantic (Filbee-Dexter et al., 10.1371/journal.pbio.3001702). 
# A ratio or quotient distribution is best described by the Cauchy distribution 
# if both original distributions were normal, but both original distributions
# are confined to [0, 100], so the ratio cannot be zero or negative. The closest
# liklihood is the F-distribution (ratio of two non-negtaive X^2 distributions),
# but this may not be the correct likelihood. Gamma would work but may not properly
# represent the ratio of two percentages. Since I know the underlying C and N
# likelihoods much better and have priors and posteriors for both, I'll just compute
# the prior and posterior ratio distributions and plot them over the data.

# 4.2 Prediction ####
# 4.2.1 Combine C and N priors and posteriors ####
C_N_prior_posterior <- C_prior_posterior %>%
  rename(C_mu = mu, C_obs = obs, C_mu_new = mu_new, C_obs_new = obs_new) %>%
  select(-c(alpha_t, sigma_s, sigma_i, nu)) %>%
  left_join(
    N_prior_posterior %>%
      rename(N_mu = mu, N_obs = obs, N_mu_new = mu_new, N_obs_new = obs_new) %>%
      select(-c(alpha_t, sigma_s, sigma_i, theta)),
    by = c("Treatment", ".chain", ".iteration", ".draw")
  ) %>%
  mutate(mu = C_mu / N_mu, # Calculate ratio distributions!
         obs = C_obs / N_obs,
         mu_new = C_mu_new / N_mu_new,
         obs_new = C_obs_new / N_obs_new) %>%
  select(Treatment, .chain, .iteration, .draw, mu, obs, mu_new, obs_new)
  
# 4.2.2 Plot predictions ####
C_N_prior_posterior %>%
  pivot_longer(cols = c(mu, obs, mu_new, obs_new), 
               values_to = "C_N", names_to = "Level") %>%
  filter(Level %in% c("mu_new", "obs_new")) %>%
  ggplot(aes(C_N, Treatment, alpha = Level)) +
    geom_density_ridges(from = 0, to = 100) +
    scale_alpha_manual(values = c(0.8, 0.2)) +
    scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks fairly gamma distributed.

# 4.2.3 Calculate difference ####
C_N_diff <- C_N_prior_posterior %>%
  filter(Treatment != "Prior") %>%
  droplevels() %>%
  pivot_wider(names_from = Treatment, values_from = c(mu, obs, mu_new, obs_new)) %>%
  mutate(mu = mu_Kelp - mu_Faeces, # calculate differences
         obs = obs_Kelp - obs_Faeces,
         mu_new = mu_new_Kelp - mu_new_Faeces,
         obs_new = obs_new_Kelp - obs_new_Faeces) %>%
  select(.chain, .iteration, .draw, mu, obs, mu_new, obs_new) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = "Parameter",
               values_to = "Difference")

# 4.2.4 Plot difference ####
C_N_diff %>%
  ggplot(aes(Difference, Parameter)) +
    geom_density_ridges(from = -50, to = 100) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-50, 100), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 4.2.5 Summarise difference ####
C_N_diff_summary <- C_N_diff %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P = mean(Difference > 0),
            n = length(Difference)) %T>%
  print()

# 4.2.6 Add labels to C_N_diff ####
C_N_diff %<>%
  left_join(C_N_diff_summary %>%
              select(Parameter, P), 
            by = "Parameter") %>%
  mutate(label_Kelp = ( P * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"),
         label_Faeces = ( (1 - P) * 100 ) %>% 
           signif(digits = 2) %>% 
           str_c("%"))

# 5. Visualisation ####
# No custom densities have to be calculated here because the data don't come
# with observational error, so can simply be plotted with geom_point.

# 5.1 Define custom theme ####
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

# 5.2 Carbon ####
Fig_2a_left_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = C, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.5, size = 5, height = 0.4) +
  stat_density_ridges(data = C_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 60, rel_min_height = 0.001, 
                      bandwidth = 1, scale = 2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 20),
                     oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#c3b300"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15), # circle, triangle, square
                     # override grey fill legend shapes because they can be confused with prior
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Carbon content (%)") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme
Fig_2a_left_top

require(geomtextpath)
Fig_2a_left_bottom <- C_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter,
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 1,
                      from = -50, to = 50, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 1,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 2,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 1,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 2,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -50, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = 0, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-50, 50), oob = scales::oob_keep,
                     breaks = seq(-50, 50, 25),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#c3b300", 0.6)),
                    guide = "none") +
  xlab("Difference (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme
Fig_2a_left_bottom

# 5.3 Nitrogen ####
Fig_2a_middle_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = N, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.5, size = 5, height = 0.4) +
  stat_density_ridges(data = N_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 3, rel_min_height = 0.001, 
                      bandwidth = 0.05, scale = 2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 3), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15),
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Nitrogen content (%)") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme
Fig_2a_middle_top

Fig_2a_middle_bottom <- N_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.05,
                      from = -2, to = 2, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 1,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 0.05, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 2,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 0.05, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 1,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.22, vjust = 0,
                   n = 2^10, bw = 0.05, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 2,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.22, vjust = 0,
                   n = 2^10, bw = 0.05, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -2, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = 0, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#c3b300", 0.6)),
                    guide = "none") +
  xlab("Difference (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme
Fig_2a_middle_bottom

# 5.4 Carbon-nitrogen ####
# Add C:N to data.
C_N %<>% mutate(C_N = C / N)

Fig_2a_right_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = C_N, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.5, size = 5, height = 0.4) +
  stat_density_ridges(data = C_N_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 120, rel_min_height = 0.001, 
                      bandwidth = 1.5, scale = 1.2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 120), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15),
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Carbon-nitrogen ratio") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme
Fig_2a_right_top

Fig_2a_right_bottom <- C_N_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 2,
                      from = -150, to = 150, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 1,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 2, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 2,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 2, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 1,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 2, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 2,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 2, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -150, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = 0, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-150, 150), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#c3b300", 0.6)),
                    guide = "none") +
  xlab("Difference") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme
Fig_2a_right_bottom

# 6. Save relevant data ####
C_N %>% write_rds(here("Biochemistry", "C_N", "RDS", "C_N.rds"))
C_prior_posterior %>% write_rds(here("Biochemistry", "C_N", "RDS", "C_prior_posterior.rds"))
C_diff %>% write_rds(here("Biochemistry", "C_N", "RDS", "C_diff.rds"))
N_prior_posterior %>% write_rds(here("Biochemistry", "C_N", "RDS", "N_prior_posterior.rds"))
N_diff %>% write_rds(here("Biochemistry", "C_N", "RDS", "N_diff.rds"))
C_N_prior_posterior %>% write_rds(here("Biochemistry", "C_N", "RDS", "C_N_prior_posterior.rds"))
C_N_diff %>% write_rds(here("Biochemistry", "C_N", "RDS", "C_N_diff.rds"))