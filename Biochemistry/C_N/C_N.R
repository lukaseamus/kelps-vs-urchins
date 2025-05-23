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
    by = c("ID", "Method")
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

tibble(n = 1:1e5,
       mu_logit = rnorm( 1e5 , log( 0.3 / (1 - 0.3) ) , 0.3 ),
       sigma = rexp( 1e5, 20 ),
       mu = 1 / ( 1 + exp(-mu_logit) ),
       nu = ( mu * (1 - mu) / sigma^2 - 1 ),
       C = rbeta( 1e5 , mu * nu , (1 - mu) * nu )) %>%
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
C_c_model <- here("Biochemistry", "C_N", "Stan", "C_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

C_nc_model <- here("Biochemistry", "C_N", "Stan", "C_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

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
# Most rhat above 1.001.

C_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# ~1% rhat above 1.001.

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
# Chains are better, but still not great. Chain 8 was
# dropped in both cases.

# 2.3.3 Pairs ####
C_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]", 
                      "sigma"))

C_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "sigma"))
# Correlation between alpha_t and alpha_s for both treatments.

C_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]",
                      "sigma"))

C_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "sigma"))
# Same as above.

# 4.6 Prior-posterior comparison ####
# 4.6.1 Sample priors ####
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

# 4.6.2 Plot prior-posterior comparison ####
C_c_prior %>% 
  prior_posterior_draws(
    posterior_samples = C_c_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "sigma"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# alpha_s and alpha_i are all clustered around zero. Generally looks fine.

C_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = C_nc_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "z_s[Treatment, Season]", 
                   "z_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "sigma"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# Same as above.

# 2.7 Stan models ####
# I can do better by using sum_to_zero_vector for alpha_s and alpha_i because
# I know a priori that they shouldn't contribute to the Treatment effect, so
# should be identifiable.

C_c_stz_model <- here("Biochemistry", "C_N", "Stan", "C_c_stz.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

C_nc_stz_model <- here("Biochemistry", "C_N", "Stan", "C_nc_stz.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

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
# Most rhat above 1.001.

C_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# ~1% rhat above 1.001.

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
# Chains are better, but still not great. Chain 8 was
# dropped in both cases.

# 2.3.3 Pairs ####
C_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]", 
                      "sigma"))

C_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "sigma"))
# Correlation between alpha_t and alpha_s for both treatments.

C_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]",
                      "sigma"))

C_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[2]",
                      "alpha_s[2,1]", "sigma_s[2]",
                      "alpha_i[2,1]", "sigma_i[2]",
                      "sigma"))
# Same as above.

# 4.6 Prior-posterior comparison ####
# 4.6.1 Sample priors ####
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

# 4.6.2 Plot prior-posterior comparison ####
C_c_prior %>% 
  prior_posterior_draws(
    posterior_samples = C_c_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "sigma"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# alpha_s and alpha_i are all clustered around zero. Generally looks fine.

C_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = C_nc_samples,
    group = C_N %>%
      select(Treatment, Season, Individual),
    parameters = c("alpha_t[Treatment]", "alpha_s[Treatment, Season]", 
                   "alpha_i[Treatment, Individual]", "z_s[Treatment, Season]", 
                   "z_i[Treatment, Individual]", "sigma_s[Treatment]", 
                   "sigma_i[Treatment]", "sigma"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Treatment", ridges = FALSE) %>%
      print() } %T>%
  { prior_posterior_plot(., group_name = "Season", ridges = FALSE) %>%
      print() } %>%
  prior_posterior_plot(group_name = "Individual", ridges = FALSE)
# Same as above.

