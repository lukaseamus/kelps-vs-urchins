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
          # Separate distinct standard curves
          group_by(Treatment, Quenched) %>%
          mutate(Group = cur_group_id()) %>%
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
              facet_grid(~ Group) +
              ggtitle(Name) +
              theme_minimal()
        )
    ) %>%
  ungroup()

enzymes %$% 
  wrap_plots(Standard_Plot_Time) %>%
  ggsave(filename = "Standard_Plot_Time.pdf",
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
  ggsave(filename = "Standard_Plot_Concentration.pdf",
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
# around 20 to 100 thousand. beta is expected to be near Fmax / 150 µM, so around
# 133 to 667 fluorescence units µM^-1. Both need a gamma distribution. 
# A truncated normal likelihood ensures positive predictions but doesn't mess 
# with the nonlinear model like a gamma likelihood would.

require(truncnorm) # R doesn't have a built-in truncated normal distribution.
tibble(n = 1:1e3,
       F0 = rexp( 1e3 , 0.01 ),
       Fmax = rgamma( 1e3 , 8e4^2 / 4e4^2 , 8e4 / 4e4^2 ),
       beta = rgamma( 1e3 , 600^2 / 300^2 , 600 / 300^2 ),
       sigma = rexp( 1e3, 1e-3 )) %>%
  expand_grid(c = seq(0, 150, 2)) %>%
  mutate(mu_rh = F0 + Fmax * beta * c / ( Fmax + beta * c ),
         mu_es = F0 + Fmax * ( 1 - exp( -beta * c / Fmax ) ),
         mu_ht = F0 + Fmax * tanh( beta * c / Fmax ),
         mu_pl = if_else(c <= Fmax / beta, F0 + beta * c, F0 + Fmax),
         F_rh = rtruncnorm( n = n() , mean = mu_rh , sd = sigma , a = 0 ),
         F_es = rtruncnorm( n = n() , mean = mu_es , sd = sigma , a = 0 ),
         F_ht = rtruncnorm( n = n() , mean = mu_ht , sd = sigma , a = 0 ),
         F_pl = rtruncnorm( n = n() , mean = mu_pl , sd = sigma , a = 0 ),) %>%
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

# 2.3 Stan models ####
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

# 2.3.2 Run models ####
require(tidybayes)
enzymes %<>%
  mutate(
    Standard_rh_1_Samples = Standard_Data %>%
      map(
        ~ standard_rh_model$sample(
          data = .x %>%
            filter(Group == 1) %>%
            select(Concentration, Fluorescence) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_rh_2_Samples = Standard_Data %>%
      map(
        ~ standard_rh_model$sample(
          data = .x %>%
            filter(Group == 2) %>%
            select(Concentration, Fluorescence) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_es_1_Samples = Standard_Data %>%
      map(
        ~ standard_es_model$sample(
          data = .x %>%
            filter(Group == 1) %>%
            select(Concentration, Fluorescence) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_es_2_Samples = Standard_Data %>%
      map(
        ~ standard_es_model$sample(
          data = .x %>%
            filter(Group == 2) %>%
            select(Concentration, Fluorescence) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_ht_1_Samples = Standard_Data %>%
      map(
        ~ standard_ht_model$sample(
          data = .x %>%
            filter(Group == 1) %>%
            select(Concentration, Fluorescence) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_ht_2_Samples = Standard_Data %>%
      map(
        ~ standard_ht_model$sample(
          data = .x %>%
            filter(Group == 2) %>%
            select(Concentration, Fluorescence) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_pl_1_Samples = Standard_Data %>%
      map(
        ~ standard_pl_model$sample(
          data = .x %>%
            filter(Group == 1) %>%
            select(Concentration, Fluorescence) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_pl_2_Samples = Standard_Data %>%
      map(
        ~ standard_pl_model$sample(
          data = .x %>%
            filter(Group == 2) %>%
            select(Concentration, Fluorescence) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      )
  ) %T>%
  print(n = 33)

# 2.4 Model checks ####
# 2.4.1 Rhat ####
enzymes %<>%
  mutate(
    summary_rh_1 = Standard_rh_1_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_rh_1 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_rh_1 = mean(rhat),
                    rhat_sd_rh_1 = sd(rhat))
      ),
    summary_rh_2 = Standard_rh_2_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_rh_2 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_rh_2 = mean(rhat),
                    rhat_sd_rh_2 = sd(rhat))
      ),
    summary_es_1 = Standard_es_1_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_es_1 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_es_1 = mean(rhat),
                    rhat_sd_es_1 = sd(rhat))
      ),
    summary_es_2 = Standard_es_2_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_es_2 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_es_2 = mean(rhat),
                    rhat_sd_es_2 = sd(rhat))
      ),
    summary_ht_1 = Standard_ht_1_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_ht_1 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_ht_1 = mean(rhat),
                    rhat_sd_ht_1 = sd(rhat))
      ),
    summary_ht_2 = Standard_ht_2_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_ht_2 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_ht_2 = mean(rhat),
                    rhat_sd_ht_2 = sd(rhat))
      ),
    summary_pl_1 = Standard_pl_1_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_pl_1 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_pl_1 = mean(rhat),
                    rhat_sd_pl_1 = sd(rhat))
      ),
    summary_pl_2 = Standard_pl_2_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_pl_2 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_pl_2 = mean(rhat),
                    rhat_sd_pl_2 = sd(rhat))
      )
  ) %>%
  unnest(cols = c(summary_rh_1, summary_es_1, summary_ht_1, summary_pl_1,
                  summary_rh_2, summary_es_2, summary_ht_2, summary_pl_2))

enzymes %>% select(Name, ends_with("rh_1"))
enzymes %>% select(Name, ends_with("es_1"))
enzymes %>% select(Name, ends_with("ht_1"))
enzymes %>% select(Name, ends_with("pl_1"))
# No rhat above 1.001.
