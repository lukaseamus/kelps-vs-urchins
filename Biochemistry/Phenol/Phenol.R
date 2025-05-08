# 1. Load data ####
# 1.1 Load raw data ####
require(tidyverse)
require(here)

phenol <- here("Biochemistry", "Phenol", "Raw") %>%
  list.files(pattern = "\\.csv$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Name = Path %>% basename() %>% 
      str_remove("\\..*$"),
    Data = Path %>% 
      map(~ read_csv(.x, skip = 10, col_select = 1:3, col_types = list("f", "f")) %>%
              rename(Absorbance = "Raw Data (765)")),
    Date = Name %>% str_remove("^X") %>% str_split_i("_", 1) %>% ymd(),
    Plate = Name %>% str_split_i("_", 3) %>% as.numeric()
  ) %>%
  select(-Path)

phenol

# 1.2 Load metadata ####
meta <- here("Biochemistry", "Phenol", "Phenol_Meta.csv") %>% 
  read_csv(col_types = list("f", "f")) %>%
  group_by(Plate) %>%
  nest(.key = "Metadata")

meta

# 1.3 Join metadata to data ####
require(magrittr)
phenol %<>%
  full_join(meta %>% rename(Name = Plate), by = "Name")
phenol

phenol %<>% 
  mutate(Data = map2(Data, Metadata, 
                     ~ full_join(.x, .y, by = "Content"))) %>%
  select(-Metadata)
phenol
phenol$Data[[3]]

# 1.4 Separate standard and samples ####
phenol %<>% 
  mutate(
    Standard_Data = Data %>% 
      map(
          ~ .x %>%
            filter(str_detect(Annotation, "^[0-9.]+$")) %>% 
            rename(Concentration = Annotation) %>%
            mutate(Concentration = as.numeric(Concentration)) %>%
            select(-Mass) # no sample mass for standards
          ), 
    Samples_Data = Data %>%
      map(
          ~ .x %>% 
            filter(!str_detect(Annotation, "^[0-9.]+$"))
          )
    ) %>%
  select(-Data)
phenol
phenol$Standard[[3]]
phenol$Samples[[2]]

# 2. Standard curve models ####
# 2.1 Visualise data ####
phenol %<>%
  mutate(
    Standard_Plot = Standard_Data %>%
      map(~ .x %>% 
            ggplot(aes(Concentration, Absorbance)) +
              geom_point() +
              theme_minimal()
          )
    )

require(patchwork)
phenol %$% wrap_plots(Standard_Plot)
# Pretty linear. In reality there is a lower bound at 0 and an upper bound
# around the maximal absorbance (2.5 a.u.). But since we know the limit of 
# linearity to be around 1 mg mL^-1, we can approximate with a linear model.

# 2.2 Prior simulation ####
# alpha (absorbance) and beta (relationship between concentration and absorbance)
# are necessarily positive, so exponential and gamma distributions are best. 
# alpha is expected to be very near zero but no other information is available. 
# beta is expected to be near the absorbance maximum (2.5 a.u.) divided by the 
# maximal concentration (1 mg mL^1), so 2.5.

tibble(n = 1:1e3,  
       alpha = rexp(n = 1e3, rate = 10),
       beta = rgamma(n = 1e3, shape = 2.5^2 / 1^2, rate = 2.5 / 1^2)) %>%
  expand_grid(c = seq(0, 1)) %>%
  mutate(A = alpha + beta * c) %>%
  ggplot(aes(c, A, group = n)) +
    geom_hline(yintercept = c(0, 2.5)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()
# Covers all reasonable possibilities.

# 2.3 Run model ####
standard_stan <- "
    data{ 
      int n;
      vector<lower=0>[n] Absorbance;
      vector<lower=0>[n] Concentration;
      }

    parameters{
      real<lower=0> alpha;
      real<lower=0> beta;
      real<lower=0> sigma;
      }

    model{
      // Priors
      alpha ~ exponential( 10 );
      beta ~ gamma( 2.5^2 / 1^2, 2.5 / 1^2 );
      sigma ~ exponential( 1 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = alpha + beta * Concentration[i];
      }

      // Truncated likelihood
      Absorbance ~ normal( mu , sigma ) T[0,];
      }
"


require(cmdstanr)
standard_model <- standard_stan %>%
  write_stan_file() %>%
  cmdstan_model()
  
require(tidybayes)
phenol %<>%
  mutate(
    Standard_Samples = Standard_Data %>%
      map(
        ~ standard_model$sample(
          data = .x %>%
            select(Absorbance, Concentration) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
        )
  )

# 2.4 Model checks ####
# 2.4.1 Rhat ####
phenol %<>%
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
  unnest(summary)

phenol
# No rhat above 1.001.

# 2.4.2 Chains ####
require(bayesplot)
phenol %<>%
  mutate(
    Standard_Chains = Standard_Samples %>%
      map(
        ~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()
      )
    )

phenol %$% wrap_plots(Standard_Chains)
# Chains look good.

# 2.4.3 Pairs ####
phenol %<>%
  mutate(
    Standard_Pairs = Standard_Samples %>%
      map(
        ~ .x$draws(format = "df") %>%
          mcmc_pairs(pars = c("alpha", "beta"))
      )
  )

phenol %$% wrap_plots(Standard_Pairs)
# Some correlation between alpha and beta, but not concerning.

# 2.5 Prior-posterior comparison ####
source("functions.R")
# 2.5.1 Sample prior ####
phenol %<>%
  mutate(
    Standard_Prior = Standard_Data %>%
      map(
        ~ prior_samples(model = standard_model,
                        data = .x %>%
                          select(Absorbance, Concentration) %>%
                          compose_data())
      )
  )

# 2.5.2 Combine prior and posterior ####
phenol %<>%
  mutate(
    Standard_Prior_Posterior = Standard_Prior %>%
      map2(Standard_Samples,
        ~ prior_posterior_draws(prior_samples = .x,
                                posterior_samples = .y,
                                parameters = c("alpha", "beta", "sigma"),
                                format = "short")
      )
  )

phenol$Standard_Prior_Posterior[[3]]

# 2.5.3 Plot comparison ####
phenol %<>%
  mutate(
    Standard_Prior_Posterior_Plot = Standard_Prior %>%
      map2(
        Standard_Samples,
        ~ prior_posterior_draws(prior_samples = .x,
                                posterior_samples = .y,
                                parameters = c("alpha", "beta", "sigma"),
                                format = "long") %>%
          prior_posterior_plot()
      )
  )

phenol %$% 
  wrap_plots(Standard_Prior_Posterior_Plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
  
# 2.6 Prediction ####
# 2.6.1 Calculate prediction ####
require(truncnorm) # R doesn't have a native truncated normal function
phenol %<>%
  mutate(
    Standard_Prediction = Standard_Prior_Posterior %>%
      map2(
        Standard_Data,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration") %>%
          mutate(mu = alpha + beta * Concentration,
                 obs = rtruncnorm( n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration) %>%
          reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                  obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
          unnest(c(mu, obs), names_sep = "_")
      )
  )

# 2.6.2 Plot prediction ####
phenol %<>%
  rowwise() %>% # rowwise allows much easier plotting syntax
  mutate(
    Standard_Prediction_Plot =
      list(
        Standard_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, aes(Concentration, Absorbance)) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
                            alpha = factor(obs_.width))) +
            geom_ribbon(data = . %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  ) 

phenol %$% wrap_plots(Standard_Prediction_Plot)

