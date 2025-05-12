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
            filter(!str_detect(Annotation, "^[0-9.]+$")) %>% 
            rename(ID = Annotation) %>%
            mutate(ID = ID %>% fct())
          )
    ) %>%
  select(-Data)
phenol
phenol$Standard_Data[[3]]
phenol$Samples_Data[[2]]

# 2. Technical triplicate models ####
# 2.1 Prepare data ####
# I no longer want the plate-nested structure of phenol because my
# first set of models looks at each sample on an individual basis.
# Therefore I need nesting by ID, each of which contains the
# technical triplicate for one sample.

# Re-nest phenol.
phenol %<>%
  unnest(cols = Samples_Data) %>%
  group_by(Name, Date, Plate, Standard_Data, ID, Mass) %>%
  nest(.key = "Technical_Data")

phenol
phenol$Technical_Data[[30]]

# Calculate mean of each triplicate.
phenol %<>%
  mutate(
    Absorbance_mean = Technical_Data %>%
      map(
        ~ .x %$% mean(Absorbance)
      )
  )

# 3.2 Visualise data ####
phenol %<>%
  rowwise() %>%
  mutate(
    Technical_Plot = 
      list(
        Technical_Data %>%
          ggplot(aes(Well, Absorbance)) +
            geom_point() +
            geom_hline(yintercept = Absorbance_mean) +
            ggtitle(ID) +
            theme_minimal() +
            theme(panel.grid = element_blank())
        )
    ) %>%
  ungroup()

require(patchwork)
phenol %$% 
  wrap_plots(Technical_Plot) %>%
  ggsave(filename = "Technical_Data.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 60, width = 60, units = "cm")
# This technical variation needs to be modelled.

# 2.3 Prior simulation ####
# One thing I know to be consistent across samples is that they need
# to be positive and they are likely between 0 and 2.6, which is roughly
# the maximal absorbance for this chromophore. The intercept alpha
# can then be normally distributed around the log mean and the log link
# function together with the gamma likelihood ensure positivity. Due
# to the low number of technical replicates there are divergences
# if I set the prior on the grand mean.
phenol$Technical_Data %>% 
  bind_rows() %$% 
  mean(Absorbance) %>% 
  signif(1)

tibble(alpha = rnorm( 1e5 , log(0.8) , 1 ),
       mu = alpha %>% exp(),
       sigma = rexp( 1e5 , 5 ),
       A = rgamma( 1e5 , mu^2 / sigma^2 , mu / sigma^2 )) %>%
  ggplot(aes(y = A)) +
    geom_hline(yintercept = c(0, 2.6)) + # absorbances are between 0 and 2.6
    geom_density(orientation = "y", colour = NA, fill = "black", alpha = 0.5) +
    scale_y_continuous(limits = c(0, 5),
                       oob = scales::oob_keep) +
    coord_cartesian(expand = F, clip = "off",
                    ylim = c(0, 5)) +
    theme_minimal()

# Therefore, I'll include the calculated mean in the Stan model,
# and centre the prior on the log of that.

# 2.4 Run model ####
technical_stan <- "
    data{
      int n;
      vector<lower=0>[n] Absorbance;
      real<lower=0> Absorbance_mean;
      }

    parameters{
      real alpha; // Intercept in the log space
      real<lower=0> sigma; // Likelihood uncertainty
      }

    model{
      // Priors
      alpha ~ normal( log(Absorbance_mean) , 0.2 );
      sigma ~ exponential( 5 );

      // Model
      real mu;
      mu = exp( alpha ); // Log link function

      // Gamma likelihood
      Absorbance ~ gamma( square(mu) / square(sigma) , mu / square(sigma) );
      }
"

require(cmdstanr)
technical_model <- technical_stan %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
phenol %<>%
  mutate(
    Technical_Samples = Technical_Data %>%
      map2(
        Absorbance_mean, # list containing calculated means
        ~ technical_model$sample(
          data = .x %>%
            select(Absorbance) %>%
            compose_data() %>%
            list_modify(Absorbance_mean = .y), # here the mean is added
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
        )
  )
# Some divergent transitions, likely due to n = 3.

# 2.5 Model checks ####
# 2.5.1 Rhat ####
phenol %<>%
  mutate(
    summary = Technical_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean = mean(rhat),
                    rhat_sd = sd(rhat))
        )
  ) %>%
  unnest(cols = summary)

phenol %>% 
  select(Name, ID, rhat_1.001, rhat_mean, rhat_sd) %>%
  print(n = 91)
# Almost no rhat above 1.001, and even the few problematic cases have a mean
# rhat of 1.00 with small sd.

# 2.5.2 Chains ####
require(bayesplot)
phenol %<>%
  rowwise() %>%
  mutate(
    Technical_Chains = 
      list(
        Technical_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(Name)
        )
      ) %>%
  ungroup()

( phenol %$% 
  wrap_plots(Technical_Chains) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") ) %>%
  ggsave(filename = "Technical_Chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 80, width = 80, units = "cm")
# Chains look fine.

# 2.6 Prior-posterior comparison ####
# 2.6.1 Sample priors ####
phenol %<>%
  mutate(
    Technical_Prior = Absorbance_mean %>%
      map(
        ~ tibble(.chain = 1:8 %>% rep(each = 1e4),
                 .iteration = 1:1e4 %>% rep(times = 8),
                 .draw = 1:8e4,
                 alpha = rnorm( 8e4 , log(.x) , 0.2 ), # variable log mean
                 sigma = rexp( 8e4 , 5 ))
      )
  )

# 2.6.2 Extract posteriors ####
phenol %<>%
  mutate(
    Technical_Posterior = Technical_Samples %>%
      map(
        ~ .x %>% spread_draws(alpha, sigma)
      )
  )

# 2.6.3 Plot comparison ####
source("functions.R")
phenol %<>%
  rowwise() %>%
  mutate(
    Technical_Prior_Posterior =
      list(
          Technical_Prior %>%
            bind_rows(Technical_Posterior) %>%
            mutate(distribution = c("prior", "posterior") %>%
                     rep(each = 8e4) %>% fct()) %>%
            pivot_longer(cols = c(alpha, sigma),
                         names_to = ".variable",
                         values_to = ".value") %>%
            prior_posterior_plot() +
            ggtitle(Name)
      )
  ) %>%
  ungroup()

( phenol %$% 
  wrap_plots(Technical_Prior_Posterior) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") ) %>%
  ggsave(filename = "Technical_Prior_Posterior.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 80, width = 80, units = "cm")


# 2.7 Prediction ####
# 2.7.1 Calculate prediction ####
phenol %<>%
  mutate(
    Technical_Prediction = Technical_Posterior %>%
      map(
        ~ .x %>%
          mutate(mu = exp( alpha ),
                 obs = rgamma( n() , mu^2 / sigma^2 , mu / sigma^2 ))
      )
  )

# 2.7.2 Plot prediction ####
require(ggdist)
phenol %<>%
  rowwise() %>%
  mutate(
    Technical_Prediction_Plot = 
      list(
          ggplot() +
            geom_point(data = Technical_Data,
                       aes(Well, Absorbance)) +
            geom_hline(yintercept = Absorbance_mean) +
            stat_eye(data = Technical_Prediction,
                     aes(2, obs), alpha = 0.5,
                     point_interval = NULL) +
            scale_y_continuous(limits = Technical_Data %$% 
                                 c(min(Absorbance)*0.8, max(Absorbance)*1.2)) +
            ggtitle(ID) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  )

phenol %$% 
  wrap_plots(Technical_Prediction_Plot) %>%
  ggsave(filename = "Technical_Prediction.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 60, width = 60, units = "cm")

# 2. Standard curve models ####
# 2.1 Visualise data ####
phenol %<>%
  rowwise() %>% # rowwise allows much easier plotting syntax
  mutate(
    Standard_Plot = 
      list(
        Standard_Data %>%
            ggplot(aes(Concentration, Absorbance)) +
              geom_point() +
              geom_smooth(se = FALSE, colour = "grey") +
              ggtitle(Name) +
              theme_minimal()
        )
    ) %>%
  ungroup() # undo rowwise

require(patchwork)
phenol %$% wrap_plots(Standard_Plot)
# The limit of linearity for this assay is assumed to be 1 mg mL^-1, 
# but this suggests that saturation already happens earlier.

# 2.2 Prior simulation ####
# There are three saturating models that can be tested alongside the linear:
# rectangular hyperbola, exponential saturation, and hyperbolic tangent.
# A0 (absorbance intercept), Amax (absorbance maximum) and beta (linear 
# relationship between concentration and absorbance) are necessarily positive. 
# A0 is expected to be very near zero but no other information is available, 
# so an exponential distribution is best. Amax is expected above 2.5 a.u. and
# beta is expected to be near 2.5 a.u. divided by 1 mg mL^1, so 2.5.
# Both need a gamma distribution. A truncated normal likelihood
# ensures positive predictions but doesn't mess with the nonlinear model like 
# a different positive-values-only likelihood such as gamma would.

require(truncnorm) # R doesn't have a built-in truncated normal distribution.
tibble(n = 1:1e3,
       A0 = rexp( 1e3 , 10 ),
       Amax = rgamma( 1e3 , 3^2 / 2^2 , 3 / 2^2 ),
       beta = rgamma( 1e3 , 2.5^2 / 1.5^2 , 2.5 / 1.5^2 ),
       sigma = rexp( 1e3, 3 )) %>%
  expand_grid(c = seq(0, 1, 0.05)) %>%
  mutate(mu_lm = A0 + beta * c,
         mu_rh = A0 + Amax * beta * c / ( Amax + beta * c ),
         mu_es = A0 + Amax * ( 1 - exp( -beta * c / Amax ) ),
         mu_ht = A0 + Amax * tanh( beta * c / Amax ),
         A_lm = rtruncnorm( n = n() , mean = mu_lm , sd = sigma , a = 0 ),
         A_rh = rtruncnorm( n = n() , mean = mu_rh , sd = sigma , a = 0 ),
         A_es = rtruncnorm( n = n() , mean = mu_es , sd = sigma , a = 0 ),
         A_ht = rtruncnorm( n = n() , mean = mu_ht , sd = sigma , a = 0 )) %>%
  pivot_longer(cols = c(starts_with("mu"), starts_with("A_")),
               names_to = "par_mod",
               values_to = "value") %>%
  separate(par_mod, into = c("par", "mod")) %>%
  mutate(mod = mod %>% fct_relevel("lm", "rh", "es"),
         par = par %>% fct_relevel("mu")) %>%
  ggplot(aes(c, value, group = n)) +
    geom_hline(yintercept = c(0, 2.5)) +
    geom_line(alpha = 0.01) +
    facet_grid(par ~ mod) +
    coord_cartesian(expand = F,
                    ylim = c(-1, 5)) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Covers all reasonable possibilities.

# 2.3 Stan models ####
# 2.3.1 Linear model ####
standard_lm_stan <- "
    data{ 
      int n;
      vector<lower=0>[n] Absorbance;
      vector<lower=0>[n] Concentration;
      }

    parameters{
      real<lower=0> A0;
      real<lower=0> beta;
      real<lower=0> sigma;
      }

    model{
      // Priors
      A0 ~ exponential( 10 );
      beta ~ gamma( square(2.5) / square(1.5) , 2.5 / square(1.5) );
      sigma ~ exponential( 3 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + beta * Concentration[i];
      }

      // Truncated normal likelihood
      Absorbance ~ normal( mu , sigma ) T[0,];
      }
"

require(cmdstanr)
standard_lm_model <- standard_lm_stan %>%
  write_stan_file() %>%
  cmdstan_model()

# 2.3.2 Rectangular hyperbola ####
standard_rh_stan <- "
    data{ 
      int n;
      vector<lower=0>[n] Absorbance;
      vector<lower=0>[n] Concentration;
      }

    parameters{
      real<lower=0> A0;
      real<lower=0> Amax;
      real<lower=0> beta;
      real<lower=0> sigma;
      }

    model{
      // Priors
      A0 ~ exponential( 10 );
      Amax ~ gamma( square(3) / square(2) , 3 / square(2) );
      beta ~ gamma( square(2.5) / square(1.5) , 2.5 / square(1.5) );
      sigma ~ exponential( 3 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + Amax * beta * Concentration[i] / 
                ( Amax + beta * Concentration[i] );
      }

      // Truncated normal likelihood
      Absorbance ~ normal( mu , sigma ) T[0,];
      }
"

standard_rh_model <- standard_rh_stan %>%
  write_stan_file() %>%
  cmdstan_model()

# 2.3.3 Exponential saturation ####
standard_es_stan <- "
    data{ 
      int n;
      vector<lower=0>[n] Absorbance;
      vector<lower=0>[n] Concentration;
      }

    parameters{
      real<lower=0> A0;
      real<lower=0> Amax;
      real<lower=0> beta;
      real<lower=0> sigma;
      }

    model{
      // Priors
      A0 ~ exponential( 10 );
      Amax ~ gamma( square(3) / square(2) , 3 / square(2) );
      beta ~ gamma( square(2.5) / square(1.5) , 2.5 / square(1.5) );
      sigma ~ exponential( 3 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + Amax * ( 1 - exp( -beta * Concentration[i] / Amax ) );
      }

      // Truncated normal likelihood
      Absorbance ~ normal( mu , sigma ) T[0,];
      }
"

standard_es_model <- standard_es_stan %>%
  write_stan_file() %>%
  cmdstan_model()

# 2.3.4 Hyperbolic tangent ####
standard_ht_stan <- "
    data{ 
      int n;
      vector<lower=0>[n] Absorbance;
      vector<lower=0>[n] Concentration;
      }

    parameters{
      real<lower=0> A0;
      real<lower=0> Amax;
      real<lower=0> beta;
      real<lower=0> sigma;
      }

    model{
      // Priors
      A0 ~ exponential( 10 );
      Amax ~ gamma( square(3) / square(2) , 3 / square(2) );
      beta ~ gamma( square(2.5) / square(1.5) , 2.5 / square(1.5) );
      sigma ~ exponential( 3 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + Amax * tanh( beta * Concentration[i] / Amax );
      }

      // Truncated normal likelihood
      Absorbance ~ normal( mu , sigma ) T[0,];
      }
"

standard_ht_model <- standard_ht_stan %>%
  write_stan_file() %>%
  cmdstan_model()

# 2.3.5 Run all models ####
require(tidybayes)
phenol %<>%
  mutate(
    Standard_lm_Samples = Standard_Data %>%
      map(
        ~ standard_lm_model$sample(
          data = .x %>%
            select(Absorbance, Concentration) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
        ),
    Standard_rh_Samples = Standard_Data %>%
      map(
        ~ standard_rh_model$sample(
          data = .x %>%
            select(Absorbance, Concentration) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_es_Samples = Standard_Data %>%
      map(
        ~ standard_es_model$sample(
          data = .x %>%
            select(Absorbance, Concentration) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
      ),
    Standard_ht_Samples = Standard_Data %>%
      map(
        ~ standard_ht_model$sample(
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

phenol

# 2.4 Model checks ####
# 2.4.1 Rhat ####
phenol %<>%
  mutate(
    summary_lm = Standard_lm_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_lm = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_lm = mean(rhat),
                    rhat_sd_lm = sd(rhat))
        ),
    summary_rh = Standard_rh_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_rh = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_rh = mean(rhat),
                    rhat_sd_rh = sd(rhat))
      ),
    summary_es = Standard_es_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_es = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_es = mean(rhat),
                    rhat_sd_es = sd(rhat))
      ),
    summary_ht = Standard_ht_Samples %>%
      map(
        ~ .x$summary() %>%
          mutate(rhat_check = rhat > 1.001) %>%
          summarise(rhat_1.001_ht = sum(rhat_check) / length(rhat), # proportion > 1.001
                    rhat_mean_ht = mean(rhat),
                    rhat_sd_ht = sd(rhat))
      )
  ) %>%
  unnest(cols = c(summary_lm, summary_rh, summary_es, summary_ht))

phenol %>% select(Name, ends_with("lm"))
phenol %>% select(Name, ends_with("rh"))
phenol %>% select(Name, ends_with("_es"))
phenol %>% select(Name, ends_with("ht"))
# No rhat above 1.001.

# 2.4.2 Chains ####
require(bayesplot)
phenol %<>%
  rowwise() %>%
  mutate(
    Standard_lm_Chains = 
      list(
        Standard_lm_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(Name) # and adding titles
      ),
    Standard_rh_Chains =
      list(
        Standard_rh_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(Name) # and adding titles
      ),
    Standard_es_Chains =
      list(
        Standard_es_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(Name) # and adding titles
      ),
    Standard_ht_Chains =
      list(
        Standard_ht_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(Name) # and adding titles
      )
    ) %>%
  ungroup()

phenol %$% wrap_plots(Standard_lm_Chains)
phenol %$% wrap_plots(Standard_rh_Chains)
phenol %$% wrap_plots(Standard_es_Chains)
phenol %$% wrap_plots(Standard_ht_Chains)
# Chains look good.

# 2.4.3 Pairs ####
phenol %<>%
  mutate(
    Standard_lm_Pairs = Standard_lm_Samples %>%
      map(
        ~ .x$draws(format = "df") %>% # mcmc_pairs does not allow adding titles
          mcmc_pairs(pars = c("A0", "beta"))
      ),
    Standard_rh_Pairs = Standard_rh_Samples %>%
      map(
        ~ .x$draws(format = "df") %>% # mcmc_pairs does not allow adding titles
          mcmc_pairs(pars = c("A0", "Amax", "beta"))
      ),
    Standard_es_Pairs = Standard_es_Samples %>%
      map(
        ~ .x$draws(format = "df") %>% # mcmc_pairs does not allow adding titles
          mcmc_pairs(pars = c("A0", "Amax", "beta"))
      ),
    Standard_ht_Pairs = Standard_ht_Samples %>%
      map(
        ~ .x$draws(format = "df") %>% # mcmc_pairs does not allow adding titles
          mcmc_pairs(pars = c("A0", "Amax", "beta"))
      )
  )

phenol %$% wrap_plots(Standard_lm_Pairs)
phenol %$% wrap_plots(Standard_rh_Pairs)
phenol %$% wrap_plots(Standard_es_Pairs)
phenol %$% wrap_plots(Standard_ht_Pairs)
# Some correlation between A0 and beta in the linear model
# and Amax and beta in the other models, but not concerning.

# 2.5 Prior-posterior comparison ####
source("functions.R")
# 2.5.1 Sample prior ####
phenol %<>%
  mutate(
    Standard_lm_Prior = Standard_Data %>%
      map(
        ~ prior_samples(model = standard_lm_model,
                        data = .x %>%
                          select(Absorbance, Concentration) %>%
                          compose_data())
      ),
    Standard_rh_Prior = Standard_Data %>%
      map(
        ~ prior_samples(model = standard_rh_model,
                        data = .x %>%
                          select(Absorbance, Concentration) %>%
                          compose_data())
      ),
    Standard_es_Prior = Standard_Data %>%
      map(
        ~ prior_samples(model = standard_es_model,
                        data = .x %>%
                          select(Absorbance, Concentration) %>%
                          compose_data())
      ),
    Standard_ht_Prior = Standard_Data %>%
      map(
        ~ prior_samples(model = standard_ht_model,
                        data = .x %>%
                          select(Absorbance, Concentration) %>%
                          compose_data())
      )
  )

# 2.5.2 Combine prior and posterior ####
phenol %<>%
  mutate(
    Standard_lm_Prior_Posterior = Standard_lm_Prior %>%
      map2(Standard_lm_Samples,
        ~ prior_posterior_draws(prior_samples = .x,
                                posterior_samples = .y,
                                parameters = c("A0", "beta", "sigma"),
                                format = "short")
      ),
    Standard_rh_Prior_Posterior = Standard_rh_Prior %>%
      map2(Standard_rh_Samples,
           ~ prior_posterior_draws(prior_samples = .x,
                                   posterior_samples = .y,
                                   parameters = c("A0", "Amax", "beta", "sigma"),
                                   format = "short")
      ),
    Standard_es_Prior_Posterior = Standard_es_Prior %>%
      map2(Standard_es_Samples,
           ~ prior_posterior_draws(prior_samples = .x,
                                   posterior_samples = .y,
                                   parameters = c("A0", "Amax", "beta", "sigma"),
                                   format = "short")
      ),
    Standard_ht_Prior_Posterior = Standard_ht_Prior %>%
      map2(Standard_ht_Samples,
           ~ prior_posterior_draws(prior_samples = .x,
                                   posterior_samples = .y,
                                   parameters = c("A0", "Amax", "beta", "sigma"),
                                   format = "short")
      )
  )

# 2.5.3 Plot comparison ####
phenol %<>%
  rowwise() %>%
  mutate(
    Standard_lm_Prior_Posterior_Plot =
      list(
          prior_posterior_draws(prior_samples = Standard_lm_Prior,
                                posterior_samples = Standard_lm_Samples,
                                parameters = c("A0", "beta", "sigma"),
                                format = "long") %>%
          prior_posterior_plot() +
          ggtitle(Name)
      ),
    Standard_rh_Prior_Posterior_Plot =
      list(
        prior_posterior_draws(prior_samples = Standard_rh_Prior,
                              posterior_samples = Standard_rh_Samples,
                              parameters = c("A0", "Amax", "beta", "sigma"),
                              format = "long") %>%
          prior_posterior_plot() +
          ggtitle(Name)
      ),
    Standard_es_Prior_Posterior_Plot =
      list(
        prior_posterior_draws(prior_samples = Standard_es_Prior,
                              posterior_samples = Standard_es_Samples,
                              parameters = c("A0", "Amax", "beta", "sigma"),
                              format = "long") %>%
          prior_posterior_plot() +
          ggtitle(Name)
      ),
    Standard_ht_Prior_Posterior_Plot =
      list(
        prior_posterior_draws(prior_samples = Standard_ht_Prior,
                              posterior_samples = Standard_ht_Samples,
                              parameters = c("A0", "Amax", "beta", "sigma"),
                              format = "long") %>%
          prior_posterior_plot() +
          ggtitle(Name)
      )
  ) %>%
  ungroup()

phenol %$% 
  wrap_plots(Standard_lm_Prior_Posterior_Plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

phenol %$% 
  wrap_plots(Standard_rh_Prior_Posterior_Plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

phenol %$% 
  wrap_plots(Standard_es_Prior_Posterior_Plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

phenol %$% 
  wrap_plots(Standard_ht_Prior_Posterior_Plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# 2.6 Prediction ####
# 2.6.1 Calculate prediction ####
phenol %<>%
  mutate(
    Standard_lm_Prediction = Standard_lm_Prior_Posterior %>%
      map2(
        Standard_Data,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration") %>%
          mutate(mu = A0 + beta * Concentration,
                 obs = rtruncnorm( n = n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration) %>%
          reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                  obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
          unnest(c(mu, obs), names_sep = "_")
      ),
    Standard_rh_Prediction = Standard_rh_Prior_Posterior %>%
      map2(
        Standard_Data,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration") %>%
          mutate(mu = A0 + Amax * beta * Concentration / 
                      ( Amax + beta * Concentration ),
                 obs = rtruncnorm( n = n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration) %>%
          reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                  obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
          unnest(c(mu, obs), names_sep = "_")
      ),
    Standard_es_Prediction = Standard_es_Prior_Posterior %>%
      map2(
        Standard_Data,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration") %>%
          mutate(mu = A0 + Amax * ( 1 - exp( -beta * Concentration / Amax ) ),
                 obs = rtruncnorm( n = n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration) %>%
          reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                  obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
          unnest(c(mu, obs), names_sep = "_")
      ),
    Standard_ht_Prediction = Standard_ht_Prior_Posterior %>%
      map2(
        Standard_Data,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration") %>%
          mutate(mu = A0 + Amax * tanh( beta * Concentration / Amax ),
                 obs = rtruncnorm( n = n() , mean = mu , sd = sigma , a = 0 )) %>%
          group_by(distribution, Concentration) %>%
          reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                  obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
          unnest(c(mu, obs), names_sep = "_")
      )
  )

# 2.6.2 Plot prediction ####
phenol %<>%
  rowwise() %>% 
  mutate(
    Standard_lm_Prediction_Plot =
      list(
        Standard_lm_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, aes(Concentration, Absorbance)) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = . %>% filter(distribution == "posterior"),
            #             aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
            #                 alpha = factor(obs_.width))) + # unhash to see observations
            geom_ribbon(data = . %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      ),
    Standard_rh_Prediction_Plot =
      list(
        Standard_rh_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, aes(Concentration, Absorbance)) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = . %>% filter(distribution == "posterior"),
            #             aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
            #                 alpha = factor(obs_.width))) + # unhash to see observations
            geom_ribbon(data = . %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      ),
    Standard_es_Prediction_Plot =
      list(
        Standard_es_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, aes(Concentration, Absorbance)) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = . %>% filter(distribution == "posterior"),
            #             aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
            #                 alpha = factor(obs_.width))) + # unhash to see observations
            geom_ribbon(data = . %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      ),
    Standard_ht_Prediction_Plot =
      list(
        Standard_ht_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, aes(Concentration, Absorbance)) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = . %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = . %>% filter(distribution == "posterior"),
            #             aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
            #                 alpha = factor(obs_.width))) + # unhash to see observations
            geom_ribbon(data = . %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  ) %>%
  ungroup()

phenol %$% wrap_plots(Standard_lm_Prediction_Plot)
phenol %$% wrap_plots(Standard_rh_Prediction_Plot)
phenol %$% wrap_plots(Standard_es_Prediction_Plot)
phenol %$% wrap_plots(Standard_ht_Prediction_Plot)
# The linear models generally overestimate absorbance somewhat for 
# concentrations near 1 mg mL^1. The saturating models all seem to 
# fit better but it is unclear which fits best.

# Plot only the mean of mu for all models.
phenol %<>%
  rowwise() %>% 
  mutate(
    Standard_Prediction = 
      list(
        bind_rows(Standard_lm_Prediction,
                  Standard_rh_Prediction,
                  Standard_es_Prediction,
                  Standard_ht_Prediction) %>%
          mutate(model = rep(c("lm", "rh", "es", "ht"), each = n() / 4))
        ),
    Standard_Prediction_Plot =
      list(
        Standard_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, aes(Concentration, Absorbance)) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y, colour = model)) +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  )

phenol %$% wrap_plots(Standard_Prediction_Plot)

# 2.6.3 Inverse prediction ####
# To predict concentrations from absorbance I need inverse prediction,
# i.e. solving for x. Given y = alpha + beta * x this becomes 
# x = ( y - alpha ) / beta, i.e. Concentration = ( Absorbance - alpha ) / beta.
# Simply solving for x with alpha and beta posteriors is equivalent to an
# inverse prediction of the mean, but I want to predict observations in x.
# Since the likelihood, an therefore also the observational error, is normal
# this should be possible. 

# Here are are some example values:
phenol$Standard_Samples[[2]]
alpha <- rnorm( 1e5 , 0.14 , 0.05 )
beta <- rnorm( 1e5 , 2.46 , 0.09 )
sigma <- rnorm( 1e5 , 0.15 , 0.03 )
Concentration <- 0.5
Absorbance <- 1

# Re-expressing Absorbance ~ normal( alpha + beta * Concentration , sigma )
# as Absorbance ~ alpha + beta * Concentration + normal( 0 , sigma ), allows
# for easier rearrangement of terms when solving for Concentration.

# Here's proof that these expressions are equivalent:
set.seed(1)
rnorm( 1e5 , alpha + beta * Concentration , sigma ) %>% hist()
set.seed(1)
( alpha + beta * Concentration + rnorm( 1e5 , 0 , sigma ) ) %>% hist()

# Now that I know I can use the rearranged form 
# Absorbance ~ alpha + beta * Concentration + normal( 0 , sigma ),
# I can rearrange the observational error along with the parameters
# and get 
# Concentration ~ ( Absorbance - alpha - normal( 0 , sigma ) ) / beta.
# This adds the additional uncertainty to our Concentration estimates:
( ( Absorbance - alpha ) / beta ) %>% hist()
set.seed(1)
( ( Absorbance - alpha - rnorm( 1e5 , 0 , sigma ) ) / beta ) %>% hist()

# Here are some equivalent expressions:
set.seed(1)
( ( Absorbance - rnorm( 1e5 , alpha , sigma ) ) / beta ) %>% hist()
set.seed(1)
( ( Absorbance - alpha ) / beta - rnorm( 1e5 , 0 , sigma ) / beta ) %>% hist()
set.seed(1)
( ( Absorbance - alpha ) / beta + rnorm( 1e5 , 0 , sigma ) / - beta ) %>% hist()

# But reversing the first step and enclosing everything in the likelihood
# function is not exactly equivalent in practice because it requires a change 
# of sign:
set.seed(1)
( ( Absorbance - alpha ) / beta + rnorm( 1e5 , 0 , sigma ) / beta ) %>% hist()
# which allows
set.seed(1)
rnorm( 1e5 , ( Absorbance - alpha ) / beta , sigma / beta ) %>% hist()
# So even though adding rather than subtracting a normal with mean
# zero is equivalent in theory, it is not exactly equivalent in practice, 
# and we'd best go with one of the former expressions. I'll go with the
# least derived: ( Absorbance - alpha - rnorm( 1e5 , 0 , sigma ) ) / beta

# 2.6.4 Calculating sample concentration ####
# I need a new list of tibbles where the prior estimates and sample
# data are joined.
phenol$Standard_Prior_Posterior[[1]]
phenol$Samples_Data[[1]]

phenol %<>%
  mutate(
    Samples_Data_Converted = Samples_Data %>%
      map2(
        Standard_Prior_Posterior,
        ~ .x %>%
          cross_join(
            .y %>% 
              filter(distribution == "posterior") %>%
              select(-distribution)
            ) %>% 
          mutate( # Here's where the magic happens!
            Concentration = ( Absorbance - alpha - rnorm( n() , 0 , sigma ) ) / beta
            )
      )
  )

phenol$Samples_Data_Converted[[3]]
# Success! It's not a problem that there are lost of negative values in this
# estimated distribution Concentration because this will be constrained to
# positive values only in the next model. 



# Calculate summary.
phenol %<>%
  mutate(
    Technical_Data_Summary = Technical_Data %>%
      map(
        ~ .x %>% 
          group_by(Well, Content, Absorbance) %>%
          summarise(Concentration_mean = mean(Concentration),
                    Concentration_sd = sd(Concentration)) %>%
          ungroup()
      )
  )

"
    data{ 
      int n;
      vector[n] Concentration_mean;
      vector<lower=0>[n] Concentration_sd;
      }

    parameters{
      vector[n] Concentration; // True estimates for measurement error
      real alpha; // Intercept in the log space
      real<lower=0> sigma; // Likelihood uncertainty
      }

    model{
      // Priors
      alpha ~ normal( log(0.5) , 0.8 );
      sigma ~ exponential( 1 );

      // Model
      real mu;
      mu = exp( alpha ); // Log link function

      // Gamma likelihood with normal measurement error (note elementwise division)
      Concentration ~ gamma( square(mu) / square(sigma) , mu / square(sigma) );
      Concentration_mean ~ gamma( square(Concentration) ./ square(Concentration_sd) , 
                                  Concentration ./ square(Concentration_sd) );
      }
"

# Model fails because non-negativity in Concentration_mean is not allowed.