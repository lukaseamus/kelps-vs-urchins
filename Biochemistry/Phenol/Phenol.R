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

# 1.5 Filter out samples that exceed the standard curve ####
# Some replicates exceed the standard curve, which cannot be extrapolated.
# In this case the whole sample was diluted and run again. I need to 
# remove all samples which exceed the standard curve, even if only in
# one technical replicate.

# Example:
phenol$Standard_Data[[5]] %>%
  print(n = 21) # Maximum standard curve Absorbance is 2.6 a.u., 
# and mean max is
phenol$Standard_Data[[5]] %>%
  group_by(Concentration) %>%
  summarise(Absorbance = mean(Absorbance)) %$%
  max(Absorbance) # 2.502 a.u.

phenol$Samples_Data[[5]] %>%
  print(n = 36) # Several samples exceed this (e.g. U9_A1)

phenol %<>%
  mutate(
    Samples_Data = Samples_Data %>%
      map2(
        Standard_Data,
        ~ .x %>%
          group_by(ID) %>%
          filter( 
              mean(Absorbance) <= .y %>%
                group_by(Concentration) %>%
                summarise(Absorbance = mean(Absorbance)) %$% 
                max(Absorbance) 
            ) %>%
          ungroup()
      )
  )

# Test:
phenol$Samples_Data[[5]] %>%
  print(n = 36) # All remaining samples are within standard curve.

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
  rowwise() %>% # rowwise allows much easier plotting syntax
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
  ungroup() # undo rowwise

require(patchwork)
phenol %$% 
  wrap_plots(Technical_Plot) %>%
  ggsave(filename = "Technical_Data.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 60, width = 60, units = "cm")
# This technical variation needs to be modelled.

# 2.3 Prior simulation ####
# One thing I know to be consistent across samples is that they need
# to be positive and they are likely between 0 and 2.8, which is roughly
# the maximal absorbance for this chromophore. Previously I modelled
# absorbance with a gamma likelihood to ensure positivity but this causes
# long distribution tails when propagating uncertainty and is not strictly
# necessary because even a normal likelihood will be positive given little
# variation. Also it is helpful to use truncation above 2.8 which is easiest
# to do with a normal likelihood, truncated at 0 and 2.8. 

# Due to the low number of technical replicates (n = 3) there are divergences 
# if I set the prior on the grand mean. Therefore, I'll include the calculated 
# mean in the Stan model, and centre the prior for mu on that with 8% relative sd. 
# sigma has the usual exponential prior but the rate is inversely related to the 
# calculated mean.

# 2.4 Stan model ####
technical_stan <- "
    data{
      int n;
      vector<lower=0,upper=2.8>[n] Absorbance;
      real Absorbance_mean;
      }

    parameters{
      real<lower=0,upper=2.8> mu;
      real<lower=0> sigma;
      }

    model{
      // Priors
      mu ~ normal( Absorbance_mean , 0.08 * Absorbance_mean ) T[ 0 , 2.8 ];
      sigma ~ exponential( 15 / ( Absorbance_mean + 0.4 ) );

      // Likelihood
      Absorbance ~ normal( mu , sigma ) T[ 0 , 2.8 ];
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
  print(n = 77)
# No rhat above 1.001.

# 2.5.2 Chains ####
require(bayesplot)
phenol %<>%
  rowwise() %>%
  mutate(
    Technical_Chains = 
      list(
        Technical_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(ID)
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
require(truncnorm) # R doesn't have a native turncated normal.
phenol %<>%
  mutate(
    Technical_Prior = Absorbance_mean %>%
      map(
        ~ tibble(.chain = 1:8 %>% rep(each = 1e4),
                 .iteration = 1:1e4 %>% rep(times = 8),
                 .draw = 1:8e4, # variable mean and sd
                 mu = rtruncnorm( n = 8e4 , mean = .x , sd = 0.08 * .x , a = 0 , b = 2.8 ), 
                 sigma = rexp( 8e4 , 15 / ( .x + 0.4 ) ))
      )
  )

# 2.6.2 Extract posteriors ####
phenol %<>%
  mutate(
    Technical_Posterior = Technical_Samples %>%
      map(
        ~ .x %>% spread_draws(mu, sigma)
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
            pivot_longer(cols = c(mu, sigma),
                         names_to = ".variable",
                         values_to = ".value") %>%
            prior_posterior_plot() +
            ggtitle(ID)
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
          mutate(Absorbance = rtruncnorm( n = n() , mean = mu , sd = sigma , a = 0 , b = 2.8 ))
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
                     aes(2, mu), alpha = 0.5,
                     point_interval = NULL) +
            stat_eye(data = Technical_Prediction,
                     aes(2, Absorbance), alpha = 0.5,
                     point_interval = NULL) +
            coord_cartesian(ylim = Technical_Data %$%
                              c( (1 - 0.4) * Absorbance_mean, 
                                 (1 + 0.4) * Absorbance_mean )) +
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

# 3. Standard curve models ####
# 3.1 Prepare data ####
# No unnesting required because standard curves are hierarhcically higher
# than samples.
phenol %<>%
  group_by(Name, Date, Plate, Standard_Data) %>%
  nest(.key = "Technical")

phenol

# 3.2 Visualise data ####
phenol %<>%
  rowwise() %>% 
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
  ungroup() 

phenol %$% wrap_plots(Standard_Plot)
# The limit of linearity for this assay is assumed to be 1 mg mL^-1, 
# but this suggests that saturation already happens earlier.

# 3.3 Prior simulation ####
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
       beta = rgamma( 1e3 , 2.5^2 / 2^2 , 2.5 / 2^2 ),
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

# 3.4 Stan models ####
# 3.4.1 Linear model ####
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
      beta ~ gamma( square(2.5) / square(2) , 2.5 / square(2) );
      sigma ~ exponential( 5 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + beta * Concentration[i];
      }

      // Truncated normal likelihood
      Absorbance ~ normal( mu , sigma ) T[0,];
      }
    
    generated quantities{
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + beta * Concentration[i];
      }
      
      vector[n] log_lik;
      for ( i in 1:n ) {
        log_lik[i] = normal_lpdf( Absorbance[i] | mu[i] , sigma ) - 
                     normal_lccdf( 0 | mu[i] , sigma );
      }
      }
"

require(cmdstanr)
standard_lm_model <- standard_lm_stan %>%
  write_stan_file() %>%
  cmdstan_model()

# 3.4.2 Rectangular hyperbola ####
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
      beta ~ gamma( square(2.5) / square(2) , 2.5 / square(2) );
      sigma ~ exponential( 5 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + Amax * beta * Concentration[i] / 
                ( Amax + beta * Concentration[i] );
      }

      // Truncated normal likelihood
      Absorbance ~ normal( mu , sigma ) T[0,];
      }
      
      generated quantities{
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + Amax * beta * Concentration[i] / 
                ( Amax + beta * Concentration[i] );
      }
      
      vector[n] log_lik;
      for ( i in 1:n ) {
        log_lik[i] = normal_lpdf( Absorbance[i] | mu[i] , sigma ) - 
                     normal_lccdf( 0 | mu[i] , sigma );
      }
      }
"

standard_rh_model <- standard_rh_stan %>%
  write_stan_file() %>%
  cmdstan_model()

# 3.4.3 Exponential saturation ####
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
      beta ~ gamma( square(2.5) / square(2) , 2.5 / square(2) );
      sigma ~ exponential( 5 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + Amax * ( 1 - exp( -beta * Concentration[i] / Amax ) );
      }

      // Truncated normal likelihood
      Absorbance ~ normal( mu , sigma ) T[0,];
      }
      
      generated quantities{
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + Amax * ( 1 - exp( -beta * Concentration[i] / Amax ) );
      }
      
      vector[n] log_lik;
      for ( i in 1:n ) {
        log_lik[i] = normal_lpdf( Absorbance[i] | mu[i] , sigma ) - 
                     normal_lccdf( 0 | mu[i] , sigma );
      }
      }
"

standard_es_model <- standard_es_stan %>%
  write_stan_file() %>%
  cmdstan_model()

# 3.4.4 Hyperbolic tangent ####
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
      beta ~ gamma( square(2.5) / square(2) , 2.5 / square(2) );
      sigma ~ exponential( 5 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + Amax * tanh( beta * Concentration[i] / Amax );
      }

      // Truncated normal likelihood
      Absorbance ~ normal( mu , sigma ) T[0,];
      }
      
      generated quantities{
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = A0 + Amax * tanh( beta * Concentration[i] / Amax );
      }
      
      vector[n] log_lik;
      for ( i in 1:n ) {
        log_lik[i] = normal_lpdf( Absorbance[i] | mu[i] , sigma ) - 
                     normal_lccdf( 0 | mu[i] , sigma );
      }
      }
"

standard_ht_model <- standard_ht_stan %>%
  write_stan_file() %>%
  cmdstan_model()

# 3.4.5 Run all models ####
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

# 3.5 Model checks ####
# 3.5.1 Rhat ####
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

# 3.5.2 Chains ####
require(bayesplot)
phenol %<>%
  rowwise() %>%
  mutate(
    Standard_lm_Chains = 
      list(
        Standard_lm_Samples$draws(format = "df") %>%
          mcmc_rank_overlay(pars = c("A0", "beta", "sigma")) +
          ggtitle(Name) # and adding titles
      ),
    Standard_rh_Chains =
      list(
        Standard_rh_Samples$draws(format = "df") %>%
          mcmc_rank_overlay(pars = c("A0", "Amax", "beta", "sigma")) +
          ggtitle(Name) # and adding titles
      ),
    Standard_es_Chains =
      list(
        Standard_es_Samples$draws(format = "df") %>%
          mcmc_rank_overlay(pars = c("A0", "Amax", "beta", "sigma")) +
          ggtitle(Name) # and adding titles
      ),
    Standard_ht_Chains =
      list(
        Standard_ht_Samples$draws(format = "df") %>%
          mcmc_rank_overlay(pars = c("A0", "Amax", "beta", "sigma")) +
          ggtitle(Name) # and adding titles
      )
    ) %>%
  ungroup()

phenol %$% wrap_plots(Standard_lm_Chains)
phenol %$% wrap_plots(Standard_rh_Chains)
phenol %$% wrap_plots(Standard_es_Chains)
phenol %$% wrap_plots(Standard_ht_Chains)
# Chains look good.

# 3.5.3 Pairs ####
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

# 3.6 Prior-posterior comparison ####
# 3.6.1 Sample prior ####
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

# 3.6.2 Combine prior and posterior ####
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

# 3.6.3 Plot comparison ####
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
# Amax seems to be more realistically parameterised in the hyperbolic tangent 
# because absorbance saturation is expected around 3 a.u. but the Amax posterior
# clearly escapes the prior bounds to much larger values in the rectangular
# hyperbola and exponential saturation parameterisations. So one point for ht.

# 3.7 Prediction ####
# 3.7.1 Calculate prediction ####
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

# 3.7.2 Plot prediction ####
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
  ) %>%
  ungroup()

phenol %$% wrap_plots(Standard_Prediction_Plot)
# Saturating models clearly fit much better than the linear model but it is still
# hard to tell which fits best. I need finer tools. Enter loo!

# 3.8 Model selection ####
require(loo)
phenol %<>%
  rowwise() %>%
  mutate(
    loo_comparison =
      list(
        loo_compare(
          list(
            lm = Standard_lm_Samples$loo(cores = parallel::detectCores()),
            rh = Standard_rh_Samples$loo(cores = parallel::detectCores()),
            es = Standard_es_Samples$loo(cores = parallel::detectCores()),
            ht = Standard_ht_Samples$loo(cores = parallel::detectCores())
          )
        )        
      )
    ) %>%
  ungroup()

phenol$loo_comparison
# The hyperbolic tangent is best in three out of five cases. Together with the more
# realistic parameterisation of Amax, hyperbolic tangent is the clear winner.

# 3.9 Inverse prediction ####
# 3.9.1 Rearrange parameters ####
# The optimal standard curve model takes the form Absorbance = A0 + Amax * 
# tanh( beta * Concentration / Amax ), so it predicts Absorbance with
# Concentration. I want to do the inverse, that is predict Concentration
# with Absorbance. Solving for Concentration, I get Concentration = Amax / beta *
# atanh( ( Absorbance - A0 ) / Amax ).

# 3.9.2 Visualise ####
# Reminder of original prediction
phenol %$% wrap_plots(Standard_ht_Prediction_Plot) 

# Calculate inverse prediction
phenol %<>%
  mutate(
    Standard_ht_Inverse_Prediction = Standard_ht_Prior_Posterior %>%
      map2(
        Standard_Data,
        ~ spread_continuous(prior_posterior_draws_short = .x %>%
                              filter(distribution == "posterior"),
                            data = .y,
                            predictor_name = "Absorbance") %>%
          mutate(mu = Amax / beta * atanh( ( Absorbance - A0 ) / Amax )) %>%
          group_by(distribution, Absorbance) %>%
          mean_qi(mu, .width = c(.5, .8, .9))
        )
    )
# Some NaNs produced, likely by atanh(), which is defined between -1 and 1,
# so if A0 = 0 and a low sample of Amax = 2.4, then some values for Absorbance
# will be higher than Amax and cause values above 1 to be passed to atanh().

# Plot inverse prediction
phenol %<>%
  rowwise() %>% 
  mutate(
    Standard_ht_Inverse_Prediction_Plot =
      list(
        Standard_ht_Inverse_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, aes(Absorbance, Concentration)) +
            geom_line(aes(Absorbance, mu)) +
            geom_ribbon(aes(Absorbance, ymin = .lower, ymax = .upper,
                            alpha = factor(.width))) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            ggtitle(Name) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  ) %>%
  ungroup()

phenol %$% wrap_plots(Standard_ht_Inverse_Prediction_Plot)
# Looks fine.

# What about those NaNs?
phenol$Standard_ht_Inverse_Prediction %>% 
  map(~ .x %>% filter(is.nan(mu)))
# Only predictions for Absorbance >= 2.46 a.u. in one case.

phenol$Standard_Data[[2]] %>%
  print(n = 21)
# Should still be fine to predict sample concentrations
# because that 2.49 a.u. was clearly and outlier, so
# there is no need to predict >= 2.46 a.u.

# 3.9.3 Re-nest data ####
phenol %<>% unnest(cols = Technical)

# Now that standard curve and technical triplicate estimates are 
# on the same nesting level, the list variables of interest are
phenol$Technical_Prediction
# and
phenol$Standard_ht_Prior_Posterior %>%
  map(~ .x %>% filter(distribution == "posterior"))

# 3.9.4 Calculate sample concentration ####
# I need a new list variable where the standard and technical
# estimates are joined. I would use cross_join() if I had to match
# estimates with each of a vector of point observations, but since 
# both are posteriors based on 8e4 samples, I can use full_join()
# to match by .chain, .iteration and .draw, which dramatically
# reduces computational cost.

phenol %<>%
  mutate(
    Samples_Data = Technical_Prediction %>%
      map2(
        Standard_ht_Prior_Posterior,
        ~ .x %>%
          select(-c(mu, sigma)) %>%
          full_join(
            .y %>% 
              filter(distribution == "posterior") %>%
              select(-c(distribution, sigma)),
            by = c(".chain", ".iteration", ".draw")
            ) %>% 
          mutate( # Here's where the magic happens!
            Concentration = Amax / beta * atanh( ( Absorbance - A0 ) / Amax )
            )
      )
  )

phenol$Samples_Data %>% 
  map(~ .x %>% filter(is.nan(Concentration)))
# Very few NaNs produced (maximum 2 per sample). That's good enough.
# Still enough samples and it adds further regularisation.

# Remove NaNs.
phenol %<>%
  mutate(
    Samples_Data = Samples_Data %>%
      map(
        ~ .x %>% filter(!is.nan(Concentration))
      )
  )

phenol$Samples_Data %>% 
  map_lgl(~ .x %$% any(is.nan(Concentration))) %>%
  any()
# No more NaNs.

# Multiply 50% diluted samples by 2.
phenol %<>%
  mutate(
    Samples_Data = if_else(ID %>% str_detect("50%"),
                           Samples_Data %>%
                             map(
                               ~ .x %>%
                                 mutate(Concentration = Concentration * 2)
                             ),
                           Samples_Data),
    ID = if_else(ID %>% str_detect("50%"),
                 ID %>% str_remove("_50%"),
                 ID) %>% fct()
    )

# Calculate summary.
phenol %<>%
  mutate(
    Samples_Data_Summary = Samples_Data %>%
      map(
        ~ .x %>%
          summarise(Concentration_mean = mean(Concentration),
                    Concentration_sd = sd(Concentration))
      )
  )

# 4. Phenol model ####
# 4.1 Visualise ####
( phenol %>%
    select(Name, ID, Samples_Data) %>%
    unnest(cols = Samples_Data) %>%
    ggplot(aes(Concentration, ID)) +
      geom_vline(xintercept = 0) +
      stat_slab(n = 2e3, height = 10, 
                colour = "black", linewidth = 0.1) +
      facet_grid(rows = vars(Name), scales = "free") +
      coord_cartesian(xlim = c(0, 2)) +
      theme_minimal() +
      theme(panel.grid = element_blank()) ) %>%
  ggsave(filename = "Samples_Data.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 40, width = 20, units = "cm")

( phenol %>%
    select(Name, ID, Samples_Data_Summary) %>%
    unnest(cols = Samples_Data_Summary) %>%
    ggplot(aes(Concentration_mean, ID)) +
      geom_vline(xintercept = 0) +
      geom_pointrange(aes(xmin = Concentration_mean - Concentration_sd,
                          xmax = Concentration_mean + Concentration_sd)) +
      facet_grid(rows = vars(Name), scales = "free") +
      coord_cartesian(xlim = c(0, 2)) +
      theme_minimal() +
      theme(panel.grid = element_blank()) ) %>%
  ggsave(filename = "Samples_Data_Summary.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 40, width = 20, units = "cm")

# 4.2 Prepare data ####
# I need to extract several other variables from ID and Date for the model.
phenol %<>%
  mutate(Urchin = if_else(ID %>% str_detect("f") | ID %>% str_detect("F"),
                          "Faeces", "Kelp") %>% fct(),
         Experiment = case_when(
                        Date %>% year() == 2023 ~ 1,
                        Date %>% year() == 2025 &
                          ID %>% str_split_i(pattern = "_",
                                             i = 2) %>%
                                 str_detect("1") ~ 2,
                        Date %>% year() == 2025 &
                          ID %>% str_split_i(pattern = "_",
                                             i = 2) %>%
                                 str_detect("2") ~ 3
                        ) %>% str_c() %>% fct(),
         Tank = if_else(Date %>% year() == 2023,
                        ID %>% str_extract("\\d+"),
                        ID %>% str_split_i(pattern = "_",
                                           i = 1) %>%
                               str_extract("\\d+")) %>%
                fct()
         )

# Here are the data I'll pass to the model:
phenol %>%
  select(Urchin, Experiment, Tank, Samples_Data_Summary) %>%
  unnest(cols = Samples_Data_Summary) %>%
  print(n = 79)

# 4.3 Prior simulation ####
# The experimental model will compare phenolic content of urchin faeces and
# kelp while accounting for season. I will model these data with a gamma
# likelihood, which ensures positivity and accounts fo the increase in 
# uncertainty with the increase in mean seen in the visualisation. In truth
# the data ultimately follow a beta distribution because phenolic content
# is given in mg mL^-1 of phloroglucinol equivalents. All extracts have
# a 10% (w/v) concentration, i.e. 100 mg mL^-1, so mass-based concentration,
# which I am interested in is mg 100 mg^-1 or %, which is capped at 100.
# However, the beta distribution is hard to meaningfully parameterise and
# the maximal phenolic content measured here is below 2%, so the gamma
# distribution is an excellent approximation. The expected mean polyphenolic
# content for Laminaria hyperborea is 1.07 ± 0.07% according to Wright et al.
# 2022 (doi: 10.1111/gcb.16299). Based on the litertaure alone, my prior 
# expectation is that the mean for both treatments must fall near 1.07%.

tibble(n = 1:1e5,
       mu_log = rnorm( 1e5 , log(1.07) , 0.35 ),
       sigma = rexp( 1e5, 5 ),
       mu = exp(mu_log),
       P = rgamma( 1e5 , mu^2 / sigma^2 , mu / sigma^2 )) %>%
  pivot_longer(cols = c(mu, P), 
               names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, parameter)) +
    geom_vline(xintercept = c(0, 2)) +
    stat_slab(alpha = 0.5, height = 2) +
    coord_cartesian(expand = F,
                    xlim = c(-1, 3)) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 4.4 Stan model ####
phenol_stan <- "
    data{ 
      int n;
      vector[n] Concentration_mean;
      vector[n] Concentration_sd;
      array[n] int Urchin; // Urchin effect
      int n_Urchin;
      //array[n] int Experiment; // Repeat
      //int n_Experiment;
      array[n] int Tank; // Tanks across repeats
      int n_Tank;
      }

    parameters{
      // True estimates for measurement error
      vector<lower=0>[n] Concentration; 
      
      // Intercepts in log space
      vector[n_Urchin] alpha_u; 
      // vector[n_Experiment] z_e; // z-scores
      vector[n_Tank] z_t; 
      
      // Experiment and tank uncertainties
      // real<lower=0> sigma_e;
      real<lower=0> sigma_t;
      
      // Likelihood uncertainty
      real<lower=0> sigma; 
      }

    model{
      // Hyperpriors
      // sigma_e ~ exponential( 1 );
      sigma_t ~ exponential( 1 );
      
      // Priors
      alpha_u ~ normal( log(1.07) , 0.35 );
      // z_e ~ normal( 0 , 1 );
      z_t ~ normal( 0 , 1 );
      sigma ~ exponential( 5 );
      
      // Convert z-scores
      //vector[n_Experiment] alpha_e;
      //alpha_e = z_e * sigma_e + 0;
      vector[n_Tank] alpha_t;
      alpha_t = z_t * sigma_t + 0;
      
      // Model with link function
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = exp( alpha_u[ Urchin[i] ] + 
                     //alpha_e[ Experiment[i] ] + 
                     alpha_t[ Tank[i] ] );
      }

      // Likelihood with measurement error
      Concentration ~ gamma( square(mu) / square(sigma) , mu / square(sigma) );
      Concentration_mean ~ normal( Concentration , Concentration_sd );
      }
    
    generated quantities{
      // Save converted z-scores
      //vector[n_Experiment] alpha_e;
      //alpha_e = z_e * sigma_e + 0;
      vector[n_Tank] alpha_t;
      alpha_t = z_t * sigma_t + 0;
      }
"

phenol_model <- phenol_stan %>%
  write_stan_file() %>%
  cmdstan_model()

phenol_samples <- phenol_model$sample(
          data = phenol %>%
            select(Urchin, Experiment, Tank, Samples_Data_Summary) %>%
            unnest(cols = Samples_Data_Summary) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
# Some divergences despite non-centred parameterisation.

# 4.5 Model checks ####
# 4.5.1 Rhat ####
phenol_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001.

# 4.5.2 Chains ####
phenol_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "Phenol_Chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 40, width = 40, units = "cm")

# 4.6 Prior-posterior comparison ####
# 4.6.1 Sample prior ####
phenol_prior <- prior_samples(
  model = phenol_model,
  data = phenol %>%
    select(Urchin, Experiment, Tank, Samples_Data_Summary) %>%
    unnest(cols = Samples_Data_Summary) %>%
    compose_data()
  )

# 4.6.2 Plot prior-posterior comparison ####
phenol_prior %>% 
  prior_posterior_draws(
    posterior_samples = phenol_samples,
    group = phenol %>%
      select(Urchin, Experiment, Tank, Samples_Data_Summary) %>%
      unnest(cols = Samples_Data_Summary) %>%
      select(Urchin, Tank),
    parameters = c("alpha_t[Tank]", 
                   #"alpha_e[Experiment]", 
                   "alpha_u[Urchin]", 
                   "z_t[Tank]", 
                   #"z_e[Experiment]",
                   "sigma", 
                   #"sigma_e", 
                   "sigma_t"),
    format = "long"
    ) %T>%
  { prior_posterior_plot(., group_name = "Urchin") %>%
      print() } %>%
  # { prior_posterior_plot(., group_name = "Experiment") %>%
  #     print() } %>%
  prior_posterior_plot(group_name = "Tank")
# alpha_e shows a consitent shift below zero so seems to be compensating for
# something. This is because of the imbalance between treatments in experiment 1.
# Try modelling experimental and tank uncertainty within rather than across treatments.


# This isn't necessarily concerning since all alphas will be added.
# alpha_u is escaping the prior space for both groups, but that is part of the 
# story, so is not concerning.

# 4.7 Prediction ####
# 4.7.1 Combine relevant priors and posteriors ####
phenol_prior_posterior <- phenol_prior %>% 
  prior_posterior_draws(
    posterior_samples = phenol_samples,
    group = phenol %>%
      select(Urchin, Experiment, Tank, Samples_Data_Summary) %>%
      unnest(cols = Samples_Data_Summary) %>%
      select(Urchin),
    parameters = c("alpha_u[Urchin]", "sigma", "sigma_e", "sigma_t"),
    format = "short"
  )

# 4.7.2 Calculate predictions for new experiments and tanks ####
phenol_prior_posterior %<>%
  mutate(mu = exp( alpha_u ),
         obs = rgamma( n() , mu^2 / sigma^2 , mu / sigma^2 ))

# 4.7.3 Plot predictions ####
phenol_prior_posterior %>% # priors are identical for both treatments
  filter(!(Urchin == "Faeces" & distribution == "prior")) %>% # remove one
  mutate(Urchin = if_else(distribution == "prior",
                          "Prior", Urchin) %>% fct()) %>%
  select(-distribution) %>%
  pivot_longer(cols = c(mu, obs), values_to = "Concentration", names_to = "Level") %>%
  ggplot(aes(Concentration, Urchin, alpha = Level)) +
    stat_slab(height = 10, n = 2e3) +
    scale_alpha_manual(values = c(0.6, 0.3)) +
    coord_cartesian(xlim = c(0, 2)) +
    theme_minimal() +
    theme(panel.grid = element_blank())
