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
# Pretty linear. In reality there is a lower bound at 0 and an upper bound
# around the maximal absorbance (~2.5 a.u.). But since I know the limit of 
# linearity to be around 1 mg mL^-1, I can approximate with a linear model.

# 2.2 Model choice ####
# The typical linear model takes the form y ~ normal( alpha + beta * x , sigma ).
# alpha (absorbance) and beta (relationship between concentration and absorbance)
# are necessarily positive. Initially I went with an exponential distribution
# for alpha, because it is expected to be very near zero but no other information 
# is available, and a gamma distributions for beta, which is expected to be near 
# the absorbance maximum (2.5 a.u.) divided by the maximal concentration (1 mg mL^1), 
# so 2.5, as well as a normal likelihood truncated at 0. 

# But now I realise these constraints at this early stage are useless because after 
# converting fixed absorbance measurements, negative concentration estimates are still 
# possible (e.g. when sample absorbance falls within the posterior of alpha). Working 
# with fewer constrains also frees up other uncertainty propagation possibilities such 
# as incorporating uncertainty in sigma. Positivity can still be enforced at a later 
# stage when I will run the model to estimate the mean of the technical triplicate. 
# So I will now go with a normal distribution centred on zero for alpha and a normal 
# distribution centred on 2.5 for beta as well as an unconstrained normal likelihood.

# 2.3 Prior simulation ####
tibble(n = 1:1e3,  
       alpha = rnorm( 1e3 , 0 , 0.2 ),
       beta = rnorm( 1e3 , 2.5 , 0.5 )) %>%
  expand_grid(c = seq(0, 1)) %>%
  mutate(A = alpha + beta * c) %>%
  ggplot(aes(c, A, group = n)) +
    geom_hline(yintercept = c(0, 2.5)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()
# Covers all reasonable possibilities.

# 2.4 Run model ####
standard_stan <- "
    data{ 
      int n;
      vector[n] Absorbance;
      vector[n] Concentration;
      }

    parameters{
      real alpha;
      real beta;
      real<lower=0> sigma;
      }

    model{
      // Priors
      alpha ~ normal( 0 , 0.2 );
      beta ~ normal( 2.5 , 0.5 );
      sigma ~ exponential( 1 );

      // Model
      vector[n] mu;
      for ( i in 1:n ) {
        mu[i] = alpha + beta * Concentration[i];
      }

      // Likelihood
      Absorbance ~ normal( mu , sigma );
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

# 2.5 Model checks ####
# 2.5.1 Rhat ####
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

# 2.5.2 Chains ####
require(bayesplot)
phenol %<>%
  rowwise() %>%
  mutate(
    Standard_Chains = 
      list(
          Standard_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(Name) # and adding titles
      )
    ) %>%
  ungroup()

phenol %$% wrap_plots(Standard_Chains)
# Chains look good.

# 2.5.3 Pairs ####
phenol %<>%
  mutate(
    Standard_Pairs = Standard_Samples %>%
      map(
        ~ .x$draws(format = "df") %>% # mcmc_pairs does not allow adding titles
          mcmc_pairs(pars = c("alpha", "beta"))
      )
  )

phenol %$% wrap_plots(Standard_Pairs)
# Some correlation between alpha and beta, but not concerning.

# 2.6 Prior-posterior comparison ####
source("functions.R")
# 2.6.1 Sample prior ####
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

# 2.6.2 Combine prior and posterior ####
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

# 2.6.3 Plot comparison ####
phenol %<>%
  rowwise() %>%
  mutate(
    Standard_Prior_Posterior_Plot =
      list(
          prior_posterior_draws(prior_samples = Standard_Prior,
                                posterior_samples = Standard_Samples,
                                parameters = c("alpha", "beta", "sigma"),
                                format = "long") %>%
          prior_posterior_plot() +
          ggtitle(Name)
      )
  ) %>%
  ungroup()

phenol %$% 
  wrap_plots(Standard_Prior_Posterior_Plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
  
# 2.7 Prediction ####
# 2.7.1 Calculate prediction ####
phenol %<>%
  mutate(
    Standard_Prediction = Standard_Prior_Posterior %>%
      map2(
        Standard_Data,
        ~ spread_continuous(prior_posterior_draws_short = .x,
                            data = .y,
                            predictor_name = "Concentration") %>%
          mutate(mu = alpha + beta * Concentration,
                 obs = rnorm( n() , mu , sigma )) %>%
          group_by(distribution, Concentration) %>%
          reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                  obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
          unnest(c(mu, obs), names_sep = "_")
      )
  )

# 2.7.2 Plot prediction ####
phenol %<>%
  rowwise() %>% 
  mutate(
    Standard_Prediction_Plot =
      list(
        Standard_Prediction %>%
          ggplot() +
            geom_point(data = Standard_Data, aes(Concentration, Absorbance)) +
            geom_line(data = . %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            # geom_ribbon(data = . %>% filter(distribution == "posterior"),
            #             aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
            #                 alpha = factor(mu_.width))) + # unhash to see mu predictions
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
  ) %>%
  ungroup()

phenol %$% wrap_plots(Standard_Prediction_Plot)
# The linear models generally fit well but overestimate absorbance
# somewhat at concentrations near 1 mg mL^1. When observational
# uncertainty is incorporated this is good enough!

# 2.7.3 Inverse prediction ####
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

# 2.7.4 Calculating sample concentration ####
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

# 3. Technical triplicate model ####
# 3.1 Prepare data ####
# I no longer want the plate-nested structure of phenol because
# the plate-specific data (standard curve) can now be disregarded.
# Instead a nesting by ID is useful because each ID is a biological
# sample and contains a technical triplicate.

# Save data associated with plate structure to separate tibble.
plates <- phenol

# Re-nest phenol.
phenol %<>%
  select(Name, Date, Plate, Samples_Data_Converted) %>%
  unnest(cols = Samples_Data_Converted) %>%
  group_by(Name, Date, Plate, ID, Mass) %>%
  nest(.key = "Technical_Data")

phenol
phenol$Technical_Data[[30]]

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

phenol
phenol$Technical_Data_Summary[[20]]

# 3.2 Visualise data ####
require(ggdist)
phenol %<>%
  rowwise() %>%
  mutate(
    Technical_Plot = 
      list(
        Technical_Data %>%
          ggplot(aes(Concentration, Well)) +
            stat_slab(height = 2, colour = "black") +
            ggtitle(ID) +
            theme_minimal()
        )
    ) %>%
  ungroup()

phenol %$% 
  wrap_plots(Technical_Plot) %>%
  ggsave(filename = "Technical_Plot.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 80, width = 80, units = "cm")

# 3.3 Prior simulation ####
# The only thing I know to be consistent across samples is that they need
# to be positive, so I will use a gamma likelihood. The intercept alpha
# can then be normally distributed.

tibble(alpha = rnorm( 1e5 , log(0.5) , 0.8 ),
       mu = alpha %>% exp()) %>%
  ggplot(aes(y = mu)) +
    geom_hline(yintercept = c(0, 1)) + # concentrations are between 0 and 1
    geom_density(orientation = "y", colour = NA, fill = "black", alpha = 0.5) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()
# Covers all reasonable possibilities.

# 3.4 Run model ####
technical_stan <- "
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

technical_model <- technical_stan %>%
  write_stan_file() %>%
  cmdstan_model()
  
phenol %<>%
  mutate(
    Technical_Samples = Technical_Data_Summary %>%
      map(
        ~ technical_model$sample(
          data = .x %>%
            select(Concentration_mean, Concentration_sd) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        )
        )
  )

# Model fails because non-negativity in Concentration_mean is not allowed.