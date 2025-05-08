# 1. Load data ####
# 1.1 Load raw data ####
require(tidyverse)
require(here)

phenol <- here("Biochemistry", "Phenol", "Raw") %>%
  list.files(pattern = "\\.csv$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Name = Path %>% basename() %>% 
      str_remove("\\..*$") %>% str_remove("^X") %>% 
      make.names(),
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
  read_csv() %>%
  group_by(Plate) %>%
  nest(.key = "Metadata")

meta

# 1.3 Join metadata to data ####
require(magrittr)
phenol %<>%
  left_join(meta %>% rename(Name = Plate), by = "Name")
phenol

phenol %<>% 
  imap(~ .x %>% bind_cols(
                  meta %>% 
                    rename(Name = Plate) %>%
                    filter(Name == str_remove(.y, "^X"))
                  )) 
phenol

# 1.4 Check success of join ####
phenol %>% 
  map(~ .x %$% 
        identical(
          Date, 
          Name %>%
            str_split_i("_", 1) %>%
            str_remove("^X") %>% 
            ymd()
          )
      )

phenol %>% 
  map(~ .x %$% 
        identical(
          Plate, 
          Name %>%
            str_split_i("_", 3) %>%
            as.numeric()
        )
  )


phenol %>% enframe()


# 1.5 Separate standard and samples ####
phenol %<>% map(~ list(
  standard = .x %>% # the = instead of <- is important to retain names
    filter(str_detect(Annotation, "^[0-9.]+$")) %>% 
    rename(Concentration = Annotation) %>%
    mutate(Concentration = as.numeric(Concentration)) %>%
    select(-Mass), # no sample mass for standards
  samples = .x %>% filter(!str_detect(Annotation, "^[0-9.]+$"))
))

phenol

# 2. Standard curve models ####
# 2.1 Visualise data ####
require(patchwork)
phenol %>%
  map(~ .x$standard %>%
        ggplot(aes(Concentration, Absorbance)) +
          geom_point() +
          geom_line() +
          # coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 2e3)) + # unhash to see intercept
          theme_minimal()
      ) %>%
  wrap_plots()
# Pretty linear (we know the limit of linearity to be around 1 mg mL^-1).

# 2.2 Prior simulation ####
# alpha (absorbance) and beta (relationship between concentration and absorbance)
# are necessarily positive, so a gamma distribution is best. alpha is expected to be
# very near zero. beta is expected to be near the absorbance maximum (2.5 a.u.) divided
# by the maximal concentration (1 mg mL^1), so 2.5.

tibble(n = 1:1e3,  # gamma is reparameterised in terms of mean and s.d.
       alpha = rgamma(n = 1e3, shape = 0.05^2 / 0.04^2, rate = 0.05 / 0.04^2),
       beta = rgamma(n = 1e3, shape = 2.5^2 / 1^2, rate = 2.5 / 1^2)) %>%
  expand_grid(c = seq(0, 1)) %>%
  mutate(A = alpha + beta * c) %>%
  ggplot(aes(c, A, group = n)) +
    geom_hline(yintercept = c(0, 2.5)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off",
                    # xlim = c(0, 1), ylim = c(0, 2e3) # unhash to check F0
                    ) +
    theme_minimal()
# covers all reasonable possibilities

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
  alpha ~ gamma( 0.05^2 / 0.04^2, 0.05 / 0.04^2 );
  beta ~ gamma( 2.5^2 / 1^2, 2.5 / 1^2 );
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
standard_mod <- standard_stan %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
standard_samples <- phenol %>%
  map(~ standard_mod$sample(data = .x$standard %>%
                                            select(Fluorescence, Concentration) %>%
                                            compose_data(),
                                        chains = 8,
                                        parallel_chains = parallel::detectCores(),
                                        iter_warmup = 1e4,
                                        iter_sampling = 1e4)
      )

# 2.1.3 Model checks ####
# check Rhat, effective sample size and chains
kinetics_standard_ht_samples %>%
  map(~ .x$summary()) %>%
  bind_rows() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

require(bayesplot)
kinetics_standard_ht_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_rank_overlay()) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# chains look good

kinetics_standard_ht_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_pairs(pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots()
# some correlation between Fmax and beta, indicating some interdependence

# 2.1.4 Prior-posterior comparison ####
source("functions.R")
# sample prior
kinetics_standard_ht_prior <- kinetics %>%
  map(~ prior_samples(model = kinetics_standard_ht_mod,
                      data = .x$standard %>%
                        select(Fluorescence, Concentration) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
kinetics_standard_ht_prior %>%
  map2(kinetics_standard_ht_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA), # no groups so this has to be an empty list or tibble
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "long")) %>%
  map(~ .x %>% prior_posterior_plot()) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# some posteriors for F0 broke out of the expected prior probability space





