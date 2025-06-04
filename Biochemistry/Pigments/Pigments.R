# 1. Load data ####
# 1.1 Load raw data ####
require(tidyverse)
require(here)

pigments <- here("Biochemistry", "Pigments", "Raw") %>%
  list.files(pattern = "Pigments.*\\.csv$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Name = Path %>% basename() %>% 
      str_remove("\\..*$"),
    Spectra = Path %>% 
      map(~ read_csv(.x, skip = 11, col_types = list("f", "f")) %>%
            pivot_longer(cols = -c(1:2), names_to = "Wavelength", values_to = "Absorbance") %>%
            rename(Well = 1, Content = 2) %>%
            mutate(Wavelength = Wavelength %>% as.numeric())),
    Date = Name %>% str_remove("^X") %>% str_split_i("_", 1) %>% ymd(),
    Plate = Name %>% str_split_i("_", 3) %>% as.numeric(),
    Concentration = Name %>% str_split_i("_", 5) %>%
                      { if_else(str_detect(., "%"),
                                str_extract(., "\\d+%"),
                                "100%") }
  ) %>%
  select(-Path)

pigments
pigments$Spectra[[1]]

# 1.2 Load metadata ####
meta <- here("Biochemistry", "Pigments", "Pigments_Meta.csv") %>% 
  read_csv(col_types = list("f", "f", "f")) %>%
  rename(Name = Plate, ID = Annotation) %>%
  group_by(Name) %>%
  nest(.key = "Metadata")

meta

# 1.3 Join metadata to data ####
require(magrittr)
pigments %<>% full_join(meta, by = "Name")

pigments
pigments$Metadata[[1]]
rm(meta)

pigments %<>% 
  mutate(Spectra = map2(Spectra, Metadata, 
                        ~ full_join(.x, .y, by = "Content"))) %>%
  select(-Metadata)

pigments$Spectra[[1]]

# 1.4 Visualise spectra ####
pigments %<>%
  rowwise() %>% # rowwise allows much easier plotting syntax
  mutate(
    Spectra_Plot = 
      list(
        Spectra %>%
          ggplot(aes(Wavelength, Absorbance)) +
            geom_point(shape = 16, alpha = 0.01) +
            ggtitle(Date, Concentration) + # max 2 titles allowed
            theme_minimal() +
            theme(panel.grid = element_blank())
        )
    ) %>%
  ungroup() # undo rowwise

require(patchwork)
pigments %$% 
  wrap_plots(Spectra_Plot, nrow = 2) %>%
  ggsave(filename = "spectra.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 20, width = 60, units = "cm")
# Clearly several plates hit the absorbance ceiling of ~3.5 a.u.
# I therefore ran dilutions. Now these plates can be removed.
# For the 2023 plates 50% dilution was enough, but the 2025
# plates needed 25% dilution despite the larger elution volume.

# 1.5 Remove unwanted plates ####
pigments %<>%
  filter(Date %>% year() == 2023 & Concentration == "50%" |
           Date %>% year() == 2025 & Concentration == "25%")

# 1.6 Separate sample and blank spectra ####
pigments %<>%
  mutate(
    Sample_Spectra = Spectra %>% 
      map(
        ~ .x %>%
          filter(ID != "96% ethanol")
      ), 
    Blank_Spectra = Spectra %>%
      map(
        ~ .x %>% 
          filter(ID == "96% ethanol") %>% 
          select(-c(ID, Mass))
      )
  ) %>%
  select(-Spectra)

pigments

# 1.7 Average the blank spectra ####
pigments %<>%
  mutate(Blank_Spectra = Blank_Spectra %>%
           map(~ .x %>%
                 group_by(Wavelength) %>%
                 summarise(Absorbance = mean(Absorbance))))
pigments

# 1.8 Join blank to sample spectra ####
pigments %<>%
  mutate(Spectra = Sample_Spectra %>%
           map2(
             Blank_Spectra,
             ~ .x %>% 
               full_join(.y %>%
                           rename(Blank = Absorbance), 
                         by = "Wavelength")
             )) %>%
  select(-c(Sample_Spectra, Blank_Spectra))

pigments$Spectra

# 1.9 Re-nest data ####
pigments %<>%
  select(-Spectra_Plot) %>%
  unnest(cols = Spectra) %>%
  group_by(Name, Date, Plate, Concentration, 
           ID, Mass, Well, Content) %>%
  nest(.key = "Spectrum")

pigments

# 1.10 Visualise spectra ####
pigments %<>%
  rowwise() %>%
  mutate(
    Spectrum_Plot = 
      list(
        Spectrum %>%
          ggplot() +
            geom_line(aes(Wavelength, Blank), colour = "grey") +
            geom_line(aes(Wavelength, Absorbance)) +
            ylab("Absorbance") +
            ggtitle(ID, Well) +
            theme_minimal() +
            theme(panel.grid = element_blank())
        )
    ) %>%
  ungroup()

pigments %$% 
  wrap_plots(Spectrum_Plot) %>%
  ggsave(filename = "spectrum.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 60, width = 80, units = "cm")

# Solvent (96% ethanol) triplicate blanks were measured for each plate, 
# but the spectral deconvolution pipeline I use does not require manual 
# subtraction of blank spectra since these are modelled as part of the 
# microplate background spectrum component of the sample spectrum and 
# thus accounted for:

# "Note that the absorbance of a blank (plate and solvent absorbance) was not 
# subtracted. Instead, this absorbance was modeled as a part of the general 
# background..." — Thrane et al. 2015, doi: 10.1371/journal.pone.0137645.

# This background modelling method has been verified for several microplates 
# including the Greiner Bio-One model I used (#655201, Fig. S2 in the paper).
# I will also prove it later by plotting the measured blank alongside the
# modelled background.

# 2. Spectral deconvolution ####
# 2.1 Get pigment function pipeline ####
# The functions need to be downloaded from File S1 in Thrane et al. 2015 at 
# https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0137645.s003.
# I have stored all required files in the folder "Deconvolution". The only 
# edit I made is to comment out the object "core" in pigment.function.R, 
# because I want to define my own "core". The package nnls for non-negative 
# least squares fitting needs to be installed. 
require(withr) # Give function script a temporary working directory.
with_dir(here("Biochemistry", "Pigments", "Deconvolution"),
         source("pigment.function.R"))

# 2.2 Define component pigments ####
# Chl b does not occur in Laminaria hyperborea, so cannot be part of
# the spectrum. The same applies to phytoplankton-specific pigments.
# Here is my list of pigments that can be found in L. hyperborea:
core <- c("Chl.a", "Chl.c1", "Chl.c2", "Phe.a", "Pheide.a",
          "Fuco", "bb.Car", "Viola", "Anth", "Zea")

# 2.3 Gaussian peak spectra model ####
# 2.3.1 Run model ####
# Thrane et al. use 400-700 nm and that's what I've used in the past 
# (Wright & Foggo 2021, doi: 10.3354/meps13886; Wright & Kregting 2023, 
# doi: 10.1007/s00227-023-04289-y). But I found that background is 
# modelled more accurately if I use the full 400-750-nm spectrum.

pigments %<>%
  mutate(GPS_Model = Spectrum %>%
           map(~ .x %$% 
                 pigment.fit(w = Wavelength, 
                             y = Absorbance)))

# 2.3.2 Visualise residuals ####
pigments %<>%
  mutate(GPS_Residuals = GPS_Model %>%
           map(~ tibble(Wavelength = .x$w,
                        Prediction = .x$m$fitted %>% c(),
                        Residuals = .x$m$residuals %>% c())))

pigments %<>%
  rowwise() %>%
  mutate(
    Residuals_Plot = 
      list(
        GPS_Residuals %>%
          ggplot() +
            geom_hline(yintercept = 0) +
            geom_point(aes(Wavelength, Residuals, size = Prediction),
                       shape = 16, alpha = 0.05) +
            scale_size_continuous(range = c(0.5, 3),
                                  guide = "none") +
            ggtitle(ID, Well) +
            coord_cartesian(clip = "off") +
            theme_minimal() +
            theme(panel.grid = element_blank())
        )
    ) %>%
  ungroup()

pigments %$% 
  wrap_plots(Residuals_Plot) %>%
  ggsave(filename = "residuals.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 60, width = 80, units = "cm")
# Residuals are fairly equally balanced on both sides of the zero line 
# with a slight skew towards positive residuals. Larger predictions 
# usually lead to larger residuals. There is definitely variability
# in uncertainty across spectra which should be propagated.

require(ggridges)
pigments %<>%
  rowwise() %>%
  mutate(
    Residuals_Density_Plot = 
      list(
        GPS_Residuals %>%
          ggplot() +
            geom_density_ridges(aes(Residuals, y = 0),
                                position = position_raincloud(height = 30),
                                jittered_points = TRUE,
                                point_alpha = 0.1, point_shape = 16) +
            geom_vline(xintercept = 0) +
            ggtitle(ID, Well) +
            theme_minimal() +
            theme(panel.grid = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank())
        )
    ) %>%
  ungroup()

pigments %$% 
  wrap_plots(Residuals_Density_Plot) %>%
  ggsave(filename = "residuals_density.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 60, width = 80, units = "cm")
# Despite the slight right skew, residuals are fairly normally distributed, 
# although a gamma likelihood clearly would be better if this could be 
# implemented in nnls, which it can't. A project for another time. I will
# approximate the GPS model uncertainty with a normal distribution.

# 2.3.3 Compute summary statistics ####
# Compute sd of residuals
pigments %<>%
  mutate(SD = GPS_Residuals %>%
           map_dbl(~ .x %$% sd(Residuals)))

# Compute root mean square deviation
pigments %<>%
  mutate(RMSD = GPS_Model %>%
           map_dbl(~ sqrt( .x$m$deviance / length(.x$y) )))

# Compare SD and RMSD
pigments %$% (SD - RMSD) %>% mean() 
# They are practically identical. I'll use RMSD because
# this is directly provided by the model and also generally
# used for nnls.

# Re-express RMSD as coefficient of variation
pigments %<>%
  rowwise() %>%
  mutate(CV = RMSD / Spectrum %$% mean(Absorbance)) %>%
  ungroup()
# This is necessary because pigment concentrations are on
# a different scale to absorbance.

# Summarise
pigments %>%
  summarise(CV_mean = mean(CV),
            CV_sd = sd(CV))
# CV = 4.02% ± 1.46% across samples. Since this is a frequentist
# metric it is not immediately intuitive, but basically indicates that 
# at each wavelength the modelled absorbance falls within ~4% of the 
# observed absorbance 68% of times on average.

# 2.3.4 Visualise modelled spectra ####
# Extract individual and summed fits.
pigments %<>%
  mutate(GPS_Fit = GPS_Model %>% 
           map(
               ~ tibble(
                   Wavelength = .x$w,
                   Total = fitted.spectrum(.x) %>% c(),
                   Background = background.spectrum(.x) %>% c(),
                   Pigments = pigment.spectrum(.x) %>% c()
                 ) %>%
                 bind_cols(
                   .x %$% (
                     X[,(k + 2):ncol(X)] %*% 
                       diag( m$x[(k + 2):length(m$x)] ) %>%
                     as_tibble() %>%
                     set_names( colnames(X)[(k + 2):ncol(X)] ) 
                     )
                 ) %>%
                 pivot_longer(cols = -Wavelength,
                              names_to = "Spectrum",
                              values_to = "Absorbance") %>%
                 mutate(Spectrum = Spectrum %>% fct())
               ))

pigments$GPS_Fit %>%
  bind_rows() %>% 
  group_by(Spectrum) %>%
  filter(all(Absorbance == 0)) %>%
  ungroup() %$%
  unique(Spectrum)
# Only Violaxanthin is absent from all samples.

# Visualise for all samples.
pigments %<>%
  rowwise() %>%
  mutate(
    GPS_Plot = 
      list(
        ggplot() +
          geom_point(data = Spectrum,
                     aes(Wavelength, Blank),
                     shape = 16, size = 0.5,
                     alpha = 0.2, colour = "grey") +
          geom_point(data = Spectrum,
                     aes(Wavelength, Absorbance),
                     shape = 16, size = 0.5,
                     alpha = 0.2) +
          geom_line(data = GPS_Fit %>%
                      mutate(Spectrum = Spectrum %>%
                               fct_relevel("Total", "Background", "Pigments",
                                           "Chl.a", "Phe.a", "Pheide.a",
                                           "Chl.c1", "Chl.c2", "Fuco",
                                           "bb.Car", "Viola", "Anth", "Zea")),
                    aes(Wavelength, Absorbance, 
                        colour = Spectrum)) +
          scale_colour_manual(values = c("black", "grey", "forestgreen",
                                         "cyan3", "darkseagreen", "olivedrab",
                                         rep("cyan4", 2), "darkorange",
                                         "tomato1", rep("goldenrod1", 3))) +
          ylab("Absorbance") +
          ggtitle(ID, Well) +
          theme_minimal() +
          theme(panel.grid = element_blank())
        )
    ) %>%
  ungroup()

pigments %$% 
  ( wrap_plots(GPS_Plot) +
    plot_layout(guides = "collect") &
      theme(legend.position = "bottom")) %>%
  ggsave(filename = "GPS.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 100, width = 100, units = "cm")

# The first thing to notice is that the modelled background is always better
# (accounts for more background and fits the sample spectrum better) than the
# measured blank. There is no benefit gained by first subtracting the blank 
# manually. On the contrary, it would cause negative values in some cases.

# We can already see that faecal samples contain very little Chl a relative to its
# degradation products (mostly Pheide a) while kelp samples contain mostly Chl a.

# Select two spectra for supplement.
spectra <- pigments %>%
  filter( (ID == "U2_A1" | ID == "U2_F1") & Well == "A08" ) %>%
  mutate(
    GPS_Fit = GPS_Fit %>%
      map(~ .x %>% 
            pivot_wider(values_from = Absorbance,
                        names_from = Spectrum)),
    Spectra = Spectrum %>%
      map2(
        GPS_Fit,
        ~ .x %>%
          full_join(.y, by = "Wavelength")
      )
  ) %>%
  select(-c(starts_with("Spectrum"), starts_with("GPS"), starts_with("Residuals"))) %>%
  unnest(Spectra) %>%
  mutate(Treatment = if_else(ID %>% str_detect("A"),
                             "Kelp", "Faeces") %>% fct_relevel("Kelp")) %>%
  select(where(~ !all(. == 0))) %>%
  mutate(Chl.c = Chl.c1 + Chl.c2) %>%
  select(-c(Chl.c1, Chl.c2))

# Visualise selected spectra.
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
                 legend.key.width = unit(.5, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "inside",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing = unit(1, "cm"),
                 text = element_text(family = "Futura"))

require(ggspectra)
Fig_S3a <- spectra %>%
  select(Treatment, Wavelength, Absorbance, Blank) %>%
  pivot_longer(cols = c(Absorbance, Blank),
               names_to = "Spectrum",
               values_to = "Absorbance") %>%
  mutate(Spectrum = Spectrum %>% 
           fct_recode(Extract = "Absorbance") %>%
           fct_relevel("Blank")) %>%
  ggplot() +
    geom_point(aes(Wavelength, Absorbance, shape = Spectrum),
               size = 2, alpha = 0.2) +
    scale_shape_manual(values = c(17, 16)) +
    geom_line(data = spectra %>%
                select(c(Treatment, Wavelength, 
                         Total, Background, Pigments)) %>%
                pivot_longer(cols = -c(Treatment, Wavelength),
                             names_to = "Spectrum",
                             values_to = "Absorbance") %>%
                mutate(Spectrum = Spectrum %>% fct()),
              aes(Wavelength, Absorbance, colour = Spectrum)) +
    scale_colour_manual(values = c("#000000", "#b5b8ba", "#8fcb24")) +
    scale_x_continuous(breaks = seq(400, 750, 50)) +
    scale_y_continuous(labels = scales::label_number(accuracy = c(1, rep(0.1, 4)))) +
    wl_guide(aes(Wavelength), ymin = 1.57, ymax = 1.6035) +
    facet_grid(~ Treatment) +
    labs(x = "Wavelength (nm)", y = "Absorbance") +
    coord_cartesian(xlim = c(400, 750), ylim = c(0, 1.6),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.position.inside = c(0.8, 0.75))

Fig_S3b <- spectra %>% # Subtract modelled background from measurements.
  mutate(Absorbance = Absorbance - Background) %>%
  ggplot() +
    geom_point(aes(Wavelength, Absorbance),
               shape = 16, size = 2, alpha = 0.2) +
    geom_line(data = spectra %>%
                select(-c(1:11, Absorbance, Blank, Total, Background)) %>%
                pivot_longer(cols = -c(Wavelength, Treatment),
                             names_to = "Spectrum",
                             values_to = "Absorbance") %>%
                mutate(Spectrum = Spectrum %>% 
                         fct_relevel("Pigments", "Chl.a", "Pheide.a", 
                                     "Chl.c", "Fuco", "bb.Car", "Zea")),
              aes(Wavelength, Absorbance, colour = Spectrum)) +
    scale_colour_manual(values = c("#8fcb24", "#009e8f", "#a29400", "#a1d5cf",
                                   "#ef8407", "#be3819", "#f1c700"),
                        labels = c("Pigments", expression("Chlorophyll "*italic("a")),
                                   expression("Pheophorbide "*italic("a")),
                                   expression("Chlorophyll "*italic("c")),
                                   "Fucoxanthin", expression(beta*"-carotene"),
                                   "Zeaxanthin")) +
    scale_x_continuous(breaks = seq(400, 750, 50)) +
    scale_y_continuous(labels = scales::label_number(accuracy = c(1, rep(0.1, 4)))) +
    # wl_guide(aes(Wavelength), ymin = 1.57, ymax = 1.6035) +
    facet_grid(~ Treatment) +
    labs(x = "Wavelength (nm)", y = "Absorbance") +
    coord_cartesian(xlim = c(400, 750), ylim = c(0, 1.6),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.position.inside = c(0.8, 0.72))

Fig_S3 <- ( ( Fig_S3a +
                theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank()) ) / 
              Fig_S3b +
                theme(strip.text = element_blank()) ) +
          plot_annotation(tag_levels = "a") &
          theme(plot.tag = element_text(family = "Futura", size = 15, face = "bold"),
                plot.tag.position = c(0.005, 0.97))

Fig_S3 %>%
  ggsave(filename = "Fig_S3.pdf", device = cairo_pdf, path = "Figures", 
         height = 21, width = 21, units = "cm")

# 2.4 Pigment concentrations ####
# 2.4.1 Pathlength ####
# Empirical pathlength
pathlength <- here("Biochemistry", "Pigments", "Raw") %>%
  list.files(pattern = "Pathlength.*\\.csv$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Name = Path %>% basename() %>% 
      str_remove("\\..*$"),
    Data = Path %>% 
      map(~ read_csv(.x, skip = 10, col_types = list("f", "f")) %>%
            rename(A997 = 3, A900 = 4)),
    Date = Name %>% str_remove("^X") %>% str_split_i("_", 1) %>% ymd(),
    Pipetting = if_else(Name %>% str_detect("epMotion"),
                        "epMotion", "Manual") %>% fct()
  ) %>%
  select(-Path)

pathlength
pathlength$Data

# Calculate pathlength in cm
pathlength %<>%
  mutate(Data = Data %>%
           map(~ .x %>% # Water-peak pathlength equation:
                 mutate(Pathlength = (A997 - A900) / 0.18)))

# Unnest data
pathlength %<>% unnest(cols = Data)

# Calculate theoretical pathlength
# According to the product data sheet for Greiner Bio-One #655201
# each well has a top diameter of 6.96 mm, a bottom diameter of
# 6.58 mm and a height of 10.9 mm. Converting diameters to radii:
r_b <- 6.58 / 2 # 3.29 mm
r_t <- 6.96 / 2 # 3.48 mm 
# The equation for the volume of frustum is V = 1/3 * pi * h * (r_b^2 +
# r_b * r_t + r_t^2), so total well volume is
1/3 * pi * 10.9 * (r_b^2 + r_b * r_t + r_t^2) # 392.4711 mm^3
# But since I only filled the well with 200 µL, I don't know height
# or r_t. The increase in radius with height, i.e. the slope of taper
#  of the well frustum is
d_r <- (r_t - r_b) / 10.9 # 0.01743119 mm mm^-1
# Now I can calculate the radius at any given height as r_h = r_b + d_r * h.
r_b + d_r * 10.9 # This gives me r_t because I used the total height.
# So I can re-express r_t in the volume equation in terms of r_b, d_r and h 
# to create a volume function that only takes height as a variable, and
# r_b and d_r as constants.
V <- function(h){
  r_b <- 6.58 / 2
  r_t <- 6.96 / 2
  d_r <- (r_t - r_b) / 10.9
  r_h <- r_b + d_r * h
  V <- 1/3 * pi * h * ( r_b^2 + r_b * r_h + r_h^2 )
  return(V)
}
V(10.9) # This returns the same volume as before: 392.4711 mm^3.
# Now my unknown isn't V, but h, but since I cannot simply solve for
# h because this is a polynomial I'll use stats::uniroot().
pathlength_theoretical <- 
  uniroot(f = function(h) V(h) - 200, # 200 µL is my known V.
          interval = c(0, 10.9))$root / 10
# The theoretical pathlength is 5.71 mm, or 0.57 cm.

# Visualise pathlength estimates
pathlength %>%
  ggplot(aes(Pathlength, Pipetting, fill = Pipetting, colour = Pipetting)) +
    geom_density_ridges(scale = 0.5, alpha = 0.6, position = "raincloud", jittered_points = TRUE, 
                        point_alpha = 0.2, point_shape = 16, point_size = 2, bandwidth = 0.003) +
            geom_vline(xintercept = pathlength_theoretical) +
            theme_minimal() +
            theme(panel.grid = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  legend.position = "top")
# The theoretical estimate is higher than the empirical ones.
# As expected, manual pipetting has more variability than the 
# pipetting machine, but the difference is actually not large.
# The manual distribution looks cleaner and more normal.

pathlength %>%
  group_by(Pipetting) %>%
  summarise(Pathlength_mean = mean(Pathlength),
            Pathlength_sd = sd(Pathlength))
# Again the difference in sd is only 0.00069 cm, but the mean
# estimate for manual pipetting is lower.

pathlength_empirical <- pathlength %>%
  filter(Pipetting == "Manual") %$%
  mean(Pathlength)
pathlength_theoretical - pathlength_empirical
# The theoretical pathlength is 0.039 cm or 0.39 mm greater than
# the empirical pathlength. This is likely because it assumes a
# flat surface but water and 96% ethanol both have a concave meniscus.
# I know that ethanol has less of a meniscus than water, so it is
# reasonable to assume that 96% ethanol is at a halfway point
# between the theoretical (no meniscus) and empirical (water meniscus) 
# estimates.
pathlength_mean <- (pathlength_theoretical + pathlength_empirical) / 2
pathlength_mean # 0.55 cm

# 2.4.2 Volume-based concentration ####
pigments %<>%
  mutate(Concentration_mg_L = GPS_Model %>%
           map(~ pigment.concentration(.x, pathl = pathlength_mean) %>%
                 as_tibble_row()))

# 2.4.3 Correct concentrations ####
# Elution volume and mass varied across extractions and of course
# solutions are diluted to varying extents. Starting with my mg L^-1
# concentration, this is converted to mg by multiplying by the elution
# volume in L, then dividing by the mass in g to get mg g^-1. But 
# dilution increases the effective elution volume, so I need to multiply 
# by the dilution factor.
pigments %<>%
  mutate(Dilution = Concentration %>% 
           str_extract("\\d+") %>%
           as.numeric() %>%
           { (./100)^-1 }, # Dilution factor
         Volume = if_else(Date %>% year() == 2023,
                          1.3, 1.8)) %>% # mL
  rowwise() %>% # mL / mg = L / g = L g^-1 -> mg L^-1 * L g^-1 = mg g^-1
  mutate(Concentration_mg_g = list( 
    tibble( Concentration_mg_L * Volume * Dilution / Mass ) 
    )) %>%
  ungroup()
# The original method used sweep to convert concentrations and mutate(across())
# is also a tidyverse-friendly option:  
# Concentration_mg_g = list(
# Concentration_mg_L %>%
#   mutate(across(
#     everything(),
#     ~ .x * Volume * Dilution / Mass
#   ))
# )

# 2.4.4 Clean up ####
rm(spectra, Fig_S3, Fig_S3a, Fig_S3b, gaussian.peaks, sac.table,
   pathlength, sac, core, d_r, r_b, r_t, pathlength_empirical,
   pathlength_mean, pathlength_theoretical)

# 3. Summary ####
# 3.1 Re-nest ####
pigments %<>%
  select(-c(Concentration_mg_L, Spectrum, Spectrum_Plot, GPS_Model, 
            GPS_Residuals, Residuals_Plot, Residuals_Density_Plot,
            GPS_Fit, GPS_Plot)) %>%
  unnest(Concentration_mg_g) %>%
  group_by(Name, Date, Plate, ID, Concentration, Dilution, Volume, Mass) %>%
  nest(.key = "Technical_Data") %>%
  ungroup()

# 3.2 Sum concentration ####
pigments %<>%
  mutate(
    Technical_Data = Technical_Data %>%
      map(
        ~ .x %>%
          rowwise() %>%
          mutate(Total = sum( c_across(Anth:Zea) ),
                 Chl.c = Chl.c1 + Chl.c2) %>%
          ungroup() %>%
          select(-c(Chl.c1, Chl.c2))
      )
  )

# 3.3 Calculate technical mean ####
pigments %<>%
  mutate(
    Technical_Mean = Technical_Data %>%
      map(
        ~ .x %>%
          summarise( across(Anth:Chl.c, mean) )
      )
  )

# 3.4 Add predictors ####
pigments %<>%
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
         ) %T>%
  print(n = 79)

# 3.5 Summarise by treatment and season ####
pigments_summary <- pigments %>%
  unnest(Technical_Mean) %>%
  mutate(Treatment = Treatment %>% fct_relevel("Kelp"),
         Season = Season %>% fct_relevel("Spring", "Summer")) %>%
  group_by(Season, Treatment) %>%
  summarise( across( Anth:Chl.c, 
                     ~ str_c( mean(.) %>% 
                                signif(2) , 
                              " ± " , 
                              sd(.) %>% 
                                signif(2) ) ),
             n = n() ) %>%
  ungroup() %>%
  mutate(Treatment = str_glue("{Treatment} (*n* = {n})")) %>%
  select(-n) %>%
  select(Season, Treatment, Chl.a, Phe.a, Pheide.a, Chl.c, Fuco, bb.Car, Viola, Anth, Zea, Total)
pigments_summary

# 3.6 Build table for supplement ####
require(gt)
pigments_gt <- pigments_summary %>%
  gt(groupname_col = "Season", 
     rowname_col = "Treatment") %>%
  tab_header(
    title = md("**Table S3**. Transformation of kelp (*Laminaria hyperborea*) 
               photosyntehtic pigments (mean ± s.d. mg g<sup>−1</sup>) by sea 
               urchins (*Strongylocentrotus droebachiensis*).")
  ) %>%
  tab_spanner(
    label = "Chlorophylls",
    columns = c(Chl.a, Phe.a, Pheide.a, Chl.c)
  ) %>%
  tab_spanner(
    label = "Major carotenoids",
    columns = c(Fuco, bb.Car)
  ) %>%
  tab_spanner(
    label = "Xanthophyll cycle pigments",
    columns = c(Viola, Anth, Zea)
  ) %>%
  cols_label(
    Chl.a = md("Chlorophyll *a*"),
    Phe.a = md("Pheophytin *a*<sup>*</sup>"),
    Pheide.a = md("Pheophorbide *a*<sup>*</sup>"),
    Chl.c = md("Chlorophyll *c*"),
    Fuco = "Fucoxanthin",
    bb.Car = "β-carotene",
    Viola = "Violaxanthin",
    Anth = "Antheraxanthin",
    Zea = "Zeaxanthin",
    Total = "Total pigment"
  ) %>%
  tab_source_note(
    source_note = md("Each *n* is the mean of a technical triplicate and 
                     represents a mature sporophyte or the sea urchin faeces 
                     derived from it. The *nnls* modelling error in each techincal 
                     replicate was not propagated for these estimates.")
  ) %>%
  tab_source_note(
    source_note = md("<sup>*</sup>Chlorophyll *a* degradation products, 
                     collectively called pheopigments.")
  ) %>%
  opt_table_font(font = "Futura") %>%
  tab_options(
    table.align = "left",
    table.font.size = px(10),
    source_notes.font.size = px(10),
    table.font.color = "black",
    data_row.padding = px(2),
    row_group.padding = px(2),
    column_labels.padding = px(2),
    source_notes.padding = px(2),
    heading.padding = px(4),
    stub.border.width = px(0),
    table.border.top.width = px(0),
    table.border.bottom.width = px(0),
    table_body.hlines.width = px(0),
    table_body.vlines.width = px(0),
    table_body.border.top.width = px(1),
    table_body.border.top.color = "black",
    table_body.border.bottom.width = px(1),
    table_body.border.bottom.color = "black",
    heading.border.bottom.width = px(1),
    heading.border.bottom.color = "black",
    column_labels.border.top.width = px(1),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    row_group.border.top.width = px(0),
    row_group.border.top.color = "black",
    row_group.border.bottom.width = px(0),
    row_group.border.bottom.color = "black"
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_stub()
  ) %>%
  tab_style(
    style = cell_text(size = px(12), 
                      align = "left"),
    locations = cells_title()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold",
                      align = "left"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_spanners()
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      weight = px(0)
    ),
    locations = cells_column_labels()
  ) %>%
  text_transform(
    fn = function(x) md(x),
    locations = cells_stub()
  ) %>%
  fmt_markdown(
    columns = Treatment
  )
  
pigments_gt %>%
  gtsave(filename = "Tab_S3.pdf", path = "Figures")
# Doesn't fit on portrait and Can't save to landscape.

pigments_gt %>%
  gtsave(filename = "Tab_S3.html", path = "Figures")
require(pagedown)
chrome_print(
  input = here("Figures", "Tab_S3.html"),
  output = here("Figures", "Tab_S3.pdf"),
  options = list(paperWidth = 29.7 / 2.54,
                 paperHeight = 21 / 2.54)
  )

# 3.7 Multivariate visualisation ####



# 4. Technical triplicate models ####
# 4.1 Sum pheopigments ####
pigments %<>% 
  mutate(Technical_Data = Technical_Data %>%
           map(~ .x %>% 
                 mutate(Pheopigments = Phe.a + Pheide.a) %>%
                 rename(Chlorophyll = Chl.a)))

# 4.2 Update technical mean ####
pigments %<>%
  mutate(
    Technical_Mean = Technical_Data %>%
      map(~ .x %>%
            summarise( across(Anth:Pheopigments, mean) ))
    )

# 4.3 Visualise data ####
pigments %<>%
  rowwise() %>%
  mutate(
    Technical_Chl_Plot = 
      list(
        Technical_Data %>%
          ggplot(aes(Well, Chlorophyll)) +
            geom_pointrange(aes(ymin = Chlorophyll - CV * Chlorophyll,
                                ymax = Chlorophyll + CV * Chlorophyll)) +
            geom_hline(data = Technical_Mean,
                       aes(yintercept = Chlorophyll)) +
            ggtitle(ID) +
            theme_minimal() +
            theme(panel.grid = element_blank())
        ),
    Technical_Pheo_Plot = 
      list(
        Technical_Data %>%
          ggplot(aes(Well, Pheopigments)) +
            geom_pointrange(aes(ymin = Pheopigments - CV * Pheopigments,
                                ymax = Pheopigments + CV * Pheopigments)) +
            geom_hline(data = Technical_Mean,
                       aes(yintercept = Pheopigments)) +
            ggtitle(ID) +
            theme_minimal() +
            theme(panel.grid = element_blank())
        ),
    Technical_Total_Plot = 
      list(
        Technical_Data %>%
          ggplot(aes(Well, Total)) +
            geom_pointrange(aes(ymin = Total - CV * Total,
                                ymax = Total + CV * Total)) +
            geom_hline(data = Technical_Mean,
                       aes(yintercept = Total)) +
            ggtitle(ID) +
            theme_minimal() +
            theme(panel.grid = element_blank())
        )
    ) %>%
  ungroup()

pigments %$% 
  wrap_plots(Technical_Chl_Plot) %>%
  ggsave(filename = "technical_chl_data.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 60, width = 60, units = "cm")

pigments %$% 
  wrap_plots(Technical_Pheo_Plot) %>%
  ggsave(filename = "technical_pheo_data.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 60, width = 60, units = "cm")

pigments %$% 
  wrap_plots(Technical_Total_Plot) %>%
  ggsave(filename = "technical_total_data.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 60, width = 60, units = "cm")

# 4.4 Stan model ####
require(cmdstanr)
technical_chl_model <- here("Biochemistry", "Pigments", "Stan", "technical_chl.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

technical_pheo_model <- here("Biochemistry", "Pigments", "Stan", "technical_pheo.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

technical_total_model <- here("Biochemistry", "Pigments", "Stan", "technical_total.stan") %>%
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
pigments %<>%
  mutate(
    Technical_Chl_Samples = Technical_Data %>%
      map2(
        Technical_Mean,
        ~ technical_chl_model$sample(
          data = .x %>%
            select(Chlorophyll, CV) %>%
            compose_data() %>%
            list_modify(Chlorophyll_mean = .y %>% 
                          pull(Chlorophyll)),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
          adapt_delta = 0.99, # force sampler to slow down
          max_treedepth = 12
        )
      ),
    Technical_Pheo_Samples = Technical_Data %>%
      map2(
        Technical_Mean,
        ~ technical_pheo_model$sample(
          data = .x %>%
            select(Pheopigments, CV) %>%
            compose_data() %>%
            list_modify(Pheopigments_mean = .y %>% 
                          pull(Pheopigments)),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
          adapt_delta = 0.99,
          max_treedepth = 12
        )
      ),
    Technical_Total_Samples = Technical_Data %>%
      map2(
        Technical_Mean,
        ~ technical_total_model$sample(
          data = .x %>%
            select(Total, CV) %>%
            compose_data() %>%
            list_modify(Total_mean = .y %>% 
                          pull(Total)),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
          adapt_delta = 0.99,
          max_treedepth = 12
        )
      )
    )
# Some divergent transitions, likely due to n = 3.

# 4.5 Model checks ####
# 4.5.1 Rhat ####
pigments$Technical_Chl_Samples %>%
  map(~ .x$summary() %>%
        mutate(rhat_check = rhat > 1.001)) %>%
  bind_rows() %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Less than 30% of rhat above 1.001. rhat = 1.00 ± 0.00494.

pigments$Technical_Pheo_Samples %>%
  map(~ .x$summary() %>%
        mutate(rhat_check = rhat > 1.001)) %>%
  bind_rows() %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Less than 30% of rhat above 1.001. rhat = 1.00 ± 0.00612.

pigments$Technical_Total_Samples %>%
  map(~ .x$summary() %>%
        mutate(rhat_check = rhat > 1.001)) %>%
  bind_rows() %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Less than 30% of rhat above 1.001. rhat = 1.00 ± 0.00715.

# 4.5.2 Chains ####
require(bayesplot)
pigments %<>%
  rowwise() %>%
  mutate(
    Technical_Chl_Chains = 
      list(
        Technical_Chl_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(ID)
        ),
    Technical_Pheo_Chains = 
      list(
        Technical_Pheo_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(ID)
      ),
    Technical_Total_Chains = 
      list(
        Technical_Total_Samples$draws(format = "df") %>%
          mcmc_rank_overlay() +
          ggtitle(ID)
      )
    ) %>%
  ungroup()

pigments %$% 
  ( wrap_plots(Technical_Chl_Chains) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom") ) %>%
  ggsave(filename = "technical_chl_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 80, width = 80, units = "cm")

pigments %$% 
  ( wrap_plots(Technical_Pheo_Chains) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom") ) %>%
  ggsave(filename = "technical_pheo_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 80, width = 80, units = "cm")

pigments %$% 
  ( wrap_plots(Technical_Total_Chains) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom") ) %>%
  ggsave(filename = "technical_total_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 80, width = 80, units = "cm")
# Chains look ok (a few convergence issues here and there).

# 4.6 Prior-posterior comparison ####
# 4.6.1 Sample priors ####
require(truncnorm) # R doesn't have a native turncated normal.
pigments %<>%
  mutate(
    Technical_Chl_Prior = Technical_Mean %>%
      map(
        ~ tibble(.chain = 1:8 %>% rep(each = 1e4),
                 .iteration = 1:1e4 %>% rep(times = 8),
                 .draw = 1:8e4,
                 mu = rtruncnorm( n = 8e4 , 
                                  mean = .x %>% pull(Chlorophyll) , 
                                  sd = 0.08 * .x %>% pull(Chlorophyll) , 
                                  a = 0 ), 
                 sigma = rexp( 8e4 , 15 / ( .x %>% pull(Chlorophyll) + 0.4 ) ))
      ),
    Technical_Pheo_Prior = Technical_Mean %>%
      map(
        ~ tibble(.chain = 1:8 %>% rep(each = 1e4),
                 .iteration = 1:1e4 %>% rep(times = 8),
                 .draw = 1:8e4,
                 mu = rtruncnorm( n = 8e4 , 
                                  mean = .x %>% pull(Pheopigments) , 
                                  sd = 0.08 * .x %>% pull(Pheopigments) , 
                                  a = 0 ), 
                 sigma = rexp( 8e4 , 15 / ( .x %>% pull(Pheopigments) + 0.4 ) ))
      ),
    Technical_Total_Prior = Technical_Mean %>%
      map(
        ~ tibble(.chain = 1:8 %>% rep(each = 1e4),
                 .iteration = 1:1e4 %>% rep(times = 8),
                 .draw = 1:8e4,
                 mu = rtruncnorm( n = 8e4 , 
                                  mean = .x %>% pull(Total) , 
                                  sd = 0.08 * .x %>% pull(Total) , 
                                  a = 0 ), 
                 sigma = rexp( 8e4 , 15 / ( .x %>% pull(Total) + 0.4 ) ))
      )
  )

# 4.6.2 Extract posteriors ####
pigments %<>%
  mutate(
    Technical_Chl_Posterior = Technical_Chl_Samples %>%
      map(
        ~ .x %>% spread_draws(mu, sigma)
      ),
    Technical_Pheo_Posterior = Technical_Pheo_Samples %>%
      map(
        ~ .x %>% spread_draws(mu, sigma)
      ),
    Technical_Total_Posterior = Technical_Total_Samples %>%
      map(
        ~ .x %>% spread_draws(mu, sigma)
      )
  )

# 4.6.3 Plot comparison ####
source("functions.R")
pigments %<>%
  rowwise() %>%
  mutate(
    Technical_Chl_Prior_Posterior =
      list(
          Technical_Chl_Prior %>%
            bind_rows(Technical_Chl_Posterior) %>%
            mutate(distribution = c("prior", "posterior") %>%
                     rep(each = 8e4) %>% fct()) %>%
            pivot_longer(cols = c(mu, sigma),
                         names_to = ".variable",
                         values_to = ".value") %>%
            prior_posterior_plot() +
            ggtitle(ID)
      ),
    Technical_Pheo_Prior_Posterior =
      list(
          Technical_Pheo_Prior %>%
            bind_rows(Technical_Pheo_Posterior) %>%
            mutate(distribution = c("prior", "posterior") %>%
                     rep(each = 8e4) %>% fct()) %>%
            pivot_longer(cols = c(mu, sigma),
                         names_to = ".variable",
                         values_to = ".value") %>%
            prior_posterior_plot() +
            ggtitle(ID)
      ),
    Technical_Total_Prior_Posterior =
      list(
          Technical_Total_Prior %>%
            bind_rows(Technical_Total_Posterior) %>%
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

pigments %$% 
  ( wrap_plots(Technical_Chl_Prior_Posterior) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom") ) %>%
  ggsave(filename = "technical_chl_prior_posterior.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 80, width = 80, units = "cm")

pigments %$% 
  ( wrap_plots(Technical_Pheo_Prior_Posterior) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom") ) %>%
  ggsave(filename = "technical_pheo_prior_posterior.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 80, width = 80, units = "cm")

pigments %$% 
  ( wrap_plots(Technical_Total_Prior_Posterior) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom") ) %>%
  ggsave(filename = "technical_total_prior_posterior.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 80, width = 80, units = "cm")
# Posteriors are a bit noisy/spiky in some places, but generally fine.

# 4.7 Prediction ####
# 4.7.1 Calculate predictions ####
pigments %<>%
  mutate(
    Technical_Chl_Posterior = Technical_Chl_Posterior %>%
      map(
        ~ .x %>%
          mutate(Chlorophyll = rtruncnorm( n = n() , mean = mu , sd = sigma , a = 0 ))
      ),
    Technical_Pheo_Posterior = Technical_Pheo_Posterior %>%
      map(
        ~ .x %>%
          mutate(Pheopigments = rtruncnorm( n = n() , mean = mu , sd = sigma , a = 0 ))
      ),
    Technical_Total_Posterior = Technical_Total_Posterior %>%
      map(
        ~ .x %>%
          mutate(Total = rtruncnorm( n = n() , mean = mu , sd = sigma , a = 0 ))
      ),
    Samples_Data = Technical_Chl_Posterior %>%
      map2(
        Technical_Pheo_Posterior,
        ~ .x %>% 
            select(-c(mu, sigma)) %>%
            full_join(.y %>%
                        select(-c(mu, sigma)),
                      by = c(".chain", ".iteration", ".draw"))
      ) %>%
      map2(
        Technical_Total_Posterior,
        ~ .x %>%
            full_join(.y %>%
                        select(-c(mu, sigma)),
                      by = c(".chain", ".iteration", ".draw"))
      )
  )

# 4.7.2 Plot prediction ####
require(ggdist)
pigments %<>%
  rowwise() %>%
  mutate(
    Technical_Chl_Prediction = 
      list(
          ggplot() +
            geom_pointrange(data = Technical_Data,
                            aes(Well, Chlorophyll, 
                                ymin = Chlorophyll - CV * Chlorophyll,
                                ymax = Chlorophyll + CV * Chlorophyll)) +
            geom_hline(data = Technical_Mean,
                       aes(yintercept = Chlorophyll)) +
            stat_eye(data = Technical_Chl_Posterior,
                     aes(2, mu), alpha = 0.5,
                     point_interval = NULL) +
            stat_eye(data = Technical_Chl_Posterior,
                     aes(2, Chlorophyll), alpha = 0.5,
                     point_interval = NULL) +
            coord_cartesian(ylim = Technical_Mean %$%
                              c( (1 - 0.4) * Chlorophyll,
                                 (1 + 0.4) * Chlorophyll )) +
            ggtitle(ID) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      ),
    Technical_Pheo_Prediction = 
      list(
          ggplot() +
            geom_pointrange(data = Technical_Data,
                            aes(Well, Pheopigments, 
                                ymin = Pheopigments - CV * Pheopigments,
                                ymax = Pheopigments + CV * Pheopigments)) +
            geom_hline(data = Technical_Mean,
                       aes(yintercept = Pheopigments)) +
            stat_eye(data = Technical_Pheo_Posterior,
                     aes(2, mu), alpha = 0.5,
                     point_interval = NULL) +
            stat_eye(data = Technical_Pheo_Posterior,
                     aes(2, Pheopigments), alpha = 0.5,
                     point_interval = NULL) +
            coord_cartesian(ylim = Technical_Mean %$%
                              c( (1 - 0.4) * Pheopigments,
                                 (1 + 0.4) * Pheopigments )) +
            ggtitle(ID) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      ),
    Technical_Total_Prediction = 
      list(
          ggplot() +
            geom_pointrange(data = Technical_Data,
                            aes(Well, Total, 
                                ymin = Total - CV * Total,
                                ymax = Total + CV * Total)) +
            geom_hline(data = Technical_Mean,
                       aes(yintercept = Total)) +
            stat_eye(data = Technical_Total_Posterior,
                     aes(2, mu), alpha = 0.5,
                     point_interval = NULL) +
            stat_eye(data = Technical_Total_Posterior,
                     aes(2, Total), alpha = 0.5,
                     point_interval = NULL) +
            coord_cartesian(ylim = Technical_Mean %$%
                              c( (1 - 0.4) * Total,
                                 (1 + 0.4) * Total )) +
            ggtitle(ID) +
            theme_minimal() +
            theme(panel.grid = element_blank())
      )
  ) %>%
  ungroup()

pigments %$% 
  wrap_plots(Technical_Chl_Prediction) %>%
  ggsave(filename = "technical_chl_prediction.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 60, width = 60, units = "cm")

pigments %$% 
  wrap_plots(Technical_Pheo_Prediction) %>%
  ggsave(filename = "technical_pheo_prediction.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 60, width = 60, units = "cm")

pigments %$% 
  wrap_plots(Technical_Total_Prediction) %>%
  ggsave(filename = "technical_total_prediction.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 60, width = 60, units = "cm")

# 4.7.3 Summarise for models ####
pigments %<>%
  mutate(
    Samples_Data_Summary = Samples_Data %>%
      map(
        ~ .x %>%
          summarise(Chlorophyll_mean = mean(Chlorophyll),
                    Chlorophyll_sd = sd(Chlorophyll),
                    Pheopigments_mean = mean(Pheopigments),
                    Pheopigments_sd = sd(Pheopigments),
                    Total_mean = mean(Total),
                    Total_sd = sd(Total))
      )
  )

# 5. Pigment models ####
# 5.1 Total pigment ####
# 5.1.1 Visualise ####
( pigments %>%
    select(Name, ID, Samples_Data) %>%
    unnest(cols = Samples_Data) %>%
    ggplot(aes(Total, ID)) +
      geom_vline(xintercept = 0) +
      stat_slab(n = 2e3, height = 10, 
                colour = "black", linewidth = 0.1) +
      facet_grid(rows = vars(Name), scales = "free") +
      coord_cartesian(xlim = c(0, 2)) +
      theme_minimal() +
      theme(panel.grid = element_blank()) ) %>%
  ggsave(filename = "samples_data.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Pigments", "Plots"),
         height = 40, width = 20, units = "cm")

( pigments %>%
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
  ggsave(filename = "samples_data_summary.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 40, width = 20, units = "cm")

# Here are the data I'll pass to the model:
pigments %>%
  select(Treatment, Season, Individual, Samples_Data_Summary) %>%
  unnest(cols = Samples_Data_Summary) %>%
  print(n = 79)

# 5.1.2 Prior simulation ####
# Total photosynthetic pigment content for Laminaria hyperborea is
# exprected between 1.95963 ± 0.10126 mg g^-1 (Wright & Foggo 2021,
# doi: 10.3354/meps13886) and 1.27 mg g^-1 (Wright & Kregting 2023,
# doi: 10.1007/s00227-023-04289-y), so around 1.6 mg g^-1.

tibble(n = 1:1e5,
       mu_log = rnorm( 1e5 , log(1.6) , 0.4 ),
       theta = rexp( 1e5, 5 ),
       mu = exp(mu_log),
       P = rgamma( 1e5 , mu / theta , 1 / theta )) %>%
  pivot_longer(cols = c(mu, P), 
               names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, parameter)) +
    geom_vline(xintercept = c(0, 2)) +
    stat_slab(alpha = 0.5, height = 2, n = 3e3) +
    coord_cartesian(expand = F,
                    xlim = c(-1, 3)) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 4.4 Stan models ####
phenol_c_model <- here("Biochemistry", "Phenol", "Stan", "phenol_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

phenol_nc_model <- here("Biochemistry", "Phenol", "Stan", "phenol_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

phenol_c_samples <- phenol_c_model$sample(
          data = phenol %>%
            select(Treatment, Season, Individual, Samples_Data_Summary) %>%
            unnest(cols = Samples_Data_Summary) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

phenol_nc_samples <- phenol_nc_model$sample(
          data = phenol %>%
            select(Treatment, Season, Individual, Samples_Data_Summary) %>%
            unnest(cols = Samples_Data_Summary) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
        )

# 4.5 Model checks ####
# 4.5.1 Rhat ####
phenol_c_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# All rhat above 1.001. rhat = 1.03 ± 0.0207.

phenol_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Most rhat above 1.001. rhat = 1.00 ± 0.00224.

# Plot comparison between centred and non-centred parameterisation.
phenol_nc_samples$summary() %>%
  left_join(phenol_c_samples$summary(),
            by = "variable") %>%
  rename(rhat_nc = rhat.x, rhat_c = rhat.y) %>%
  ggplot(aes(rhat_c, rhat_nc)) +
    geom_abline(slope = 1) +
    geom_point() +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Warning because z-scores are dropped as they have no equivalent in
# the centred parameteristion. The non-centred model is better.

# 4.5.2 Chains ####
phenol_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "phenol_c_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains look ok.

phenol_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "phenol_nc_chains.pdf", device = cairo_pdf, 
         path = here("Biochemistry", "Phenol", "Plots"),
         height = 40, width = 40, units = "cm")
# Chains look good.

# 4.5.3 Pairs ####
phenol_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_t[1]", 
                      "alpha_s[1,1]", "sigma_s[1]",
                      "alpha_i[1,1]", "sigma_i[1]", 
                      "theta[1]"))
# Correlation between alpha_t and alpha_s for Treatment 1 (Faeces).
# alpha_s is being pulled far away from its zero mean and there is
# a strong inflation of sigma_s.
# 5.2 Chlorophyll vs. Pheopigments ####

