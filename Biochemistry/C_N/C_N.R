# 1. Load data ####
# 1.1 Load raw data ####
require(tidyverse)
require(here)

cn <- here("Biochemistry", "C_N", "C_N.csv") %>%
  read_csv() %>%
  mutate(Date = Date %>% dmy(),
         ID = ID %>% fct(),
         Method = Method %>% fct())
cn

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

Fig_S1a <- cn %>%
  filter(year(Date) != 2025) %>%
  mutate(Oven = if_else(ID %>% str_detect("_oven"),
                        "Yes", "No") %>% fct(),
         ID = ID %>% str_remove("_oven") %>% fct()) %>%
  ggplot(aes(Mass, C, colour = Method, 
             shape = Oven, group = ID)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_path(alpha = 0.6, linewidth = 1) +
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


Fig_S1b <- cn %>%
  filter(year(Date) != 2025) %>%
  mutate(Oven = if_else(ID %>% str_detect("_oven"),
                        "Yes", "No") %>% fct(),
         ID = ID %>% str_remove("_oven") %>% fct()) %>%
  ggplot(aes(Mass, N, colour = Method, 
             shape = Oven, group = ID)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_path(alpha = 0.6, linewidth = 1) +
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

require(patchwork)
Fig_S1 <- ( ( Fig_S1a + 
                theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank()) ) / 
              Fig_S1b ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top",
        legend.justification = 0)

Fig_S1 %>%
  ggsave(filename = "Fig_S1.pdf", device = cairo_pdf, path = "Figures", 
         height = 21, width = 21, units = "cm")

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

cn %>%
  drop_na(d15N, d13C) %>%
  mutate(Treatment = if_else(ID %>% str_detect("f"),
                             "Faeces", "Kelp") %>% fct()) %>%
  ggplot(aes(d15N, d13C, colour = Treatment)) +
    geom_point(size = 3, shape = 16, alpha = 0.5) +
    mytheme

# Faeces and Kelp could not be teased apart with stable isotope analysis.

# 1.3 Prepare data ####
cn %>%
  filter(year(Date) != 2024) %>% # filter out third run
  left_join( # then join it by column with the other runs
    cn %>% 
      filter(year(Date) == 2024) %>%
      select(ID, N, C) %>%
      rename(N_3 = N, C_3 = C),
    by = "ID"
  ) %>% # re-dried estimates are removed by left_join
  mutate(N_mean = if_else(Method == "IRMS",
                          (N + N_3)/2, N),
         C_mean = if_else(Method == "IRMS",
                          (C + C_3)/2, C)) %>% pull(C_mean)

  
  
  
  group_by(ID) %>%
  mutate(N = if_else(Method == "IRMS", mean(N), N),
         C = if_else(Method == "IRMS", mean(C), C)) %>%
  ungroup()


