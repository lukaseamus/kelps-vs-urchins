# 1. Prepare data ####
require(tidyverse)
require(here)
# 1.1 Load carbon and nitrogen data ####
C_N <- 
  here("Biochemistry", "C_N", "RDS", "C_N.rds") %>%
  read_rds()
C_prior_posterior <- 
  here("Biochemistry", "C_N", "RDS", "C_prior_posterior.rds") %>%
  read_rds()
C_diff  <- 
  here("Biochemistry", "C_N", "RDS", "C_diff.rds") %>%
  read_rds()
N_prior_posterior <- 
  here("Biochemistry", "C_N", "RDS", "N_prior_posterior.rds") %>%
  read_rds()
N_diff <- 
  here("Biochemistry", "C_N", "RDS", "N_diff.rds") %>%
  read_rds()
C_N_prior_posterior <- 
  here("Biochemistry", "C_N", "RDS", "C_N_prior_posterior.rds") %>%
  read_rds()
C_N_diff <- 
  here("Biochemistry", "C_N", "RDS", "C_N_diff.rds") %>%
  read_rds()

# 1.2 Load phenol data ####
phenol <- 
  here("Biochemistry", "Phenol", "RDS", "phenol.rds") %>% 
  read_rds()
phenol_ID_dens <- 
  here("Biochemistry", "Phenol", "RDS", "phenol_ID_dens.rds") %>% 
  read_rds()
phenol_prior_posterior <- 
  here("Biochemistry", "Phenol", "RDS", "phenol_prior_posterior.rds") %>% 
  read_rds()
phenol_diff <- 
  here("Biochemistry", "Phenol", "RDS", "phenol_diff.rds") %>%
  read_rds()

# 1.3 Load pigments data ####
pigments <- 
  here("Biochemistry", "Pigments", "RDS", "pigments.rds") %>% 
  read_rds()
total_ID_dens <- 
  here("Biochemistry", "Pigments", "RDS", "total_ID_dens.rds") %>% 
  read_rds()
chlvspheo_ID_dens <- 
  here("Biochemistry", "Pigments", "RDS", "chlvspheo_ID_dens.rds") %>% 
  read_rds()
total_prior_posterior <- 
  here("Biochemistry", "Pigments", "RDS", "total_prior_posterior.rds") %>% 
  read_rds()
chlvspheo_prior_posterior <- 
  here("Biochemistry", "Pigments", "RDS", "chlvspheo_prior_posterior.rds") %>% 
  read_rds()
total_diff <- 
  here("Biochemistry", "Pigments", "RDS", "total_diff.rds") %>% 
  read_rds()
chlvspheo_diff <- 
  here("Biochemistry", "Pigments", "RDS", "chlvspheo_diff.rds") %>% 
  read_rds()

# 1.4 Rescale ID densities ####
phenol_ID_dens <- phenol %>%
  select(Treatment, Season, Individual, ID, Samples_Data) %>%
  unnest(cols = Samples_Data) %>%
  group_by(Treatment, Season, ID) %>% 
  # Starts from negative because some densities are not strictly positive.
  reframe(x = density(Concentration, n = 2^10, from = -0.1, to = 2)$x,
          y = density(Concentration, n = 2^10, from = -0.1, to = 2)$y) %>%
  group_by(ID) %>%
  mutate(y_area = y * 0.014 / ( sum(y) * ( x[2] - x[1] ) )) %>%
  ungroup() %>%
  filter(y > 0.1) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = c(x, x %>% rev()),
          y_area = c(y_area, -y_area %>% rev())) %>%
  ungroup()

total_ID_dens <- pigments %>%
  select(Treatment, Season, Individual, ID, Samples_Data) %>%
  unnest(cols = Samples_Data) %>%
  select(-c(Chlorophyll, Pheopigments)) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = density(Total, n = 2^10, from = 0, to = 3)$x,
          y = density(Total, n = 2^10, from = 0, to = 3)$y) %>%
  group_by(ID) %>%
  mutate(y_area = y * 0.02 / ( sum(y) * ( x[2] - x[1] ) )) %>%
  ungroup() %>%
  filter(y > 0.095) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = c(x, x %>% rev()),
          y_area = c(y_area, -y_area %>% rev())) %>%
  ungroup()

chlvspheo_ID_dens <- pigments %>%
  select(Treatment, Season, Individual, ID, Samples_Data) %>%
  unnest(cols = Samples_Data) %>%
  mutate(Proportion = Chlorophyll / ( Chlorophyll + Pheopigments ) * 100) %>%
  select(-c(Total, Chlorophyll, Pheopigments)) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = density(Proportion, n = 2^10, from = 0, to = 100)$x,
          y = density(Proportion, n = 2^10, from = 0, to = 100)$y) %>%
  group_by(ID) %>%
  mutate(y_area = y * 0.7 / ( sum(y) * ( x[2] - x[1] ) )) %>%
  ungroup() %>%
  filter(y > 0.003) %>%
  group_by(Treatment, Season, ID) %>%
  reframe(x = c(x, x %>% rev()),
          y_area = c(y_area, -y_area %>% rev())) %>%
  ungroup()

# 2. Visualise ####
# 2.1 Define theme ####
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0, 0.5, 0, 0, unit = "cm"),
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

# 2.2 Carbon ####
# 2.2.1 Prediction ####
Fig_2a_left_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = C, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.4, size = 2.5, height = 0.36) +
  stat_density_ridges(data = C_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 60, rel_min_height = 0.001, 
                      bandwidth = 0.6, scale = 1.2, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 20),
                     oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE, order = 1)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15), # circle, triangle, square
                     # override grey fill legend shapes because they can be confused with prior
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Carbon content (%)") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

# 2.2.2 Difference ####
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
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 2,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.75, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 1,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.35, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 2,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.32, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -50, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-50, 50), oob = scales::oob_keep,
                     breaks = seq(-50, 50, 25),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ carbon content (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 2.3 Nitrogen ####
# 2.3.1 Prediction ####
Fig_2a_middle_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = N, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.4, size = 2.5, height = 0.36) +
  stat_density_ridges(data = N_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 3, rel_min_height = 0.001, 
                      bandwidth = 0.03, scale = 1.2, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 3), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE, order = 1)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15),
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Nitrogen content (%)") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

# 2.3.2 Difference ####
Fig_2a_middle_bottom <- N_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.04,
                      from = -2, to = 2, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 1,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.7, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 2,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.65, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 1,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.3, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 2,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.35, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -2, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ nitrogen content (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 2.4 Carbon-nitrogen ####
# 2.4.1 Prediction ####
Fig_2a_right_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = C_N, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.4, size = 2.5, height = 0.36) +
  stat_density_ridges(data = C_N_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 120, rel_min_height = 0.001, 
                      bandwidth = 1.2, scale = 1.2, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 120), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE, order = 1)) +
  scale_colour_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15),
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Carbon-nitrogen ratio") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

# 2.4.2 Difference ####
Fig_2a_right_bottom <- C_N_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 2.4,
                      from = -120, to = 120, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 1,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 2.4, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 2,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.7, vjust = 0,
                   n = 2^10, bw = 2.4, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 1,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.41, vjust = 0,
                   n = 2^10, bw = 2.4, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 2,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.43, vjust = 0,
                   n = 2^10, bw = 2.4, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -120, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-120, 120), breaks = seq(-120, 120, 60),
                     oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ carbon-nitrogen ratio") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 2.5 Phenol ####
# 2.5.1 Prediction ####
Fig_2b_top <- ggplot() +
  geom_polygon(data = phenol_ID_dens %>% # Stratify by Treatment
                 mutate(y_area = y_area + if_else(Treatment == "Faeces", 0.5, 1.5)) %>%
                 group_by(ID) %>% # Jitter
                 mutate(y_area = y_area + runif( 1 , -0.36 , 0.36 )),
               aes(x = x, y = y_area, group = ID, 
                   fill = Treatment), 
               alpha = 0.4) +
  stat_density_ridges(data = phenol_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), colour = NA, n = 2^10,
                      from = 0, to = 2, rel_min_height = 0.001, 
                      bandwidth = 0.02, scale = 4, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5),
                     labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1)),
                     oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = "none") +
  xlab("Phenolic content (%)") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

# 2.5.2 Difference ####
Fig_2b_bottom <- phenol_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.04,
                      from = -2, to = 2, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.515 + 1,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura",
                   size = 3.5, hjust = 0.84, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.515 + 2,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura",
                   size = 3.5, hjust = 0.805, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.515 + 1, 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.32, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.515 + 2, 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -2, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                     breaks = seq(-2, 2, 1),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ phenolic content (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 2.6 Total pigment ####
# 2.6.1 Prediction ####
Fig_2c_left_top <- ggplot() +
  geom_polygon(data = total_ID_dens %>% # Stratify by Treatment
                 mutate(y_area = y_area + if_else(Treatment == "Faeces", 0.5, 1.5)) %>%
                 group_by(ID) %>% # Jitter
                 mutate(y_area = y_area + runif( 1 , -0.36 , 0.36 )),
               aes(x = x, y = y_area, group = ID, 
                   fill = Treatment), 
               alpha = 0.4) +
  stat_density_ridges(data = total_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), colour = NA, n = 2^10,
                      from = 0, to = 3, rel_min_height = 0.001, 
                      bandwidth = 0.03, scale = 1.2, alpha = 0.7) +
  geom_textpath(data = pigments %>%
                  select(ID, Treatment, Season, Samples_Data_Summary) %>%
                  unnest(cols = Samples_Data_Summary) %>%
                  filter(Season == "Spring" & Treatment == "Kelp") %>%
                  summarise(min = min(Total_mean - Total_sd),
                            max = max(Total_mean + Total_sd)) %$%
                  tibble(x = c(rep(min, 2), rep(max, 2)),
                         y = c(2.1, 2.2, 2.2, 2.1)),
                aes(x = x, y = y, label = "Spring"),
                colour = "#dabc23", linejoin = "mitre", lineend = "square",
                family = "Futura", size = 3.5, vjust = -0.1) +
  scale_x_continuous(limits = c(0, 3), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = "none") +
  xlab(expression("Total pigment (mg g"^-1*")")) +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme + # Adjust the title margin to counteract the superscript.
  theme(axis.title.x = element_text(margin = margin(t = 0)))

# 2.6.2 Difference ####
Fig_2c_left_bottom <- total_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.04,
                      from = -2, to = 2, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 1.07 + 1,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura",
                   size = 3.5, hjust = 0.71, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 1.07 + 2,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura",
                   size = 3.5, hjust = 0.68, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 1.07 + 1, 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.35, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 1.07 + 2, 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.38, vjust = 0,
                   n = 2^10, bw = 0.04, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -2, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                     breaks = seq(-2, 2, 1),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab(expression("Δ total pigment (mg g"^-1*")")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme + # The previous plot also affected this one in the same way.
  theme(axis.title.x = element_text(margin = margin(t = 0)))

# 2.7 Chlorophyll vs. pheopigments ####
# 2.7.1 Prediction ####
Fig_2c_right_top <- ggplot() +
  geom_polygon(data = chlvspheo_ID_dens %>% # Stratify by Treatment
                 mutate(y_area = y_area + if_else(Treatment == "Faeces", 0.5, 1.5)) %>%
                 group_by(ID) %>% # Jitter
                 mutate(y_area = y_area + runif( 1 , -0.36 , 0.36 )),
               aes(x = x, y = y_area, group = ID, 
                   fill = Treatment), 
               alpha = 0.4) +
  stat_density_ridges(data = chlvspheo_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), colour = NA, n = 2^10,
                      from = 0, to = 100, rel_min_height = 0.001, 
                      bandwidth = 1, scale = 1.2, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#dabc23", "#b5b8ba"),
                    guide = "none") +
  xlab(expression("Intact chlorophyll (%)")) +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

# 2.7.2 Difference ####
Fig_2c_right_bottom <- chlvspheo_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 1.6,
                      from = -80, to = 80, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 42.5 + 1,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura",
                   size = 3.5, hjust = 0.9, vjust = 0,
                   n = 2^10, bw = 1.6, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 42.5 + 2,
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura",
                   size = 3.5, hjust = 0.82, vjust = 0,
                   n = 2^10, bw = 1.6, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 42.5 + 1, 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.37, vjust = 0,
                   n = 2^10, bw = 1.6, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 42.5 + 2, 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.36, vjust = 0,
                   n = 2^10, bw = 1.6, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -80, y = c(1, 2),
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = -0.2, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-80, 80), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.7), alpha("#dabc23", 0.7)),
                    guide = "none") +
  xlab("Δ intact chlorophyll (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 2.8 Combine ####
Fig_2 <- ( ( Fig_2a_left_top | Fig_2a_middle_top | Fig_2a_right_top ) /
           ( ( Fig_2a_left_bottom | Fig_2a_middle_bottom | Fig_2a_right_bottom ) +
                # margin allows finer adjustment than plot_spacer
                theme(plot.margin = margin(0.2, 0.5, 0, 0, unit = "cm")) ) /
           plot_spacer() / # create space between the two grouped rows
           ( Fig_2b_top | Fig_2c_left_top | Fig_2c_right_top ) / 
           ( ( Fig_2b_bottom | Fig_2c_left_bottom | Fig_2c_right_bottom ) +
                theme(plot.margin = margin(0.06, 0.5, 0, 0, unit = "cm")) ) ) +
  plot_annotation(tag_levels = list(c("a", rep("", 5), "b", "c", rep("", 4))),
                  theme = theme(plot.margin = margin(0.75, 0.45, 0.5, 0.5, unit = "cm"))) +
  plot_layout(heights = c(1, 0.5, 0.2, 1, 0.5), 
              guides = "collect") &
  theme(plot.tag = element_text(family = "Futura",
                                size = 15, face = "bold"),
        plot.tag.position = c(-0.04, 1.06),
        legend.position = "top",
        legend.justification = 0)

Fig_2 %>%
  ggsave(filename = "Fig_2.pdf", device = cairo_pdf, path = "Figures", 
         height = 21, width = 21, units = "cm")