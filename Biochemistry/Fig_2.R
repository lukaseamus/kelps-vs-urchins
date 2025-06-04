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

Fig_2b_top <- ggplot() +
  geom_polygon(data = ID_dens %>% # Stratify by Treatment
                 mutate(y_area = y_area + if_else(Treatment == "Faeces", 0.6, 1.6)) %>%
                 group_by(ID) %>% # Jitter
                 mutate(y_area = y_area + runif( 1 , -0.35 , 0.35 )),
               aes(x = x, y = y_area, group = ID, 
                   fill = Treatment), 
               alpha = 0.2) +
  stat_density_ridges(data = phenol_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), colour = NA, n = 2^10,
                      from = 0, to = 2, rel_min_height = 0.001, 
                      bandwidth = 0.03, scale = 3, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5),
                     labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1)),
                     oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  xlab("Phenolic content (%)") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

require(geomtextpath)
Fig_2b_bottom <- phenol_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.03,
                      from = -2, to = 2, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.515 + 1,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura",
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 0.03, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.515 + 2,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura",
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 0.03, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.515 + 1, 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.35, vjust = 0,
                   n = 2^10, bw = 0.03, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.515 + 2, 
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura",
                   size = 3.5, hjust = 0.32, vjust = 0,
                   n = 2^10, bw = 0.03, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -2, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = 0, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                     breaks = seq(-2, 2, 1),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#c3b300", 0.6)),
                    guide = "none") +
  xlab("Difference (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

Fig_2b <- ( Fig_2b_top / Fig_2b_bottom ) +
  plot_layout(heights = c(1, 2/3))

Fig_2b %>%
  ggsave(filename = "Phenol.pdf", device = cairo_pdf, path = "Figures", 
         height = 12, width = 7, units = "cm")

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

# 5.2 Carbon ####
Fig_2a_left_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = C, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.5, size = 5, height = 0.4) +
  stat_density_ridges(data = C_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 60, rel_min_height = 0.001, 
                      bandwidth = 1, scale = 2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 20),
                     oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#c3b300"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15), # circle, triangle, square
                     # override grey fill legend shapes because they can be confused with prior
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Carbon content (%)") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

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
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 2,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 1,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 17 + 2,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 1, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -50, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = 0, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-50, 50), oob = scales::oob_keep,
                     breaks = seq(-50, 50, 25),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#c3b300", 0.6)),
                    guide = "none") +
  xlab("Difference (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 5.3 Nitrogen ####
Fig_2a_middle_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = N, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.5, size = 5, height = 0.4) +
  stat_density_ridges(data = N_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 3, rel_min_height = 0.001, 
                      bandwidth = 0.05, scale = 2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 3), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15),
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Nitrogen content (%)") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

Fig_2a_middle_bottom <- N_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 0.05,
                      from = -2, to = 2, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 1,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 0.05, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 2,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 0.05, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 1,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.22, vjust = 0,
                   n = 2^10, bw = 0.05, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 0.92 + 2,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.22, vjust = 0,
                   n = 2^10, bw = 0.05, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -2, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = 0, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#c3b300", 0.6)),
                    guide = "none") +
  xlab("Difference (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme

# 5.3 Carbon-nitrogen ####
# Add C:N to data.
C_N %<>% mutate(C_N = C / N)

Fig_2a_right_top <- ggplot() +
  geom_jitter(data = C_N %>%
                mutate(Season = Season %>% fct_relevel("Spring", "Summer")),
              aes(x = C_N, y = Treatment %>% as.numeric() - 0.5, 
                  colour = Treatment, shape = Season), 
              alpha = 0.5, size = 5, height = 0.4) +
  stat_density_ridges(data = C_N_prior_posterior %>%
                        mutate(Treatment = Treatment %>% fct_relevel("Faeces", "Kelp")),
                      aes(x = obs_new, y = Treatment %>% as.numeric(), 
                          fill = Treatment), 
                      colour = NA, n = 2^10,
                      from = 0, to = 120, rel_min_height = 0.001, 
                      bandwidth = 1.5, scale = 1.2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 120), oob = scales::oob_keep) +
  scale_fill_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = c("#7030a5", "#c3b300", "#b5b8ba"),
                      guide = "none") +
  scale_shape_manual(values = c(16, 17, 15),
                     guide = guide_legend(override.aes = list(shape = c(1, 2, 0)))) +
  xlab("Carbon-nitrogen ratio") +
  coord_cartesian(ylim = c(0, 4), expand = FALSE, clip = "off") +
  mytheme

Fig_2a_right_bottom <- C_N_diff %>% 
  filter(Parameter %in% c("mu_new", "obs_new")) %>%
  mutate(Parameter = Parameter %>% fct_relevel("obs_new")) %>%
  ggplot() +
  stat_density_ridges(aes(x = Difference, y = Parameter, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 2,
                      from = -150, to = 150, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 1,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 2, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 2,
                       label = label_Kelp),
                   colour = "#c3b300", family = "Futura", 
                   size = 3.5, hjust = 0.8, vjust = 0,
                   n = 2^10, bw = 2, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "obs_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 1,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 2, text_only = TRUE) +
  geom_textdensity(data = . %>% filter(Parameter == "mu_new"),
                   aes(x = Difference, y = after_stat(density) * 50 + 2,
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.25, vjust = 0,
                   n = 2^10, bw = 2, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("text", x = -150, y = c(1, 2), 
           label = c("italic(tilde('y'))", "italic('µ')"),
           hjust = 0, vjust = 0, family = "Futura", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(limits = c(-150, 150), oob = scales::oob_keep,
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c(alpha("#7030a5", 0.6), alpha("#c3b300", 0.6)),
                    guide = "none") +
  xlab("Difference (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme


Fig_2a <- ( ( Fig_2a_left_top / Fig_2a_left_bottom ) +
              plot_layout(heights = c(1, 0.3333333)) | 
            ( Fig_2a_middle_top / Fig_2a_middle_bottom ) +
              plot_layout(heights = c(1, 0.3333333)) |
            ( Fig_2a_right_top / Fig_2a_right_bottom ) +
              plot_layout(heights = c(1, 0.3333333)) ) 

Fig_2a %>%
  ggsave(filename = "Fig_2a.pdf", device = cairo_pdf, path = "Figures", 
         height = 12, width = 21, units = "cm")


