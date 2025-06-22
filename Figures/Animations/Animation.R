# 1. Sinking speed ####
# 1.1 Load data ####
require(tidyverse)
require(here)
speed <- 
  read_rds(here("Sinking", "RDS", "speed.rds"))
speed_prior_posterior <- 
  read_rds(here("Sinking", "RDS", "speed_prior_posterior.rds"))
speed_diff <- 
  read_rds(here("Sinking", "RDS", "speed_diff.rds"))

# 1.2 Add animation aesthetics to data points ####
require(magrittr)
speed %<>%
  mutate(y = -0.5 + runif( n() , -0.36 , 0.36 ), # Manual jitter
         alpha = 0.4, # Baseline transparency
         colour = if_else( # Colours
           Tissue == "Kelp",
           "#dabc23", "#7030a5"
         )) %>%
  rename(x = Speed) %>% # Rename for consistency
  select(-c(Mass, Area)) %T>% # Remove irrelevant variables
  print()

# 1.3 Build density polygons ####
speed_dens <- speed_prior_posterior %>%
  group_by(Tissue) %>% # Close polygons at zero and max with c(0, density, max).
  reframe(x = c(0, density(obs, n = 2^10, from = 0, to = 25, bw = 25 * 0.01)$x, 25),
          y = c(0, density(obs, n = 2^10, from = 0, to = 25, bw = 25 * 0.01)$y, 0)) %>%
  group_by(Tissue) %>% # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate(y = y * 2.6 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

# 1.4 Add animation aesthetics to polygons ####
speed_dens %<>%
  mutate(fill = case_when(
    Tissue == "Kelp" ~ "#dabc23", 
    Tissue == "Faeces" ~ "#7030a5",
    Tissue == "Prior" ~ "#b5b8ba"
  )) %T>%
  print()

# 1.5 Static plot ####
# Define custom theme
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

# Plot
ggplot() +
  geom_point(data = speed,
             aes(x = x, y = y, colour = colour, alpha = alpha),
             shape = 16, size = 2.5) +
  geom_polygon(data = speed_dens,
               aes(x = x, y = y, fill = fill)) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 25), ylim = c(-1, 2), 
                  expand = FALSE, clip = "off") +
  xlab(expression("Sinking speed (cm s"^-1*")")) +
  mytheme
# Looks fine

# 1.6 Tween points ####
# As per usual there is an uneven sample size, so not all points can be paired.
# Define enter/exit functions
enter <- function(data){
  data$alpha <- 0
  data
}

exit <- function(data){
  data$alpha <- 0
  data
}

# Tween
require(tweenr)
speed_ani <- bind_rows( # Filter for non-existent prior data to add prior
  tween_state(speed %>% filter(Tissue == "Prior"), 
              speed %>% filter(Tissue == "Kelp"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50),
  tween_state(speed %>% filter(Tissue == "Kelp"),
              speed %>% filter(Tissue == "Faeces"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150), # Ensure uniqueness of .frame
  tween_state(speed %>% filter(Tissue == "Faeces"),
              speed %>% filter(Tissue == "Prior"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300) 
) %T>%
  print() 
# .id is automatically assigned based on row number, i.e. arbitrary pairings.

# 1.7 Tween polygons ####
require(transformr)
speed_dens_ani <- bind_rows(
  tween_polygon(speed_dens %>% filter(Tissue == "Prior"),
                speed_dens %>% filter(Tissue == "Kelp"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(speed_dens %>% filter(Tissue == "Kelp"),
                speed_dens %>% filter(Tissue == "Faeces"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(speed_dens %>% filter(Tissue == "Faeces"),
                speed_dens %>% filter(Tissue == "Prior"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()
# Note that tween_polygon doesn't always work for densities,
# especially when they have uneven lengths. In those cases
# tween_path can be more stable. No enter and exit required here.

# 1.8 Dynamic plot ####
require(gganimate) # Needed to animate frames
require(ggnewscale) # Needed to add legend
require(ggtext) # Needed to use Markdown syntax in gif
( ggplot() +
  geom_point(data = speed_ani,
             aes(x = x, y = y, colour = colour, alpha = alpha),
             shape = 16, size = 2.5) +
  geom_polygon(data = speed_dens_ani,
               aes(x = x, y = y, fill = fill)) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_fill_identity() +
  new_scale_fill() +
  geom_polygon(data = speed_dens %>%
                 mutate(x = Inf), # Evict static plot from plotting space
               aes(x = x, y = y, fill = Tissue)) +
  scale_fill_manual(values = c("#b5b8ba", "#dabc23", "#7030a5")) +
  coord_cartesian(xlim = c(0, 25), ylim = c(-1, 2), 
                  expand = FALSE, clip = "off") +
  # Superscript for 12 pt font defaults to 8.4 pt font in expression().
  xlab("Sinking speed (cm s<sup><span style='font-size:8.4pt'>−1</span></sup>)") +
  transition_manual(.frame) +
  mytheme +
  theme(axis.title = element_markdown()) ) %>%
  animate(nframes = 450, duration = 10, # 10 s, so 45 fps
          width = 21, height = 10 * 2/3,
          units = "cm", res = 500, renderer = gifski_renderer()) %>%
  anim_save(filename = "speed.gif", path = here("Figures", "Animations"))

# 1.9 Matching static plot of difference ####
require(ggridges)
require(geomtextpath)
( speed_diff %>% 
  filter(parameter == "obs") %>%
  ggplot() +
  stat_density_ridges(aes(x = diff, y = 0, 
                          fill = if_else(after_stat(x) < 0,
                                         "Faeces", "Kelp")), 
                      geom = "density_ridges_gradient", n = 2^10,
                      colour = NA, linewidth = 0, bandwidth = 16 * 2 * 0.01,
                      from = -16, to = 16, rel_min_height = 0.001,
                      scale = 1) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Kelp),
                   colour = "#dabc23", family = "Futura", 
                   size = 3.5, hjust = 0.7, vjust = 0,
                   n = 2^10, bw = 16 * 2 * 0.01, text_only = TRUE) +
  geom_textdensity(aes(x = diff, y = after_stat(density),
                       label = label_Faeces),
                   colour = "#7030a5", family = "Futura", 
                   size = 3.5, hjust = 0.47, vjust = 0,
                   n = 2^10, bw = 16 * 2 * 0.01, text_only = TRUE) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-16, 16), oob = scales::oob_keep,
                     breaks = seq(-16, 16, 8),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_fill_manual(values = c("#7030a5", "#dabc23"),
                    guide = "none") +
  xlab("Δ sinking speed (cm s<sup><span style='font-size:8.4pt'>−1</span></sup>)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme +
  theme(axis.title = element_markdown()) ) %>%
  ggsave(filename = "speed.png", path = here("Figures", "Animations"),
         width = 21, height = 10 * 1/3, units = "cm", dpi = 500)


# 2. Distance ####
# 2.1 Load data ####
distance <- 
  read_rds(here("Sinking", "RDS", "distance.rds"))
distance_prior_posterior <- 
  read_rds(here("Sinking", "RDS", "distance_prior_posterior.rds"))
distance_diff <- 
  read_rds(here("Sinking", "RDS", "distance_diff.rds"))

# 2.2 Add animation aesthetics to data points ####
distance %<>%
  mutate(y = -0.5 + runif( n() , -0.36 , 0.36 ), # Manual jitter
         alpha = 0.4, # Baseline transparency
         colour = if_else( # Colours
           Tissue == "Kelp",
           "#dabc23", "#7030a5"
         )) %>%
  rename(x = Distance) %T>% # Rename for consistency
  print()

# 2.3 Build density polygons ####
distance_dens <- distance_prior_posterior %>%
  group_by(Tissue) %>%
  reframe(x = c(0, density(obs, n = 2^10, from = 0, to = 400, bw = 400 * 0.01)$x, 400),
          y = c(0, density(obs, n = 2^10, from = 0, to = 400, bw = 400 * 0.01)$y, 0)) %>%
  group_by(Tissue) %>% # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate(y = y * 45 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

# 2.4 Add animation aesthetics to polygons ####
distance_dens %<>%
  mutate(fill = case_when(
    Tissue == "Kelp" ~ "#dabc23", 
    Tissue == "Faeces" ~ "#7030a5",
    Tissue == "Prior" ~ "#b5b8ba"
  )) %T>%
  print()

# 2.5 Static plot ####
# Plot
ggplot() +
  geom_point(data = distance,
             aes(x = x, y = y, colour = colour, alpha = alpha),
             shape = 16, size = 2.5) +
  geom_polygon(data = distance_dens,
               aes(x = x, y = y, fill = fill)) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 400), ylim = c(-1, 2), 
                  expand = FALSE, clip = "off") +
  xlab("Export distance (km)") +
  mytheme
# Looks fine

# 2.6 Tween points ####
# Tween
distance_ani <- bind_rows(
  tween_state(distance %>% filter(Tissue == "Prior"), 
              distance %>% filter(Tissue == "Kelp"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50),
  tween_state(distance %>% filter(Tissue == "Kelp"),
              distance %>% filter(Tissue == "Faeces"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_state(distance %>% filter(Tissue == "Faeces"),
              distance %>% filter(Tissue == "Prior"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300) 
) %T>%
  print()


# 2.7 Tween polygons ####
distance_dens_ani <- bind_rows(
  tween_polygon(distance_dens %>% filter(Tissue == "Prior"),
                distance_dens %>% filter(Tissue == "Kelp"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(distance_dens %>% filter(Tissue == "Kelp"),
                distance_dens %>% filter(Tissue == "Faeces"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(distance_dens %>% filter(Tissue == "Faeces"),
                distance_dens %>% filter(Tissue == "Prior"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

# 2.8 Dynamic plot ####
( ggplot() +
    geom_point(data = distance_ani,
               aes(x = x, y = y, colour = colour, alpha = alpha),
               shape = 16, size = 2.5) +
    geom_polygon(data = distance_dens_ani,
                 aes(x = x, y = y, fill = fill)) +
    scale_colour_identity() +
    scale_alpha_identity() +
    scale_fill_identity() +
    new_scale_fill() +
    geom_polygon(data = distance_dens %>%
                   mutate(x = Inf),
                 aes(x = x, y = y, fill = Tissue %>% fct_relevel("Prior", "Kelp"))) +
    scale_fill_manual(values = c("#b5b8ba", "#dabc23", "#7030a5")) +
    coord_cartesian(xlim = c(0, 400), ylim = c(-1, 2), 
                    expand = FALSE, clip = "off") +
    xlab("Export distance (km)") +
    transition_manual(.frame) +
    mytheme ) %>%
  animate(nframes = 450, duration = 10, # 10 s, so 45 fps
          width = 21, height = 10 * 2/3,
          units = "cm", res = 500, renderer = gifski_renderer()) %>%
  anim_save(filename = "distance.gif", path = here("Figures", "Animations"))

# 2.9 Matching static plot of difference ####
( distance_diff %>% 
    filter(parameter == "obs") %>%
    ggplot() +
      stat_density_ridges(aes(x = diff, y = 0, 
                            fill = if_else(after_stat(x) < 0,
                                           "Faeces", "Kelp")), 
                        geom = "density_ridges_gradient", n = 2^10,
                        colour = NA, linewidth = 0, bandwidth = 400 * 2 * 0.01,
                        from = -400, to = 400, rel_min_height = 0.001,
                        scale = 1) +
      geom_textdensity(aes(x = diff, y = after_stat(density),
                         label = label_Kelp),
                     colour = "#dabc23", family = "Futura", 
                     size = 3.5, hjust = 0.6, vjust = 0,
                     n = 2^10, bw = 400 * 2 * 0.01, text_only = TRUE) +
      geom_textdensity(aes(x = diff, y = after_stat(density),
                         label = label_Faeces),
                     colour = "#7030a5", family = "Futura", 
                     size = 3.5, hjust = 0.4, vjust = 0,
                     n = 2^10, bw = 400 * 2 * 0.01, text_only = TRUE) +
      geom_vline(xintercept = 0) +
      scale_x_continuous(limits = c(-400, 400), oob = scales::oob_keep,
                       breaks = seq(-400, 400, 200),
                       labels = scales::label_number(style_negative = "minus")) +
      scale_fill_manual(values = c("#7030a5", "#dabc23"),
                      guide = "none") +
      xlab("Δ export distance (km)") +
      coord_cartesian(expand = FALSE, clip = "off") +
      mytheme ) %>%
  ggsave(filename = "distance.png", path = here("Figures", "Animations"),
         width = 21, height = 10 * 1/3, units = "cm", dpi = 500)


# 3. Carbon ####
# 3.1 Load data ####
carbon <- 
  read_rds(here("Biochemistry", "C_N", "RDS", "C_N.rds"))
carbon_prior_posterior <- 
  read_rds(here("Biochemistry", "C_N", "RDS", "C_prior_posterior.rds"))
carbon_diff <- 
  read_rds(here("Biochemistry", "C_N", "RDS", "C_diff.rds"))

# 3.2 Add animation aesthetics to data points ####
carbon %<>%
  mutate(y = -0.5 + runif( n() , -0.36 , 0.36 ), # Manual jitter
         alpha = 0.4, # Baseline transparency
         colour = if_else( # Colours
           Treatment == "Kelp",
           "#dabc23", "#7030a5"
         )) %>%
  rename(x = C) %>% 
  select(-c(N, C_N, Season)) %T>%
  print()

# 3.3 Build density polygons ####
carbon_dens <- carbon_prior_posterior %>%
  group_by(Treatment) %>%
  reframe(x = c(0, density(obs_new, n = 2^10, from = 0, to = 60, bw = 60 * 0.015)$x, 60),
          y = c(0, density(obs_new, n = 2^10, from = 0, to = 60, bw = 60 * 0.015)$y, 0)) %>%
  group_by(Treatment) %>%
  mutate(y = y * 30 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

# 3.4 Add animation aesthetics to polygons ####
carbon_dens %<>%
  mutate(fill = case_when(
    Treatment == "Kelp" ~ "#dabc23", 
    Treatment == "Faeces" ~ "#7030a5",
    Treatment == "Prior" ~ "#b5b8ba"
  )) %T>%
  print()

# 3.5 Static plot ####
# Plot
ggplot() +
  geom_point(data = carbon,
             aes(x = x, y = y, colour = colour, alpha = alpha),
             shape = 16, size = 2.5) +
  geom_polygon(data = carbon_dens,
               aes(x = x, y = y, fill = fill)) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 60), ylim = c(-1, 2), 
                  expand = FALSE, clip = "off") +
  xlab("Carbon content (%)") +
  mytheme
# Looks fine

# 3.6 Tween points ####
# Tween
carbon_ani <- bind_rows(
  tween_state(carbon %>% filter(Treatment == "Prior"), 
              carbon %>% filter(Treatment == "Kelp"),
              ease = "cubic-in-out", nframes = 100,
              id = Individual, # Here I pair actual sample pairs
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50),
  tween_state(carbon %>% filter(Treatment == "Kelp"),
              carbon %>% filter(Treatment == "Faeces"),
              ease = "cubic-in-out", nframes = 100,
              id = Individual,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_state(carbon %>% filter(Treatment == "Faeces"),
              carbon %>% filter(Treatment == "Prior"),
              ease = "cubic-in-out", nframes = 100,
              id = Individual,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300) 
) %T>%
  print()

# 3.7 Tween polygons ####
carbon_dens_ani <- bind_rows(
  tween_polygon(carbon_dens %>% filter(Treatment == "Prior"),
                carbon_dens %>% filter(Treatment == "Kelp"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(carbon_dens %>% filter(Treatment == "Kelp"),
                carbon_dens %>% filter(Treatment == "Faeces"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(carbon_dens %>% filter(Treatment == "Faeces"),
                carbon_dens %>% filter(Treatment == "Prior"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

# 3.8 Dynamic plot ####
( ggplot() +
    geom_point(data = carbon_ani,
             aes(x = x, y = y, colour = colour, alpha = alpha),
             shape = 16, size = 2.5) +
    geom_polygon(data = carbon_dens_ani,
               aes(x = x, y = y, fill = fill)) +
    scale_colour_identity() +
    scale_alpha_identity() +
    scale_fill_identity() +
    new_scale_fill() +
    geom_polygon(data = carbon_dens %>%
                 mutate(x = Inf),
               aes(x = x, y = y, fill = Treatment %>% fct_relevel("Prior", "Kelp"))) +
    scale_fill_manual(values = c("#b5b8ba", "#dabc23", "#7030a5")) +
    coord_cartesian(xlim = c(0, 60), ylim = c(-1, 2), 
                  expand = FALSE, clip = "off") +
    xlab("Carbon content (%)") +
    transition_manual(.frame) +
    mytheme ) %>%
  animate(nframes = 450, duration = 10, # 10 s, so 45 fps
          width = 21, height = 10 * 2/3,
          units = "cm", res = 500, renderer = gifski_renderer()) %>%
  anim_save(filename = "carbon.gif", path = here("Figures", "Animations"))

# 3.9 Matching static plot of difference ####
( carbon_diff %>% 
    filter(parameter == "obs_new") %>%
    ggplot() +
    stat_density_ridges(aes(x = diff, y = 0, 
                            fill = if_else(after_stat(x) < 0,
                                           "Faeces", "Kelp")), 
                        geom = "density_ridges_gradient", n = 2^10,
                        colour = NA, linewidth = 0, bandwidth = 50 * 2 * 0.01,
                        from = -50, to = 50, rel_min_height = 0.001,
                        scale = 1) +
    geom_textdensity(aes(x = diff, y = after_stat(density),
                         label = label_Kelp),
                     colour = "#dabc23", family = "Futura", 
                     size = 3.5, hjust = 0.72, vjust = 0,
                     n = 2^10, bw = 50 * 2 * 0.01, text_only = TRUE) +
    geom_textdensity(aes(x = diff, y = after_stat(density),
                         label = label_Faeces),
                     colour = "#7030a5", family = "Futura", 
                     size = 3.5, hjust = 0.42, vjust = 0,
                     n = 2^10, bw = 50 * 2 * 0.01, text_only = TRUE) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-50, 50), oob = scales::oob_keep,
                       breaks = seq(-50, 50, 25),
                       labels = scales::label_number(style_negative = "minus")) +
    scale_fill_manual(values = c("#7030a5", "#dabc23"),
                      guide = "none") +
    xlab("Δ carbon content (%)") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme ) %>%
  ggsave(filename = "carbon.png", path = here("Figures", "Animations"),
         width = 21, height = 10 * 1/3, units = "cm", dpi = 500)


# 4. Phenol ####
# 4.1 Load data ####
phenol <- 
  read_rds(here("Biochemistry", "Phenol", "RDS", "phenol.rds"))
phenol_prior_posterior <- 
  read_rds(here("Biochemistry", "Phenol", "RDS", "phenol_prior_posterior.rds"))
phenol_diff <- 
  read_rds(here("Biochemistry", "Phenol", "RDS", "phenol_diff.rds"))

# 4.2 Add animation aesthetics to data points ####
phenol %<>%
  select(Treatment, Individual, Samples_Data_Summary) %>%
  unnest(cols = Samples_Data_Summary) %>%
  mutate(y = -0.5 + runif( n() , -0.36 , 0.36 ), # Manual jitter
         alpha = 0.4, # Baseline transparency
         colour = if_else( # Colours
           Treatment == "Kelp",
           "#dabc23", "#7030a5"
         )) %>%
  rename(x = Concentration_mean) %>% 
  select(-Concentration_sd) %T>%
  print()

# 4.3 Build density polygons ####
phenol_dens <- phenol_prior_posterior %>%
  group_by(Treatment) %>% # No reasonable joint bandwidth.
  reframe(x = c(0, density(obs_new, n = 2^10, from = 0, to = 2)$x, 2),
          y = c(0, density(obs_new, n = 2^10, from = 0, to = 2)$y, 0)) %>%
  group_by(Treatment) %>%
  mutate(y = y * 0.6 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

# 4.4 Add animation aesthetics to polygons ####
phenol_dens %<>%
  mutate(fill = case_when(
    Treatment == "Kelp" ~ "#dabc23", 
    Treatment == "Faeces" ~ "#7030a5",
    Treatment == "Prior" ~ "#b5b8ba"
  )) %T>%
  print()

# 4.5 Static plot ####
# Plot
ggplot() +
  geom_point(data = phenol,
             aes(x = x, y = y, colour = colour, alpha = alpha),
             shape = 16, size = 2.5) +
  geom_polygon(data = phenol_dens,
               aes(x = x, y = y, fill = fill)) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_fill_identity() +
  scale_x_continuous(labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1))) +
  coord_cartesian(xlim = c(0, 2), ylim = c(-1, 2), 
                  expand = FALSE, clip = "off") +
  xlab("Phenolic content (%)") +
  mytheme
# Looks fine. I'll let the probability density for Faeces shoot through the roof
# to see more of the Prior and Kelp probability densities.

# 4.6 Tween points ####
phenol_ani <- bind_rows(
  tween_state(phenol %>% filter(Treatment == "Prior"), 
              phenol %>% filter(Treatment == "Kelp"),
              ease = "cubic-in-out", nframes = 100,
              id = Individual, # Here I also have meaningful pairs
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50),
  tween_state(phenol %>% filter(Treatment == "Kelp"),
              phenol %>% filter(Treatment == "Faeces"),
              ease = "cubic-in-out", nframes = 100,
              id = Individual,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_state(phenol %>% filter(Treatment == "Faeces"),
              phenol %>% filter(Treatment == "Prior"),
              ease = "cubic-in-out", nframes = 100,
              id = Individual,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300) 
) %T>%
  print()

# 4.7 Tween polygons ####
phenol_dens_ani <- bind_rows(
  tween_polygon(phenol_dens %>% filter(Treatment == "Prior"),
                phenol_dens %>% filter(Treatment == "Kelp"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(phenol_dens %>% filter(Treatment == "Kelp"),
                phenol_dens %>% filter(Treatment == "Faeces"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(phenol_dens %>% filter(Treatment == "Faeces"),
                phenol_dens %>% filter(Treatment == "Prior"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

# 4.8 Dynamic plot ####
( ggplot() +
    geom_point(data = phenol_ani,
               aes(x = x, y = y, colour = colour, alpha = alpha),
               shape = 16, size = 2.5) +
    geom_polygon(data = phenol_dens_ani,
                 aes(x = x, y = y, fill = fill)) +
    scale_colour_identity() +
    scale_alpha_identity() +
    scale_fill_identity() +
    new_scale_fill() +
    geom_polygon(data = phenol_dens %>%
                   mutate(x = Inf),
                 aes(x = x, y = y, fill = Treatment %>% fct_relevel("Prior", "Kelp"))) +
    scale_fill_manual(values = c("#b5b8ba", "#dabc23", "#7030a5")) +
    scale_x_continuous(labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1))) +
    coord_cartesian(xlim = c(0, 2), ylim = c(-1, 2), 
                    expand = FALSE, clip = "off") +
    xlab("Phenolic content (%)") +
    transition_manual(.frame) +
    mytheme ) %>%
  animate(nframes = 450, duration = 10, # 10 s, so 45 fps
          width = 21, height = 10 * 2/3,
          units = "cm", res = 500, renderer = gifski_renderer()) %>%
  anim_save(filename = "phenol.gif", path = here("Figures", "Animations"))

# 4.9 Matching static plot of difference ####
( phenol_diff %>% 
    filter(parameter == "obs_new") %>%
    ggplot() +
    stat_density_ridges(aes(x = diff, y = 0, 
                            fill = if_else(after_stat(x) < 0,
                                           "Faeces", "Kelp")), 
                        geom = "density_ridges_gradient", n = 2^10,
                        colour = NA, linewidth = 0, bandwidth = 2 * 2 * 0.01,
                        from = -2, to = 2, rel_min_height = 0.001,
                        scale = 1) +
    geom_textdensity(aes(x = diff, y = after_stat(density),
                         label = label_Kelp),
                     colour = "#dabc23", family = "Futura", 
                     size = 3.5, hjust = 0.8, vjust = 0,
                     n = 2^10, bw = 2 * 2 * 0.01, text_only = TRUE) +
    geom_textdensity(aes(x = diff, y = after_stat(density),
                         label = label_Faeces),
                     colour = "#7030a5", family = "Futura", 
                     size = 3.5, hjust = 0.44, vjust = 0,
                     n = 2^10, bw = 2 * 2 * 0.01, text_only = TRUE) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-2, 2), oob = scales::oob_keep,
                       breaks = seq(-2, 2, 1),
                       labels = scales::label_number(style_negative = "minus")) +
    scale_fill_manual(values = c("#7030a5", "#dabc23"),
                      guide = "none") +
    xlab("Δ phenolic content (%)") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme ) %>%
  ggsave(filename = "phenol.png", path = here("Figures", "Animations"),
         width = 21, height = 10 * 1/3, units = "cm", dpi = 500)


# 5. Chlorophyll ####
# 5.1 Load data ####
chlorophyll <- 
  read_rds(here("Biochemistry", "Pigments", "RDS", "pigments.rds"))
chlorophyll_prior_posterior <- 
  read_rds(here("Biochemistry", "Pigments", "RDS", "chlvspheo_prior_posterior.rds"))
chlorophyll_diff <- 
  read_rds(here("Biochemistry", "Pigments", "RDS", "chlvspheo_diff.rds"))

# 5.2 Add animation aesthetics to data points ####
chlorophyll %<>%
  select(Treatment, Individual, Samples_Data_Summary) %>%
  unnest(cols = Samples_Data_Summary) %>% # Calculate percentage of intact chlorophyll
  mutate(x = Chlorophyll_mean / (Chlorophyll_mean + Pheopigments_mean) * 100) %>%
  mutate(y = -0.5 + runif( n() , -0.36 , 0.36 ), # Manual jitter
         alpha = 0.4, # Baseline transparency
         colour = if_else( # Colours
           Treatment == "Kelp",
           "#dabc23", "#7030a5"
         )) %>%
  select(Treatment, Individual, x, y, alpha, colour) %T>%
  print()

# 5.3 Build density polygons ####
chlorophyll_dens <- chlorophyll_prior_posterior %>%
  group_by(Treatment) %>%
  reframe(x = c(0, density(obs_new, n = 2^10, from = 0, to = 100, bw = 100 * 0.025)$x, 100),
          y = c(0, density(obs_new, n = 2^10, from = 0, to = 100, bw = 100 * 0.025)$y, 0)) %>%
  group_by(Treatment) %>%
  mutate(y = y * 66 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

# 5.4 Add animation aesthetics to polygons ####
chlorophyll_dens %<>%
  mutate(fill = case_when(
    Treatment == "Kelp" ~ "#dabc23", 
    Treatment == "Faeces" ~ "#7030a5",
    Treatment == "Prior" ~ "#b5b8ba"
  )) %T>%
  print()

# 5.5 Static plot ####
# Plot
ggplot() +
  geom_point(data = chlorophyll,
             aes(x = x, y = y, colour = colour, alpha = alpha),
             shape = 16, size = 2.5) +
  geom_polygon(data = chlorophyll_dens,
               aes(x = x, y = y, fill = fill)) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 100), ylim = c(-1, 2), 
                  expand = FALSE, clip = "off") +
  xlab("Intact chlorophyll (%)") +
  mytheme
# Looks fine.

# 5.6 Tween points ####
chlorophyll_ani <- bind_rows(
  tween_state(chlorophyll %>% filter(Treatment == "Prior"), 
              chlorophyll %>% filter(Treatment == "Kelp"),
              ease = "cubic-in-out", nframes = 100,
              id = Individual,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50),
  tween_state(chlorophyll %>% filter(Treatment == "Kelp"),
              chlorophyll %>% filter(Treatment == "Faeces"),
              ease = "cubic-in-out", nframes = 100,
              id = Individual,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_state(chlorophyll %>% filter(Treatment == "Faeces"),
              chlorophyll %>% filter(Treatment == "Prior"),
              ease = "cubic-in-out", nframes = 100,
              id = Individual,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300) 
) %T>%
  print()

# 5.7 Tween polygons ####
chlorophyll_dens_ani <- bind_rows(
  tween_polygon(chlorophyll_dens %>% filter(Treatment == "Prior"),
                chlorophyll_dens %>% filter(Treatment == "Kelp"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(chlorophyll_dens %>% filter(Treatment == "Kelp"),
                chlorophyll_dens %>% filter(Treatment == "Faeces"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(chlorophyll_dens %>% filter(Treatment == "Faeces"),
                chlorophyll_dens %>% filter(Treatment == "Prior"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

# 5.8 Dynamic plot ####
( ggplot() +
    geom_point(data = chlorophyll_ani,
               aes(x = x, y = y, colour = colour, alpha = alpha),
               shape = 16, size = 2.5) +
    geom_polygon(data = chlorophyll_dens_ani,
                 aes(x = x, y = y, fill = fill)) +
    scale_colour_identity() +
    scale_alpha_identity() +
    scale_fill_identity() +
    new_scale_fill() +
    geom_polygon(data = chlorophyll_dens %>%
                   mutate(x = Inf),
                 aes(x = x, y = y, fill = Treatment %>% fct_relevel("Prior", "Kelp"))) +
    scale_fill_manual(values = c("#b5b8ba", "#dabc23", "#7030a5")) +
    coord_cartesian(xlim = c(0, 100), ylim = c(-1, 2), 
                    expand = FALSE, clip = "off") +
    xlab("Intact chlorophyll (%)") +
    transition_manual(.frame) +
    mytheme ) %>%
  animate(nframes = 450, duration = 10, # 10 s, so 45 fps
          width = 21, height = 10 * 2/3,
          units = "cm", res = 500, renderer = gifski_renderer()) %>%
  anim_save(filename = "chlorophyll.gif", path = here("Figures", "Animations"))

# 5.9 Matching static plot of difference ####
( chlorophyll_diff %>% 
    filter(parameter == "obs_new") %>%
    ggplot() +
    stat_density_ridges(aes(x = diff, y = 0, 
                            fill = if_else(after_stat(x) < 0,
                                           "Faeces", "Kelp")), 
                        geom = "density_ridges_gradient", n = 2^10,
                        colour = NA, linewidth = 0, bandwidth = 80 * 2 * 0.01,
                        from = -80, to = 80, rel_min_height = 0.001,
                        scale = 1) +
    geom_textdensity(aes(x = diff, y = after_stat(density),
                         label = label_Kelp),
                     colour = "#dabc23", family = "Futura", 
                     size = 3.5, hjust = 0.8, vjust = 0,
                     n = 2^10, bw = 80 * 2 * 0.01, text_only = TRUE) +
    geom_textdensity(aes(x = diff, y = after_stat(density),
                         label = label_Faeces),
                     colour = "#7030a5", family = "Futura", 
                     size = 3.5, hjust = 0.44, vjust = 0,
                     n = 2^10, bw = 80 * 2 * 0.01, text_only = TRUE) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-80, 80), oob = scales::oob_keep,
                       labels = scales::label_number(style_negative = "minus")) +
    scale_fill_manual(values = c("#7030a5", "#dabc23"),
                      guide = "none") +
    xlab("Δ intact chlorophyll (%)") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme ) %>%
  ggsave(filename = "chlorophyll.png", path = here("Figures", "Animations"),
         width = 21, height = 10 * 1/3, units = "cm", dpi = 500)


# 6. Grazing ####
# 6.1 Load data ####
grazing <- 
  read_rds(here("Urchins", "RDS", "grazing.rds"))
grazing_prior_posterior <- 
  read_rds(here("Urchins", "RDS", "grazing_prior_posterior.rds"))
grazing_prediction_summary <- 
  read_rds(here("Urchins", "RDS", "grazing_prediction_summary.rds"))
consumption_prior_posterior <- 
  read_rds(here("Urchins", "RDS", "consumption_prior_posterior.rds"))
defecation_prior_posterior <- 
  read_rds(here("Urchins", "RDS", "defecation_prior_posterior.rds"))

# 6.2 Prepare points and lines ####
grazing_points <- grazing %>%
  select(Season, Consumption, Defecation) %>%
  mutate(alpha = 0.4) %>%
  rename(x = Consumption, y = Defecation) %T>%
  print()

grazing_jitter <- grazing %>% 
  mutate(x = Defecation / Consumption * 100, # Calculate percent defecation
         y = -0.5 + runif( n() , -0.36 , 0.36 ), # Manual jitter
         alpha = 0.4) %>% 
  select(Season, x, y, alpha) %T>%
  print()

grazing_line <- grazing_prediction_summary %>%
  filter(!Season %in% c("Prior", "Annual")) %>%
  droplevels() %>%
  select(Season, Consumption, mu) %>%
  rename(x = Consumption, y = mu) %T>%
  print()

# 6.3 Build ribbon and density polygons ####
grazing_ribbon <- grazing_prediction_summary %>%
  filter(!Season %in% c("Prior", "Annual")) %>%
  droplevels() %>%
  group_by(Season, .width) %>%
  reframe(x = c(Consumption, rev(Consumption)),
          y = c(mu.upper, rev(mu.lower))) %>%
  mutate(alpha = case_when(
    .width == 0.9 ~ 0.5, 
    .width == 0.8 ~ 0.4, 
    .width == 0.5 ~ 0.3
  )) %>%
  select(-.width) %T>%
  print()

grazing_dens <- grazing_prior_posterior %>%
  filter(!Season %in% c("Prior", "Annual")) %>%
  droplevels() %>%
  group_by(Season) %>%
  reframe(x = c(0, density(beta, n = 2^10, from = 0, to = 100, bw = 100 * 0.01)$x, 100),
          y = c(0, density(beta, n = 2^10, from = 0, to = 100, bw = 100 * 0.01)$y, 0)) %>%
  group_by(Season) %>%
  mutate(y = y * 16 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

consumption_dens <- consumption_prior_posterior %>%
  filter(!Season %in% c("Prior", "Annual")) %>%
  droplevels() %>%
  group_by(Season) %>%
  reframe(x = c(0, density(mu, n = 2^10, from = 0, to = 10, bw = 10 * 0.01)$x, 10),
          y = c(0, density(mu, n = 2^10, from = 0, to = 10, bw = 10 * 0.01)$y, 0)) %>%
  group_by(Season) %>%
  mutate(y = y * 0.75 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

defecation_dens <- defecation_prior_posterior %>%
  filter(!Season %in% c("Prior", "Annual")) %>%
  droplevels() %>%
  group_by(Season) %>% # Reversed x and y positions because this density is on the y axis
  reframe(x = c(0, density(mu, n = 2^10, from = 0, to = 5, bw = 5 * 0.01)$y, 0),
          y = c(0, density(mu, n = 2^10, from = 0, to = 5, bw = 5 * 0.01)$x, 5)) %>%
  group_by(Season) %>%
  mutate(x = x * 0.5 / ( sum(x) * ( y[3] - y[2] ) )) %>%
  ungroup() %T>%
  print()

# 6.4 Static plots ####
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

ggplot() +
  geom_textline(data = tibble(x = c(0, 5), y = c(0, 5)), aes(x, y),
                label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
  geom_point(data = grazing_points,
             aes(x = x, y = y, alpha = alpha),
             shape = 16, size = 2.5, colour = "#7030a5") +
  geom_line(data = grazing_line,
            aes(x = x, y = y, group = Season),
            colour = "#7030a5") +
  geom_polygon(data = grazing_ribbon,
               aes(x = x, y = y, alpha = alpha, 
                   # Interaction grouping is only needed in the static version
                   group = interaction(alpha, Season)),
               fill = "#7030a5") +
  geom_polygon(data = consumption_dens,
               aes(x = x, y = y),
               fill = "#7030a5") +
  geom_polygon(data = defecation_dens,
               aes(x = x, y = y),
               fill = "#7030a5") +
  scale_alpha_identity() +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 5),
                  expand = FALSE, clip = "off") +
  labs(x = expression("Consumption (mg g"^-1*" d"^-1*")"),
       y = expression("Defecation (mg g"^-1*" d"^-1*")")) +
  mytheme

ggplot() +
  geom_point(data = grazing_jitter,
             aes(x = x, y = y, alpha = alpha),
             shape = 16, size = 2.5, colour = "#7030a5") +
  geom_polygon(data = grazing_dens,
               aes(x = x, y = y), fill = "#7030a5") +
  scale_alpha_identity() +
  coord_cartesian(xlim = c(0, 100), ylim = c(-1, 2), 
                  expand = FALSE, clip = "off") +
  xlab("Defecation (%)") +
  mytheme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())
# All looks to be in order.

# 6.5 Tween points ####
grazing_points_ani <- bind_rows( # Here it's arbitrary pairings again
  tween_state(grazing_points %>% filter(Season == "Spring"), 
              grazing_points %>% filter(Season == "Summer"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50),
  tween_state(grazing_points %>% filter(Season == "Summer"), 
              grazing_points %>% filter(Season == "Autumn"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_state(grazing_points %>% filter(Season == "Autumn"), 
              grazing_points %>% filter(Season == "Spring"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300) 
) %T>%
  print()

grazing_jitter_ani <- bind_rows(
  tween_state(grazing_jitter %>% filter(Season == "Spring"), 
              grazing_jitter %>% filter(Season == "Summer"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50),
  tween_state(grazing_jitter %>% filter(Season == "Summer"), 
              grazing_jitter %>% filter(Season == "Autumn"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_state(grazing_jitter %>% filter(Season == "Autumn"), 
              grazing_jitter %>% filter(Season == "Spring"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300) 
) %T>%
  print()

# 6.6 Tween lines and polygons ####
grazing_line_ani <- bind_rows(
  tween_path(grazing_line %>% filter(Season == "Spring"), 
             grazing_line %>% filter(Season == "Summer"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_path(grazing_line %>% filter(Season == "Summer"), 
             grazing_line %>% filter(Season == "Autumn"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_path(grazing_line %>% filter(Season == "Autumn"), 
             grazing_line %>% filter(Season == "Spring"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300) 
) %T>%
  print()
# tween_state also works for lines of equal length, but 
# I believe it is more stable to use tween_path.

grazing_dens_ani <- bind_rows(
  tween_polygon(grazing_dens %>% filter(Season == "Spring"),
                grazing_dens %>% filter(Season == "Summer"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(grazing_dens %>% filter(Season == "Summer"),
                grazing_dens %>% filter(Season == "Autumn"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(grazing_dens %>% filter(Season == "Autumn"),
                grazing_dens %>% filter(Season == "Spring"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

grazing_ribbon_ani <- bind_rows(
  tween_polygon(grazing_ribbon %>% filter(Season == "Spring"),
                grazing_ribbon %>% filter(Season == "Summer"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(grazing_ribbon %>% filter(Season == "Summer"),
                grazing_ribbon %>% filter(Season == "Autumn"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(grazing_ribbon %>% filter(Season == "Autumn"),
                grazing_ribbon %>% filter(Season == "Spring"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

consumption_dens_ani <- bind_rows(
  tween_polygon(consumption_dens %>% filter(Season == "Spring"),
             consumption_dens %>% filter(Season == "Summer"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(consumption_dens %>% filter(Season == "Summer"),
             consumption_dens %>% filter(Season == "Autumn"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(consumption_dens %>% filter(Season == "Autumn"),
             consumption_dens %>% filter(Season == "Spring"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

defecation_dens_ani <- bind_rows(
  tween_polygon(defecation_dens %>% filter(Season == "Spring"),
                defecation_dens %>% filter(Season == "Summer"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(defecation_dens %>% filter(Season == "Summer"),
                defecation_dens %>% filter(Season == "Autumn"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(defecation_dens %>% filter(Season == "Autumn"),
                defecation_dens %>% filter(Season == "Spring"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

# 6.7 Dynamic plots ####
( ggplot() +
    geom_textline(data = tibble(x = c(0, 5), y = c(0, 5)), aes(x, y),
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_point(data = grazing_points_ani,
               aes(x = x, y = y, alpha = alpha),
               shape = 16, size = 2.5, colour = "#7030a5") +
    geom_line(data = grazing_line_ani,
              aes(x = x, y = y, group = Season),
              colour = "#7030a5") +
    geom_polygon(data = grazing_ribbon_ani,
                 aes(x = x, y = y, alpha = alpha,
                     group = alpha),
                 fill = "#7030a5") +
    geom_polygon(data = consumption_dens_ani,
                 aes(x = x, y = y),
                 fill = "#7030a5") +
    geom_polygon(data = defecation_dens_ani,
                 aes(x = x, y = y),
                 fill = "#7030a5") +
    geom_text(data = grazing_line_ani %>%
                distinct(Season, .frame),
              aes(x = 0.3, y = 4.8, label = Season),
              hjust = 0, size = 5.3, family = "Futura",
              colour = "#7030a5", fontface = "bold") +
    scale_alpha_identity() +
    scale_x_continuous(breaks = seq(0, 10, 2)) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 5),
                    expand = FALSE, clip = "off") +
    labs(x = "Consumption (mg g<sup><span style='font-size:8.4pt'>−1</span></sup> 
         d<sup><span style='font-size:8.4pt'>−1</span></sup>)",
         y = "Defecation (mg g<sup><span style='font-size:8.4pt'>−1</span></sup> 
         d<sup><span style='font-size:8.4pt'>−1</span></sup>)") +
    transition_manual(.frame) +
    mytheme +
    theme(axis.title = element_markdown()) ) %>%
  animate(nframes = 450, duration = 10, # 10 s, so 45 fps
          width = 21 * 1/2, height = 10,
          units = "cm", res = 500, renderer = gifski_renderer()) %>%
  anim_save(filename = "grazing_left.gif", path = here("Figures", "Animations"))

( ggplot() +
    geom_point(data = grazing_jitter_ani,
               aes(x = x, y = y, alpha = alpha),
               shape = 16, size = 2.5, colour = "#7030a5") +
    geom_polygon(data = grazing_dens_ani,
                 aes(x = x, y = y), fill = "#7030a5") +
    scale_alpha_identity() +
    coord_cartesian(xlim = c(0, 100), ylim = c(-1, 2), 
                    expand = FALSE, clip = "off") +
    xlab("Defecation (%)<sup><span style='color:white;font-size:8.4pt'>−1</span></sup>") +
    transition_manual(.frame) +
    mytheme +
    theme(axis.title = element_markdown(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) ) %>%
  animate(nframes = 450, duration = 10, # 10 s, so 45 fps
          width = 21 * 1/2, height = 10,
          units = "cm", res = 500, renderer = gifski_renderer()) %>%
  anim_save(filename = "grazing_right.gif", path = here("Figures", "Animations"))

# # Adjust theme for higher resolution plotting.
# theme(axis.title = element_markdown(size = 12 * 4),
#       axis.text.x = element_text(size = 10 * 4, margin = margin(t = 10)),
#       axis.line = element_line(linewidth = 0.5 * 4),
#       axis.ticks = element_line(linewidth = 0.5 * 4),
#       axis.ticks.length = unit(0.25 * 4, "cm"),
#       legend.text = element_text(size = 12 * 4, margin = margin(l = 20)),
#       legend.key.width = unit(0.25 * 4, "cm"),
#       legend.key.spacing.x = unit(0.5 * 4, "cm"),
#       plot.margin = margin(c(0.2, 0.5, 0.2, 0.2) * 4, unit = "cm")) ) %>%


# # 4.2 Build ID density polygons 
# phenol_ID_dens <- phenol %>%
#   select(Treatment, Individual, ID, Samples_Data) %>%
#   unnest(cols = Samples_Data) %>%
#   group_by(Treatment, Individual, ID) %>% 
#   # Starts from negative because some densities are not strictly positive.
#   reframe(x = density(Concentration, n = 2^10, from = -0.1, to = 2)$x,
#           y = density(Concentration, n = 2^10, from = -0.1, to = 2)$y) %>%
#   group_by(ID) %>%
#   mutate(y = y * 0.006 / ( sum(y) * ( x[2] - x[1] ) )) %>%
#   ungroup() %>%
#   #filter(y > 0.001) %>%
#   group_by(Treatment, Individual, ID) %>%
#   reframe(x = c(x, x %>% rev()),
#           y = c(y, -y %>% rev())) %>%
#   ungroup()
# 
# # 4.3 Add animation aesthetics to ID densities
# phenol_ID_dens %<>%
#   group_by(ID) %>% # Grouping is important so there is no jitter within group
#   mutate(y = y - 0.5 + runif( 1 , -0.36 , 0.36 ), # Also specify n = 1
#          alpha = 0.4, # Baseline transparency
#          fill = if_else( # Fills
#            Treatment == "Kelp",
#            "#dabc23", "#7030a5"
#          )) %>%
#   ungroup()
# phenol_ID_dens