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

# 1.2 Add animation aesthetics to points ####
speed %<>%
  mutate(y = -0.5 + runif( n() , -0.36 , 0.36 ), # Manual jitter
         alpha = 0.5, # Baseline transparency
         colour = if_else( # Colours
           Tissue == "Kelp",
           "#dabc23", "#7030a5"
         )) %>%
  rename(x = Speed) %>% # Rename for consistency
  select(-c(Mass, Area)) # Remove irrelevant variables
speed

# 1.3 Build desnity polygons ####
speed_dens <- speed_prior_posterior %>%
  group_by(Tissue) %>% # Close polygons at zero with c(0,).
  reframe(x = c(0, density(y_new, n = 2^10, from = 0, to = 0.25)$x),
          y = c(0, density(y_new, n = 2^10, from = 0, to = 0.25)$y)) %>%
  group_by(Tissue) %>% # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate(y = y * 0.02 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup()

# 1.4 Add animation aesthetics to polygons ####
speed_dens %<>%
  mutate(fill = case_when(
    Tissue == "Kelp" ~ "#dabc23", 
    Tissue == "Faeces" ~ "#7030a5",
    TRUE ~ "#b5b8ba"
  ))
speed_dens

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
             shape = 16, size = 3) +
  geom_polygon(data = speed_dens,
               aes(x = x, y = y, fill = fill)) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 0.25),
                     labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01, 0.1, 0.01))) +
  coord_cartesian(ylim = c(-1, 2), expand = FALSE, clip = "off") +
  xlab(expression("Sinking speed (m s"^-1*")")) +
  mytheme
# Looks fine

# 1.6 Tween points ####
require(tweenr)
# ani_data <- tween_state(sinking %>% filter(Tissue == "Kelp"), 
#                         sinking %>% filter(Tissue == "Faeces"),
#                         ease = "cubic-in-out", nframes = 100,
#                         enter = function(df){
#                                   df$alpha <- 0 
#                                   df
#                                   },
#                         exit = function(df){
#                                 df$alpha <- 0 
#                                 df
#                                 }) %>% 
#   #keep_state(20) %>% 
#   tween_state(sinking %>% filter(Tissue == "Kelp"), 
#               ease = "cubic-in-out", nframes = 100,
#               enter = function(df){
#                 df$alpha <- 0 
#                 df
#               },
#               exit = function(df){
#                 df$alpha <- 0 
#                 df
#               })

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
speed_ani <- bind_rows( # Filter for non-existent prior data
  tween_state(speed %>% filter(Tissue == "Prior"), 
              speed %>% filter(Tissue == "Kelp"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit),
  tween_state(speed %>% filter(Tissue == "Kelp"),
              speed %>% filter(Tissue == "Faeces"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    mutate(.frame = .frame + 100),
  tween_state(speed %>% filter(Tissue == "Faeces"),
              speed %>% filter(Tissue == "Prior"),
              ease = "cubic-in-out", nframes = 100,
              enter = enter, exit = exit) %>%
    mutate(.frame = .frame + 200)
)
  
speed_ani

# 1.7 Tween polygons ####
require(transformr)
# For some reason tween_polygon doesn't work for densities,
# but tween_path does.
speed_dens_ani <- bind_rows(
  tween_path(speed_dens %>% filter(Tissue == "Prior"),
             speed_dens %>% filter(Tissue == "Kelp"),
             ease = "cubic-in-out", nframes = 100),
  tween_path(speed_dens %>% filter(Tissue == "Kelp"),
             speed_dens %>% filter(Tissue == "Faeces"),
             ease = "cubic-in-out", nframes = 100) %>%
    mutate(.frame = .frame + 100),
  tween_path(speed_dens %>% filter(Tissue == "Faeces"),
             speed_dens %>% filter(Tissue == "Prior"),
             ease = "cubic-in-out", nframes = 100) %>%
    mutate(.frame = .frame + 200)
)

speed_dens_ani


# 1.8 Dynamic plot ####
require(gganimate)
require(ggtext)
speed_plot <- ggplot() +
  geom_point(data = speed_ani,
             aes(x = x, y = y, colour = colour, alpha = alpha),
             shape = 16, size = 3) +
  geom_polygon(data = speed_dens_ani,
               aes(x = x, y = y, fill = fill)) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 0.25),
                     labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01, 0.1, 0.01))) +
  coord_cartesian(ylim = c(-1, 2), expand = FALSE, clip = "off") +
  xlab("Sinking speed (m s<sup><span style='font-size:8.4pt'>−1</span></sup>)") +
  transition_manual(.frame) +
  mytheme +
  theme(axis.title = element_markdown())
speed_plot
# Looks good

# 1.9 Save as .gif ####
speed_plot %>%
  animate(nframes = 100, duration = 6,
          width = 1000, height = 500,
          renderer = gifski_renderer()) %>%
  anim_save(filename = "speed.gif", path = here("Figures", "Animations"))





ggplot(data = sinking_dens,aes(group = Group)) +
  geom_polygon(data = sinking_dens,
               aes(x = x, y = y, fill = Tissue))



sinking %>%
  ggplot(aes(Speed)) +
  geom_boxplot() +
  transition_states(Tissue) +
  # enter_fade() + 
  # exit_shrink() +
  ease_aes() +
  theme_minimal()

sinking %>%
  group_by(Tissue) %>%
  summarise(length(Speed))




# sinking_dens <- sinking %>%
#   group_by(Tissue) %>% 
#   reframe(x = c(0, density(Speed, n = 2^10, from = 0, to = 0.3)$x),
#           y = c(0, density(Speed, n = 2^10, from = 0, to = 0.3)$y)) %>%
#   group_by(Tissue) %>%
#   mutate(.id = row_number()) %>%
#   ungroup()
# ani_poly <- tween_polygon(sinking_dens %>% filter(Tissue == "Kelp"),
#                           sinking_dens %>% filter(Tissue == "Faeces"),
#                           ease = "cubic-in-out", id = .id, nframes = 100) %>%
#             keep_state(20)

# Make arbitrary pairings
# sinking %<>%
#   group_by(Tissue) %>%
#   mutate(Group = row_number() %>% str_c() %>% fct()) %>%
#   ungroup()