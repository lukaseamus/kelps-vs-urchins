# 1. Load data ####
require(tidyverse)
require(magrittr)
require(here)

# The first experiment in autumn ran from 
# 14th September 2023
start <- dmy_hm("14.09.23 21:00") # t0
# to 26th September 2023
end <- dmy_hm("26.09.23 15:00") # t5
# Temperature loggers were retrieved at intervals
t2 <- dmy_hm("20.09.23 15:00")
t4 <- dmy_hm("24.09.23 14:00")

light <- here("Urchins", "Light.csv") %>%
  read_csv() %>%
  mutate(
    Date = Date %>% dmy(),
    Day = date(start) %--% Date / ddays(),
    Treatment = case_when(
      Tank %>% str_detect("U") ~ "Grazed",
      Tank %>% str_detect("C") ~ "Control",
      Tank %>% str_detect("M") ~ "Mechanical"
    ),
    Season = if_else(Initials == "LSW", "Autumn", "Spring")
  ) %T>%
  print()

temp <- here("Urchins", "Temperature") %>%
  list.files(pattern = "\\.csv$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    ID = Path %>% basename() %>% str_remove("\\..*$"),
    Tank = case_when(
      ID == 10640480 ~ "U5",
      ID == 20884709 ~ "C5",
      ID == 20884712 ~ "M5",
      ID == 20884701 ~ "U10",
      ID == 20884710 ~ "C10",
      ID == 20884714 ~ "M10",
      ID == 20884707 ~ "U15",
      ID == 20884711 ~ "C15",
      ID == 21508210 ~ "M15"
    ),
    Treatment = case_when(
      Tank %>% str_detect("U") ~ "Grazed",
      Tank %>% str_detect("C") ~ "Control",
      Tank %>% str_detect("M") ~ "Mechanical"
    ),
    Data = Path %>% map(
        ~ .x %>% 
          read_csv(skip = 1, col_select = 2:3) %>%
          rename(Datetime = 1, Temperature = 2) %>%
          mutate(
            Datetime = Datetime %>% mdy_hms(),
            Day = start %--% Datetime / ddays()
          ) %>%
          drop_na(Temperature)
      )
  ) %>%
  select(-Path) %T>%
  print()

# 2. Explore data ####
# 2.1 Light ####
light %>%
  ggplot(aes(0, PAR)) +
    geom_jitter() +
    facet_grid(~ Season + Treatment)
# PAR is consistent across treatments and seasons. 
# A different PAR meter or setting was used in spring 
# because accuracy is much lower. 

light %>%
  filter(Season == "Autumn") %>%
  ggplot(aes(Day, PAR)) +
    geom_point() +
    facet_grid(~Treatment)
# Slight tendency of increase in the grazed treatment
# but reasonably consistent over time.

# 2.2 Temperature ####
require(glue)
temp %<>%
  rowwise() %>%
  mutate(
    Plot = list(
      Data %>%
        ggplot(aes(Datetime, Temperature)) +
          geom_point() +
          geom_vline(xintercept = c(start, t2, t4, end)) +
          labs(title = glue("{Treatment} ({Tank})"))
    )
  ) %>%
  ungroup() %T>%
  print()

require(patchwork)
temp %$% wrap_plots(Plot)
# All data before start and after end should be filtered out
# because the logger could be out of water as can be seen by
# the temperature spikes. Additionally, for loggers retrieved
# earlier, data after t2 and t4 need to be filtered out.

# Filter
temp %<>%
  mutate(
    Number = Tank %>% str_extract("\\d+") %>% as.numeric(),
    Data = Data %>%
      map2(
        Number,
        ~ if(.y == 5) {
          .x %>% filter(Datetime >= start & Datetime < t2)
        } else if(.y == 10) {
          .x %>% filter(Datetime >= start & Datetime < t4)
        } else {
          .x %>% filter(Datetime >= start & Datetime < end)
        }
      )
  ) %T>%
  print()

# Re-plot
temp %<>%
  rowwise() %>%
  mutate(
    Plot = list(
      Data %>%
        ggplot(aes(Datetime, Temperature)) +
          geom_point() +
          geom_vline(xintercept = c(start, t2, t4, end)) +
          labs(title = glue("{Treatment} ({Tank})"))
    )
  ) %>%
  ungroup() %T>%
  print()

temp %$% wrap_plots(Plot)
# Good

# Unnested
temp %>%
  unnest(Data) %>%
  ggplot(aes(Day, Temperature)) +
    geom_point(shape = 16, alpha = 0.2) +
    facet_grid(~ Treatment)
# Except for that spike in one tank of the mechanical treatment,
# the seasonal decline in ambient water temperature is consistent.
# The drop in temperature is especially prominent around day 6.

# 3. Calculate summaries ####
# 3.1 Light ####
light %>%
  summarise(
    mean = mean(PAR),
    sd = sd(PAR),
    n = n(),
    .by = c(Season, Treatment)
  ) %>%
  mutate(
    across(
      c(mean, sd),
      ~ .x %>% signif(2)
    ),
    PAR = glue("{mean} ± {sd}")
  ) %>%
  select(Season, Treatment, PAR, n)
# Season Treatment  PAR            n
# Autumn Grazed     9.4 ± 1      120
# Autumn Control    9.8 ± 0.93   120
# Autumn Mechanical 9.4 ± 1      120
# Spring Control    9.3 ± 1.9     60
# Spring Grazed     9.5 ± 1.5     60

light %>%
  summarise(
    mean = mean(PAR),
    sd = sd(PAR),
    n = n(),
    .by = Treatment
  ) %>%
  mutate(
    across(
      c(mean, sd),
      ~ .x %>% signif(2)
    ),
    PAR = glue("{mean} ± {sd}")
  ) %>%
  select(Treatment, PAR, n)
# Treatment  PAR           n
# Grazed     9.5 ± 1.2   180
# Control    9.7 ± 1.3   180
# Mechanical 9.4 ± 1     120

light %>%
  summarise(
    mean = mean(PAR),
    sd = sd(PAR),
    n = n()
  ) %>%
  mutate(
    across(
      c(mean, sd),
      ~ .x %>% signif(2)
    ),
    PAR = glue("{mean} ± {sd}")
  ) %>%
  select(PAR, n)
# PAR           n
# 9.5 ± 1.2   480

# 3.2 Temperature ####
temp %>%
  unnest(Data) %>%
  summarise(
    mean = mean(Temperature),
    sd = sd(Temperature),
    n = n(),
    .by = Treatment
  ) %>%
  mutate(
    mean = mean %>% round(1),
    sd = sd %>% signif(2),
    Temperature = glue("{mean} ± {sd}")
  ) %>%
  select(Treatment, Temperature, n)
# Treatment  Temperature     n
# Grazed     15 ± 0.73    1306
# Control    14.9 ± 0.72    1306
# Mechanical 15 ± 0.78    1306

temp %>%
  unnest(Data) %>%
  summarise(
    mean = mean(Temperature),
    sd = sd(Temperature),
    n = n()
  ) %>%
  mutate(
    mean = mean %>% round(1),
    sd = sd %>% signif(2),
    Temperature = glue("{mean} ± {sd}")
  ) %>%
  select(Temperature, n)
# Temperature     n
# 15 ± 0.74    3918

temp %>%
  unnest(Data) %>%
  mutate(Phase = if_else(Day < 6, "Early", "Late")) %>%
  summarise(
    mean = mean(Temperature),
    sd = sd(Temperature),
    n = n(),
    .by = Phase
  ) %>%
  mutate(
    mean = mean %>% round(1),
    sd = sd %>% signif(2),
    Temperature = glue("{mean} ± {sd}")
  ) %>%
  select(Phase, Temperature, n)
# Phase Temperature     n
# Early 15.5 ± 0.31  2556
# Late  14 ± 0.2     1362