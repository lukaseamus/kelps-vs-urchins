# 1. Prepare data ####
# 1.1 Load raw data ####
require(tidyverse)
require(magrittr)
require(here)

enzymes <- here("Microbes", "Enzymes", "Raw") %>%
  list.files(pattern = "\\.csv$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Name = Path %>% basename() %>% 
      str_remove("\\..*$"),
    Data = Path %>% 
      map(~ read_csv(.x, skip = 30, col_select = -1) %>%
            rename(Temperature = starts_with("T°"))),
    Date = Name %>% str_remove("^X") %>% str_split_i("_", 1) %>% ymd(),
    Enzyme = Name %>% str_split_i("_", 2),
    Plate = Name %>% str_split_i("_", 3) %>% 
      as.numeric() %>%
      replace_na(1)
  ) %>%
  select(-Path) %T>%
  print(n = 25)

enzymes$Data[[1]] %>%
  print(n = 451)

enzymes$Data[[1]] %>%
  ggplot(aes(Time, Temperature)) +
    geom_point(alpha = 0.3, shape = 16) +
    theme_minimal()

enzymes$Data[[1]] %>%
  ggplot(aes(Time, A5)) +
    geom_point(alpha = 0.3, shape = 16) +
    theme_minimal()

enzymes$Data[[1]] %>%
  ggplot(aes(Time, G8)) +
    geom_point(alpha = 0.3, shape = 16) +
    theme_minimal()




############

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


