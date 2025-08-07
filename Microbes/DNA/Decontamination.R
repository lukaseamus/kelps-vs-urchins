# 1. Prepare data ####
# 1.1 Load ASVs ####
require(tidyverse)
require(here)
ASV <- here("Microbes", "DNA", "RDS", "ASV.rds") %>%
  read_rds()

# 1.2 Remove contamination ####
# 1.2.1 Lab- and cross-contamination ####
# decontam's frequency method requires DNA concentration, 
# so I'll read in the Nanodrop and sequencing yield results.
require(magrittr)
conc <- here("Microbes", "DNA", "Concentration.csv") %>%
  read_csv() %>%
  rename(Sample = `Sample Name`, Concentration = `Nucleic Acid(ng/uL)`) %>%
  select(Sample, Concentration) %T>%
  print(n = 21)
# Concentration is negative for sample 21. I'll add a constant
# to make all concentrations positive.

conc %<>% 
  mutate(Concentration = Concentration + 0.03) %T>%
  print(n = 21)

# I also have PCR yield estimates.
yield <- here("Microbes", "DNA", "Yield.csv") %>%
  read_csv() %T>%
  print(n = 21)

meta <- yield %>%
  bind_cols(conc %>%
              select(Concentration)) %>%
  mutate(Sample = Sample %>% 
           str_split_i(pattern = "-", i = 1)) %T>%
  print(n = 21)

# Add meta to sample data
sample_data(ASV) %<>%
  as_tibble() %>%
  # rownames_to_column("File") somehow doesn't work
  mutate(File = sample_data(ASV) %>%
           as.data.frame() %>%
           rownames()) %>%
  left_join(meta, by = "Sample") %>%
  column_to_rownames("File") %>%
  sample_data() %T>%
  print()

sample_data(ASV) %>%
  as_tibble() %>%
  ggplot(aes(Concentration, Yield, 
             colour = Sample == "blank")) +
    geom_point() +
    theme_minimal()
# Blank has lowest DNA concentration and second-lowest
# PCR yield.

sample_data(ASV) %>%
  as_tibble() %>%
  mutate(Abundance_sum = sample_sums(ASV)) %>%
  arrange(Abundance_sum) %>%
  mutate(Index = row_number()) %>%
  ggplot(aes(Index, Abundance_sum, 
             colour = Sample == "blank")) +
    geom_point() +
    theme_minimal()
# All samples have more reads than the blank except for one.

# decontam frequency method
require(decontam)
cont_freq <- isContaminant(
  ASV, method = "frequency", 
  conc = "Concentration"
)

cont_freq %>%
  count(contaminant)
# Only identified 22 contaminants

cont_freq %>%
  filter(contaminant == TRUE)

(
  plot_frequency(
    ASV, cont_freq %>%
      filter(contaminant == TRUE) %>% 
      rownames(),
    conc = "Concentration"
  ) +
    xlab(expression("DNA concentration (ng µl"^-1*")"))
) %>%
  ggsave(filename = "cont_freq.pdf", 
         path = here("Microbes", "DNA", "Plots"),
         device = cairo_pdf, width = 40, height = 20,
         units = "cm")
# This is not convincing at all. There are very few points
# to base the regression on, both models seem to fit equally
# well aand in all cases the putative contaminant is not 
# present in the blank.

# decontam prevalence method
sample_data(ASV)$blank <- sample_data(ASV)$Sample == "blank"
sample_data(ASV)

cont_prev <- isContaminant(
  ASV, method = "prevalence",
  neg = "blank"
)

cont_prev %>%
  count(contaminant)
# Identified 58 contaminants

cont_prev %>%
  filter(contaminant == TRUE)

blank <- ASV %>%
  # transform_sample_counts(fun = function(x) 1 * (x > 0)) %>%
  transform_sample_counts(fun = function(x) x / sum(x)) %>%
  prune_samples(samples = sample_data(ASV)$blank == TRUE) %T>%
  print()

samples <- ASV %>%
  # transform_sample_counts(fun = function(x) 1 * (x > 0)) %>%
  transform_sample_counts(fun = function(x) x / sum(x)) %>%
  prune_samples(samples = sample_data(ASV)$blank == FALSE) %T>%
  print()

tibble(
  Blank = taxa_sums(blank),
  Samples = taxa_sums(samples),
  Contaminant = cont_prev$contaminant
) %>%
  ggplot(aes(Blank %>% log10(), Samples %>% log10(), 
             colour = Contaminant)) +
    geom_abline(slope = 1) +
    geom_point(shape = 16, alpha = 0.3, size = 2) +
    labs(y = "Relative abundance across samples",
         x = "Relative abundance in blank") +
    theme_minimal()
# Also seems a bit random. ASVs that are only present in
# the blank are not flagged as contaminants and those
# that are present in both blank and samples are sometimes
# randomly interspersed with those that weren't flagged.
# This is based on p-values from Fisher's exact test so...
# Also, here the comparsion is between sums but on the
# one hand I have a single blank and on the other 20 samples.

# microDecon method
# Generally disappointed with decontam. I'll try microDecon.
# Format data fro microDecon
ASV_df <- otu_table(ASV) %>%
  as.data.frame() %>%
  rownames_to_column("File") %>%
  as_tibble() %>%
  mutate(Sample = File %>% 
           str_split_i(pattern = "_", i = 1)) %>%
  select(-File) %>%
  pivot_longer(cols = -Sample,
               names_to = "ASV",
               values_to = "Abundance") %>%
  pivot_wider(names_from = Sample,
              values_from = Abundance) %>%
  # full_join(
  #   tax_table(ASV) %>% 
  #     as.data.frame() %>%
  #     rownames_to_column("ASV") %>%
  #     as_tibble() %>%
  #     nest(taxa = -ASV), 
  #   by = "ASV"
  # ) %>%
  as.data.frame() %T>%
  print()

ASV_decon <- ASV_df %>%
  decon(numb.blanks = 1,
        numb.ind = 20, # treat samples as one population
        taxa = FALSE,
        runs = 2)

ASV_decon$decon.table
ASV_decon$reads.removed
ASV_decon$OTUs.removed

ASV_decon_df <- ASV_df %>%
  as_tibble() %>%
  full_join(
    ASV_decon$decon.table %>%
      as_tibble() %>%
      rename_with(.cols = -c(ASV, blank),
                  .fn = ~ .x %>% str_c("_decon")),
    by = c("ASV", "blank")
  ) %>%
  pivot_longer(
    cols = -c(ASV, blank),
    names_to = c("sample", "type"),
    names_sep = "_"
  ) %>% 
  mutate(type = type %>% replace_na("original")) %>%
  pivot_wider(names_from = type,
              values_from = value) %>%
  mutate(difference = blank - original, # calculate absolute difference
         removed = original - decon) %T>% # calculate removed reads
  print()


ASV_decon_df %>%
  filter(blank > 0 & original > 0) %>%
  ggplot(aes(difference, removed)) +
    geom_point() +
    theme_minimal()
# No clear trend beetween absolute read difference between blank and
# sample and the number of reads removed from sample, i.e. some ASVs
# had more reads removed by decon() even though they had more reads 
# in the sample than in the blank (negative difference).

ASV_decon_df %>%
  filter(blank > 0 & original > 0) %>%
  ggplot(aes(blank, removed)) +
    geom_abline(slope = 1) +
    geom_point() +
    theme_minimal()
# Removed reads never exceed blank reads (this makes sense), but are
# concentrated towards lower removals, which only makes sense if the
# sample has fewer reads than the blank. Let's look at those cases.

ASV_decon_df %>%
  filter(blank > original) %T>%
  {
    (
      ggplot(., aes(difference, removed)) +
        geom_point() +
        theme_minimal()
    ) %>%
      print()
  } %>%
  print(n = 100)
# Here specifically I'd expect the number of reads the blank has over
# the sample to relate to the number of reads removed but it doesn't.
# Take for example ASV2. It has more reads than any sample but seemingly
# arbitrarily different numbers of reads are removed from each sample:
# e.g. 18 from C15 with originally 1616, 262 from M13P with originally
# 1593 and 303 from M15 with originally only 438.

# Generally disappointed with microDecon too. I'll devise a manual method.
# Create subset phyloseq for blank ASVs
blank <- ASV %>%
  subset_samples(blank == TRUE) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %T>%
  print()

# Mark putative contaminant ASVs
tax_table(ASV) %<>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  mutate(Contaminant = ASV %in% taxa_names(blank)) %>%
  column_to_rownames("ASV") %>%
  as.matrix() %>%
  tax_table()

# Create abundance and identity tibble
ASV_tibble <- otu_table(ASV) %>%
  as.data.frame() %>%
  rownames_to_column("File") %>%
  as_tibble() %>%
  mutate(Sample = File %>% 
           str_split_i(pattern = "_", i = 1)) %>%
  pivot_longer(cols = -c(File, Sample),
               names_to = "ASV",
               values_to = "Abundance") %>%
  full_join(
    tax_table(ASV) %>% 
      as.data.frame() %>%
      rownames_to_column("ASV") %>%
      as_tibble(), by = "ASV"
  ) %T>%
  print()

# Calculate relative abundance
ASV_tibble %<>%
  group_by(Sample) %>%
  mutate(Relative = Abundance / sum(Abundance)) %>%
  ungroup()

# Compare raw reads for blank with mean sample
ASV_tibble %>%
  mutate(Type = if_else(Sample == "blank",
                        "Blank", "Sample")) %>%
  group_by(Type, ASV, Contaminant) %>%
  summarise(Abundance_mean = mean(Abundance)) %>%
  ungroup() %>%
  pivot_wider(names_from = Type, 
              values_from = Abundance_mean) %>%
   ggplot(aes(Blank %>% log10(), 
              Sample %>% log10(),
              colour = Contaminant)) +
    geom_point(shape = 16, size = 2, alpha = 0.3) +
    geom_abline(slope = 1) +
    labs(x = "Blank reads",
         y = "Mean sample reads") +
    theme_minimal()
# Some cases are clear. ASVs that are only present in the
# blank are obviously to be excluded. But many putative 
# contaminant ASVs are present in blank and samples. About 
# half of these are more abundant in the average sample.

# Compare raw reads for blank with each sample
ASV_tibble %>%
  select(Sample, ASV, Contaminant, Abundance) %>%
  pivot_wider(names_from = Sample, 
              values_from = Abundance) %>%
  pivot_longer(cols = -c(ASV, Contaminant, blank),
               names_to = "Sample",
               values_to = "Abundance") %>%
  {
    ggplot(., aes(blank %>% log10(), 
                  Abundance %>% log10(),
                  colour = Contaminant)) +
      geom_point(shape = 16, size = 2, alpha = 0.3) +
      geom_abline(slope = 1) +
      facet_wrap(~ Sample) +
      labs(x = "Blank reads",
           y = "Sample reads") +
      theme_minimal()
  } %>%
  ggsave(filename = "cont_scatter.pdf", 
         path = here("Microbes", "DNA", "Plots"),
         device = cairo_pdf, width = 40, height = 20,
         units = "cm")
# Patterns are quite variable across samples with some cases
# where most putative contaminants are more abundant in the
# sample than in the blank (e.g. M13C). The variability also
# suggests that these ASVs are important indicators of inter-
# sample differences and need to be retained. It's hard to
# distinguish between lab- and cross-contamination but a
# general rule of thumb would be that if any sample has more
# of an ASV than the blank, cross-contamination is more likely.
# Conversely, if the blank has the highest amount of an ASV
# out of all the samples, lab contamination is more likely.
# This can be visualised as the maximum reads out of all 
# samples compared to the blank reads for each ASV.

ASV_tibble %>%
  mutate(Type = if_else(Sample == "blank",
                        "Blank", "Sample")) %>%
  group_by(Type, ASV, Contaminant) %>%
  summarise(Abundance_max = max(Abundance)) %>%
  ungroup() %>%
  pivot_wider(names_from = Type, 
              values_from = Abundance_max) %>%
   ggplot(aes(Blank %>% log10(), 
              Sample %>% log10(),
              colour = Contaminant)) +
    geom_point(shape = 16, size = 2, alpha = 0.3) +
    geom_abline(slope = 1) +
    labs(x = "Blank reads",
         y = "Max sample reads") +
    theme_minimal()
# Now it's clear that most of the shared ASVs likely come 
# from cross-contamination and only some are most abundant 
# in the blank. Let's have a look at relative abundance.

ASV_tibble %>%
  mutate(Type = if_else(Sample == "blank",
                        "Blank", "Sample")) %>%
  group_by(Type, ASV, Contaminant) %>%
  summarise(Relative_max = max(Relative),
            Relative_mean = mean(Relative)) %>%
  ungroup() %>%
  pivot_wider(names_from = Type, 
              values_from = c(Relative_max, Relative_mean)) %T>%
  { (
      ggplot(., aes(Relative_mean_Blank %>% log10(), 
                    Relative_mean_Sample %>% log10(),
                    colour = Contaminant)) +
        geom_point(shape = 16, size = 2, alpha = 0.3) +
        geom_abline(slope = 1) +
        labs(x = "Blank relative abundance",
             y = "Mean sample relative abundance") +
        theme_minimal()
    ) %>% print() } %>%
  ggplot(aes(Relative_max_Blank %>% log10(), 
             Relative_max_Sample %>% log10(),
             colour = Contaminant)) +
    geom_point(shape = 16, size = 2, alpha = 0.3) +
    geom_abline(slope = 1) +
    labs(x = "Blank relative abundance",
         y = "Max sample relative abundance") +
    theme_minimal()
# Perhaps unexpectedly, relative abundance found more
# putative contaminant ASVs to be more abundant in the
# blank because the blank has less total reads. It seems
# best to go with raw reads, assuming that they are 
# similarly amplified in blank and samples. Alternatively,
# relative abundance can be calculated only for the subset
# of ASVs hat occur in the blank and sample, thereby
# removing the bias introduced by the other ASVs in the
# sample (this is one of the initial steps in microDecon).

# Mark exclusive contaminant ASVs
ASV_tibble %<>%
  group_by(ASV) %>%
  mutate(Exclusive = all(Abundance == 0 | Sample == "blank")) %>%
  ungroup()

ASV_tibble %>%
  filter(Exclusive == TRUE) %>%
  group_by(ASV) %>%
  summarise(across(.cols = c(Class, Order, Family, Genus, Species),
                   .fns = unique)) %>%
  print(n = 76)
# 76 exclusive contaminant ASVs, only few of which seem to be of marine
# origin.

# Mark ASVs that are more abundant in blank
ASV_tibble %>%
  filter(Abundance > 0) %>%
  group_by(Sample) %>%
  summarise(n_reads = mean(Abundance),
            n_ASVs = n()) %>%
  ggplot(aes(n_ASVs, n_reads, colour = Sample == "blank")) +
    geom_point() +
    # geom_label(aes(label = Sample)) +
    theme_minimal()
# At low numbers of ASVs, the mean reads per ASV exponentially increase.
# Even though the blank mean reads are comparable to other mid to high
# diversity samples, they are much lower than the low diversity samples.
# Therefore the assumption that ASVs are similarly amplified in blank 
# and samples does not hold and prevalence must be determined using
# the subset method mentioned above.

ASV_shared <- ASV_tibble %>%
  # Filter ASVs that are shared by blank and sample(s).
  filter(Contaminant == TRUE & Exclusive == FALSE) %>% 
  select(Sample, ASV, Abundance) %>%
  pivot_wider(names_from = Sample, 
              values_from = Abundance) %>%
  pivot_longer(cols = -c(ASV, blank),
               names_to = "Sample",
               values_to = "Abundance") %>%
  group_by(Sample) %>%
  mutate(Relative_blank = blank / sum(blank),
         Relative_sample = Abundance / sum(Abundance)) %>%
  ungroup() %T>%
  print()

ASV_shared %>%
  {
    ggplot(., aes(Relative_blank %>% log10(), 
                  Relative_sample %>% log10())) +
      geom_point(shape = 16, size = 2, alpha = 0.3) +
      geom_abline(slope = 1) +
      facet_wrap(~ Sample) +
      labs(x = "Blank reads (proportion of subset)",
           y = "Sample reads (proportion of subset)") +
      theme_minimal()
  } %>%
  ggsave(filename = "cont_scatter_prop.pdf", 
         path = here("Microbes", "DNA", "Plots"),
         device = cairo_pdf, width = 40, height = 20,
         units = "cm")
# Sometimes very similar to raw reads, sometimes flags more
# ASVs as more abundant in blank, sometimes less.

ASV_shared %>%
  group_by(ASV, Relative_blank) %>%
  summarise(Relative_sample_max = max(Relative_sample)) %>%
  ungroup() %>%
  ggplot(aes(Relative_blank %>% log10(), 
             Relative_sample_max %>% log10())) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  geom_abline(slope = 1) +
  labs(x = "Blank reads (proportion of subset)",
       y = "Max sample reads (proportion of subset)") +
  theme_minimal()
# Looks like more ASVs are flagged as most abundant in the
# blank than any of the samples using this method rather
# than absolute reads, i.e. more ASVs to be completely 
# removed rather than just deducted from some samples.

ASV_shared %<>%
  mutate(Prevalent = Relative_blank >= Relative_sample) %T>%
  print(n = 50)

# Update plot to check if correct ASVs are marked
ASV_shared %>%
  {
    ggplot(., aes(Relative_blank %>% log10(), 
                  Relative_sample %>% log10(),
                  colour = Prevalent)) +
      geom_point(shape = 16, size = 2, alpha = 0.3) +
      geom_abline(slope = 1) +
      facet_wrap(~ Sample) +
      labs(x = "Blank reads (proportion of subset)",
           y = "Sample reads (proportion of subset)") +
      theme_minimal()
  } %>%
  ggsave(filename = "cont_scatter_prop.pdf", 
         path = here("Microbes", "DNA", "Plots"),
         device = cairo_pdf, width = 40, height = 20,
         units = "cm")
# Looks fine

# Join prevalence column to main tibble
ASV_tibble %<>%
  left_join(
    ASV_shared %>%
      select(ASV, Sample, Prevalent),
    by = c("ASV", "Sample")
  ) %T>%
  print()

ASV_tibble %<>%
  mutate(Abundance_clean = case_when(
    Prevalent == TRUE ~ 0,
    Prevalent == FALSE ~ Abundance,
    is.na(Prevalent) ~ Abundance
  )) %T>%
  print()

ASV_tibble %>%
  filter(Prevalent) %>%
  select(Sample, ASV, Abundance, Abundance_clean) %>%
  print(n = 50)

ASV_tibble %>%
  filter(Prevalent == FALSE) %>%
  select(Sample, ASV, Abundance, Abundance_clean) %>%
  print(n = 50)

# CLean phyloseq object
# 1. Remove blank and exclusive ASVs
ASV %<>% 
  subset_samples(blank == FALSE) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %T>%
  print()

# 2. Replace otu_table
otu_table(ASV) <- ASV_tibble %>%
  filter(Sample != "blank" & Exclusive == FALSE) %>%
  select(File, ASV, Abundance_clean) %>%
  pivot_wider(names_from = ASV, values_from = Abundance_clean) %>%
  column_to_rownames("File") %>%
  otu_table(taxa_are_rows = FALSE)

# 3. Remove ASVs that are most abundant in blank
ASV %<>%
  prune_taxa(taxa = taxa_sums(.) > 0) %T>%
  print()

# 1.2.2 Chloroplasts (host contamination) ####
ASV %<>%
  subset_taxa(Order != "Chloroplast") %T>%
  print()

ASV %>%
  write_rds(here("Microbes", "DNA", "RDS", "ASV_clean.rds"))
