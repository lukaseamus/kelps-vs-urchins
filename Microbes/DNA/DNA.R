# 1. Check example sequences ####
require(tidyverse)
require(magrittr)
require(here)
require(Biostrings)

# Forward example
example_forward <- here("Microbes", "DNA", "Sequences", 
     "C13-16S_V3-V4_AACLLM5M5_AGTCAGACGA-ACCAGCGACA_L001_R1.fastq.gz") %>%
  readDNAStringSet(format = "fastq") %T>%
  print()

# Reverse example
example_reverse <- here("Microbes", "DNA", "Sequences", 
     "C13-16S_V3-V4_AACLLM5M5_AGTCAGACGA-ACCAGCGACA_L001_R2.fastq.gz") %>%
  readDNAStringSet(format = "fastq") %T>%
  print()

# 2. Load all sequences ####
require(dada2)
sequences <- here("Microbes", "DNA", "Sequences") %>%
  list.files(pattern = "\\.fastq.gz$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Sample = Path %>% basename() %>% str_split_i("-", 1) %>% fct(),
    Direction = if_else(
      Path %>% basename() %>% str_detect("R1"),
      "Forward", "Reverse"
    ) %>% fct()
  ) %T>%
  print()

# 3. Check for primers ####
primer_341F <- DNAString("CCTAYGGGRBGCASCAG") %T>% print()
primer_806R <- DNAString("GGACTACNNGGGTATCTAAT") %T>% print()

# Check forward sequence with forward primer
example_forward %>%
  sample(1) %>% .[[1]] %>%
  matchPattern(pattern = primer_341F, fixed = FALSE) 
# Only one match that starts at the beginning of the DNA sequence.
# Aligned this looks like this:
example_forward %>%
  sample(1) %>% .[[1]] %>%
  pairwiseAlignment(pattern = primer_341F)

# Check forward sequence with reverse primer
example_forward %>%
  sample(1) %>% .[[1]] %>%
  matchPattern(pattern = primer_806R, fixed = FALSE) 
# As expected, no match.

# Check reverse sequence with reverse primer
example_reverse %>%
  sample(1) %>% .[[1]] %>%
  matchPattern(pattern = primer_806R, fixed = FALSE) 
# Only one match that starts at the beginning of the DNA sequence.
example_reverse %>%
  sample(1) %>% .[[1]] %>%
  pairwiseAlignment(pattern = primer_806R)

# Check reverse sequence with forward primer
example_reverse %>%
  sample(1) %>% .[[1]] %>%
  matchPattern(pattern = primer_341F, fixed = FALSE) 
# No match.

# So I can simply remove the first 17 bases from the forward primer
# and the first 20 bases from the reverse primer.

# 4. Quality control ####
require(glue)
sequences %<>%
  rowwise() %>% # rowwise makes plotting from many variables easier
  mutate( 
    # Here I include the sample name and read direction along with
    # a conditional guessed trim for reference.
    Quality_Plot = 
      list(
        Path %>%
          plotQualityProfile() +
            geom_vline(xintercept = if (Direction == "Forward") 
                                    290 else 220) +
            ggtitle(glue("{Sample} ({Direction})")) 
            # Alternatively, ggtitle(Sample, Direction) creates 
            # title and subtitle.
      )
  ) %>%
  ungroup() # rowwise needs to be undone

require(patchwork)
sequences %$%
  wrap_plots(Quality_Plot) %>%
  ggsave(filename = "Quality_Plot.pdf", device = cairo_pdf, 
         path = here("Microbes", "DNA", "Plots"),
         height = 60, width = 60, units = "cm")
# Trim positions 290 and 220 generally look good.

# 5. Filter and trim ####
sequences %<>%
  pivot_wider(names_from = Direction, values_from = c(Path, Quality_Plot)) %>%
  mutate(
    Path_Forward_Filtered = here("Microbes", "DNA", "Sequences", "Filtered") %>%
      file.path(str_c(Sample, "_F_filt.fastq.gz")),
    Path_Reverse_Filtered = here("Microbes", "DNA", "Sequences", "Filtered") %>%
      file.path(str_c(Sample, "_R_filt.fastq.gz"))
  ) %T>%
  print()

filter <- sequences %$%
  filterAndTrim( # The paths need to be named to be matched correctly.
    fwd = Path_Forward %>% setNames(Sample), filt = Path_Forward_Filtered %>% setNames(Sample),
    rev = Path_Reverse %>% setNames(Sample), filt.rev = Path_Reverse_Filtered %>% setNames(Sample), 
    truncLen = c(290, 220), # Trim forward to 290 bases and reverse to 220
    trimLeft = c(length(primer_341F), length(primer_806R)), # Remove primers
    maxN = 0, maxEE = c(2, 2), truncQ = 2, 
    rm.phix = TRUE, compress = TRUE, 
    matchIDs = TRUE, multithread = TRUE
  )

sequences %<>%
  mutate(Reads = filter[,1],
         Reads_Filtered = filter[,2],
         Prop_Filtered = Reads_Filtered / Reads) %T>%
  print()

rm(filter)

sequences %$% mean(Prop_Filtered)
# Around 81% of reads survived filtering.

sequences %>%
  ggplot(aes(Reads, Reads_Filtered)) +
      geom_abline(slope = c(1, sequences %$% mean(Prop_Filtered)),
                  colour = c("black", "grey")) +
      geom_label(aes(label = Sample)) +
      theme_minimal()

# 6. Quality control ####
sequences %<>%
  rowwise() %>%
  mutate(
    Quality_Plot_Forward_Filtered = 
      list(
        Path_Forward_Filtered %>%
          plotQualityProfile() +
            ggtitle(glue("{Sample} (Forward)"))
      ),
    Quality_Plot_Reverse_Filtered = 
      list(
        Path_Reverse_Filtered %>%
          plotQualityProfile() +
            ggtitle(glue("{Sample} (Reverse)"))
      )
  ) %>%
  ungroup()

sequences %$%
  wrap_plots(c(Quality_Plot_Forward_Filtered, 
               Quality_Plot_Reverse_Filtered)) %>%
  ggsave(filename = "Quality_Plot_Filtered.pdf", device = cairo_pdf, 
         path = here("Microbes", "DNA", "Plots"),
         height = 60, width = 60, units = "cm")
# Looks good.

# 7. Learn error rates ####
error_forward <- sequences %$% 
  learnErrors(Path_Forward_Filtered,
              multithread = TRUE,
              randomize = TRUE)
error_reverse <- sequences %$% 
  learnErrors(Path_Reverse_Filtered,
              multithread = TRUE,
              randomize = TRUE)

plotErrors(error_forward, nominalQ = TRUE)
plotErrors(error_reverse, nominalQ = TRUE)

# 8. Infer ASVs ####
# dada_forward <- sequences %$%
#   dada(
#     derep = Path_Forward_Filtered,
#     err = error_forward,
#     pool = "pseudo",
#     multithread = TRUE
#   )
# 
# dada_reverse <- sequences %$% 
#   dada(
#     derep = Path_Reverse_Filtered,
#     err = error_reverse,
#     pool = "pseudo",
#     multithread = TRUE
#   )
# 
# dada_forward[[1]]
# dada_reverse[[1]]

sequences %<>%
  mutate(
    DADA_Forward = dada(
      derep = Path_Forward_Filtered,
      err = error_forward,
      pool = "pseudo",
      multithread = TRUE
    ),
    DADA_Reverse = dada(
      derep = Path_Reverse_Filtered,
      err = error_reverse,
      pool = "pseudo",
      multithread = TRUE
    )
  ) %T>%
  print()
  
# 9. Merge forward and reverse ####
# merge <- sequences %$%
#   mergePairs(
#   dadaF = dada_forward, derepF = Path_Forward_Filtered %>% setNames(Sample),
#   dadaR = dada_reverse, derepR = Path_Reverse_Filtered %>% setNames(Sample),
#   verbose = TRUE
# )

sequences %<>%
  mutate(
    DADA_Merged = mergePairs(
      dadaF = DADA_Forward, derepF = Path_Forward_Filtered,
      dadaR = DADA_Reverse, derepR = Path_Reverse_Filtered,
      verbose = TRUE
    )
  ) %T>%
  print()

# 10. Sequence table ####
# Create naive sequence table
sequence_table <- sequences %$% 
  makeSequenceTable(DADA_Merged)

sequence_table %>%
  getSequences() %>%
  str_length() %>%
  tibble(Basepairs = .) %>%
  count(Basepairs, name = "ASVs") %T>%
  print(n = 63) %>%
  plot()

sequence_table %>%
  as_tibble() %>%
  print(n = 21)

# Remove chimeras
sequence_table_clean <- removeBimeraDenovo(
  sequence_table, method = "consensus", 
  multithread = TRUE, verbose = TRUE
)

sequence_table_clean %>%
  getSequences() %>%
  str_length() %>%
  tibble(Basepairs = .) %>%
  count(Basepairs, name = "ASVs") %T>%
  print(n = 63) %>%
  plot()

sequence_table_clean %>%
  as_tibble() %>%
  print(n = 21)

1 - ncol(sequence_table_clean) / ncol(sequence_table)
# 96% of inferred ASVs are chimeras
1 - sum(sequence_table_clean) / sum(sequence_table)
# But chimeras only account for 19% of merged reads

# 11. Track reads and ASVs ####
# Count reads at various stages after filtering
sequences %<>%
  mutate(Reads_Forward_Denoised = DADA_Forward %>% 
           map_dbl(~ .x %>% getUniques() %>% sum()),
         Reads_Reverse_Denoised = DADA_Reverse %>% 
           map_dbl(~ .x %>% getUniques() %>% sum()),
         Reads_Merged = DADA_Merged %>% 
           map_dbl(~ .x %>% getUniques() %>% sum()),
         Reads_Sequence_Table = sequence_table %>%
           rowSums(),
         Reads_Sequence_Table_Clean = sequence_table_clean %>%
           rowSums()) %T>%
  print()

sequences %>%
  select(Sample, starts_with("Reads")) %>%
  print(n = 21)

sequences %<>%
  mutate(Prop_Chimera = 1 - Reads_Sequence_Table_Clean / Reads_Sequence_Table,
         Prop_Retained = Reads_Sequence_Table_Clean / Reads) %T>%
  print()

sequences %>%
  select(Sample, starts_with("Prop")) %>%
  print(n = 21)

# Count ASVs per sample before and after chimera removal
sequences %<>%
  mutate(ASVs = rowSums(sequence_table > 0),
         ASVs_Clean = rowSums(sequence_table_clean > 0),
         Prop_ASVs = ASVs_Clean / ASVs) %T>%
  print()

sequences %>%
  select(Sample, contains("ASV")) %>%
  print(n = 21)

# 12. Assign taxonomy ####
# 12.1 Genus-level naive Bayesian classification ####
taxa_table <- assignTaxonomy(
  seqs = sequence_table_clean,
  refFasta = here("Microbes", "DNA", "SILVA", 
                  "silva_nr99_v138.2_toGenus_trainset.fa.gz"),
  multithread = TRUE
)

taxa_table %>% as_tibble()

# 12.2 Species exact matching ####
taxa_table %<>%
  addSpecies(refFasta = here("Microbes", "DNA", "SILVA", 
                             "silva_v138.2_assignSpecies.fa.gz"),
             allowMultiple = TRUE)

taxa_table %>% as_tibble()

species_table <- assignSpecies(
  seqs = sequence_table_clean,
  refFasta = here("Microbes", "DNA", "SILVA", 
                  "silva_v138.2_assignSpecies.fa.gz"),
  allowMultiple = TRUE
)

species_table %>% as_tibble()

taxa_table %>%
  as_tibble() %>%
  select(Genus, Species) %>%
  bind_cols(species_table %>% 
              as_tibble()) %>%
  print(n = 1000)

# Clearly addSpecies() is not independent because mistakes made by
# assignTaxonomy() at the genus level are carried through. It is now
# possible to run naive Bayesian classification to species level:

# 12.3 Species-level naive Bayesian classification ####
taxa_table_species <- assignTaxonomy(
  seqs = sequence_table_clean,
  refFasta = here("Microbes", "DNA", "SILVA", 
                  "silva_nr99_v138.2_toSpecies_trainset.fa.gz"),
  multithread = TRUE
)

taxa_table_species %>% as_tibble()

taxa_table_species %>%
  as_tibble() %>%
  select(Genus, Species) %>%
  bind_cols(species_table %>% 
              as_tibble()) %>%
  print(n = 1000)

# 12.4 Resolve species-level classification ####
# The best native dada2 taxonomy assignment seems to be to run the 
# naive Bayesian classifier to species level and independently run 
# exact species matching. The two separate estiamtes can then be 
# cross-validated to get the most complete ASV table. When only
# one approach yields species, that species is accepted. Similarly,
# when both appraoches agree, there is no conflict. However,
# when both approaches yield different species, the choice is made 
# according to these rules:
# 1. If genus and/or species clash, exact matching is trusted
#    over naive Bayesian classification.
# 2. If genus aligns and exact matching identifies several species
#    including the one identified by naive Bayesian classification
#    the single species identified by the latter method is trusted.
# 3. If none of the species identified by exact matching match the
#    species identified by naive Bayesian classification, the list 
#    of species is favoured over the single species.

taxa_table_resolved <- 
  taxa_table_species %>%
  as_tibble() %>%
  bind_cols(
    species_table %>%
      as_tibble() %>% # rename also doesn't work here
      mutate(Genus_exact = Genus, Species_exact = Species) %>%
      select(Genus_exact, Species_exact)
  ) %>%
  mutate(
    Genus_resolved =
      case_when(
        !is.na(Genus) & is.na(Genus_exact) ~ Genus,
        is.na(Genus) & !is.na(Genus_exact) ~ Genus_exact,
        Genus == Genus_exact ~ Genus,
        Genus != Genus_exact ~ Genus_exact
      ),
    Species_resolved =
      case_when(
        !is.na(Species) & is.na(Species_exact) ~ Species,
        is.na(Species) & !is.na(Species_exact) ~ Species_exact,
        Species == Species_exact ~ Species,
        Species != Species_exact & 
          !Species_exact %>% str_detect("/") ~ Species_exact,
        Species_exact %>% str_detect("/") &
          Species_exact %>% str_detect(Species) ~ Species,
        Species_exact %>% str_detect("/") &
          !Species_exact %>% str_detect(Species) ~ Species_exact
      )
  ) %T>%
  print()

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  print(n = 500)
# Looks fine but better double-check specific cases.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(Genus == Genus_exact & Species == Species_exact) %>%
  print(n = 100)
# 68 exact genus and species matches.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(Genus == Genus_exact & Species != Species_exact) %>%
  print(n = 100)
# 61 cases where genus matches but species clashes. Most cases
# have several options from exact matching that include the naive
# Bayesian classifier species. In few cases there is a complete
# species mismatch despite genus match. Everything os resolved
# properly.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(Genus != Genus_exact) %>%
  print(n = 100)
# Few cases where genus is actually mismatched. Mostly different
# spelling. The only case I'll manually correct is Clostridium
# innocuum, which is identified by the naive Bayesian classifier
# but not properly split into the binomial. Exact matching could
# not differentiate between Clostridium aff. and C. innocuum.

taxa_table_resolved %<>%
  mutate(Species_resolved = if_else(
    Species_resolved == "aff./innocuum",
    "innocuum", Species_resolved
    ))

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(Genus != Genus_exact) %>%
  print(n = 100)
# Resolved.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(is.na(Species) & !is.na(Species_exact)) %>%
  print(n = 100)
# 79 cases where naive Bayesian classifier didn't get the species
# but exact matching did. Note that sometimes exact matching didn't
# get the genus despite getting the species.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(is.na(Genus_exact) & !is.na(Species_exact)) %>%
  print(n = 100)
# 6 cases where exact matching got species but not genus. Based
# on two cases where binomials match between methods, the genus
# was taken from the naive Bayesian classifier in all cases.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(!is.na(Species) & is.na(Species_exact)) %>%
  print(n = 500)
# 389 cases where naive Bayesian classifier got the species but
# exact matching didn't.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(is.na(Genus) & !is.na(Species_exact)) %>%
  print(n = 100)
# Only 4 cases where naive Bayesian classifier got neither genus nor
# species but exact matching did.

# 12.5 New DECIPHER classification algorithm ####
require(DECIPHER)
# The pipe doesn't work with load().
load(here("Microbes", "DNA", "SILVA", 
          "SILVA_SSU_r138_2_2024.RData"))

taxa_table_decipher <- IdTaxa(
  test = sequence_table_clean %>%
    getSequences() %>%
    DNAStringSet(),
  trainingSet = trainingSet,
  strand = "top", processors = NULL
)

# Reformat to match taxa_table.
taxa_table_decipher %<>%
  map(
    ~ {
      m <- match(c("domain", "phylum", "class", "order", 
                   "family", "genus", "species"), .x$rank)
      taxa <- .x$taxon[m]
      taxa[startsWith(taxa, "unclassified_")] <- NA
      taxa
    }
  ) %>% 
  simplify2array() %>%
  t() %>%
  `colnames<-`( c("domain", "phylum", "class", "order", 
                  "family", "genus", "species") ) %>%
  `rownames<-`( getSequences(sequence_table_clean) )

taxa_table_decipher %>% as_tibble()
  
taxa_table_resolved %>%
  select(Genus_resolved, Species_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus, species)
  ) %T>%
  View()
# Generally DECIPHER seems to be worse at classifying taxa, no 
# matter what the paper (doi 10.1186/s40168-018-0521-5) says.
# Let's have a closer look.

taxa_table_decipher %>%
  as_tibble() %>%
  drop_na(species)
# It didn't get a single species.

taxa_table_resolved %>%
  select(Genus_resolved, Species_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus, species)
  ) %>%
  filter(is.na(Genus_resolved) & !is.na(genus)) %>%
  print(n = 100)
# 264 times DECIPHER got the genus where the older method didn't,
# but most are Incertae Sedis, so placeholders.

taxa_table_resolved %>%
  select(Genus_resolved, Species_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus, species)
  ) %>%
  filter(!is.na(Genus_resolved) & is.na(genus)) %>%
  print(n = 100)
# 913 times the older method got the genus but DECIPHER didn't.

taxa_table_resolved %>%
  select(Genus_resolved, Species_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus, species)
  ) %>%
  filter(Genus_resolved == genus) %>%
  print(n = 100)
# Genus matched 403 times.

taxa_table_resolved %>%
  select(Genus_resolved, Species_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus, species)
  ) %>%
  filter(Genus_resolved != genus) %>%
  print(n = 100)
# 1057 genus mismatches, but mostly due to different spelling.

taxa_table_resolved %>%
  select(Genus_resolved, Species_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus, species)
  ) %>%
  filter(Genus_resolved != genus &
           genus %>% str_detect(Genus_resolved %>% 
                                  fixed())) %>%
  print(n = 100)
# 1020 cases (96%) it's due to different spelling 
# (the old name is contained within the new one).

taxa_table_resolved %>%
  select(Genus_resolved, Species_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus, species)
  ) %>%
  filter(Genus_resolved != genus &
           !genus %>% str_detect(Genus_resolved %>% 
                                   fixed())) %>%
  print(n = 100)
# In 37 cases it's actually different, but I trust exact matching
# cross-validated with the older method more.

# 12.6 Resolve all methods ####
# The only part that would improve my current taxa table would be 
# adding the 264 missing genera, but only if the other taxa match
# or are absent in the older method.
taxa_table_resolved %<>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>% # rename also doesn't work here
      mutate(Phylum_decipher = phylum, Class_decipher = class, 
             Order_decipher = order, Family_decipher = family, 
             Genus_decipher = genus) %>%
      select(Phylum_decipher, Class_decipher, Order_decipher, 
             Family_decipher, Genus_decipher)
  ) %>%
  mutate(
    Genus_resolved_decipher =
      if_else(
        is.na(Genus_resolved) &
          ( Phylum_decipher == Phylum | 
              Phylum_decipher %>%
              str_detect(Phylum %>% fixed()) |
              is.na(Phylum) ) &
          ( Class_decipher == Class | 
              Class_decipher %>%
                str_detect(Class %>% fixed()) |
              is.na(Class) ) &
          ( Order_decipher == Order | 
              Order_decipher %>%
                str_detect(Order %>% fixed()) |
              is.na(Order) ) &
          ( Family_decipher == Family | 
              Family_decipher %>%
                str_detect(Family %>% fixed()) |
              is.na(Family) ),
        Genus_decipher, Genus_resolved
      )
  ) %T>%
  print()

taxa_table_resolved %>%
  filter(is.na(Genus_resolved) & !is.na(Genus_decipher)) %>%
  View()
# Only in two cases the older method had a different order and/or family
# assigned to an ASV so the NA was not replaced with the genus name.

# Resolve other taxa where they genus was added
taxa_table_resolved %<>%
  mutate(
    Phylum = if_else(is.na(Genus_resolved) & !is.na(Genus_resolved_decipher),
                     coalesce(Phylum, Phylum_decipher), Phylum),
    Class = if_else(is.na(Genus_resolved) & !is.na(Genus_resolved_decipher),
                    coalesce(Class, Class_decipher), Class),
    Order = if_else(is.na(Genus_resolved) & !is.na(Genus_resolved_decipher),
                    coalesce(Order, Order_decipher), Order),
    Family = if_else(is.na(Genus_resolved) & !is.na(Genus_resolved_decipher),
                     coalesce(Family, Family_decipher), Family)
  )
  
# Clean taxa_table_resolved  
taxa_table_resolved %<>%
  select(Kingdom, Phylum, Class, Order, Family, 
         Genus_resolved_decipher, Species_resolved) %>%
  mutate(Genus = Genus_resolved_decipher, Species = Species_resolved) %>%
  select(-c(Genus_resolved_decipher, Species_resolved)) %T>%
  print()

# Calculate final resolution
taxa_table_resolved %>%
  filter(!is.na(Species)) %>%
  count() /
  taxa_table_resolved %>%
  count()
# 18% of ASVs are classified to species.

taxa_table_resolved %>%
  filter(!is.na(Genus)) %>%
  count() /
  taxa_table_resolved %>%
  count()
# 80% of ASVs are classified to genus.

# Convert back to matrix
taxa_table_species %>% str()
taxa_table_decipher %>% str()
# This is the matrix structure I need.

taxa_table_resolved %<>%
  as.matrix() %>%
  `rownames<-`( rownames(taxa_table_species) ) %T>%
  # Alternatively, get sequences from sequence_table_clean.
  str()
# Looks fine.

# 13. Save ASVs as phyloseq object ####
# Add metadata
time <- here("Urchins", "Biomass.csv") %>%
  read_csv() %>%
  filter(Season == "Autumn") %>%
  mutate(Deployment = Deployment %>% dmy_hm(),
         Retrieval = Retrieval %>% dmy_hm(),
         Days = Deployment %--% Retrieval / ddays()) %>%
  select(Tank, Days) %T>%
  print(n = 45)

sequences %<>%
  mutate(
    Tank = if_else(Sample == "blank",
                   NA, Sample %>% str_sub(end = 3)),
    Treatment = case_when(
      Tank %>% str_detect("C") ~ "Control",
      Tank %>% str_detect("M") ~ "Mechanical",
      Tank %>% str_detect("U") ~ "Urchin"
    ),
    Antibiotic = case_when(
      Sample %>% str_sub(start = 4) == "C" ~ "Control",
      Sample %>% str_sub(start = 4) == "P" ~ "Penicillin",
      Sample %>% str_sub(start = 4) == "N" ~ "Nystatin",
      Sample %>% str_sub(start = 4) == "PN" ~ "Penicillin + Nystatin"
    ),
    Tissue = case_when(
      Sample %>% str_detect("f") ~ "Faeces",
      Sample == "blank" ~ NA,
      TRUE ~ "Kelp"
    )
  ) %>%
  left_join(time, by = "Tank") %T>%
  print(n = 21)

# Build phyloseq object
require(phyloseq)
ASV <- phyloseq(
  otu_table(
    sequence_table_clean, 
    taxa_are_rows = FALSE
  ),
  sample_data(
    sequences %>% 
      select(Sample, Tank, Treatment,
             Antibiotic, Tissue) %>%
      as.data.frame() %>%
      `rownames<-`( rownames(sequence_table_clean) )
  ),
  tax_table(
    taxa_table_resolved
  )
) %T>%
  print()

# Simplify ASV names and save DNA in refseq
ASV %<>%
  taxa_names() %>%
  DNAStringSet() %>%
  `names<-`( ASV %>% taxa_names() ) %>%
  merge_phyloseq(ASV) %>%
  `taxa_names<-`( str_c("ASV", ASV %>% ntaxa() %>% seq()) ) %T>%
  print()

# Simple test
ASV %>%
  plot_richness(x = "Treatment", color = "Antibiotic",
                measures = "Simpson")
# Seems to work.

# Save phyloseq object
ASV %>%
  write_rds(file = here("Microbes", "DNA", "RDS", "ASV.rds"))