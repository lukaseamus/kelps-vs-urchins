# 1. Prepare data ####
# 1.1 Check example sequence ####
require(tidyverse)
require(here)
require(Biostrings)

# Forward example
here("Microbes", "DNA", "Sequences", 
     "C13-16S_V3-V4_AACLLM5M5_AGTCAGACGA-ACCAGCGACA_L001_R1.fastq.gz") %>%
  readDNAStringSet(format = "fastq")

# Reverse example
here("Microbes", "DNA", "Sequences", 
     "C13-16S_V3-V4_AACLLM5M5_AGTCAGACGA-ACCAGCGACA_L001_R2.fastq.gz") %>%
  readDNAStringSet(format = "fastq")

# 1.2 Load sequences ####
require(magrittr)
require(dada2)
sequences <- here("Microbes", "DNA", "Sequences") %>%
  list.files(pattern = "\\.fastq.gz$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Sample = Path %>% basename() %>% str_split_i("-", 1),
    Direction = if_else(
      Path %>% basename() %>% str_detect("R1"),
      "Forward", "Reverse"
    ) %>% fct()
  ) %T>%
  print()

# 1.3 Quality control ####
sequences %<>%
  mutate(Quality_Plot = Path %>% map(plotQualityProfile))

require(patchwork)
sequences %>%
  rowwise() %>%
  mutate(Quality_Plot_Labelled = 
           list(
             Quality_Plot +
               ggtitle(Sample, Direction)
             )) %$%
  wrap_plots(Quality_Plot_Labelled) %>%
  ggsave(filename = "Quality_Plot.pdf", device = cairo_pdf, 
         path = here("Microbes", "DNA", "Plots"),
         height = 60, width = 60, units = "cm")

# 1.4 Filter and trim ####
sequences %<>%
  pivot_wider(names_from = Direction, values_from = c(Path, Quality_Plot)) %>%
  mutate(
    Path_Forward_Filtered = here("Microbes", "DNA", "Sequences") %>%
      file.path(str_c(Sample, "_F_filt.fastq.gz")),
    Path_Reverse_Filtered = here("Microbes", "DNA", "Sequences") %>%
      file.path(str_c(Sample, "_R_filt.fastq.gz"))
  ) %T>%
  print()

sequences %<>%
  rowwise() %>%
  mutate(Filter = filterAndTrim(
    fwd = Path_Forward, filt = Path_Forward_Filtered,
    rev = Path_Reverse, filt.rev = Path_Reverse_Filtered, 
    truncLen = c(290, 200), maxN = 0, maxEE = c(2, 2), 
    truncQ = 2, rm.phix = TRUE, compress = TRUE, 
    multithread = TRUE
    )) %>%
  ungroup()

sequences %<>%
  mutate(Reads = Filter[,1],
         Reads_Filtered = Filter[,2]) %>%
  select(-Filter)

sequences %>%
  ggplot(aes(Reads, Reads_Filtered)) +
    geom_label(aes(label = Sample)) +
    geom_abline(slope = 1) +
    theme_minimal()

# 1.5 Learn error rates ####
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

# 1.6 Infer ASVs ####
sequences %<>%
  mutate(DADA_Forward = dada(
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
         )) %T>%
  print()

sequences$DADA_Forward[[1]]
sequences$DADA_Reverse[[1]]

# 1.7 Merge forward and reverse ####
sequences %<>%
  mutate(Merge = mergePairs(
    dadaF = DADA_Forward, derepF = Path_Forward_Filtered,
    dadaR = DADA_Reverse, derepR = Path_Reverse_Filtered,
    verbose = TRUE
  )) %T>%
  print()

sequences$Merge[[1]]

# 1.8 Make sequence table ####
sequences %<>%
  mutate(Sequence_Table = makeSequenceTable(Merge)) %T>%
  print()

sequences$Sequence_Table[[1]]

seqtab <- sequences %$% makeSequenceTable(Merge)
dim(seqtab)

seqtab %>%
  getSequences() %>%
  str_length() %>%
  tibble(Basepairs = .) %>%
  count(Basepairs, name = "ASVs") %>%
  print(n = 55)

seqtab_clean <- removeBimeraDenovo(seqtab, method = "consensus", 
                                   multithread = TRUE, verbose = TRUE)
dim(seqtab_clean)
