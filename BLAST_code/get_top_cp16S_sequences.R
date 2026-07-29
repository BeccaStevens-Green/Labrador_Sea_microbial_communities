library(phyloseq)
library(dplyr)
library(tibble)
library(Biostrings)

# ============================================================
# 1. Top 20 cp16S ASVs per cruise (by mean relative abundance)
# ============================================================
N <- 20   # change here to adjust the cutoff everywhere

chloro     <- prune_samples(sample_sums(chloro) > 0, chloro)
chloro_rel <- transform_sample_counts(chloro, function(x) x / sum(x))
meta       <- as.data.frame(sample_data(chloro_rel))   # refresh after pruning
cruises    <- sort(unique(na.omit(meta$cruise_short)))

top20_list <- lapply(cruises, function(cr) {
  keep  <- rownames(meta)[which(meta$cruise_short == cr)]
  ps_cr <- prune_samples(keep, chloro_rel)
  ps_cr <- prune_taxa(taxa_sums(ps_cr) > 0, ps_cr)
  
  abund <- sort(taxa_sums(ps_cr) / nsamples(ps_cr), decreasing = TRUE)
  n     <- min(N, length(abund))
  
  data.frame(
    cruise         = cr,
    ASV            = names(abund)[1:n],
    rank           = 1:n,
    mean_rel_abund = round(unname(abund[1:n]) * 100, 3),
    row.names      = NULL
  )
})

top20_by_cruise <- bind_rows(top20_list)
print(top20_by_cruise)

# ============================================================
# 2. Compile: union of all top 20s, with taxonomy
# ============================================================
tax_df <- as.data.frame(tax_table(chloro)) %>% rownames_to_column("ASV")

compiled <- top20_by_cruise %>%
  group_by(ASV) %>%
  summarise(
    n_cruises  = n(),
    cruises    = paste(cruise, collapse = "; "),
    best_rank  = min(rank),
    mean_abund = round(mean(mean_rel_abund), 3),
    .groups    = "drop"
  ) %>%
  left_join(tax_df, by = "ASV") %>%
  arrange(desc(n_cruises), best_rank)

cat("Unique ASVs across all cruise top", N, ":", nrow(compiled), "\n")
cat("In top", N, "of ALL", length(cruises), "cruises:",
    sum(compiled$n_cruises == length(cruises)), "\n")
print(compiled)

write.csv(top20_by_cruise, "cp16S_top20_by_cruise.csv", row.names = FALSE)
write.csv(compiled,        "cp16S_top20_compiled.csv",  row.names = FALSE)

# ============================================================
# 3. Pull sequences from the FASTA
# ============================================================
fasta <- readDNAStringSet("../16S-fuhrman/250610-1850.labsea.16S.dna-sequences.fasta")
cat("Total sequences in fasta:", length(fasta), "\n")

top_subset <- fasta[names(fasta) %in% compiled$ASV]
cat("Top", N, "ASVs found in fasta:", length(top_subset), "/", nrow(compiled), "\n")

missing <- setdiff(compiled$ASV, names(fasta))
if (length(missing) > 0) {
  cat("Missing from fasta:\n"); print(missing)
}

writeXStringSet(top_subset, "cp16S_top20_compiled.fa")