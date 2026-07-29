library(phyloseq)
library(dplyr)
library(tibble)
library(Biostrings)

# ============================================================
# 0. Settings
# ============================================================
# PS is the object everything below runs on. `euks` is the full 659-sample
# set with raw cruise codes (AZ19/AZ22/CE22/CE24). Switch to `euk2` if these
# sequences need to correspond to the taxa in Figure 3 — that's the
# station / D1 / lat-long-boxed subset with recoded cruise names.
PS <- euks

N          <- 20   # ASVs per cruise
MIN_PREV   <- 1    # ASV must appear in at least this many samples of a cruise
FASTA_FILE <- "labsea_18S_all_ASVs.fasta"   # built by build_combined_fasta.R
PREFIX     <- "euk18S"

# ============================================================
# 1. Top N ASVs per cruise (by mean relative abundance)
# ============================================================
PS      <- prune_samples(sample_sums(PS) > 0, PS)
ps_rel  <- transform_sample_counts(PS, function(x) x / sum(x))

meta    <- as.data.frame(sample_data(ps_rel))
cruises <- sort(unique(na.omit(meta$cruise_short)))

topN_list <- lapply(cruises, function(cr) {
  keep  <- rownames(meta)[which(meta$cruise_short == cr)]
  ps_cr <- prune_samples(keep, ps_rel)
  ps_cr <- prune_taxa(taxa_sums(ps_cr) > 0, ps_cr)
  
  # prevalence filter: drops singletons that can otherwise tie their way
  # into the bottom of the list on the strength of one sample
  prev  <- rowSums(otu_table(ps_cr) > 0)
  if (!taxa_are_rows(ps_cr)) prev <- colSums(otu_table(ps_cr) > 0)
  ps_cr <- prune_taxa(prev >= MIN_PREV, ps_cr)
  
  abund <- sort(taxa_sums(ps_cr) / nsamples(ps_cr), decreasing = TRUE)
  n     <- min(N, length(abund))
  vals  <- unname(abund[1:n])
  
  data.frame(
    cruise         = cr,
    ASV            = names(abund)[1:n],
    rank           = 1:n,
    mean_rel_abund = round(vals * 100, 3),
    n_samples      = nsamples(ps_cr),
    # TRUE where the value is shared with another ASV in this cruise, i.e.
    # the ordering at that position is arbitrary
    tied           = vals %in% vals[duplicated(vals)],
    row.names      = NULL
  )
})

topN_by_cruise <- bind_rows(topN_list)
print(topN_by_cruise)

if (any(topN_by_cruise$tied)) {
  cat("\nTied ranks (arbitrary ordering) by cruise:\n")
  topN_by_cruise %>%
    filter(tied) %>%
    count(cruise, mean_rel_abund, name = "n_ASVs") %>%
    print()
}

# ============================================================
# 2. Compile: union of all top Ns, with taxonomy
# ============================================================
tax_df <- as.data.frame(tax_table(PS)) %>% rownames_to_column("ASV")

compiled <- topN_by_cruise %>%
  group_by(ASV) %>%
  summarise(
    n_cruises  = n(),
    cruises    = paste(cruise, collapse = "; "),
    best_rank  = min(rank),
    mean_abund = round(mean(mean_rel_abund), 3),
    any_tied   = any(tied),
    .groups    = "drop"
  ) %>%
  left_join(tax_df, by = "ASV") %>%
  arrange(desc(n_cruises), best_rank)

cat("\nUnique ASVs across all cruise top", N, ":", nrow(compiled), "\n")
cat("In top", N, "of ALL", length(cruises), "cruises:",
    sum(compiled$n_cruises == length(cruises)), "\n")
print(compiled, n = Inf)

write.csv(topN_by_cruise, paste0(PREFIX, "_top", N, "_by_cruise.csv"), row.names = FALSE)
write.csv(compiled,       paste0(PREFIX, "_top", N, "_compiled.csv"),  row.names = FALSE)

# ============================================================
# 3. Pull sequences from the combined FASTA
# ============================================================
stopifnot(file.exists(FASTA_FILE))

fasta <- readDNAStringSet(FASTA_FILE)
cat("\nTotal sequences in fasta:", length(fasta), "\n")

missing <- setdiff(compiled$ASV, names(fasta))
if (length(missing)) {
  cat("Missing from fasta:\n"); print(missing)
  stop("FASTA does not cover every selected ASV — rebuild it before continuing.")
}

# order sequences to match the compiled table rather than FASTA order
top_subset <- fasta[compiled$ASV]
cat("Top", N, "ASVs found in fasta:", length(top_subset), "/", nrow(compiled), "\n")

writeXStringSet(top_subset, paste0(PREFIX, "_top", N, "_compiled.fa"))

# Second copy with taxonomy in the header — easier to read in a tree viewer
# or BLAST result than a bare hash. Original IDs stay in the file above.
annot <- compiled %>%
  mutate(header = paste(ASV,
                        ifelse(is.na(Genus) | Genus == "" |
                                 grepl("unidentified", Genus, ignore.case = TRUE),
                               Phylum, Genus),
                        paste0("nCruise", n_cruises),
                        sep = "_"))

named_subset <- top_subset
names(named_subset) <- annot$header
writeXStringSet(named_subset, paste0(PREFIX, "_top", N, "_compiled_named.fa"))