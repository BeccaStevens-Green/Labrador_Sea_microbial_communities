library(Biostrings)
library(phyloseq)

# ============================================================
# Build one combined FASTA per marker, then check it covers
# every ASV in the phyloseq object.
# Run once; downstream scripts just read the combined file.
# ============================================================

BASE <- "/Users/beccastevens-green/Desktop/R master/projects/R_CES22_DNA/MicrobiomeAnalyses/amplicon_2025"

FASTA_18S <- c(
  file.path(BASE, "18S-fuhrman/250609-2016.labsea.18S.dna-sequences.fasta"),
  file.path(BASE, "ce24_addition/18S-fuhrman/250818-1749.ce24.18S.dna-sequences.fasta")
)

OUT_18S <- "labsea_18S_all_ASVs.fasta"   # written to getwd() (project_LabSea_2)

# ------------------------------------------------------------
# Merge, with a hash-collision check
# ------------------------------------------------------------
# ASV IDs are MD5 hashes of the sequence, so the same ID appearing in two
# files must carry the same sequence. If it doesn't, the two runs used
# different denoising or trimming and the objects shouldn't be merged
# on ID at all — hence the stop() rather than a warning.

build_fasta <- function(paths, out_path) {
  
  stopifnot(all(file.exists(paths)))
  
  parts <- lapply(paths, function(p) {
    s <- readDNAStringSet(p)
    names(s) <- sub("\\s.*$", "", names(s))
    cat(sprintf("%-55s %6d seqs\n", basename(p), length(s)))
    s
  })
  
  all_seqs <- do.call(c, parts)
  
  dup_ids <- unique(names(all_seqs)[duplicated(names(all_seqs))])
  cat("\nIDs present in more than one file:", length(dup_ids), "\n")
  
  if (length(dup_ids)) {
    conflict <- vapply(dup_ids, function(id) {
      length(unique(as.character(all_seqs[names(all_seqs) == id]))) > 1
    }, logical(1))
    
    if (any(conflict)) {
      cat("\nSame ID, different sequence:\n")
      print(dup_ids[conflict])
      stop("Hash collision — the runs are not comparable on ASV ID.")
    }
    cat("  all duplicates are byte-identical — safe to dedupe\n")
  }
  
  combined <- all_seqs[!duplicated(names(all_seqs))]
  cat("Unique sequences written:", length(combined), "\n")
  
  writeXStringSet(combined, out_path)
  combined
}

fasta_18S <- build_fasta(FASTA_18S, OUT_18S)

# ------------------------------------------------------------
# Coverage check against the phyloseq object
# ------------------------------------------------------------

check_coverage <- function(fasta, ps, label) {
  ids     <- taxa_names(ps)
  missing <- setdiff(ids, names(fasta))
  
  cat(sprintf("\n%s: %d / %d ASVs covered\n",
              label, length(ids) - length(missing), length(ids)))
  
  if (length(missing)) {
    cat("  still missing:\n")
    print(missing)
  }
  invisible(missing)
}

check_coverage(fasta_18S, euks, "euks (unfiltered)")

# Fail loudly rather than silently dropping sequences later
stopifnot(all(taxa_names(euks) %in% names(fasta_18S)))

# ------------------------------------------------------------
# Reuse downstream
# ------------------------------------------------------------
# In the top-20 script, replace the readDNAStringSet() call with:
#
#   fasta <- readDNAStringSet("labsea_18S_all_ASVs.fasta")
#
# Headers are already stripped, so no sub() needed.

# ------------------------------------------------------------
# Same treatment for 16S, if you want it
# ------------------------------------------------------------
# Your diatom extraction lost 3 of 74 cp16S ASVs for the identical reason.
# Fill in the ce24 16S FASTA path and uncomment:
#
# FASTA_16S <- c(
#   file.path(BASE, "16S-fuhrman/250610-1850.labsea.16S.dna-sequences.fasta"),
#   file.path(BASE, "ce24_addition/16S-fuhrman/<ce24 16S fasta>.fasta")
# )
# fasta_16S <- build_fasta(FASTA_16S, "labsea_16S_all_ASVs.fasta")
# check_coverage(fasta_16S, bact_chloro2, "bact_chloro2")