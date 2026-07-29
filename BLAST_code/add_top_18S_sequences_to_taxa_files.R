library(readr)
library(dplyr)

# ============================================================
# Update PR2 taxonomy with BLAST results for the top-20 euk ASVs
# Writes new files (originals untouched): *.csv -> *_blast.csv
# ============================================================

TAXA1_PATH <- "../18S-fuhrman/taxonomy/taxa_all_18S_PR2.csv"
TAXA2_PATH <- "../ce24_addition/18S-fuhrman/taxonomy/taxa_all_18S_ce24_PR2.csv"

rank_cols   <- c("Domain","Supergroup","Phylum","Class","Order","Family","Genus","Species")
UPDATE_CONF <- TRUE   # set FALSE to keep the original classifier confidence values

# ------------------------------------------------------------
# 1. Load
# ------------------------------------------------------------
new_taxa <- read_csv("euk18S_top20_compiled_BLAST.csv", show_col_types = FALSE) %>%
  select(-any_of(c("...18", "OTU")))          # empty trailing col + duplicate of ASV

taxa1 <- read_csv(TAXA1_PATH, show_col_types = FALSE)
taxa2 <- read_csv(TAXA2_PATH, show_col_types = FALSE)

stopifnot(!any(duplicated(new_taxa$ASV)))

# ------------------------------------------------------------
# 2. Update function
# ------------------------------------------------------------
update_tax <- function(tax, new_taxa, id_col = "#OTUID",
                       cols = rank_cols, update_conf = UPDATE_CONF) {
  
  keep <- c("ASV", cols, if (update_conf) "confidence", "blast_notes")
  
  new <- new_taxa %>%
    select(all_of(keep)) %>%
    rename_with(~ paste0(.x, "_new"), -ASV)
  
  out <- tax %>% left_join(new, by = setNames("ASV", id_col))
  
  for (cl in setdiff(keep, "ASV")) {
    if (cl %in% names(out)) {
      out[[cl]] <- coalesce(out[[paste0(cl, "_new")]], out[[cl]])
    } else {
      out[[cl]] <- out[[paste0(cl, "_new")]]   # blast_notes: new column
    }
  }
  
  out %>% select(-ends_with("_new"))
}

taxa1_upd <- update_tax(taxa1, new_taxa)
taxa2_upd <- update_tax(taxa2, new_taxa)

# ------------------------------------------------------------
# 3. Checks
# ------------------------------------------------------------
cat("taxa1: ", sum(taxa1$`#OTUID` %in% new_taxa$ASV), "of", nrow(taxa1), "rows updated\n")
cat("taxa2: ", sum(taxa2$`#OTUID` %in% new_taxa$ASV), "of", nrow(taxa2), "rows updated\n")
cat("ASVs matching neither file:",
    length(setdiff(new_taxa$ASV, c(taxa1$`#OTUID`, taxa2$`#OTUID`))), "\n")

stopifnot(nrow(taxa1_upd) == nrow(taxa1),
          nrow(taxa2_upd) == nrow(taxa2))

# spot-check a few updated rows
taxa1_upd %>%
  filter(`#OTUID` %in% new_taxa$ASV) %>%
  select(`#OTUID`, Class, Family, Genus, Species, blast_notes) %>%
  head(10) %>%
  print(width = Inf)

# ------------------------------------------------------------
# 4. Write
# ------------------------------------------------------------
OUT1 <- sub("\\.csv$", "_blast.csv", TAXA1_PATH)
OUT2 <- sub("\\.csv$", "_blast.csv", TAXA2_PATH)

write_csv(taxa1_upd, OUT1, na = "")
write_csv(taxa2_upd, OUT2, na = "")

cat("\nWritten:\n ", OUT1, "\n ", OUT2, "\n")
