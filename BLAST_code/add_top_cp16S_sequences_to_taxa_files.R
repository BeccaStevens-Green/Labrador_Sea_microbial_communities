new_taxa <- read_csv("cp16S_top20_compiled_BLAST.csv")
taxa1 <- read_csv("../16S-fuhrman/taxonomy/taxa_bact_chloro_16S_blast.csv")
taxa2 <- read_csv("../ce24_addition/16S-fuhrman/taxonomy/taxa_bact_chloro_16S_ce24_blast.csv")

str(new_taxa)
str(taxa1)
str(taxa2)

library(dplyr)

rank_cols <- c("Domain","Phylum","Class","Order","Family","Genus","Species","Species2_phytoref")

taxa1 %>%
  filter(`#OTUID` %in% new_taxa$ASV) %>%
  select(all_of(rank_cols)) %>%
  head(10)

update_tax <- function(tax, id_col, new_taxa, cols = rank_cols,
                       update_conf = TRUE, update_notes = TRUE) {
  new <- new_taxa %>%
    select(ASV, all_of(cols), confidence, blast_notes) %>%
    rename_with(~ paste0(.x, "_new"), -ASV)
  
  out <- tax %>% left_join(new, by = setNames("ASV", id_col))
  
  for (cl in cols) out[[cl]] <- coalesce(out[[paste0(cl, "_new")]], out[[cl]])
  if (update_conf)  out$confidence  <- coalesce(out$confidence_new,  out$confidence)
  if (update_notes) out$blast_notes <- coalesce(out$blast_notes_new, out$blast_notes)
  
  out %>% select(-ends_with("_new"))
}

taxa1_upd <- update_tax(taxa1, "#OTUID",     new_taxa)
taxa2_upd <- update_tax(taxa2, "feature_id", new_taxa)

sum(taxa1$`#OTUID`   %in% new_taxa$ASV)   # rows hit in taxa1
sum(taxa2$feature_id %in% new_taxa$ASV)   # rows hit in taxa2
setdiff(new_taxa$ASV, c(taxa1$`#OTUID`, taxa2$feature_id))  # should be character(0)

write_csv(taxa1_upd, "../16S-fuhrman/taxonomy/taxa_bact_chloro_16S_blast.csv", na = "")
write_csv(taxa2_upd, "../ce24_addition/16S-fuhrman/taxonomy/taxa_bact_chloro_16S_ce24_blast.csv", na = "")
