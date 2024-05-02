# Data load -------------------------------------------------------------------#
radboud_ps_abs <- readRDS(paste0(data_00$outdir, "/01_absolute_phyloseq.rds"))

# Data wrangling --------------------------------------------------------------#
# In paired samples, should contain an unique sample ID column, then str_extract should be removed!
radboud_ps_abs_bac <- radboud_ps_abs %>% 
  subset_taxa(Domain == "Bacteria") %>% 
  removeZeros()

# normalize by bacterial domain
radboud_ps_rel_bac_norm <- radboud_ps_abs_bac %>% 
  transform_sample_counts(function(x) x / sum(x))


# Save wrangled file ----------------------------------------------------------#
saveRDS(radboud_ps_rel_bac_norm, paste0(data_00$outdir, "02_ps_rel_bac_norm.rds"))
saveRDS(radboud_ps_abs_bac, paste0(data_00$outdir, "/01_ps_abs_bac.rds"))
