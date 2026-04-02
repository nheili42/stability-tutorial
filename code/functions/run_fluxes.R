# This function loops across biomass bootstrap to estimate fluxes


# =============================================================================
# FUNCTION 2: run_fluxes()
# Run fluxing() across all bootstraps for one site x temporal scale
# =============================================================================

run_fluxes <- function(community_agg, diet_mat) {
  
  sp_order <- rownames(diet_mat)
  boot_ids <- unique(community_agg$boot_id)
  
  flux_list <- lapply(boot_ids, function(b) {
    
    df_b <- community_agg %>%
      filter(boot_id == b) %>%
      arrange(match(taxon_id, sp_order))
    
    if (nrow(df_b) != nrow(diet_mat)) {
      warning(paste("Boot", b, "row mismatch — skipping"))
      return(NULL)
    }
    if (!all(df_b$taxon_id == sp_order)) {
      warning(paste("Boot", b, "order mismatch — skipping"))
      return(NULL)
    }
    
    tryCatch(
      fluxing(
        mat          = diet_mat,
        biomasses    = df_b$bio_g_m,
        losses       = df_b$metrate_J_day_gram,
        efficiencies = df_b$ae_temp,
        bioms.prefs  = FALSE,  # GCA diet proportions used as-is
        bioms.losses = TRUE,   # losses are per-gram, scaled internally
        ef.level     = "prey"  # ae is a property of the prey
      ),
      error = function(e) { warning(paste("Boot", b, ":", e$message)); NULL }
    )
  })
  
  names(flux_list) <- boot_ids
  Filter(Negate(is.null), flux_list)
}
