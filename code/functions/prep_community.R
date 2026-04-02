# This function prepares the community df for fluxing()


# =============================================================================
# FUNCTION 1: prep_community()
# Filter community data by site x season, aggregate boots, fix zero biomass
# =============================================================================

prep_community <- function(community_df, site_name, season_name = NULL, diet_mat) {
  
  sp_order <- rownames(diet_mat)
  
  df <- community_df %>% filter(site == site_name)
  if (!is.null(season_name)) df <- df %>% filter(season == season_name)
  if (nrow(df) == 0) return(NULL)
  
  if (!all(sp_order %in% unique(df$taxon_id))) {
    missing <- setdiff(sp_order, unique(df$taxon_id))
    warning(paste("Missing taxa for", site_name,
                  ifelse(is.null(season_name), "annual", season_name), ":",
                  paste(missing, collapse = ", ")))
    return(NULL)
  }
  
  agg <- df %>%
    group_by(boot_id, taxon_id) %>%
    summarise(
      bio_g_m            = mean(bio_g_m,            na.rm = TRUE),
      metrate_J_day_gram = mean(metrate_J_day_gram,  na.rm = TRUE),
      ae_temp            = mean(ae_temp,             na.rm = TRUE),
      .groups = "drop"
    )
  
  # global fallback so always-zero taxa are also covered ---
  global_min_bio  <- agg %>% filter(bio_g_m > 0)            %>% pull(bio_g_m)            %>% min(na.rm = TRUE)
  global_min_loss <- agg %>% filter(metrate_J_day_gram > 0) %>% pull(metrate_J_day_gram) %>% min(na.rm = TRUE)
  
  # Per-taxon minimum where available, global minimum as backstop
  placeholders <- agg %>%
    filter(bio_g_m > 0 & metrate_J_day_gram > 0) %>%
    group_by(taxon_id) %>%
    summarise(
      ph_bio  = min(bio_g_m,            na.rm = TRUE) * 0.01,
      ph_loss = min(metrate_J_day_gram,  na.rm = TRUE) * 0.01,
      .groups = "drop"
    )
  
  # Full species list so always-zero taxa get the global fallback
  all_taxa <- tibble(taxon_id = sp_order) %>%
    left_join(placeholders, by = "taxon_id") %>%
    mutate(
      ph_bio  = coalesce(ph_bio,  global_min_bio  * 0.01),
      ph_loss = coalesce(ph_loss, global_min_loss * 0.01)
    )
  
  agg %>%
    left_join(all_taxa, by = "taxon_id") %>%
    mutate(
      bio_g_m            = case_when(
        is.na(bio_g_m)  | bio_g_m  == 0 ~ ph_bio,
        TRUE                             ~ bio_g_m),
      metrate_J_day_gram = case_when(
        is.na(metrate_J_day_gram) | metrate_J_day_gram == 0 ~ ph_loss,
        TRUE                                                  ~ metrate_J_day_gram)
    ) %>%
    select(-ph_bio, -ph_loss)
}
