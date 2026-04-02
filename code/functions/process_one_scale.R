# =============================================================================
# FUNCTION 6: process_one_scale()
# Runs fluxes + Jacobian + stability + interaction strengths
# for one site x temporal scale combination
# =============================================================================

process_one_scale <- function(community_agg, diet_mat) {
  
  sp_order  <- rownames(diet_mat)
  flux_list <- run_fluxes(community_agg, diet_mat)
  if (length(flux_list) == 0) return(NULL)
  
  boot_results <- imap(flux_list, function(flux_mat, b) {
    
    df_b <- community_agg %>%
      filter(boot_id == b) %>%
      arrange(match(taxon_id, sp_order))
    
    J <- flux_to_jacobian(
      flux_mat    = flux_mat,
      biomass_vec = df_b$bio_g_m,
      loss_vec    = df_b$metrate_J_day_gram,
      ae_vec      = df_b$ae_temp
    )
    
    # Skip this boot if Jacobian failed
    if (is.null(J)) return(NULL)
    
    # Check for non-finite values even if no zero biomass was caught
    if (any(!is.finite(J))) {
      warning(paste("Non-finite values in Jacobian for boot", b, "— skipping"))
      return(NULL)
    }
    
    list(
      stability             = jacobian_stability(J)            %>% mutate(boot_id = b),
      interaction_strengths = extract_interaction_strengths(J) %>% mutate(boot_id = b),
      jacobian              = J
    )
  })
  
  # Remove failed boots
  boot_results <- Filter(Negate(is.null), boot_results)
  if (length(boot_results) == 0) return(NULL)
  
  list(
    flux_list             = flux_list,
    stability_summary     = map_dfr(boot_results, "stability"),
    interaction_strengths = map_dfr(boot_results, "interaction_strengths"),
    jacobians             = map(boot_results, "jacobian")
  )
}
