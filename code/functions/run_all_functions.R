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


# =============================================================================
# FUNCTION 3: flux_to_jacobian()
# Bazin et al. create.jacob formulation WITH self-damping diagonal
#
# Off-diagonal: identical to Bazin —
#   jacob = t(flux_mat) * ae - flux_mat
#   jacob = sweep(jacob, MARGIN=2, biomasses, '/')
#
# Diagonal (added here, absent in Bazin):
#   J[i,i] = -(total metabolic loss + total predation loss on i) / biomass[i]
#   = -(loss_per_gram[i] * biomass[i] + sum(flux_mat[i,])) / biomass[i]
#   = -loss_per_gram[i] - sum(flux_mat[i,]) / biomass[i]
# =============================================================================

flux_to_jacobian <- function(flux_mat, biomass_vec, loss_vec, ae_vec) {
  
  n        <- nrow(flux_mat)
  sp_order <- rownames(flux_mat)
  
  # Guard: check for zero/NA biomass before proceeding
  if (any(is.na(biomass_vec)) || any(biomass_vec == 0)) {
    problem_sp <- sp_order[is.na(biomass_vec) | biomass_vec == 0]
    warning(paste("Zero/NA biomass for:", paste(problem_sp, collapse = ", "),
                  "— Jacobian will have Inf. Returning NULL."))
    return(NULL)
  }
  
  # Off-diagonal: Bazin formulation
  # t(flux_mat) * ae_vec: gains (prey j -> predator i, scaled by ae of prey j)
  # - flux_mat: losses (flux out of i to all consumers)
  jacob <- t(flux_mat) * ae_vec - flux_mat
  
  # Divide each column j by biomass of species j
  jacob <- sweep(jacob, MARGIN = 2, biomass_vec, "/")
  
  # Self-damping diagonal (departure from Bazin)
  # Total loss on species i = metabolic loss + predation loss
  # Both divided by biomass[i] to get per-capita rate
  for (i in seq_len(n)) {
    if (biomass_vec[i] > 0) {
      metabolic_loss  <- loss_vec[i]                           # already per gram
      predation_loss  <- sum(flux_mat[i, ]) / biomass_vec[i]  # flux out / biomass
      jacob[i, i]     <- -(metabolic_loss + predation_loss)
    } else {
      jacob[i, i] <- 0
    }
  }
  
  rownames(jacob) <- colnames(jacob) <- sp_order
  jacob
}


# =============================================================================
# FUNCTION 4: jacobian_stability()
# Following the approach of Bazin et al. stability metrics
# Resilience  = max(Re(eigenvalues(J)))       positive = unstable
# Reactivity  = max(eigenvalues(H))           H = (J + J^T) / 2
# =============================================================================

jacobian_stability <- function(J) {
  
  eigs        <- eigen(J, only.values = TRUE)$values
  resilience  <- max(Re(eigs))
  
  J_hermitian <- (J + t(Conj(J))) / 2
  reactivity  <- max(Re(eigen(J_hermitian, only.values = TRUE)$values))
  
  tibble(
    resilience = resilience,    # λJmax: more positive = more unstable (long-term)
    reactivity = reactivity,    # λHmax: more positive = more reactive (short-term)
    is_stable  = resilience < 0
  )
}


# =============================================================================
# FUNCTION 5: extract_interaction_strengths()
# Off-diagonal Jacobian elements as tidy tibble
# =============================================================================

extract_interaction_strengths <- function(J) {
  
  sp_order <- rownames(J)
  
  expand.grid(recipient = sp_order, actor = sp_order,
              stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    filter(recipient != actor) %>%
    mutate(
      interaction_strength = map2_dbl(recipient, actor, ~ J[.x, .y]),
      direction = case_when(
        interaction_strength > 0 ~ "positive",
        interaction_strength < 0 ~ "negative",
        TRUE                     ~ "zero"
      )
    ) %>%
    filter(interaction_strength != 0)
}

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
