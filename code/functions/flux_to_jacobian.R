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
