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
