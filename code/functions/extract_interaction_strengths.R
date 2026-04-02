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
