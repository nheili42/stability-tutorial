# =============================================================================
# FOOD WEB FLUX AND JACOBIAN STABILITY ANALYSIS
# Hengill geothermal streams — stream macroinvertebrate communities
#
# Methods:
#   - Annual and seasonal energy fluxes via fluxweb::fluxing()
#   - Basal resource biomass back-calculated from:
#       B_basal = flux_out / pb_daily_scaled
#     where pb_daily_scaled uses the Junker et al. (2024) temperature scaling:
#       P:B(T) = P:B_ref * 1.07^(T - T_ref)  [~7% increase per °C]
#   - Jacobian matrix (Bazin formulation + self-damping diagonal)
#   - Stability metrics: resilience (max Re eigenvalue), reactivity
#   - Interaction strengths: all off-diagonal Jacobian elements
#   - Sensitivity analysis across P:B low / central / high scenarios
#
# Required objects in environment (from prior data prep):
#   community   : tibble — site, date_id, taxon_id, boot_id,
#                          bio_g_m, metrate_J_day_gram, ae_temp
#   diet        : nested list — diet[[site]][["annual"]]
#                               diet[[site]][["season_means"]][[season]]
#                 (diet matrices, rownames = species order)
#   site_temps  : named numeric vector — site mean annual temperatures (°C)
#                 names must match names(diet)
#                 e.g. c(IS1 = 5.3, IS5 = 9.4, IS6 = 11.2, ...)
#
# Optional (improves seasonal basal back-calculation):
#   site_temps_seasonal : named list of named vectors
#                 e.g. list(first  = c(IS1 = 4.1, ...),
#                           second = c(IS1 = 6.8, ...),
#                           third  = c(IS1 = 5.1, ...))
#
# Junker et al. (2024) Ecology 105(6):e4314 — doi:10.1002/ecy.4314
# =============================================================================



library(fluxweb)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(lubridate)

# =============================================================================
# 1.  GLOBAL CONSTANTS
# =============================================================================

# Basal resource taxa (no measured standing-stock biomass)
BASAL_TAXA <- c(
  "amorphous_detritus", "animal", "cyanobacteria",
  "diatom", "filamentous", "green_algae", "plant_material"
)

# Season assignment: month → season label
# first = Jan–Apr, second = May–Aug, third = Sep–Dec
SEASONS <- c("first", "second", "third")

season_from_month <- function(m) {
  case_when(
    m %in%  1:4  ~ "first",
    m %in%  5:8  ~ "second",
    m %in%  9:12 ~ "third"
  )
}

# Junker et al. (2024) temperature scaling for P:B
# ~7 % increase in biomass turnover per 1 °C (Fig. 1c)
PB_TEMP_COEF <- 1.07
T_REF        <- 5.4   # Reference temperature (°C) — coldest Hengill stream

# =============================================================================
# 2.  P:B LOOKUP TABLE FOR BASAL TAXA
#     Units: per day, at reference temperature T_REF = 5 °C
#     pb_lo / pb_hi used for sensitivity analysis
#
#     Sources:
#       diatom, filamentous, green_algae, cyanobacteria:
#         epilithic algae P:B from Wetzel (2001), Lamberti & Resh (1983),
#         scaled to daily; Demars et al. (2011) for Hengill system context
#       plant_material:
#         stream macrophyte / bryophyte annual P:B ~ 2–5 y⁻¹
#         (Westlake 1982; Huryn & Wallace 2000)
#       amorphous_detritus:
#         microbial decomposition rates for CPOM/FPOM
#         (Webster & Benfield 1986; Tank et al. 2010)
#       animal:
#         meiofaunal / small invertebrate material — treated as slow consumer
#         (conservative; Strayer 1991)
# =============================================================================

BASAL_PB <- tibble(
  taxon_id   = BASAL_TAXA,
  pb_lo      = c(0.002, 0.008, 0.020, 0.050, 0.030, 0.040, 0.003),
  pb_central = c(0.005, 0.015, 0.060, 0.100, 0.070, 0.080, 0.008),
  pb_hi      = c(0.010, 0.025, 0.120, 0.200, 0.150, 0.150, 0.015)
)

# Assimilation efficiencies for basal taxa
# From Junker et al. (2024) Methods / Benke & Wallace (1980, 1997)
BASAL_AE <- tibble(
  taxon_id = BASAL_TAXA,
  ae = c(0.10,  # amorphous_detritus
         0.70,  # animal
         0.10,  # cyanobacteria
         0.30,  # diatom
         0.30,  # filamentous
         0.30,  # green_algae
         0.10)  # plant_material
)

PB_SCENARIOS <- c("lo", "central", "hi")

# =============================================================================
# 3.  HELPER — temperature-scaled P:B
# =============================================================================

#' Scale reference P:B to site temperature using Junker et al. (2024)
#' @param pb_ref   Numeric. Reference daily P:B at T_REF
#' @param temp     Numeric. Site mean temperature (°C)
#' @return Numeric. Temperature-adjusted daily P:B
pb_scaled <- function(pb_ref, temp) {
  pb_ref * PB_TEMP_COEF ^ (temp - T_REF)
}

# =============================================================================
# 4.  HELPER — get site temperature
# =============================================================================

#' Retrieve site temperature for a given site and (optionally) season
#' Falls back to annual mean if seasonal not available
#' @param s    Character. Site name
#' @param seas Character or NULL. Season name
#' @return Numeric. Temperature (°C)
get_temp <- function(s, seas = NULL) {
  if (!is.null(seas) && exists("site_temps_seasonal")) {
    t <- site_temps_seasonal[[seas]][[s]]
    if (!is.null(t) && !is.na(t)) return(t)
  }
  if (exists("site_temps") && !is.null(site_temps[[s]])) {
    return(site_temps[[s]])
  }
  message("  [WARN] No temperature for site '", s, "' season '",
          ifelse(is.null(seas), "annual", seas),
          "' — using T_REF = ", T_REF)
  T_REF
}

# =============================================================================
# 5.  PREP COMMUNITY
#     Aggregate bootstraps for one site × temporal scale.
#     Basal taxa rows are added with placeholder biomass = 1;
#     these will be replaced after the first-pass flux estimate.
# =============================================================================

#' Build aggregated community data frame for one site × scale
#' @param comm_df   Full community tibble
#' @param s         Character. Site name
#' @param seas      Character or NULL. Season (NULL = annual)
#' @param diet_mat  Matrix. Diet matrix for this site × scale
#' @return Tibble or NULL if data insufficient
prep_community <- function(comm_df, s, seas = NULL, diet_mat) {

  sp_order <- rownames(diet_mat)

  # Filter to site (and season if specified)
  df <- comm_df %>%
    mutate(season = season_from_month(month(date_id))) %>%
    filter(site == s)
  if (!is.null(seas)) df <- df %>% filter(season == seas)
  if (nrow(df) == 0) return(NULL)

  # Consumer taxa (those with measured biomass)
  consumer_sp <- setdiff(sp_order, BASAL_TAXA)
  present_sp  <- intersect(consumer_sp, unique(df$taxon_id))
  if (length(present_sp) == 0) {
    message("  [SKIP] No consumer taxa found: site=", s,
            " seas=", ifelse(is.null(seas), "annual", seas))
    return(NULL)
  }

  # Aggregate across sampling dates within each bootstrap
  agg <- df %>%
    filter(taxon_id %in% sp_order) %>%
    group_by(boot_id, taxon_id) %>%
    summarise(
      bio_g_m            = mean(bio_g_m,            na.rm = TRUE),
      metrate_J_day_gram = mean(metrate_J_day_gram,  na.rm = TRUE),
      ae_temp            = mean(ae_temp,             na.rm = TRUE),
      .groups = "drop"
    )

  # Fix zero / NA biomass for consumers (replace with 1% of taxon min)
  agg <- agg %>%
    group_by(taxon_id) %>%
    mutate(
      min_bio  = min(bio_g_m[bio_g_m > 0 & !taxon_id %in% BASAL_TAXA],
                    na.rm = TRUE),
      min_loss = min(metrate_J_day_gram[metrate_J_day_gram > 0 &
                       !taxon_id %in% BASAL_TAXA], na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      bio_g_m = case_when(
        taxon_id %in% BASAL_TAXA        ~ bio_g_m,   # handled separately
        is.na(bio_g_m) | bio_g_m <= 0  ~ min_bio * 0.01,
        TRUE                            ~ bio_g_m
      ),
      metrate_J_day_gram = case_when(
        taxon_id %in% BASAL_TAXA                           ~ metrate_J_day_gram,
        is.na(metrate_J_day_gram) | metrate_J_day_gram <= 0 ~ min_loss * 0.01,
        TRUE                                               ~ metrate_J_day_gram
      )
    ) %>%
    select(-min_bio, -min_loss)

  # Ensure all sp_order present for each boot; add missing as NA rows
  boot_ids    <- unique(agg$boot_id)
  full_grid   <- expand_grid(boot_id = boot_ids, taxon_id = sp_order)
  agg         <- full_grid %>%
    left_join(agg, by = c("boot_id", "taxon_id"))

  # Add basal taxa rows with:
  #   - placeholder biomass = 1 (replaced post-flux)
  #   - metabolic loss = 0 (basal resources have no self-respiration in fluxweb)
  #   - ae from BASAL_AE lookup
  agg <- agg %>%
    left_join(BASAL_AE %>% rename(ae_basal = ae), by = "taxon_id") %>%
    mutate(
      bio_g_m = case_when(
        taxon_id %in% BASAL_TAXA          ~ 1.0,
        is.na(bio_g_m)                    ~ 1e-6,   # safety fallback
        TRUE                              ~ bio_g_m
      ),
      metrate_J_day_gram = case_when(
        taxon_id %in% BASAL_TAXA          ~ 0.0,
        is.na(metrate_J_day_gram)         ~ 1e-6,
        TRUE                              ~ metrate_J_day_gram
      ),
      ae_temp = case_when(
        taxon_id %in% BASAL_TAXA          ~ ae_basal,
        is.na(ae_temp)                    ~ 0.3,
        TRUE                              ~ ae_temp
      )
    ) %>%
    select(-ae_basal)

  agg
}

# =============================================================================
# 6.  RUN FLUXES
#     Call fluxweb::fluxing() for every bootstrap replicate.
# =============================================================================

#' Estimate flux matrices for all bootstrap replicates
#' @param agg      Aggregated community tibble (from prep_community)
#' @param diet_mat Diet matrix
#' @return Named list of flux matrices (one per boot_id); NULL entries dropped
run_fluxes <- function(agg, diet_mat) {

  sp_order <- rownames(diet_mat)
  boot_ids <- sort(unique(agg$boot_id))

  flux_list <- lapply(boot_ids, function(b) {
    df_b <- agg %>%
      filter(boot_id == b) %>%
      arrange(match(taxon_id, sp_order))

    # Guard: need exactly sp_order rows in correct order
    if (nrow(df_b) != length(sp_order) ||
        !identical(df_b$taxon_id, sp_order)) {
      return(NULL)
    }

    tryCatch(
      fluxing(
        mat          = diet_mat,
        biomasses    = df_b$bio_g_m,
        losses       = df_b$metrate_J_day_gram,
        efficiencies = df_b$ae_temp,
        bioms.prefs  = FALSE,
        bioms.losses = TRUE,
        ef.level     = "prey"
      ),
      error   = function(e) NULL,
      warning = function(w) {
        tryCatch(
          fluxing(
            mat          = diet_mat,
            biomasses    = df_b$bio_g_m,
            losses       = df_b$metrate_J_day_gram,
            efficiencies = df_b$ae_temp,
            bioms.prefs  = FALSE,
            bioms.losses = TRUE,
            ef.level     = "prey"
          ),
          error = function(e2) NULL
        )
      }
    )
  })

  names(flux_list) <- boot_ids
  Filter(Negate(is.null), flux_list)
}

# =============================================================================
# 7.  BACK-CALCULATE BASAL BIOMASS
#     B_basal[i] = sum(flux_mat[i, ]) / pb_daily_temp_scaled[i]
# =============================================================================

#' Replace placeholder basal biomasses with back-calculated estimates
#' @param agg         Community tibble with placeholder basal biomass = 1
#' @param flux_list   Named list of flux matrices from run_fluxes()
#' @param diet_mat    Diet matrix
#' @param temp        Site mean temperature (°C)
#' @param pb_scenario Character. One of "lo", "central", "hi"
#' @return Updated community tibble
back_calc_basal <- function(agg, flux_list, diet_mat, temp,
                             pb_scenario = "central") {

  sp_order <- rownames(diet_mat)
  pb_col   <- paste0("pb_", pb_scenario)

  map_dfr(names(flux_list), function(b) {

    df_b     <- agg %>%
      filter(boot_id == b) %>%
      arrange(match(taxon_id, sp_order))

    flux_mat <- flux_list[[b]]
    rownames(flux_mat) <- colnames(flux_mat) <- sp_order

    df_b %>%
      left_join(BASAL_PB %>% select(taxon_id, !!pb_col := !!sym(pb_col)),
                by = "taxon_id") %>%
      mutate(
        bio_g_m = case_when(
          taxon_id %in% BASAL_TAXA ~ {
            # Total flux leaving this basal taxon to all consumers
            flux_out   <- rowSums(flux_mat, na.rm = TRUE)[taxon_id]
            pb_ref     <- .data[[pb_col]]
            pb_adj     <- pb_scaled(pb_ref, temp)
            # If no flux out (taxon not consumed), use a small positive biomass
            ifelse(flux_out > 0 & !is.na(pb_adj) & pb_adj > 0,
                   flux_out / pb_adj,
                   1e-4)
          },
          TRUE ~ bio_g_m
        )
      ) %>%
      select(-any_of(pb_col))
  })
}

# =============================================================================
# 8.  FLUX → JACOBIAN
#
#   Off-diagonal (Bazin et al. formulation):
#     J[i,j] = (ae[i] * flux[j→i]  -  flux[i→j]) / B[j]
#            = gain to i from eating j, net of loss of i to j, per unit B[j]
#
#   Note on indexing convention used here:
#     flux_mat[prey, predator]:  rows = prey, cols = predator
#     J[recipient, donor] convention for community matrix
#
#   Diagonal (self-damping):
#     J[i,i] = -(metabolic_loss_rate[i] + total_flux_out[i] / B[i])
#     For basal taxa metrate = 0, so: J[i,i] = -total_flux_out[i] / B[i]
# =============================================================================

#' Convert flux matrix to Jacobian community matrix
#' @param flux_mat   Matrix [prey × predator]. Flux rates (J d⁻¹ m⁻²)
#' @param biomass    Named numeric vector. Biomass (g m⁻²)
#' @param metrate    Named numeric vector. Metabolic loss rate (J d⁻¹ g⁻¹)
#' @param ae         Named numeric vector. Assimilation efficiency
#' @return Square Jacobian matrix
flux_to_jacobian <- function(flux_mat, biomass, metrate, ae) {

  n        <- nrow(flux_mat)
  sp_order <- names(biomass)
  J        <- matrix(0, n, n, dimnames = list(sp_order, sp_order))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      sp_i <- sp_order[i]
      sp_j <- sp_order[j]
      Bj   <- biomass[sp_j]
      if (Bj <= 0) next

      # gain to i from eating j (i is predator, j is prey)
      gain_i_from_j <- ae[sp_i] * flux_mat[sp_j, sp_i]

      # loss of i due to being eaten by j (j is predator, i is prey)
      loss_i_to_j   <- flux_mat[sp_i, sp_j]

      J[i, j] <- (gain_i_from_j - loss_i_to_j) / Bj
    }
  }

  # Self-damping diagonal
  for (i in seq_len(n)) {
    sp_i <- sp_order[i]
    Bi   <- biomass[sp_i]
    if (Bi <= 0) next
    flux_out_rate  <- sum(flux_mat[sp_i, ], na.rm = TRUE) / Bi
    J[i, i] <- -(metrate[sp_i] + flux_out_rate)
  }

  J
}

# =============================================================================
# 9.  STABILITY METRICS
# =============================================================================

#' Compute resilience and reactivity from Jacobian
#' Resilience  = max(Re(eigenvalues(J)))    — more negative = more stable
#' Reactivity  = max(eigenvalues((J+Jᵀ)/2)) — negative = stable to transients
#' @param J  Square numeric matrix
#' @return Tibble: resilience, reactivity, is_stable
jacobian_stability <- function(J) {

  ev_real    <- Re(eigen(J, only.values = TRUE)$values)
  resilience <- max(ev_real)

  J_sym      <- (J + t(J)) / 2
  reactivity <- max(Re(eigen(J_sym, only.values = TRUE)$values))

  tibble(
    resilience = resilience,
    reactivity = reactivity,
    is_stable  = resilience < 0
  )
}

# =============================================================================
# 10. INTERACTION STRENGTHS
#     All off-diagonal Jacobian elements as tidy tibble
# =============================================================================

#' Extract off-diagonal Jacobian elements as interaction strengths
#' @param J  Jacobian matrix with named rows/cols
#' @return Tibble: recipient, donor, interaction_strength, link_type
extract_interactions <- function(J) {

  sp_order <- rownames(J)

  expand_grid(recipient = sp_order, donor = sp_order) %>%
    filter(recipient != donor) %>%
    mutate(
      interaction_strength = map2_dbl(
        recipient, donor, ~ J[.x, .y]
      ),
      recipient_is_basal = recipient %in% BASAL_TAXA,
      donor_is_basal     = donor     %in% BASAL_TAXA,
      link_type = case_when(
        !donor_is_basal & !recipient_is_basal ~ "consumer_consumer",
        donor_is_basal  & !recipient_is_basal ~ "basal_to_consumer",
        !donor_is_basal &  recipient_is_basal ~ "consumer_to_basal",
        TRUE                                   ~ "basal_basal"
      ),
      direction = case_when(
        interaction_strength > 0 ~ "positive",
        interaction_strength < 0 ~ "negative",
        TRUE                     ~ "zero"
      )
    ) %>%
    select(-recipient_is_basal, -donor_is_basal) %>%
    filter(interaction_strength != 0)
}

# =============================================================================
# 11. CORE PIPELINE — one site × temporal scale × P:B scenario
# =============================================================================

#' Full pipeline for one site × scale combination
#'
#' Two-pass approach:
#'   Pass 1: fluxes with placeholder basal biomass (B = 1)
#'   Back-calc: replace basal biomass from flux-out / temp-scaled P:B
#'   Pass 2: recompute fluxes with corrected basal biomass
#'   Jacobian + stability + interaction strengths on Pass 2 fluxes
#'
#' @param agg         Community tibble from prep_community()
#' @param diet_mat    Diet matrix
#' @param temp        Site temperature (°C)
#' @param pb_scenario "lo", "central", or "hi"
#' @return List: stability_df, interactions_df, basal_biomass_df, flux_list
run_pipeline <- function(agg, diet_mat, temp, pb_scenario = "central") {

  sp_order <- rownames(diet_mat)

  # ---- Pass 1: fluxes with B_basal = 1 ----
  fl_pass1 <- run_fluxes(agg, diet_mat)
  if (length(fl_pass1) == 0) {
    message("    [FAIL] No valid flux estimates in pass 1")
    return(NULL)
  }

  # ---- Back-calculate basal biomass ----
  agg_updated <- back_calc_basal(agg, fl_pass1, diet_mat, temp, pb_scenario)

  # ---- Pass 2: recompute fluxes with corrected basal biomass ----
  fl_pass2 <- run_fluxes(agg_updated, diet_mat)
  if (length(fl_pass2) == 0) {
    message("    [FAIL] No valid flux estimates in pass 2")
    return(NULL)
  }

  # ---- Per-bootstrap Jacobian, stability, interactions ----
  boot_results <- imap(fl_pass2, function(flux_mat, b) {

    rownames(flux_mat) <- colnames(flux_mat) <- sp_order

    df_b <- agg_updated %>%
      filter(boot_id == b) %>%
      arrange(match(taxon_id, sp_order))

    if (nrow(df_b) != length(sp_order)) return(NULL)

    biomass <- setNames(df_b$bio_g_m,            sp_order)
    metrate <- setNames(df_b$metrate_J_day_gram,  sp_order)
    ae      <- setNames(df_b$ae_temp,             sp_order)

    J <- flux_to_jacobian(flux_mat, biomass, metrate, ae)

    list(
      stability    = jacobian_stability(J)   %>% mutate(boot_id = b),
      interactions = extract_interactions(J) %>% mutate(boot_id = b),
      basal_bio    = df_b %>%
        filter(taxon_id %in% BASAL_TAXA) %>%
        select(boot_id, taxon_id, bio_g_m)
    )
  })

  boot_results <- Filter(Negate(is.null), boot_results)
  if (length(boot_results) == 0) return(NULL)

  list(
    stability_df     = map_dfr(boot_results, "stability"),
    interactions_df  = map_dfr(boot_results, "interactions"),
    basal_biomass_df = map_dfr(boot_results, "basal_bio"),
    flux_list        = fl_pass2,
    n_boots          = length(boot_results)
  )
}

# =============================================================================
# 12. MASTER LOOP
#     Iterate: sites × {annual, seasonal} × P:B scenarios
# =============================================================================

sites <- names(diet)

message("=== Starting master loop ===")
message("Sites: ", paste(sites, collapse = ", "))
message("P:B scenarios: ", paste(PB_SCENARIOS, collapse = ", "))

# Storage: nested list [pb_scenario][[site]][[scale]]
# scale = "annual" or one of SEASONS
all_results <- setNames(
  lapply(PB_SCENARIOS, function(pb_sc) {
    setNames(vector("list", length(sites)), sites)
  }),
  PB_SCENARIOS
)

for (pb_sc in PB_SCENARIOS) {
  message("\n--- P:B scenario: ", pb_sc, " ---")

  for (s in sites) {
    message("  Site: ", s)
    all_results[[pb_sc]][[s]] <- list()

    # ---- Annual ----
    diet_ann <- diet[[s]][["annual"]]
    if (!is.null(diet_ann)) {
      message("    Annual...")
      temp_ann <- get_temp(s, seas = NULL)
      comm_ann <- prep_community(community, s, seas = NULL, diet_ann)

      if (!is.null(comm_ann)) {
        res <- run_pipeline(comm_ann, diet_ann, temp_ann, pb_sc)
        if (!is.null(res)) {
          message("      OK — ", res$n_boots, " bootstraps")
          all_results[[pb_sc]][[s]][["annual"]] <- res
        }
      }
    }

    # ---- Seasonal ----
    for (seas in SEASONS) {
      message("    Season: ", seas, "...")
      diet_seas <- diet[[s]][["season_means"]][[seas]]
      if (is.null(diet_seas)) next

      temp_seas <- get_temp(s, seas = seas)
      comm_seas <- prep_community(community, s, seas = seas, diet_seas)

      if (!is.null(comm_seas)) {
        res <- run_pipeline(comm_seas, diet_seas, temp_seas, pb_sc)
        if (!is.null(res)) {
          message("      OK — ", res$n_boots, " bootstraps")
          all_results[[pb_sc]][[s]][[seas]] <- res
        }
      }
    }
  }
}

message("\n=== Master loop complete ===\n")

# =============================================================================
# 13. FLATTEN TO ANALYSIS-READY TIBBLES
# =============================================================================

#' Helper: pull one element from all_results into a flat tibble
flatten_element <- function(element_name) {
  map_dfr(PB_SCENARIOS, function(pb_sc) {
    map_dfr(sites, function(s) {
      scales <- names(all_results[[pb_sc]][[s]])
      map_dfr(scales, function(sc) {
        x <- all_results[[pb_sc]][[s]][[sc]][[element_name]]
        if (is.null(x)) return(NULL)
        x %>% mutate(
          pb_scenario = pb_sc,
          site        = s,
          scale       = ifelse(sc == "annual", "annual", "seasonal"),
          season      = ifelse(sc == "annual", NA_character_, sc)
        )
      })
    })
  })
}

stability_raw    <- flatten_element("stability_df")
interactions_raw <- flatten_element("interactions_df")
basal_bio_raw    <- flatten_element("basal_biomass_df")

# =============================================================================
# 14. SUMMARY TABLES
# =============================================================================

# ---- 14a. Stability summary: mean ± 95 % CI per site × scale × season × scenario ----
stability_summary <- stability_raw %>%
  group_by(pb_scenario, site, scale, season) %>%
  summarise(
    n_boots         = n(),
    resilience_mean = mean(resilience,            na.rm = TRUE),
    resilience_lo95 = quantile(resilience, 0.025, na.rm = TRUE),
    resilience_hi95 = quantile(resilience, 0.975, na.rm = TRUE),
    reactivity_mean = mean(reactivity,            na.rm = TRUE),
    reactivity_lo95 = quantile(reactivity, 0.025, na.rm = TRUE),
    reactivity_hi95 = quantile(reactivity, 0.975, na.rm = TRUE),
    prop_stable     = mean(is_stable,             na.rm = TRUE),
    .groups = "drop"
  )

# ---- 14b. Sensitivity: range of stability metrics across P:B scenarios ----
stability_sensitivity <- stability_summary %>%
  select(pb_scenario, site, scale, season,
         resilience_mean, reactivity_mean, prop_stable) %>%
  pivot_wider(
    names_from  = pb_scenario,
    values_from = c(resilience_mean, reactivity_mean, prop_stable)
  ) %>%
  mutate(
    resilience_range = resilience_mean_hi - resilience_mean_lo,
    reactivity_range = reactivity_mean_hi - reactivity_mean_lo
  )

# ---- 14c. Interaction strength summary (central scenario, per link_type) ----
interaction_summary <- interactions_raw %>%
  filter(pb_scenario == "central") %>%
  group_by(site, scale, season, link_type, direction) %>%
  summarise(
    n_links           = n(),
    mean_strength     = mean(interaction_strength,           na.rm = TRUE),
    median_strength   = median(interaction_strength,         na.rm = TRUE),
    sd_strength       = sd(interaction_strength,             na.rm = TRUE),
    abs_mean_strength = mean(abs(interaction_strength),      na.rm = TRUE),
    lo95_strength     = quantile(interaction_strength, 0.025, na.rm = TRUE),
    hi95_strength     = quantile(interaction_strength, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 14d. Per-link interaction strengths across all scenarios (sensitivity) ----
interaction_sensitivity <- interactions_raw %>%
  group_by(pb_scenario, site, scale, season, recipient, donor, link_type) %>%
  summarise(
    median_strength = median(interaction_strength, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = pb_scenario,
    values_from = median_strength,
    names_prefix = "is_"
  ) %>%
  mutate(
    is_range = abs(is_hi - is_lo)
  )

# ---- 14e. Basal biomass summary (back-calculated) ----
basal_bio_summary <- basal_bio_raw %>%
  group_by(pb_scenario, site, scale, season, taxon_id) %>%
  summarise(
    bio_mean = mean(bio_g_m,            na.rm = TRUE),
    bio_lo95 = quantile(bio_g_m, 0.025, na.rm = TRUE),
    bio_hi95 = quantile(bio_g_m, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# 15. OUTPUT SUMMARY
# =============================================================================

message("Results overview:")
message("  stability_raw rows:       ", nrow(stability_raw))
message("  interactions_raw rows:    ", nrow(interactions_raw))
message("  basal_bio_raw rows:       ", nrow(basal_bio_raw))
message("  stability_summary rows:   ", nrow(stability_summary))
message("  stability_sensitivity:    ", nrow(stability_sensitivity))
message("  interaction_summary rows: ", nrow(interaction_summary))
message("  interaction_sensitivity:  ", nrow(interaction_sensitivity))
message("  basal_bio_summary rows:   ", nrow(basal_bio_summary))

# =============================================================================
# 16. EXAMPLE QUERIES
# =============================================================================

# --- Stability for one site across seasons (central scenario) ---
# stability_summary %>%
#   filter(site == "IS6", pb_scenario == "central") %>%
#   arrange(scale, season)

# --- Which site × season is most sensitive to P:B choice? ---
# stability_sensitivity %>%
#   arrange(desc(resilience_range)) %>%
#   head(10)

# --- Strongest basal→consumer interaction strengths (central, annual) ---
# interaction_summary %>%
#   filter(pb_scenario == "central", scale == "annual",
#          link_type == "basal_to_consumer") %>%
#   arrange(desc(abs_mean_strength))

# --- Back-calculated basal biomass by taxon across sites ---
# basal_bio_summary %>%
#   filter(pb_scenario == "central", scale == "annual") %>%
#   ggplot(aes(x = site, y = bio_mean, fill = taxon_id)) +
#   geom_col(position = "stack") +
#   geom_errorbar(aes(ymin = bio_lo95, ymax = bio_hi95),
#                 position = position_stack(), width = 0.2)

# --- Compare resilience across seasons within sites ---
# stability_summary %>%
#   filter(pb_scenario == "central", scale == "seasonal") %>%
#   ggplot(aes(x = season, y = resilience_mean,
#              ymin = resilience_lo95, ymax = resilience_hi95,
#              colour = site, group = site)) +
#   geom_point(size = 2) +
#   geom_line() +
#   geom_errorbar(width = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(y = "Resilience (max Re eigenvalue)", x = "Season")
