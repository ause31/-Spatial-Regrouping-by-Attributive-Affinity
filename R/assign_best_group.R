#' Rattacher chaque commune candidate au groupe voisin le plus proche
#'
#' @param data Un objet sf avec colonne "candidate" et "id_row"
#' @param group_var Variable de groupe
#' @param vars_attr Variables attributaires
#' @param profiles Profils moyens par groupe (utilises en fallback si local_profile = FALSE)
#' @param nb Matrice de voisinage (spdep) — obligatoire si local_profile = TRUE
#' @param local_profile Si TRUE (defaut), compare au profil des membres adjacents du groupe
#'   cible plutot qu'au centroide global du groupe. Plus representatif pour les grands groupes.
#' @param weights Vecteur de poids pour les variables (meme longueur que vars_attr)
#' @param standardize Standardiser les variables avant calcul (defaut TRUE)
#' @return data.frame avec id_row et best_group
#' @export
assign_best_group <- function(data, group_var, vars_attr, profiles,
                              nb            = NULL,
                              local_profile = TRUE,
                              weights       = NULL,
                              standardize   = TRUE) {

  candidates     <- data |> dplyr::filter(candidate)
  non_candidates <- data |> dplyr::filter(!candidate)

  df_nc <- sf::st_drop_geometry(non_candidates)
  df_c  <- sf::st_drop_geometry(candidates)

  # ── Standardisation et ponderation (coherentes avec compute_affinity) ────────
  # On applique les memes transformations sur les donnees completes
  df_all_attr <- sf::st_drop_geometry(data)[, vars_attr, drop = FALSE]

  if (standardize) {
    # Parametres de standardisation calcules sur l'ensemble des donnees
    col_means <- colMeans(df_all_attr, na.rm = TRUE)
    col_sds   <- apply(df_all_attr, 2, sd, na.rm = TRUE)
    col_sds[col_sds == 0] <- 1  # eviter division par zero

    # Standardiser candidates et non-candidates avec les memes parametres
    df_nc_attr <- sweep(sweep(df_nc[, vars_attr, drop = FALSE],
                               2, col_means, "-"), 2, col_sds, "/")
    df_c_attr  <- sweep(sweep(df_c[, vars_attr, drop = FALSE],
                               2, col_means, "-"), 2, col_sds, "/")
  } else {
    df_nc_attr <- df_nc[, vars_attr, drop = FALSE]
    df_c_attr  <- df_c[, vars_attr, drop = FALSE]
  }

  if (!is.null(weights)) {
    if (length(weights) != length(vars_attr))
      stop("'weights' doit avoir la meme longueur que 'vars_attr'")
    w_norm     <- weights / sum(weights)
    w_sqrt     <- sqrt(w_norm)
    df_nc_attr <- sweep(df_nc_attr, 2, w_sqrt, "*")
    df_c_attr  <- sweep(df_c_attr,  2, w_sqrt, "*")
  }

  touches <- sf::st_touches(candidates, non_candidates)

  results <- lapply(seq_len(nrow(candidates)), function(i) {

    xi             <- as.numeric(df_c_attr[i, ])
    idx            <- touches[[i]]
    groupe_origine <- as.character(df_c[[group_var]][i])

    if (length(idx) == 0) {
      idx <- sf::st_nearest_feature(candidates[i, ], non_candidates)
    }

    neighbor_groups <- unique(df_nc[[group_var]][idx])
    neighbor_groups <- neighbor_groups[neighbor_groups != groupe_origine]

    if (length(neighbor_groups) == 0) {
      return(data.frame(id_row = candidates$id_row[i], best_group = groupe_origine))
    }

    best_group <- NA_character_
    min_dist   <- Inf

    for (g in neighbor_groups) {

      if (local_profile && !is.null(nb)) {
        # ── Point 3 : profil local = membres adjacents du groupe cible ─────────
        # Indices dans non_candidates des voisins appartenant au groupe g
        local_idx <- idx[df_nc[[group_var]][idx] == g]

        if (length(local_idx) > 0) {
          local_mat <- as.matrix(df_nc_attr[local_idx, , drop = FALSE])
          d <- mean(apply(local_mat, 1,
                          function(row) sqrt(sum((xi - row)^2, na.rm = TRUE))),
                    na.rm = TRUE)
        } else {
          # Fallback : centroide global si aucun voisin direct
          profile <- profiles |>
            dplyr::filter(.data[[group_var]] == g) |>
            dplyr::select(dplyr::all_of(vars_attr))
          if (nrow(profile) == 0) next

          # Appliquer la meme standardisation sur le profil global
          prof_std <- (as.numeric(profile) - col_means) / col_sds
          if (!is.null(weights)) prof_std <- prof_std * w_sqrt
          d <- sqrt(sum((xi - prof_std)^2, na.rm = TRUE))
        }

      } else {
        # ── Centroide global (comportement original) ───────────────────────────
        profile <- profiles |>
          dplyr::filter(.data[[group_var]] == g) |>
          dplyr::select(dplyr::all_of(vars_attr))
        if (nrow(profile) == 0) next
        d <- sqrt(sum((xi - as.numeric(profile))^2, na.rm = TRUE))
      }

      if (d < min_dist) {
        min_dist   <- d
        best_group <- g
      }
    }

    if (is.na(best_group)) best_group <- groupe_origine

    data.frame(id_row = candidates$id_row[i], best_group = best_group)
  })

  dplyr::bind_rows(results)
}
