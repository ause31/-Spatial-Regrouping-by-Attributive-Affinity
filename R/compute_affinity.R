#' Calculer les affinites attributaires des communes frontalières
#'
#' @param data Un objet sf
#' @param nb Matrice de voisinage (spdep)
#' @param group_var Nom de la variable de groupe (ex: "EPCI")
#' @param vars_attr Vecteur de noms de variables attributaires
#' @param method Methode de distance : "euclidean" (defaut) ou "mahalanobis"
#' @param threshold Seuil d'affinite au-dela duquel une unite est candidate (defaut 0)
#' @param weights Vecteur de poids pour les variables (meme longueur que vars_attr,
#'   defaut NULL = poids egaux). Applique apres standardisation.
#' @param standardize Standardiser les variables avant calcul (defaut TRUE)
#' @return Un data.frame avec dist_intra, dist_extra, affinite et candidat
#' @export
compute_affinity <- function(data, nb, group_var, vars_attr,
                             method      = "euclidean",
                             threshold   = 0,
                             weights     = NULL,
                             standardize = TRUE) {

  df      <- sf::st_drop_geometry(data)
  df_attr <- df[, vars_attr, drop = FALSE]

  # ── Standardisation z-score (prerequis pour que les poids soient comparables)
  if (standardize) {
    df_attr <- as.data.frame(scale(df_attr))
  }

  # ── Ponderation : chaque variable multipliee par sa racine de poids ─────────
  # On pondere la racine pour que la distance euclidienne ponderee soit coherente
  # avec d² = sum(w_i * (x_i - y_i)²)
  if (!is.null(weights)) {
    if (length(weights) != length(vars_attr))
      stop("'weights' doit avoir la meme longueur que 'vars_attr'")
    w_norm  <- weights / sum(weights)
    df_attr <- sweep(df_attr, 2, sqrt(w_norm), "*")
  }

  # ── Matrice de covariance pour Mahalanobis (calculee une seule fois) ────────
  cov_mat <- if (method == "mahalanobis") {
    tryCatch(solve(cov(df_attr)), error = function(e) NULL)
  } else NULL

  # ── Identifier les unites frontalières ──────────────────────────────────────
  frontaliers <- sapply(seq_along(nb), function(i) {
    voisins <- unlist(nb[[i]])
    voisins <- voisins[!is.na(voisins) & voisins > 0]
    if (length(voisins) == 0) return(FALSE)
    any(df[[group_var]][voisins] != df[[group_var]][i], na.rm = TRUE)
  })

  # ── Fonction de distance (euclidienne ponderee ou mahalanobis) ───────────────
  dist_fn <- function(xi, mat) {
    if (method == "mahalanobis" && !is.null(cov_mat)) {
      mean(apply(mat, 1, function(row) {
        tryCatch(
          sqrt(as.numeric(t(xi - row) %*% cov_mat %*% (xi - row))),
          error = function(e) sqrt(sum((xi - row)^2, na.rm = TRUE))
        )
      }), na.rm = TRUE)
    } else {
      mean(apply(mat, 1, function(row) {
        sqrt(sum((xi - row)^2, na.rm = TRUE))
      }), na.rm = TRUE)
    }
  }

  dist_mean <- function(xi, ids) {
    if (length(ids) == 0) return(NA_real_)
    mat <- as.matrix(df_attr[ids, , drop = FALSE])
    dist_fn(xi, mat)
  }

  affinite_list <- lapply(which(frontaliers), function(i) {
    group_i <- df[[group_var]][i]
    voisins <- unlist(nb[[i]])
    voisins <- voisins[!is.na(voisins) & voisins > 0]
    xi      <- as.numeric(df_attr[i, ])
    intra   <- voisins[df[[group_var]][voisins] == group_i]
    extra   <- voisins[df[[group_var]][voisins] != group_i]

    data.frame(
      id_row     = i,
      dist_intra = dist_mean(xi, intra),
      dist_extra = dist_mean(xi, extra),
      affinite   = dist_mean(xi, intra) - dist_mean(xi, extra)
    )
  })

  aff          <- dplyr::bind_rows(affinite_list)
  aff$candidat <- !is.na(aff$affinite) & aff$affinite > threshold

  return(aff)
}
