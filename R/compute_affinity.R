#' Calculer les affinités attributaires des communes frontalières
#' @param data Un objet sf
#' @param nb Matrice de voisinage (spdep)
#' @param group_var Nom de la variable de groupe (ex: "EPCI")
#' @param vars_attr Vecteur de noms de variables attributaires
#' @return Un data.frame avec dist_intra, dist_extra, affinite
#' @export
compute_affinity <- function(data, nb, group_var, vars_attr) {

  df      <- sf::st_drop_geometry(data)
  df_attr <- df[, vars_attr, drop = FALSE]

  frontaliers <- sapply(seq_along(nb), function(i) {
    voisins <- unlist(nb[[i]])
    voisins <- voisins[!is.na(voisins) & voisins > 0]
    if (length(voisins) == 0) return(FALSE)
    any(df[[group_var]][voisins] != df[[group_var]][i], na.rm = TRUE)
  })

  dist_mean <- function(xi, ids) {
    if (length(ids) == 0) return(NA_real_)
    mat <- as.matrix(df_attr[ids, , drop = FALSE])
    mean(apply(mat, 1, function(row) sqrt(sum((xi - row)^2, na.rm = TRUE))),
         na.rm = TRUE)
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

  dplyr::bind_rows(affinite_list)
}
