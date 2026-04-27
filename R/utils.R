#' @importFrom rlang .data :=
NULL

utils::globalVariables(c(
  "candidate",
  "best_group",
  ".data"
))

#' Calculer la coherence interne d'un groupement (eta-carre moyen)
#'
#' Mesure la part de variance inter-groupe sur la variance totale pour chaque
#' variable attributaire, puis en fait la moyenne. Valeur entre 0 et 1 :
#' plus elle est elevee, plus le groupement est coherent avec les variables.
#' Equivalent a un eta² (rapport de correlation) sans necessiter de modele mixte.
#'
#' @param data Un objet sf ou data.frame
#' @param group_var Nom de la variable de groupe
#' @param vars_attr Vecteur de noms de variables attributaires
#' @return Un scalaire entre 0 et 1 (moyenne des eta² par variable)
#' @export
compute_cohesion <- function(data, group_var, vars_attr) {
  df <- sf::st_drop_geometry(data)
  g  <- as.character(df[[group_var]])

  eta2_vals <- sapply(vars_attr, function(v) {
    x <- as.numeric(df[[v]])
    x <- x[!is.na(x)]
    g_v <- g[!is.na(df[[v]])]
    if (length(unique(g_v)) < 2) return(NA_real_)

    ss_total  <- sum((x - mean(x))^2)
    if (ss_total == 0) return(NA_real_)

    group_means <- tapply(x, g_v, mean, na.rm = TRUE)
    group_ns    <- table(g_v)
    ss_between  <- sum(as.numeric(group_ns) *
                         (group_means[names(group_ns)] - mean(x))^2)

    ss_between / ss_total
  })

  mean(eta2_vals, na.rm = TRUE)
}
