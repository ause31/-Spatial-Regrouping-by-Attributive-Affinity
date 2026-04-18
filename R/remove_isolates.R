#' Réintégrer les communes candidates isolées
#' @param data Un objet sf avec colonne "candidate"
#' @param nb Matrice de voisinage
#' @return data avec candidats isolés remis à FALSE
#' @export
remove_isolates <- function(data, nb) {
  candidate_ids <- which(data$candidate)
  isolat_idx <- sapply(candidate_ids, function(i) {
    voisins <- unlist(nb[[i]])
    voisins <- voisins[!is.na(voisins) & voisins > 0]
    if (length(voisins) == 0) return(TRUE)
    !any(data$candidate[voisins])
  })
  a_corriger <- candidate_ids[isolat_idx]
  if (length(a_corriger) > 0) {
    message(length(a_corriger), " isolat(s) réintégré(s)")
    data$candidate[a_corriger] <- FALSE
  } else {
    message("Aucun isolat détecté")
  }
  return(data)
}
