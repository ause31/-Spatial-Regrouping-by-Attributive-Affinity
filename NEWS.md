# spatialRegroup 0.2.1

## Bug fixes & improvements

* `spatialRegroup()` : correction du compteur `n_reclassed` — seules les unités
  dont l'affectation a **réellement changé** sont désormais comptées (au lieu du
  nombre de candidats). Cela permet une détection correcte de la convergence.
* `spatialRegroup()` : ajout d'une détection de cycle — l'algorithme s'arrête si
  un état de partition déjà observé réapparaît, évitant les boucles infinies en
  cas d'oscillation entre deux états.
* `spatialRegroup()` : le voisinage spatial (`build_neighbors()`) est maintenant
  calculé **une seule fois avant la boucle** (gain de ~25 % sur le temps
  d'exécution pour 6 itérations).
* Mise à jour de la documentation roxygen (`@param`, `@return`).

# spatialRegroup 0.1.0

## New features

* First release of spatialRegroup.
* `spatialRegroup()` : main function for spatial regrouping by attributive affinity.
* `compute_affinity()` : computes intra vs extra group attributive distance.
* `build_neighbors()` : builds spatial neighborhood matrix.
* `remove_isolates()` : reintegrates isolated candidates.
* `compute_profiles()` : computes mean attributive profiles per group.
* `assign_best_group()` : assigns each candidate to the most similar neighboring group.

