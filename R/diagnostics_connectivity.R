#' Check MET Connectivity (Island Detection)
#'
#' Builds the trial-by-trial matrix of shared genotypes and finds connected
#' components. Two trials are linked when they share at least
#' \code{min_shared} genotypes. Disconnected trials ("islands") make genetic
#' covariances between trial groups inestimable; trials with only a few links
#' are flagged because factor analytic covariances need far more than one
#' common genotype to be estimated well.
#'
#' @param data Data frame containing the MET data.
#' @param trial_col Character. Trial/environment column.
#' @param gen_col Character. Genotype column.
#' @param min_shared Integer. Minimum shared genotypes for two trials to be
#'   linked (default 1).
#' @param warn_below Integer. Report trials whose best link to any other trial
#'   shares fewer than this many genotypes (default 5).
#' @param stop_if_disconnected Logical. Abort if islands are found.
#' @return A list with \code{n_clusters}, \code{membership} (\code{Trial},
#'   \code{Cluster}), \code{islands} (trials outside the largest cluster),
#'   \code{shared} (matrix of shared genotype counts) and \code{weak}
#'   (trials whose best link is below \code{warn_below}). Rows with missing
#'   trial or genotype are ignored.
#' @export
check_met_connectivity <- function(data, trial_col = "Trial", gen_col = "Genotype", min_shared = 1,
                                   warn_below = 5, stop_if_disconnected = FALSE) {
  for (col in c(trial_col, gen_col)) {
    if (!col %in% names(data)) cli::cli_abort("Column {.val {col}} not found in data.")
  }
  tr <- as.character(data[[trial_col]])
  ge <- as.character(data[[gen_col]])
  keep <- !is.na(tr) & !is.na(ge)
  if (any(!keep)) cli::cli_alert_info("Ignoring {sum(!keep)} row{?s} with missing trial or genotype.")
  tr <- tr[keep]
  ge <- ge[keep]
  trials <- unique(tr)

  inc <- table(factor(ge, levels = unique(ge)), factor(tr, levels = trials)) > 0
  shared <- crossprod(inc + 0)
  dimnames(shared) <- list(trials, trials)

  adj <- shared >= min_shared
  diag(adj) <- FALSE
  cluster <- connected_components(adj)
  sizes <- table(cluster)
  main <- as.integer(names(sizes)[which.max(sizes)])
  membership <- data.frame(Trial = trials, Cluster = unname(cluster), stringsAsFactors = FALSE)
  islands <- trials[cluster != main]

  off <- shared
  diag(off) <- 0
  best_link <- if (length(trials) > 1) apply(off, 1, max) else stats::setNames(0, trials)
  weak <- trials[best_link < warn_below & length(trials) > 1]

  if (length(islands) > 0) {
    cli::cli_alert_warning("Found {length(islands)} disconnected trial{?s} (islands): {.val {utils::head(islands, 10)}}")
    if (stop_if_disconnected) cli::cli_abort("Disconnected trial network. Remove or bridge the islands.")
  } else {
    cli::cli_alert_success("All {length(trials)} trial{?s} are connected.")
  }
  if (length(weak) > 0) {
    cli::cli_alert_warning("{length(weak)} trial{?s} share fewer than {warn_below} genotypes with every other trial.")
  }

  list(
    n_clusters = length(sizes), membership = membership, islands = islands,
    shared = shared, weak = weak
  )
}

# Connected components of a symmetric logical adjacency matrix (BFS).
connected_components <- function(adj) {
  n <- nrow(adj)
  comp <- integer(n)
  id <- 0L
  for (s in seq_len(n)) {
    if (comp[s] != 0L) next
    id <- id + 1L
    comp[s] <- id
    frontier <- s
    while (length(frontier) > 0) {
      nb <- which(colSums(adj[frontier, , drop = FALSE]) > 0 & comp == 0L)
      comp[nb] <- id
      frontier <- nb
    }
  }
  stats::setNames(comp, rownames(adj))
}
