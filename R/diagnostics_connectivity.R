#' Check MET Connectivity (Island Detection)
#'
#' Scans a Multi-Environment Trial (MET) dataset to identify disconnected clusters
#' of environments. Disconnected trials ("islands") will result in singular
#' factor analytic models.
#'
#' @param data Dataframe containing the MET data.
#' @param trial_col Character. Column name for Trial/Environment.
#' @param gen_col Character. Column name for Genotype.
#' @param stop_if_disconnected Logical. If TRUE, stops execution if islands are found.
#' @return A list containing:
#'   \item{n_clusters}{Integer. Number of disconnected components.}
#'   \item{membership}{Dataframe mapping Trial to Cluster.}
#'   \item{islands}{Character vector of Trial names that are disconnected from the main group.}
#' @import cli
#' @export
check_met_connectivity <- function(data, trial_col = "Trial", gen_col = "Genotype", stop_if_disconnected = FALSE) {
    if (!trial_col %in% names(data)) cli::cli_abort("Column {.val {trial_col}} not found in data.")
    if (!gen_col %in% names(data)) cli::cli_abort("Column {.val {gen_col}} not found in data.")

    trials <- unique(as.character(data[[trial_col]]))
    n_t <- length(trials)

    if (n_t < 2) {
        return(list(n_clusters = 1, membership = data.frame(Trial = trials, Cluster = 1), islands = character(0)))
    }

    # 1. Build Adjacency using Sparse Logic -------------------------------------
    pairs <- unique(data[, c(trial_col, gen_col)])
    gen_list <- split(pairs[[gen_col]], pairs[[trial_col]])

    adj_list <- vector("list", n_t)
    names(adj_list) <- trials

    # Efficient Connectivity Check
    # We iterate through unique genotypes to see which trials they connect
    # This is O(G * T_per_G^2) which is usually better than O(T^2 * G) if set intersections are slow
    # Actually, Jaccard intersection approach is fine for T < 500

    for (i in 1:(n_t - 1)) {
        t1 <- trials[i]
        g1 <- gen_list[[t1]]
        # optimization: store lengths?
        for (j in (i + 1):n_t) {
            t2 <- trials[j]
            g2 <- gen_list[[t2]]

            # Intersection
            if (length(intersect(g1, g2)) > 0) {
                adj_list[[t1]] <- c(adj_list[[t1]], t2)
                adj_list[[t2]] <- c(adj_list[[t2]], t1)
            }
        }
    }

    # 2. BFS for Connected Components -------------------------------------------
    visited <- rep(FALSE, n_t)
    names(visited) <- trials
    cluster_id <- rep(NA, n_t)
    names(cluster_id) <- trials
    cid <- 0

    for (t in trials) {
        if (!visited[t]) {
            cid <- cid + 1
            queue <- c(t)
            visited[t] <- TRUE
            cluster_id[t] <- cid

            while (length(queue) > 0) {
                curr <- queue[1]
                queue <- queue[-1]

                neighbors <- adj_list[[curr]]
                for (nb in neighbors) {
                    if (!visited[nb]) {
                        visited[nb] <- TRUE
                        cluster_id[nb] <- cid
                        queue <- c(queue, nb)
                    }
                }
            }
        }
    }

    # 3. Analyze & Report -------------------------------------------------------
    counts <- table(cluster_id)
    main_cluster <- as.numeric(names(counts)[which.max(counts)])
    cluster_df <- data.frame(Trial = names(cluster_id), Cluster = cluster_id, stringsAsFactors = FALSE)

    islands <- cluster_df$Trial[cluster_df$Cluster != main_cluster]

    res <- list(
        n_clusters = length(unique(cluster_id)),
        membership = cluster_df,
        islands = islands
    )

    if (length(islands) > 0) {
        cli::cli_alert_warning("MET Connectivity Issue: Found {length(islands)} disconnected trials (Islands).")
        cli::cli_ul(head(islands, 10))

        if (stop_if_disconnected) {
            cli::cli_abort("Stopping analysis due to disconnected trial network. Please remove islands or bridge them.")
        }
    } else {
        cli::cli_alert_success("MET Connectivity: All trials are connected.")
    }

    return(res)
}
