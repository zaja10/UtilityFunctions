#' Check MET Connectivity (Island Detection)
#'
#' Scans a Multi-Environment Trial (MET) dataset to identify disconnected clusters
#' of environments. Disconnected trials ("islands") will result in singular
#' factor analytic models.
#'
#' @param data Dataframe containing the MET data.
#' @param trial_col Character. Column name for Trial/Environment.
#' @param gen_col Character. Column name for Genotype.
#' @return A list containing:
#'   \item{n_clusters}{Integer. Number of disconnected components.}
#'   \item{membership}{Dataframe mapping Trial to Cluster.}
#'   \item{islands}{Character vector of Trial names that are disconnected from the main group.}
#' @export
check_met_connectivity <- function(data, trial_col = "Trial", gen_col = "Genotype") {
    trials <- unique(as.character(data[[trial_col]]))
    n_t <- length(trials)

    # Adjacency Matrix
    # Two trials are connected if they share at least one genotype
    # Using matrix math for speed: X'X where X is Genotype x Trial incidence

    # Create Incidence Matrix
    # We can use xtabs
    form <- as.formula(paste("~", gen_col, "+", trial_col))
    # Sparse logical incidence
    # Use explicit tabulation to avoid huge generic matrix if many genotypes
    # Actually, simpler: just get unique pairs
    pairs <- unique(data[, c(trial_col, gen_col)])

    # This might be slow for massive data, but robust.
    # Alternative: Use split
    gen_list <- split(pairs[[gen_col]], pairs[[trial_col]])

    # Build Adjacency List (BFS Graph)
    adj_list <- vector("list", n_t)
    names(adj_list) <- trials

    # Naive O(N^2) comparison is fine for N_sites < 200
    # Jaccard check logic from before, but binary
    for (i in 1:(n_t - 1)) {
        t1 <- trials[i]
        g1 <- gen_list[[t1]]
        for (j in (i + 1):n_t) {
            t2 <- trials[j]
            g2 <- gen_list[[t2]]

            common <- intersect(g1, g2)
            if (length(common) > 0) {
                adj_list[[t1]] <- c(adj_list[[t1]], t2)
                adj_list[[t2]] <- c(adj_list[[t2]], t1)
            }
        }
    }

    # BFS for Connected Components
    visited <- rep(FALSE, n_t)
    names(visited) <- trials
    cluster_id <- rep(NA, n_t)
    names(cluster_id) <- trials
    cid <- 0

    for (t in trials) {
        if (!visited[t]) {
            cid <- cid + 1
            # BFS Queue
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

    # Analyze Clusters
    counts <- table(cluster_id)
    main_cluster <- as.numeric(names(counts)[which.max(counts)])
    cluster_df <- data.frame(Trial = names(cluster_id), Cluster = cluster_id, stringsAsFactors = FALSE)

    islands <- cluster_df$Trial[cluster_df$Cluster != main_cluster]

    if (length(islands) > 0) {
        cli::cli_alert_warning("MET Connectivity Issue: Found {length(islands)} disconnected trials (Islands).")
        cli::cli_ul(head(islands, 5))
    } else {
        cli::cli_alert_success("MET Connectivity: All trials are connected.")
    }

    return(list(
        n_clusters = length(unique(cluster_id)),
        membership = cluster_df,
        islands = islands
    ))
}
