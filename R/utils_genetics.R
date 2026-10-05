#' Match Germplasm Names with Genomic Matrix IDs
#'
#' For each germplasm name, returns the name itself if it is a valid genomic
#' ID, otherwise the first of its comma-separated synonyms that is.
#' Vectorised over \code{name} and \code{syn_string}.
#'
#' @param name Primary germplasm name(s) (character or factor).
#' @param syn_string Comma-separated synonym string(s), recycled to the length
#'   of \code{name}. \code{NA} or \code{""} means no synonyms.
#' @param valid_ids Character vector of valid IDs in the genomic relationship matrix.
#' @return Character vector the length of \code{name}: the matched ID, or
#'   \code{NA_character_} where nothing matches.
#' @export
match_germplasm <- function(name, syn_string, valid_ids) {
  name <- as.character(name)
  syn_string <- rep_len(as.character(syn_string), length(name))
  valid_ids <- as.character(valid_ids)

  out <- ifelse(name %in% valid_ids, name, NA_character_)
  todo <- which(is.na(out) & !is.na(syn_string) & nzchar(syn_string))
  for (i in todo) {
    syns <- trimws(strsplit(syn_string[i], ",", fixed = TRUE)[[1]])
    hit <- syns[syns %in% valid_ids]
    if (length(hit) > 0) out[i] <- hit[1]
  }
  out
}

#' Subset and Pad a Genomic Relationship Matrix for One Trial
#'
#' Subsets a master GRM to the genotyped lines in a trial, rescales it so the
#' mean diagonal is 1, and appends ungenotyped lines as independent
#' (identity) rows/columns so every trial line has a relationship entry.
#'
#' @param master_grm Square relationship matrix (not inverse) with dimnames.
#' @param trial_data Data frame of the trial.
#' @param genotype_col Character. Genotype column in \code{trial_data}
#'   (default \code{"germplasmName"}).
#' @param ridge Numeric. Value added to the diagonal to guarantee positive
#'   definiteness (default \code{1e-4}).
#' @return Square matrix ordered as \code{levels(factor(trial_data[[genotype_col]]))}.
#' @export
prepare_trial_grm <- function(master_grm, trial_data, genotype_col = "germplasmName", ridge = 1e-4) {
  if (!genotype_col %in% names(trial_data)) {
    cli::cli_abort("Column {.val {genotype_col}} not found in {.arg trial_data}.")
  }
  if (is.null(rownames(master_grm))) cli::cli_abort("{.arg master_grm} must have row names.")

  geno <- factor(trial_data[[genotype_col]])
  trial_lines <- levels(geno)
  genotyped <- intersect(trial_lines, rownames(master_grm))
  ungenotyped <- setdiff(trial_lines, genotyped)
  if (length(genotyped) == 0) cli::cli_abort("No trial lines were found in the GRM.")

  sub <- master_grm[genotyped, genotyped, drop = FALSE]
  sub <- sub / mean(diag(sub))

  n_g <- length(genotyped)
  n_u <- length(ungenotyped)
  out <- matrix(0, n_g + n_u, n_g + n_u, dimnames = list(c(genotyped, ungenotyped), c(genotyped, ungenotyped)))
  out[seq_len(n_g), seq_len(n_g)] <- sub
  if (n_u > 0) diag(out)[n_g + seq_len(n_u)] <- 1
  diag(out) <- diag(out) + ridge

  out[trial_lines, trial_lines, drop = FALSE]
}

#' Invert a GRM into ASReml's Sparse Inverse Format
#'
#' Blends the relationship matrix towards the identity,
#' \eqn{G_b = (1 - blend) G + blend I}, inverts it, and returns the lower
#' triangle as the three-column (row, column, value) matrix ASReml-R expects
#' for \code{vm(..., ginv)}, with \code{rowNames} and \code{INVERSE}
#' attributes. (Dense inverse matrices can crash ASReml-R 4.2.)
#'
#' @param G Square relationship matrix with dimnames.
#' @param blend Numeric in [0, 1). Weight on the identity (default 0.02).
#' @param tol Numeric. Asymmetry tolerance before symmetrising (default 1e-8).
#' @return A numeric matrix with columns \code{row}, \code{column}, \code{value}.
#' @export
prepare_asreml_grm <- function(G, blend = 0.02, tol = 1e-8) {
  G <- as.matrix(G)
  if (nrow(G) != ncol(G)) cli::cli_abort("{.arg G} must be square.")
  ids <- rownames(G) %||% colnames(G)
  if (is.null(ids)) cli::cli_abort("{.arg G} must have row names.")
  if (blend < 0 || blend >= 1) cli::cli_abort("{.arg blend} must be in [0, 1).")

  if (max(abs(G - t(G))) > tol) {
    cli::cli_inform(c("!" = "GRM is not symmetric; using (G + t(G)) / 2."))
    G <- (G + t(G)) / 2
  }

  Gb <- (1 - blend) * G + blend * diag(nrow(G))
  R <- tryCatch(chol(Gb), error = function(e) NULL)
  if (is.null(R)) {
    cli::cli_abort(c(
      "x" = "Blended GRM is not positive definite.",
      "i" = "Increase {.arg blend} (currently {blend})."
    ))
  }
  Gi <- chol2inv(R)

  idx <- which(lower.tri(Gi, diag = TRUE), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  out <- cbind(row = idx[, 1], column = idx[, 2], value = Gi[idx])
  attr(out, "rowNames") <- ids
  attr(out, "INVERSE") <- TRUE
  out
}

#' Identifiers of a relationship matrix
#'
#' Works for dense matrices (row names) and ASReml sparse matrices
#' (\code{rowNames} attribute).
#' @keywords internal
#' @noRd
grm_ids <- function(grm) {
  ids <- attr(grm, "rowNames") %||% rownames(grm)
  if (is.null(ids)) {
    cli::cli_abort("Cannot determine genotype IDs: {.arg grm} has neither row names nor a {.field rowNames} attribute.")
  }
  as.character(ids)
}
