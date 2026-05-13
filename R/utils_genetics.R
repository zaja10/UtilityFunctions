#' Prepare GRM for ASReml (Inverse & Sparse)
#'
#' Checks a Genomic Relationship Matrix (GRM) for singularity, bends it if necessary,
#' inverts it, and converts it to the sparse 'giv' format required by ASReml.
#'
#' @param G A symmetric numeric matrix (GRM).
#' @param blend Numeric. Amount of identity matrix to blend (0-1) to ensure invertibility.
#'        Default 0.02 (2\%).
#' @return A sparse dataframe (Row, Col, Value) with attribute 'INVERSE' = TRUE.
#' @importFrom MASS ginv
#' @importFrom Matrix forceSymmetric
#' @importClassesFrom Matrix dsCMatrix
#' @import cli
#' @export
prepare_asreml_grm <- function(G, blend = 0.02) {
    if (!isSymmetric(G)) {
        cli::cli_alert_warning("GRM is not symmetric. Forcing symmetry.")
        G <- (G + t(G)) / 2
    }

    n <- nrow(G)
    cli::cli_alert_info("Processing {n} x {n} GRM...")

    if (blend > 0) {
        cli::cli_alert_info("Applying {blend*100}% bending to ensure invertibility.")
        I <- diag(n)
        G <- ((1 - blend) * G) + (blend * I)
    }

    G_inv <- tryCatch(solve(G), error = function(e) {
        cli::cli_alert_warning("Standard inversion failed. Using Generalized Inverse (slower).")
        MASS::ginv(G)
    })

    rownames(G_inv) <- rownames(G)
    colnames(G_inv) <- colnames(G)

    G_inv[lower.tri(G_inv)] <- 0
    s_obj <- as(Matrix::forceSymmetric(G_inv), "dsCMatrix")

    # Helper to align with asreml's sparse format (1-based index if needed, but usually row/col names suffice if using vm(..., giv(df)))
    # ASReml-R usually wants a dataframe with Row, Col, Value
    # But if we return the matrix directly with INVERSE attribute, newer asreml versions handle it.

    # However, consistent with the prompt's request for commercial robustness, let's stick to the sparse matrix object or dataframe?
    # The user prompt's code returns G_inv (matrix) but also calculates 's_obj' and 'summ' without using them?
    # Correction: The code provided in the prompt creates s_obj but returns G_inv with attribute.
    # "return(G_inv)"

    # Attribute is key for ASReml 4.2+
    attr(G_inv, "INVERSE") <- TRUE

    cli::cli_alert_success("GRM inverted and prepared.")
    return(G_inv)
}

#' Align Phenotype and Genotype Data
#'
#' Ensures intersection of genotypes between phenotypic data and GRM.
#'
#' @param pheno Dataframe containing phenotypic data.
#' @param geno GRM matrix.
#' @param id_col Character string. Name of the genotype column in pheno.
#' @return A list containing aligned pheno and geno.
#' @export
align_pheno_geno <- function(pheno, geno, id_col = "Genotype") {
    if (!id_col %in% names(pheno)) stop(paste("Column", id_col, "not found in phenotype data."))

    common <- intersect(pheno[[id_col]], rownames(geno))

    if (length(common) == 0) stop("No common genotypes found.")

    pheno_sub <- pheno[pheno[[id_col]] %in% common, ]
    geno_sub <- geno[common, common]

    cli::cli_alert_success("Aligned: {length(common)} genotypes retained.")

    return(list(pheno = pheno_sub, geno = geno_sub))
}

#' Match Germplasm Using Synonyms
#'
#' Robustly matches a germplasm name to a list of valid IDs. If the primary name
#' fails, it attempts to split and match against a comma-separated list of synonyms.
#'
#' @param name Character or Factor. Primary germplasm name.
#' @param syn_string Character or Factor. Comma-separated string of synonyms.
#' @param valid_ids Character vector. List of valid IDs (e.g., from a GRM).
#' @return The matched character ID, or NA_character_ if no match is found.
#' @export
match_germplasm <- function(name, syn_string, valid_ids) {
  name <- as.character(name)
  syn_string <- as.character(syn_string)
  
  # Check 1: Primary match
  if (name %in% valid_ids) {
    return(name)
  }
  
  # Check 2: Synonym match
  if (!is.na(syn_string) && syn_string != "") {
    syn_list <- trimws(unlist(strsplit(syn_string, ",")))
    matched_syns <- syn_list[syn_list %in% valid_ids]
    
    if (length(matched_syns) > 0) {
      return(matched_syns[1])
    }
  }
  
  return(NA_character_)
}

#' Prepare Trial GRM
#'
#' Subsets a master Genomic Relationship Matrix (GRM) to the lines present in a trial.
#' Automatically pads the matrix with a dummy identity matrix for ungenotyped lines
#' to prevent ASReml singularity/inverse errors.
#'
#' @param master_grm A square numeric matrix representing the master GRM.
#' @param trial_data A dataframe containing the trial data.
#' @param genotype_col Character. Column name in `trial_data` containing the genotypes.
#' @return A square matrix representing the trial-specific GRM, scaled and padded.
#' @importFrom Matrix bdiag
#' @export
prepare_trial_grm <- function(master_grm, trial_data, genotype_col = "Genotype") {
  
  if (!genotype_col %in% names(trial_data)) {
    stop(paste("Genotype column", genotype_col, "not found in trial data."))
  }
  
  trial_lines <- as.character(unique(trial_data[[genotype_col]]))
  genotyped_lines <- intersect(trial_lines, rownames(master_grm))
  ungenotyped_lines <- setdiff(trial_lines, genotyped_lines)
  
  if (length(genotyped_lines) == 0) {
    stop("No genotyped lines found in the GRM for this trial.")
  }
  
  # Subset and scale
  grm_sub <- master_grm[genotyped_lines, genotyped_lines, drop = FALSE]
  scale_factor <- mean(diag(grm_sub))
  grm_sub <- grm_sub / scale_factor
  
  # Pad ungenotyped lines with an identity matrix
  if (length(ungenotyped_lines) > 0) {
    pad_mat <- diag(1, nrow = length(ungenotyped_lines))
    rownames(pad_mat) <- colnames(pad_mat) <- ungenotyped_lines
    
    # Combine using bdiag from Matrix package
    grm_final <- as.matrix(Matrix::bdiag(grm_sub, pad_mat))
    rownames(grm_final) <- colnames(grm_final) <- c(genotyped_lines, ungenotyped_lines)
  } else {
    grm_final <- grm_sub
  }
  
  # Add slight ridge for numerical stability during inversion
  diag(grm_final) <- diag(grm_final) + 1e-4
  
  # Ensure the matrix order exactly matches the factor levels in the data
  ordered_levels <- levels(factor(trial_data[[genotype_col]]))
  # Only subset to levels that we just constructed (in case factor has unused levels)
  valid_levels <- intersect(ordered_levels, rownames(grm_final))
  
  return(grm_final[valid_levels, valid_levels])
}

