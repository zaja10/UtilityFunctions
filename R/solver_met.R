# ==============================================================================
# SMART SOLVER & EXTRACTION ENGINE
# ==============================================================================

#' Fit and Extract Multi-Environment Trial (MET) Models
#'
#' @description
#' The core engine of GenomicFlow. It fits a Factor Analytic (FA) model using ASReml-R
#' based on a `design_tableau` structure, or extracts parameters from an already fitted model.
#'
#' It performs post-analysis rotation (SVD/PC or Varimax) to align loadings for biological
#' interpretation and selection.
#'
#' @param object A `design_tableau` object (for fitting) OR a fitted `asreml` object (for extraction).
#' @param k Integer. Number of factors for the FA model (default 1).
#' @param genotype Character. Column name of the Genotype variable (required if fitting).
#' @param site Character. Column name of the Site/Environment variable (required if fitting).
#' @param rotation Character. Rotation method: "pc" (Principal Component/SVD) or "varimax".
#' @param ... Additional arguments passed to \code{asreml()}.
#'
#' @return An object of class \code{fa_asreml} containing loadings, scores, and FAST indices.
#' @importFrom stats coef varimax
#' @export
fit_met_model <- function(object, k = 1, genotype = NULL, site = NULL, rotation = "pc", response = NULL, ...) {
  model <- NULL

  # 1. FIT MODE (if object is design_tableau)
  if (inherits(object, "design_tableau")) {
    if (!requireNamespace("asreml", quietly = TRUE)) {
      stop("ASReml-R is required for fitting models. Please install it.")
    }

    if (is.null(genotype) || is.null(site)) stop("Must specify 'genotype' and 'site' column names for fitting.")

    cat("--- Fitting MET Model with ASReml ---\n")

    # Construct FA Term: fa(Site, k):Genotype
    fa_term <- paste0("fa(", site, ", ", k, "):", genotype)

    # Random Formula: Tableau Random + FA Term
    rand_form <- stats::update(object$formulas$random, paste("~ . +", fa_term))

    # Fixed Formula handling
    fix_form <- object$formulas$fixed

    # If response is provided, prepend it to the formula
    if (!is.null(response)) {
      # Check if formula already has a response
      if (length(fix_form) == 3) {
        warning("Fixed formula already has a response variable. Ignoring 'response' argument.")
      } else {
        fix_form <- stats::update(fix_form, paste(response, "~ ."))
      }
    }

    # Validation: Ensure fixed formula has a response
    if (length(fix_form) < 3) {
      stop("Fixed formula must have a response variable (e.g., 'yield ~ ...'). Please provide 'response' argument.")
    }

    # Prepare arguments
    args <- list(...)
    # Default residual if not provided
    if (!"residual" %in% names(args)) {
      # A common MET residual is ~ dsum(~id(units)|Site)
    }

    call_args <- c(
      list(
        fixed = fix_form,
        random = rand_form,
        data = object$data
      ),
      args
    )

    model <- do.call(asreml::asreml, call_args)

    if (!model$converge) warning("ASReml model did not converge.")
  } else if (inherits(object, "asreml")) {
    # 2. EXTRACT MODE (legacy support)
    model <- object
  } else {
    stop("Object must be of class 'design_tableau' or 'asreml'.")
  }

  # 3. EXTRACTION LOGIC
  # We use the robust extraction logic adapted from the original fa.asreml
  extracted <- .extract_fa_components(model, k)

  lambda_mat <- extracted$loadings
  f_mat <- extracted$scores
  psi_df <- extracted$psi

  # 4. ROTATION
  cat(sprintf("-> Applying Rotation: %s\n", rotation))

  lambda_rot <- lambda_mat
  f_rot <- f_mat
  rot_mat <- diag(k)

  if (k > 1) {
    if (rotation == "pc") {
      # SVD for PC solution (Factor 1 = Max Var)
      svd_res <- svd(lambda_mat)
      rot_mat <- svd_res$v
      lambda_rot <- lambda_mat %*% rot_mat

      # Flip sign if mean loading is negative (interpretability)
      if (mean(lambda_rot[, 1]) < 0) {
        rot_mat[, 1] <- -rot_mat[, 1]
        lambda_rot[, 1] <- -lambda_rot[, 1]
      }
      f_rot <- f_mat %*% rot_mat
    } else if (rotation == "varimax") {
      # Varimax on loadings
      if (ncol(lambda_mat) < 2) {
        warning("Varimax rotation requires >= 2 factors. Skipping.")
      } else {
        vm <- stats::varimax(lambda_mat)
        rot_mat <- vm$rotmat
        lambda_rot <- lambda_mat %*% rot_mat
        f_rot <- f_mat %*% rot_mat
      }
    }
  }

  colnames(lambda_rot) <- colnames(f_rot) <- paste0("Fac", 1:k)

  # 5. CONSTRUCT OUTPUT
  # Re-calculate G and Cor based on rotated solution
  # (Psi is invariant to rotation)
  common_sites <- intersect(rownames(lambda_rot), psi_df$Site)
  if (length(common_sites) == 0 && nrow(psi_df) > 0) warning("Site mismatch between Lambda and Psi.")

  G_est <- matrix(NA, nrow(lambda_rot), nrow(lambda_rot))
  C_est <- matrix(NA, nrow(lambda_rot), nrow(lambda_rot))
  vaf_df <- NULL

  if (length(common_sites) > 0) {
    lam_ord <- lambda_rot[common_sites, , drop = FALSE]
    psi_ord <- diag(psi_df$Psi[match(common_sites, psi_df$Site)])
    G_est <- (lam_ord %*% t(lam_ord)) + psi_ord
    C_est <- tryCatch(stats::cov2cor(G_est), error = function(e) G_est)

    vaf_pct <- ifelse(diag(G_est) > 1e-8, (diag(lam_ord %*% t(lam_ord)) / diag(G_est)) * 100, 0)
    vaf_df <- data.frame(Site = common_sites, VAF = vaf_pct)
  }

  # FAST indices
  fast_df <- NULL
  if (nrow(f_rot) > 0) {
    OP <- f_rot[, 1] * mean(lambda_rot[, 1], na.rm = TRUE)
    RMSD <- if (k > 1) apply(f_rot[, 2:k, drop = F] %*% t(lambda_rot[, 2:k, drop = F]), 1, function(x) sqrt(mean(x^2))) else rep(0, nrow(f_rot))
    fast_df <- data.frame(Genotype = rownames(f_rot), OP = OP, RMSD = RMSD)
    fast_df <- fast_df[order(fast_df$OP, decreasing = TRUE), ]
  }

  out <- list(
    loadings = list(raw = lambda_mat, rotated = lambda_rot),
    scores = list(raw = f_mat, rotated = f_rot),
    var_comp = list(psi = psi_df, vaf = vaf_df),
    matrices = list(G = G_est, Cor = C_est),
    fast = fast_df,
    rotation_matrix = rot_mat,
    meta = list(k = k, rotation = rotation)
  )

  class(out) <- "fa_asreml"
  return(out)
}

#' @importFrom dplyr %>% mutate select arrange desc bind_rows distinct
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_extract
#' @importFrom tibble column_to_rownames
.extract_fa_components <- function(model, k) {
  vc <- summary(model)$varcomp
  vc_names <- rownames(vc)
  coefs <- coef(model)$random
  coef_names <- rownames(coefs)

  # --- LOADINGS ---
  lambda_df <- data.frame()
  # Regex to find fa(Site, k).SiteName!fa_i or similar
  # Matches: anything!fa_i where i is 1..k
  for (i in 1:k) {
    # This pattern looks for !fa_1 or !fa1 at end of string
    pat <- paste0("!fa[_]?", i, "$")
    rows <- grep(pat, vc_names)
    if (length(rows) > 0) {
      # Extract Site Name.
      # Valid formats: "fa(Site,2).SiteA!fa_1" -> SiteA
      # "rr(Site,2).SiteA!rr_1" -> SiteA
      # We scan for the part between the term close and the !fa part
      # Regex: [something] . [SiteName] ! [fa_i]
      # Simplified: Remove the suffix !fa_i, then take everything after the last dot?
      # Or ASReml 3 style: fa(Site,2)_SiteA!fa_1 ?

      # We use a robust split on '!' and '.'
      for (r in rows) {
        term <- vc_names[r]
        val <- vc[r, "component"]

        # Try to parse site
        # Strategy: Remove the !fa part.
        base <- sub("!.*", "", term)
        # Now base is "fa(Site,2).SiteA".
        # Site is everything after the last dot or underscore?
        # Or between ) and end?

        site_match <- regexpr("(?<=\\.|_)[^\\._]+$", base, perl = TRUE)
        if (site_match > 0) {
          site <- substring(base, site_match, site_match + attr(site_match, "match.length") - 1)
        } else {
          # Fallback: if term is simple?
          site <- base
        }

        lambda_df <- rbind(lambda_df, data.frame(Site = site, Value = val, Factor = i))
      }
    }
  }

  if (nrow(lambda_df) > 0) {
    lambda_mat <- tidyr::pivot_wider(lambda_df, names_from = Factor, values_from = Value, values_fill = 0) %>%
      tibble::column_to_rownames("Site") %>%
      as.matrix()
    # Correct order if needed
    if (ncol(lambda_mat) < k) {
      # Fill missing cols with 0
      miss <- setdiff(1:k, colnames(lambda_mat))
      for (m in miss) lambda_mat <- cbind(lambda_mat, 0)
      colnames(lambda_mat)[(ncol(lambda_mat) - length(miss) + 1):ncol(lambda_mat)] <- miss
      lambda_mat <- lambda_mat[, order(as.numeric(colnames(lambda_mat)))]
    }
  } else {
    lambda_mat <- matrix(0, 0, k)
  }

  # --- PSI ---
  # Look for !var terms
  psi_df <- data.frame()
  var_rows <- grep("!var$", vc_names)
  if (length(var_rows) > 0) {
    for (r in var_rows) {
      term <- vc_names[r]
      val <- vc[r, "component"]
      # Remove !var
      base <- sub("!var", "", term)
      # Extract site same as above
      site_match <- regexpr("(?<=\\.|_)[^\\._]+$", base, perl = TRUE)
      if (site_match > 0) {
        site <- substring(base, site_match, site_match + attr(site_match, "match.length") - 1)
        psi_df <- rbind(psi_df, data.frame(Site = site, Psi = val))
      }
    }
  }

  # --- SCORES ---
  # Look for Comp_i or Fac_i in coefs
  f_mat <- matrix(0, 0, k)
  scores_df <- data.frame()

  for (i in 1:k) {
    # Fac_1 or Comp_1
    pat <- paste0("^(Fac|Comp)[_]?", i)
    rows <- grep(pat, coef_names)
    if (length(rows) > 0) {
      for (r in rows) {
        term <- coef_names[r]
        val <- coefs[r, 1]

        # Extract Genotype (everything after last delimiter)
        # "fa(Site,2):Genotype_G1" -> G1
        # "Comp_1_G1"

        parts <- strsplit(term, "[:_]")[[1]]
        gen <- tail(parts, 1) # Naive

        scores_df <- rbind(scores_df, data.frame(Genotype = gen, Value = val, Factor = i))
      }
    }
  }

  if (nrow(scores_df) > 0) {
    f_mat <- tidyr::pivot_wider(scores_df, names_from = Factor, values_from = Value, values_fill = 0)
    if ("Genotype" %in% names(f_mat)) {
      f_mat <- tibble::column_to_rownames(f_mat, "Genotype") %>% as.matrix()
    } else {
      f_mat <- as.matrix(f_mat) # Should not happen
    }
  }

  return(list(loadings = lambda_mat, scores = f_mat, psi = psi_df))
}
