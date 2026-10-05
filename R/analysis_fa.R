#' Extract and Rotate Factor Analytic Model Parameters
#'
#' A parser for ASReml-R factor analytic (\code{fa}) and reduced rank
#' (\code{rr}) terms implementing the factor analytic selection tools (FAST)
#' of Smith & Cullis (2018).
#'
#' @param model A fitted object of class \code{asreml}.
#' @param classify Character (optional). The FA/RR term to extract, e.g.
#'   \code{"fa(Site, 2):vm(Genotype, Ginv)"}. A unique prefix such as
#'   \code{"fa(Site, 2)"} is accepted when only one term matches. If
#'   \code{NULL}, the model must contain exactly one FA/RR term.
#' @param psi_term Character (optional). For RR models, the term holding the
#'   specific variances, e.g. \code{"diag(Site):vm(Genotype, Ginv)"}.
#' @param rotate Rotation: \code{"svd"} (default; principal-component rotation
#'   with factor 1 oriented so its mean loading is positive, as required by
#'   FAST) or \code{"none"}. \code{TRUE}/\code{FALSE} map to these. The old
#'   values \code{"mean"} and \code{"varimax"} are accepted as aliases of
#'   \code{"svd"} (no varimax rotation was ever performed).
#' @param annotation Optional data frame or named list passed to
#'   \code{\link{annotate_model}}.
#' @param grm Optional relationship matrix \strong{G} (not its inverse) used
#'   in the model, or a named vector of its diagonal. Needed for correct
#'   reliabilities (\eqn{r^2_g}) with genomic (\code{vm()}) models because the variance of a
#'   genotype's factor score is \eqn{g_{ii}}, not 1.
#'
#' @return An object of class \code{fa_model} with elements \code{loadings},
#'   \code{scores}, \code{var_comp}, \code{matrices}, \code{fast} and
#'   \code{meta}.
#'
#' @details
#' Prediction error variances (PEVs) of rotated scores are computed from the
#' standard errors of the raw scores, \eqn{\sum_r v_{r1}^2 \mathrm{PEV}(f_{ir})}.
#' This ignores prediction error covariances between a genotype's raw scores,
#' which ASReml does not report through \code{summary(coef = TRUE)}; the
#' resulting PEV and reliability are approximations.
#'
#' \code{meta$var_explained} is the percentage of the total genetic variance
#' (\eqn{\mathrm{tr}(G)}) attributable to each rotated factor.
#'
#' @references Smith, A. B. & Cullis, B. R. (2018). Plant breeding selection
#'   tools built on factor analytic mixed models for multi-environment trial
#'   data. \emph{Euphytica} 214, 143.
#' @export
fa.asreml <- function(model, classify = NULL, psi_term = NULL, rotate = "svd",
                      annotation = NULL, grm = NULL) {
  if (!inherits(model, "asreml")) cli::cli_abort("Object must be of class {.cls asreml}.")
  check_asreml_availability()
  rotate <- normalise_rotation(rotate)

  vc <- asreml_components(model)
  term <- resolve_fa_term(names(vc), classify)
  spec <- parse_fa_term(term)
  cli::cli_alert_info("Extracting {.val {term}} (type = {spec$type}, k = {spec$k}).")

  # 1. Loadings and specific variances -----------------------------------
  par <- fa_parameter_table(vc, term)
  sites <- unique(par$Group)
  lambda <- matrix(NA_real_, length(sites), spec$k, dimnames = list(sites, paste0("Fac", seq_len(spec$k))))
  for (i in seq_len(spec$k)) {
    rows <- par[par$Param == paste0("fa", i), , drop = FALSE]
    lambda[rows$Group, i] <- rows$Value
  }
  if (anyNA(lambda)) {
    miss <- which(is.na(lambda), arr.ind = TRUE)
    cli::cli_abort("Missing loadings for {nrow(miss)} site/factor combination{?s} in {.val {term}}.")
  }

  psi <- extract_psi(vc, par, sites, spec, psi_term)

  # 2. Scores -------------------------------------------------------------
  sc <- fa_scores(model, term, spec)
  has_scores <- !is.null(sc)

  # 3. Rotation -----------------------------------------------------------
  V <- diag(spec$k)
  if (rotate == "svd") {
    V <- svd(lambda)$v
    if (mean((lambda %*% V)[, 1]) < 0) V[, 1] <- -V[, 1]
  }
  fac_names <- paste0("Fac", seq_len(spec$k))
  lambda_rot <- lambda %*% V
  dimnames(lambda_rot) <- list(sites, fac_names)
  f_rot <- pev_rot <- NULL
  if (has_scores) {
    f_rot <- sc$scores %*% V
    pev_rot <- sc$pev %*% (V^2)
    dimnames(f_rot) <- dimnames(pev_rot) <- list(rownames(sc$scores), fac_names)
  }

  # 4. Genetic covariance and variance accounted for ------------------------
  G_est <- lambda_rot %*% t(lambda_rot) + diag(psi, nrow = length(psi))
  dimnames(G_est) <- list(sites, sites)
  C_est <- if (all(diag(G_est) > 0)) stats::cov2cor(G_est) else NULL

  vaf <- data.frame(Group = sites, stringsAsFactors = FALSE)
  g_diag_env <- diag(G_est)
  for (i in seq_len(spec$k)) {
    vaf[[paste0("VAF_", fac_names[i])]] <- round(100 * lambda_rot[, i]^2 / g_diag_env, 2)
  }
  vaf$Total_VAF <- round(100 * rowSums(lambda_rot^2) / g_diag_env, 2)
  var_exp <- 100 * colSums(lambda_rot^2) / sum(g_diag_env)

  res <- list(
    loadings = list(raw = lambda, rotated = lambda_rot),
    scores = list(
      raw = if (has_scores) sc$scores, rotated = f_rot,
      pev_raw = if (has_scores) sc$pev, pev_rotated = pev_rot
    ),
    var_comp = list(psi = data.frame(Group = sites, Psi = psi, stringsAsFactors = FALSE), vaf = vaf),
    matrices = list(G = G_est, Cor = C_est, rotation = V),
    fast = NULL,
    meta = list(
      k = spec$k, type = spec$type, term = term, group = spec$site_var,
      genotype = spec$geno_var, rotation = rotate, var_explained = var_exp
    )
  )
  class(res) <- "fa_model"

  # 5. FAST indices --------------------------------------------------------
  if (has_scores && rotate == "svd") {
    res$fast <- calculate_fast_indices(res, g_diag = genotype_variance(grm, rownames(f_rot)))
  } else if (has_scores) {
    cli::cli_alert_info("FAST indices require {.code rotate = \"svd\"}; not computed.")
  }

  if (!is.null(annotation)) {
    res <- if (is.data.frame(annotation)) {
      annotate_model(res, df = annotation)
    } else {
      do.call(annotate_model, c(list(object = res), annotation))
    }
  }

  cli::cli_alert_success("Extraction complete (k = {spec$k}, rotation = {rotate}).")
  res
}

normalise_rotation <- function(rotate) {
  if (is.logical(rotate)) {
    return(if (isTRUE(rotate)) "svd" else "none")
  }
  if (rotate %in% c("mean", "varimax")) {
    cli::cli_inform(c("i" = "{.code rotate = \"{rotate}\"} is deprecated; using the SVD rotation ({.code \"svd\"})."))
    return("svd")
  }
  match.arg(rotate, c("svd", "none"))
}

# Find the full FA/RR term from variance component names.
resolve_fa_term <- function(vc_names, classify) {
  prefixes <- unique(sub("!.*$", "", vc_names))
  fa_terms <- prefixes[grepl("^(fa|rr)\\(", prefixes)]
  if (length(fa_terms) == 0) cli::cli_abort("No {.code fa()} or {.code rr()} term found in the model.")

  if (is.null(classify)) {
    if (length(fa_terms) > 1) {
      cli::cli_abort(c("x" = "Several FA/RR terms found; specify {.arg classify}.", "i" = "Terms: {.val {fa_terms}}"))
    }
    return(fa_terms)
  }
  strip <- function(x) gsub("\\s+", "", x)
  key <- strip(classify)
  exact <- fa_terms[strip(fa_terms) == key]
  if (length(exact) == 1) {
    return(exact)
  }
  partial <- fa_terms[startsWith(strip(fa_terms), paste0(key, ":"))]
  if (length(partial) == 1) {
    return(partial)
  }
  cli::cli_abort(c(
    "x" = "{.arg classify} = {.val {classify}} matches {length(partial)} FA/RR terms.",
    "i" = "Available: {.val {fa_terms}}"
  ))
}

parse_fa_term <- function(term) {
  m <- regmatches(term, regexec("^(fa|rr)\\(\\s*([^,]+?)\\s*,\\s*([0-9]+)\\s*\\)(:(.+))?$", term))[[1]]
  if (length(m) == 0) cli::cli_abort("Cannot parse FA/RR term {.val {term}}.")
  rest <- m[6]
  geno_var <- if (nzchar(rest)) sub("^vm\\(\\s*([^,]+),.*$", "\\1", rest) else NA_character_
  list(type = m[2], site_var = m[3], k = as.integer(m[4]), fa_part = sub(":.*$", "", term), rest = rest, geno_var = geno_var)
}

# Long table of the parameters belonging exactly to `term`.
fa_parameter_table <- function(vc, term) {
  nm <- names(vc)
  own <- startsWith(nm, paste0(term, "!"))
  rest <- substring(nm[own], nchar(term) + 2)
  ok <- grepl("^.+!(fa[0-9]+|var)$", rest)
  tab <- data.frame(
    Group = sub("!(fa[0-9]+|var)$", "", rest[ok]),
    Param = sub("^.*!", "", rest[ok]),
    Value = unname(vc[own][ok]),
    stringsAsFactors = FALSE
  )
  dup <- duplicated(tab[, c("Group", "Param")])
  if (any(dup)) cli::cli_abort("Duplicate parameters for {.val {term}}; refusing to combine them.")
  tab
}

extract_psi <- function(vc, par, sites, spec, psi_term) {
  if (!is.null(psi_term)) {
    prefix <- sub("^diag\\(([^)]+)\\)", "\\1", gsub("\\s+", " ", psi_term))
    keys <- paste0(prefix, "!", spec$site_var, "_", sites)
    hit <- match(keys, names(vc))
    if (anyNA(hit)) {
      cli::cli_abort("Specific variances for {sum(is.na(hit))} site{?s} not found under {.val {psi_term}} (expected e.g. {.val {keys[1]}}).")
    }
    return(stats::setNames(unname(vc[hit]), sites))
  }
  v <- par[par$Param == "var", , drop = FALSE]
  psi <- stats::setNames(v$Value[match(sites, v$Group)], sites)
  if (spec$type == "fa" && anyNA(psi)) {
    cli::cli_abort("Specific variances missing for {sum(is.na(psi))} site{?s} of the FA term.")
  }
  if (spec$type == "rr" && (anyNA(psi) || all(psi == 0))) {
    cli::cli_alert_info("RR term without {.arg psi_term}: specific variances taken as 0.")
    psi[is.na(psi)] <- 0
  }
  psi
}

# Scores (solutions) and PEVs (squared SEs) for the Comp levels of the term.
fa_scores <- function(model, term, spec) {
  if (is.na(spec$geno_var)) {
    cli::cli_warn("Term {.val {term}} has no genotype component; scores not extracted.")
    return(NULL)
  }
  cr <- tryCatch(summary(model, coef = TRUE)$coef.random, error = function(e) NULL)
  if (is.null(cr)) {
    cli::cli_warn("Random coefficients unavailable; scores and FAST indices not extracted.")
    return(NULL)
  }
  sol_col <- intersect(c("solution", "Solution"), colnames(cr))[1]
  se_col <- intersect(c("std.error", "std error"), colnames(cr))[1]
  nm <- rownames(cr)

  ids <- NULL
  scores <- pev <- NULL
  for (i in seq_len(spec$k)) {
    prefix <- paste0(spec$fa_part, "_Comp", i, ":", spec$rest, "_")
    geno <- strip_prefix(nm, prefix)
    rows <- which(!is.na(geno))
    if (length(rows) == 0) {
      cli::cli_warn("No score coefficients found for factor {i} (expected names starting {.val {prefix}}).")
      return(NULL)
    }
    if (anyDuplicated(geno[rows])) cli::cli_abort("Duplicate score coefficients for factor {i}.")
    if (is.null(ids)) {
      ids <- geno[rows]
      scores <- pev <- matrix(NA_real_, length(ids), spec$k, dimnames = list(ids, NULL))
    }
    pos <- match(geno[rows], ids)
    if (anyNA(pos)) cli::cli_abort("Factor {i} scores cover different genotypes from factor 1.")
    scores[pos, i] <- cr[rows, sol_col]
    if (!is.na(se_col)) pev[pos, i] <- cr[rows, se_col]^2
  }
  if (anyNA(scores)) cli::cli_abort("Incomplete score matrix for {.val {term}}.")

  # Under ASReml's sigma parameterisation (no explicit residual, sigma2 != 1)
  # the score coefficients are scaled by sqrt(sigma2) while the loadings in
  # the variance component table are not. Return scores with unit prior
  # variance (times g_ii) so that OP, RMSD and reliability are on the
  # loadings' scale. With an explicit residual sigma2 is 1 and this is a no-op.
  s2 <- model$sigma2
  if (is.numeric(s2) && length(s2) == 1 && is.finite(s2) && s2 > 0) {
    scores <- scores / sqrt(s2)
    pev <- pev / s2
  }
  list(scores = scores, pev = pev)
}

# Remainder of `x` after an exact (whitespace-tolerant) prefix, else NA.
strip_prefix <- function(x, prefix) {
  out <- rep(NA_character_, length(x))
  hit <- startsWith(x, prefix)
  if (any(hit)) {
    out[hit] <- substring(x[hit], nchar(prefix) + 1)
    return(out)
  }
  escaped <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", prefix)
  m <- regexpr(paste0("^", gsub("\\s+", "\\\\s*", escaped)), x)
  hit <- m == 1
  out[hit] <- substring(x[hit], attr(m, "match.length")[hit] + 1)
  out
}

# Diagonal of G for the requested genotypes (1 when no GRM is supplied).
genotype_variance <- function(grm, ids) {
  if (is.null(grm)) {
    return(NULL)
  }
  if (isTRUE(attr(grm, "INVERSE"))) {
    cli::cli_warn("{.arg grm} is an inverse; supply G itself for reliabilities. Assuming g_ii = 1.")
    return(NULL)
  }
  d <- if (is.matrix(grm) || inherits(grm, "Matrix")) {
    stats::setNames(diag(as.matrix(grm)), rownames(grm))
  } else {
    grm
  }
  if (is.null(names(d))) cli::cli_abort("{.arg grm} must carry genotype names.")
  out <- d[ids]
  if (anyNA(out)) cli::cli_warn("{sum(is.na(out))} genotype{?s} missing from {.arg grm}; their reliability is NA.")
  stats::setNames(unname(out), ids)
}

#' Calculate FAST Indices (OP, RMSD, PEV, Reliability)
#'
#' Overall performance \eqn{OP_i = \bar\lambda_1 f_{i1}} and stability
#' \eqn{RMSD_i = \sqrt{\frac{1}{p}\sum_j (\sum_{r \ge 2} \lambda_{jr} f_{ir})^2}}
#' from an SVD-rotated FA solution (Smith & Cullis 2018).
#'
#' @param fa_object An \code{fa_model} or a list with \code{loadings$rotated},
#'   \code{scores$rotated} and optionally \code{scores$pev_rotated}.
#' @param k Number of factors (default: all).
#' @param g_diag Optional named vector of \eqn{g_{ii}} (variance of a
#'   genotype's scores); defaults to 1 (independent genotypes).
#' @return Data frame with \code{Genotype}, \code{OP}, \code{PEV_OP},
#'   \code{Reliability} (\eqn{r^2_g = 1 - PEV(f_{i1}) / g_{ii}}, the squared
#'   accuracy of OP) and \code{RMSD},
#'   ordered by decreasing OP.
#' @export
calculate_fast_indices <- function(fa_object, k = NULL, g_diag = NULL) {
  rot <- fa_object$meta$rotation
  if (!is.null(rot) && rot != "svd") {
    cli::cli_abort("FAST indices require an SVD-rotated solution (got {.val {rot}}).")
  }
  L <- fa_object$loadings$rotated
  F_sc <- fa_object$scores$rotated
  PEV_rot <- fa_object$scores$pev_rotated
  if (is.null(L) || is.null(F_sc)) {
    return(NULL)
  }
  if (is.null(k)) k <- ncol(L)

  lbar1 <- mean(L[, 1])
  df <- data.frame(Genotype = rownames(F_sc), OP = F_sc[, 1] * lbar1, stringsAsFactors = FALSE)

  if (!is.null(PEV_rot) && !all(is.na(PEV_rot))) {
    g <- if (is.null(g_diag)) rep(1, nrow(F_sc)) else unname(g_diag[rownames(F_sc)])
    df$PEV_OP <- PEV_rot[, 1] * lbar1^2
    rel <- 1 - PEV_rot[, 1] / g
    n_neg <- sum(rel < 0, na.rm = TRUE)
    if (n_neg > 0) {
      cli::cli_warn("{n_neg} reliabilit{?y is/ies are} negative (PEV > genetic variance); set to 0. Check {.arg grm}.")
      rel[!is.na(rel) & rel < 0] <- 0
    }
    df$Reliability <- rel
  }

  df$RMSD <- 0
  if (k > 1) {
    I_mat <- F_sc[, 2:k, drop = FALSE] %*% t(L[, 2:k, drop = FALSE])
    df$RMSD <- sqrt(rowMeans(I_mat^2))
  }
  rownames(df) <- NULL
  df[order(df$OP, decreasing = TRUE), , drop = FALSE]
}

#' Interaction Classes (iClasses)
#'
#' Groups environments by the signs of their rotated loadings on the chosen
#' factors (Smith et al. 2021) and predicts each genotype's performance within
#' each iClass,
#' \eqn{\hat\beta_{ic} = \sum_{r \in \{1\} \cup R} \bar\lambda_{rc} f_{ir}},
#' where \eqn{\bar\lambda_{rc}} is the mean rotated loading in iClass \eqn{c}.
#'
#' @param fa_object An SVD-rotated \code{fa_model} with scores.
#' @param factors Integer vector of factors (excluding 1) defining the
#'   classes (default 2; up to \eqn{2^{|R|}} classes).
#' @return A list with \code{environments} (\code{Group}, \code{iClass}) and
#'   \code{genotypes} (long: \code{Genotype}, \code{iClass}, \code{n_env},
#'   \code{iClassOP}).
#' @references Smith, A., Norman, A., Kuchel, H. & Cullis, B. (2021). Plant
#'   variety selection using interaction classes derived from factor analytic
#'   linear mixed models. \emph{Frontiers in Plant Science} 12, 737462.
#' @export
calculate_i_classes <- function(fa_object, factors = 2) {
  if (!identical(fa_object$meta$rotation, "svd")) cli::cli_abort("iClasses require an SVD-rotated {.cls fa_model}.")
  L <- fa_object$loadings$rotated
  F_sc <- fa_object$scores$rotated
  if (is.null(F_sc)) cli::cli_abort("Scores are required.")
  k <- ncol(L)
  factors <- as.integer(factors)
  if (k < 2 || any(factors < 2) || any(factors > k)) {
    cli::cli_abort("{.arg factors} must lie in 2..{k} (model has k = {k}).")
  }

  sign_lab <- ifelse(L[, factors, drop = FALSE] >= 0, "p", "n")
  env_class <- apply(sign_lab, 1, paste, collapse = "")
  use <- c(1L, factors)

  geno <- do.call(rbind, lapply(sort(unique(env_class)), function(cl) {
    in_c <- env_class == cl
    lbar <- colMeans(L[in_c, use, drop = FALSE])
    data.frame(
      Genotype = rownames(F_sc), iClass = cl, n_env = sum(in_c),
      iClassOP = as.vector(F_sc[, use, drop = FALSE] %*% lbar),
      stringsAsFactors = FALSE
    )
  }))
  rownames(geno) <- NULL

  list(
    environments = data.frame(Group = rownames(L), iClass = env_class, stringsAsFactors = FALSE, row.names = NULL),
    genotypes = geno
  )
}

#' Rank Genotypes using a Selection Index
#'
#' Ranks genotypes by \eqn{I = OP - b \times RMSD}.
#'
#' @param fast_df Data frame with \code{OP} and \code{RMSD} columns, or an
#'   \code{fa_model}.
#' @param weight_stability Numeric. Penalty \eqn{b} on RMSD (default 1).
#'   OP and RMSD share the trait's units, so \eqn{b} is the yield a breeder
#'   would trade for one unit less crossover GxE.
#' @return The data frame sorted by \code{Index}, with a \code{Rank} column.
#' @export
rank_genotypes <- function(fast_df, weight_stability = 1) {
  if (inherits(fast_df, "fa_model")) fast_df <- fast_df$fast
  if (is.null(fast_df) || !all(c("OP", "RMSD") %in% names(fast_df))) {
    cli::cli_abort("Input must contain {.field OP} and {.field RMSD} columns.")
  }
  fast_df$Index <- fast_df$OP - weight_stability * fast_df$RMSD
  fast_df <- fast_df[order(fast_df$Index, decreasing = TRUE), , drop = FALSE]
  fast_df$Rank <- seq_len(nrow(fast_df))
  fast_df
}

#' Extract Factor Analytic Loadings and VAF
#'
#' @param fa_object An \code{fa_model}.
#' @return A list with \code{Loadings} (rotated), \code{Loadings_raw} and
#'   \code{VAF}.
#' @export
extract_fa_loadings <- function(fa_object) {
  if (!inherits(fa_object, "fa_model")) cli::cli_abort("Object must be of class {.cls fa_model}.")
  list(Loadings = fa_object$loadings$rotated, Loadings_raw = fa_object$loadings$raw, VAF = fa_object$var_comp$vaf)
}

#' Extract Factor Analytic Scores and FAST Indices
#'
#' @param fa_object An \code{fa_model}.
#' @return A list with \code{Scores} (rotated), \code{Scores_raw} and
#'   \code{FAST}.
#' @export
extract_fa_blups <- function(fa_object) {
  if (!inherits(fa_object, "fa_model")) cli::cli_abort("Object must be of class {.cls fa_model}.")
  list(Scores = fa_object$scores$rotated, Scores_raw = fa_object$scores$raw, FAST = fa_object$fast)
}

#' Summary Method for Factor Analytic Models
#'
#' @param object An \code{fa_model}.
#' @param n Integer. Number of top genotypes to report (default 5).
#' @param ... Ignored.
#' @return An object of class \code{summary.fa_model}.
#' @export
summary.fa_model <- function(object, n = 5, ...) {
  L <- object$loadings$rotated
  top_env <- data.frame(
    Factor = colnames(L),
    Group = rownames(L)[apply(abs(L), 2, which.max)],
    Loading = L[cbind(apply(abs(L), 2, which.max), seq_len(ncol(L)))],
    stringsAsFactors = FALSE
  )
  structure(
    list(
      meta = object$meta,
      var_explained = object$meta$var_explained,
      mean_vaf = mean(object$var_comp$vaf$Total_VAF, na.rm = TRUE),
      top_environments = top_env,
      top_genotypes = if (!is.null(object$fast)) utils::head(object$fast, n)
    ),
    class = "summary.fa_model"
  )
}

#' @export
print.summary.fa_model <- function(x, ...) {
  cli::cli_h1("Factor Analytic Model Summary")
  cli::cli_text("Term: {.val {x$meta$term}} | Type: {toupper(x$meta$type)} | k = {x$meta$k} | rotation: {x$meta$rotation}")
  cli::cli_h2("Variance explained (% of total genetic variance)")
  for (i in seq_along(x$var_explained)) cli::cli_text("Factor {i}: {round(x$var_explained[i], 1)}%")
  cli::cli_text("Mean total VAF across environments: {round(x$mean_vaf, 1)}%")
  cli::cli_h2("Largest |loading| per factor")
  print(x$top_environments, row.names = FALSE)
  cli::cli_h2("Top genotypes (OP)")
  if (is.null(x$top_genotypes)) cli::cli_alert_warning("No FAST indices available.") else print(x$top_genotypes, row.names = FALSE)
  invisible(x)
}

#' Print Method for Factor Analytic Models
#'
#' @param x An \code{fa_model}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.fa_model <- function(x, ...) {
  cli::cli_h1("Factor Analytic Model")
  cli::cli_text("Term: {.val {x$meta$term}} | Type: {toupper(x$meta$type)} | k = {x$meta$k} | rotation: {x$meta$rotation}")
  vaf <- x$var_comp$vaf$Total_VAF
  if (!is.null(vaf)) cli::cli_text("Mean total VAF across environments: {round(mean(vaf, na.rm = TRUE), 1)}%")
  n_g <- if (is.null(x$scores$rotated)) 0 else nrow(x$scores$rotated)
  cli::cli_text("{nrow(x$loadings$rotated)} environments, {n_g} genotypes with scores.")
  cli::cli_text("Use {.code summary(x)} for details or {.code plot(x)} for figures.")
  invisible(x)
}
