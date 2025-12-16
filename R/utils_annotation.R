#' Annotate FA Model with Group Metadata
#'
#' Attaches metadata (e.g., "Frost", "Irrigated", "Exclude") to the groups (Sites/Traits)
#' in a fitted Factor Analytic model. These annotations automatically influence
#' downstream plots (color/shape) and summaries.
#'
#' @param object An object of class \code{fa_model}.
#' @param ... Named arguments defining the annotation.
#'        Format: \code{"SiteName" = "Tag"}.
#'        Example: \code{SiteA = "Frost", SiteB = "Frost"}.
#' @param df Optional dataframe with columns "Group" and "Tag".
#'        Overrides \code{...} if provided.
#' @param default_tag Character. Default label for non-specified groups. Default "Standard".
#'
#' @return The modified \code{fa_model} object with updated metadata.
#' @import cli
#' @export
annotate_model <- function(object, ..., df = NULL, default_tag = "Standard") {
    if (!inherits(object, "fa_model")) cli::cli_abort("Object must be of class {.cls fa_model}.")

    # 1. Get existing groups from the model
    # Note: Loadings are the source of truth for group names
    model_groups <- rownames(object$loadings$rotated)

    # 2. Parse Inputs
    if (!is.null(df)) {
        if (!all(c("Group", "Tag") %in% names(df))) {
            cli::cli_abort("Annotation dataframe must have columns {.val Group} and {.val Tag}.")
        }
        input_tags <- df$Tag
        names(input_tags) <- df$Group
    } else {
        # Capture dots
        dots <- list(...)
        if (length(dots) == 0) cli::cli_abort("No annotations provided.")
        input_tags <- unlist(dots)
    }

    # 3. Validation
    # Check if user provided names that don't exist in the model
    unknowns <- setdiff(names(input_tags), model_groups)
    if (length(unknowns) > 0) {
        cli::cli_warn("Ignored annotations for unknown groups: {.val {unknowns}}")
    }

    # 4. Construct Metadata Table
    meta_df <- data.frame(
        Group = model_groups,
        Tag = default_tag,
        stringsAsFactors = FALSE
    )

    # Apply updates
    matches <- intersect(names(input_tags), model_groups)
    if (length(matches) > 0) {
        meta_df$Tag[match(matches, meta_df$Group)] <- input_tags[matches]
    }

    # 5. Attach to Object
    object$meta$annotation <- meta_df

    # Report
    counts <- table(meta_df$Tag)
    count_str <- paste(names(counts), counts, sep = ": ", collapse = ", ")
    cli::cli_alert_success("Model annotated: {count_str}")

    return(object)
}
