#' @keywords internal
"_PACKAGE"

#' @useDynLib UtilityFunctions, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        ".data", "Status", "Genotype", "RMSD", "OP", "Type",
        "Env1", "Env2", "Correlation",
        "Value", "Loading", "Effect", "ColorGroup", "Slope", "Factor",
        "X", "Y", "Env", "VAF", "Site", "Impact_Pct",
        "Pred_Neg", "Pred_Pos", "Class", "Value_To_Plot",
        "Method", "Index", "X_Num"
    ))
}
