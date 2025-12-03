#' Calculate Vegetation Spectral Indices
#'
#' This function calculates all possible vegetation indices from the 'Awesome
#' Spectral Indices' list based on the columns present in the input data frame.
#' It requires the 'dplyr' and 'rlang' packages.
#'
#' @param spectra_df A data frame where each row is an observation and
#'   each column is a spectral band. Column names must match the expected
#'   band names (e.g., N, R, G, B, RE1, etc.).
#'   Possible band names are: A, B, G, G1, N, N2, PAR, R, RE1, RE2, RE3, S1, S2
#' @param g Parameter for certain indices. Default is 2.5.
#' @param C1 Parameter for certain indices. Default is 6.0.
#' @param C2 Parameter for certain indices. Default is 7.5.
#' @param L Parameter for certain indices. Default is 1.0.
#' @param L_savi Parameter for certain indices. Default is 0.5.
#' @param gamma Parameter for certain indices. Default is 1.0.
#' @param alpha Parameter for certain indices. Default is 0.1.
#' @param sla Parameter for certain indices. Default is 1.0.
#' @param slb Parameter for certain indices. Default is 0.0.
#' @param n Parameter for certain indices. Default is 2.0.
#' @param lambda_N Parameter for certain indices. Default is 800.0.
#' @param lambda_R Parameter for certain indices. Default is 670.0.
#' @param lambda_G Parameter for certain indices. Default is 550.0.
#' @param k Parameter for certain indices. Default is 0.0.
#' @param c Parameter for certain indices. Default is 1.0.
#' @param epsilon Parameter for certain indices. Default is 0.1.
#' @param fdelta Parameter for certain indices. Default is 0.0.
#'
#' @return A data frame with the original spectral bands plus new columns
#'   for each successfully calculated spectral index.
#'
#' @importFrom dplyr mutate
#' @importFrom rlang parse_expr !!!
#' @export
#'
#' @examples
#' \dontrun{
#' # Install packages if you don't have them
#' # install.packages(c("dplyr", "rlang"))
#' library(dplyr)
#'
#' # Create dummy data with ONLY N, R, G, B bands
#' # This matches your common use case
#' spectra_data_simple <- data.frame(
#'   N = c(0.8, 0.75),
#'   R = c(0.1, 0.12),
#'   G = c(0.2, 0.22),
#'   B = c(0.1, 0.11)
#' )
#'
#' # This will only calculate indices that use N, R, G, and/or B
#' # (e.g., NDVI, GNDVI, GRVI, etc.)
#' # It will print a message: "Calculating 23 indices: ..."
#' simple_indices <- calculate_vegetation_indices(spectra_data_simple)
#' print(names(simple_indices))
#'
#' # Create more complex data with Red Edge and SWIR bands
#' spectra_data_complex <- data.frame(
#'   N = c(0.8, 0.75), R = c(0.1, 0.12), G = c(0.2, 0.22), B = c(0.1, 0.11),
#'   RE1 = c(0.3, 0.33), RE2 = c(0.4, 0.44), RE3 = c(0.5, 0.55),
#'   S1 = c(1.6, 1.65), S2 = c(2.1, 2.15), A = c(0.05, 0.06),
#'   PAR = c(2000, 2100), G1 = c(0.21, 0.23), N2 = c(0.81, 0.76)
#' )
#'
#' # This will calculate all 129 indices
#' # It will print a message: "Calculating 129 indices: ..."
#' complex_indices <- calculate_vegetation_indices(spectra_data_complex)
#' print(names(complex_indices))
#' }
calculate_vegetation_indices <- function(spectra_df, g = 2.5, C1 = 6.0, C2 = 7.5, L = 1.0, L_savi = 0.5, gamma = 1.0, alpha = 0.1, sla = 1.0, slb = 0.0, n = 2.0, lambda_N = 800.0, lambda_R = 670.0, lambda_G = 550.0, k = 0.0, c = 1.0, epsilon = 0.1, fdelta = 0.0) {

  # --- 1. Dependency Checks ---
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The 'dplyr' package is required. Please install it with install.packages('dplyr')")
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop("The 'rlang' package is required. Please install it with install.packages('rlang')")
  }

  # --- 2. Define Index Requirements and Formulas ---

  # A list mapping each index to the vector of band names it requires
  .index_band_requirements <- list(
    "AFRI1600" = c("N", "S1"),
    "AFRI2100" = c("N", "S2"),
    "ARI" = c("G", "RE1"),
    "ARI2" = c("G", "N", "RE1"),
    "ARVI" = c("B", "N", "R"),
    "ATSAVI" = c("N", "R"),
    "AVI" = c("N", "R"),
    "BCC" = c("B", "G", "R"),
    "BNDVI" = c("B", "N"),
    "BWDRVI" = c("B", "N"),
    "CCI" = c("G1", "R"),
    "CIG" = c("G", "N"),
    "CIRE" = c("N", "RE1"),
    "CVI" = c("G", "N", "R"),
    "DSI" = c("N", "S1"),
    "DSWI1" = c("N", "S1"),
    "DSWI2" = c("G", "S1"),
    "DSWI3" = c("R", "S1"),
    "DSWI4" = c("G", "R"),
    "DSWI5" = c("G", "N", "R", "S1"),
    "DVI" = c("N", "R"),
    "DVIplus" = c("G", "N", "R"),
    "EBI" = c("B", "G", "R"),
    "EVI" = c("B", "N", "R"),
    "EVI2" = c("N", "R"),
    "ExG" = c("B", "G", "R"),
    "ExGR" = c("B", "G", "R"),
    "ExR" = c("G", "R"),
    "FCVI" = c("B", "G", "N", "R"),
    "GARI" = c("B", "G", "N", "R"),
    "GBNDVI" = c("B", "G", "N"),
    "GCC" = c("B", "G", "R"),
    "GDVI" = c("N", "R"),
    "GEMI" = c("N", "R"),
    "GLI" = c("B", "G", "R"),
    "GM1" = c("G", "RE2"),
    "GM2" = c("RE1", "RE2"),
    "GNDVI" = c("G", "N"),
    "GOSAVI" = c("G", "N"),
    "GRNDVI" = c("G", "N", "R"),
    "GRVI" = c("G", "N"),
    "GSAVI" = c("G", "N"),
    "GVMI" = c("N", "S2"),
    "IAVI" = c("B", "N", "R"),
    "IKAW" = c("B", "R"),
    "IPVI" = c("N", "R"),
    "IRECI" = c("R", "RE1", "RE2", "RE3"),
    "MCARI" = c("G", "R", "RE1"),
    "MCARI1" = c("G", "N", "R"),
    "MCARI2" = c("G", "N", "R"),
    "MCARI705" = c("G", "RE1", "RE2"),
    "MCARIOSAVI" = c("G", "N", "R", "RE1"),
    "MCARIOSAVI705" = c("G", "RE1", "RE2"),
    "MGRVI" = c("G", "R"),
    "MNDVI" = c("N", "S2"),
    "MNLI" = c("N", "R"),
    "MRBVI" = c("B", "R"),
    "MSAVI" = c("N", "R"),
    "MSI" = c("N", "S1"),
    "MSR" = c("N", "R"),
    "MSR705" = c("RE1", "RE2"),
    "MTCI" = c("R", "RE1", "RE2"),
    "MTVI1" = c("G", "N", "R"),
    "MTVI2" = c("G", "N", "R"),
    "ND705" = c("RE1", "RE2"),
    "NDDI" = c("G", "N", "R"),
    "NDGI" = c("G", "N", "R"),
    "NDII" = c("N", "S1"),
    "NDMI" = c("N", "S1"),
    "NDPI" = c("N", "R", "S1"),
    "NDREI" = c("N", "RE1"),
    "NDVI" = c("N", "R"),
    "NDVI705" = c("RE1", "RE2"),
    "NDYI" = c("B", "G"),
    "NGRDI" = c("G", "R"),
    "NIRv" = c("N", "R"),
    "NIRvH2" = c("N", "R"),
    "NIRvP" = c("N", "PAR", "R"),
    "NLI" = c("N", "R"),
    "NMDI" = c("N", "S1", "S2"),
    "NRFIg" = c("G", "S2"),
    "NRFIr" = c("R", "S2"),
    "NormG" = c("G", "N", "R"),
    "NormNIR" = c("G", "N", "R"),
    "NormR" = c("G", "N", "R"),
    "OCVI" = c("G", "N", "R"),
    "OSAVI" = c("N", "R"),
    "PSRI" = c("B", "R", "RE2"),
    "RCC" = c("B", "G", "R"),
    "RDVI" = c("N", "R"),
    "REDSI" = c("R", "RE1", "RE3"),
    "RENDVI" = c("RE1", "RE2"),
    "RGBVI" = c("B", "G", "R"),
    "RGRI" = c("G", "R"),
    "RI" = c("G", "R"),
    "RVI" = c("R", "RE2"),
    "S2REP" = c("R", "RE1", "RE2", "RE3"),
    "SARVI" = c("B", "N", "R"),
    "SAVI" = c("N", "R"),
    "SAVI2" = c("N", "R"),
    "SEVI" = c("N", "R"),
    "SI" = c("B", "G", "R"),
    "SIPI" = c("A", "N", "R"),
    "SLAVI" = c("N", "R", "S2"),
    "SR" = c("N", "R"),
    "SR2" = c("G", "N"),
    "SR3" = c("G", "N2", "RE1"),
    "SR555" = c("G", "RE2"),
    "SR705" = c("RE1", "RE2"),
    "SeLI" = c("N2", "RE1"),
    "TCARI" = c("G", "R", "RE1"),
    "TCARIOSAVI" = c("G", "N", "R", "RE1"),
    "TCARIOSAVI705" = c("G", "RE1", "RE2"),
    "TCI" = c("G", "R", "RE1"),
    "TDVI" = c("N", "R"),
    "TGI" = c("B", "G", "R"),
    "TRRVI" = c("N", "R", "RE2"),
    "TSAVI" = c("N", "R"),
    "TTVI" = c("N2", "RE2", "RE3"),
    "TVI" = c("N", "R"),
    "TriVI" = c("G", "N", "R"),
    "VARI" = c("B", "G", "R"),
    "VARI700" = c("B", "R", "RE1"),
    "VI700" = c("R", "RE1"),
    "VIG" = c("G", "R"),
    "WDRVI" = c("N", "R"),
    "WDVI" = c("N", "R"),
    "mND705" = c("A", "RE1", "RE2"),
    "mSR705" = c("A", "RE2")
  )

  # A list mapping each index to its formula as an R expression
  .index_formulas <- list(
    "AFRI1600" = rlang::parse_expr("(N - 0.66 * S1) / (N + 0.66 * S1)"),
    "AFRI2100" = rlang::parse_expr("(N - 0.5 * S2) / (N + 0.5 * S2)"),
    "ARI" = rlang::parse_expr("(1 / G) - (1 / RE1)"),
    "ARI2" = rlang::parse_expr("N * ((1 / G) - (1 / RE1))"),
    "ARVI" = rlang::parse_expr("(N - (R - gamma * (R - B))) / (N + (R - gamma * (R - B)))"),
    "ATSAVI" = rlang::parse_expr("sla * (N - sla * R - slb) / (sla * N + R - sla * slb + 0.08 * (1 + sla^2.0))"),
    "AVI" = rlang::parse_expr("sign(N*(1.0-R)*(N-R)) * abs(N*(1.0-R)*(N-R))^(1/3)"),
    "BCC" = rlang::parse_expr("B / (R + G + B)"),
    "BNDVI" = rlang::parse_expr("(N - B) / (N + B)"),
    "BWDRVI" = rlang::parse_expr("(alpha * N - B) / (alpha * N + B)"),
    "CCI" = rlang::parse_expr("(G1 - R) / (G1 + R)"),
    "CIG" = rlang::parse_expr("(N / G) - 1.0"),
    "CIRE" = rlang::parse_expr("(N / RE1) - 1"),
    "CVI" = rlang::parse_expr("(N * R) / (G^2.0)"),
    "DSI" = rlang::parse_expr("S1 / N"),
    "DSWI1" = rlang::parse_expr("N / S1"),
    "DSWI2" = rlang::parse_expr("S1 / G"),
    "DSWI3" = rlang::parse_expr("S1 / R"),
    "DSWI4" = rlang::parse_expr("G / R"),
    "DSWI5" = rlang::parse_expr("(N + G) / (S1 + R)"),
    "DVI" = rlang::parse_expr("N - R"),
    "DVIplus" = rlang::parse_expr("((lambda_N - lambda_R) / (lambda_N - lambda_G)) * G + (1.0 - ((lambda_N - lambda_R) / (lambda_N - lambda_G))) * N - R"),
    "EBI" = rlang::parse_expr("(R + G + B) / ((G / B) * (R - B + epsilon))"),
    "EVI" = rlang::parse_expr("g * (N - R) / (N + C1 * R - C2 * B + L)"),
    "EVI2" = rlang::parse_expr("g * (N - R) / (N + 2.4 * R + L)"),
    "ExG" = rlang::parse_expr("2 * G - R - B"),
    "ExGR" = rlang::parse_expr("(2.0 * G - R - B) - (1.3 * R - G)"),
    "ExR" = rlang::parse_expr("1.3 * R - G"),
    "FCVI" = rlang::parse_expr("N - ((R + G + B) / 3.0)"),
    "GARI" = rlang::parse_expr("(N - (G - (B - R))) / (N - (G + (B - R)))"),
    "GBNDVI" = rlang::parse_expr("(N - (G + B)) / (N + (G + B))"),
    "GCC" = rlang::parse_expr("G / (R + G + B)"),
    "GDVI" = rlang::parse_expr("((N^n) - (R^n)) / ((N^n) + (R^n))"),
    "GEMI" = rlang::parse_expr("((2.0 * ((N^2.0) - (R^2.0)) + 1.5 * N + 0.5 * R) / (N + R + 0.5)) * (1.0 - 0.25 * ((2.0 * ((N^2.0) - (R^2)) + 1.5 * N + 0.5 * R) / (N + R + 0.5))) - ((R - 0.125) / (1 - R))"),
    "GLI" = rlang::parse_expr("(2.0 * G - R - B) / (2.0 * G + R + B)"),
    "GM1" = rlang::parse_expr("RE2 / G"),
    "GM2" = rlang::parse_expr("RE2 / RE1"),
    "GNDVI" = rlang::parse_expr("(N - G) / (N + G)"),
    "GOSAVI" = rlang::parse_expr("(N - G) / (N + G + 0.16)"),
    "GRNDVI" = rlang::parse_expr("(N - (G + R)) / (N + (G + R))"),
    "GRVI" = rlang::parse_expr("N / G"),
    "GSAVI" = rlang::parse_expr("(1.0 + L_savi) * (N - G) / (N + G + L_savi)"),
    "GVMI" = rlang::parse_expr("((N + 0.1) - (S2 + 0.02)) / ((N + 0.1) + (S2 + 0.02))"),
    "IAVI" = rlang::parse_expr("(N - (R - gamma * (B - R))) / (N + (R - gamma * (B - R)))"),
    "IKAW" = rlang::parse_expr("(R - B) / (R + B)"),
    "IPVI" = rlang::parse_expr("N / (N + R)"),
    "IRECI" = rlang::parse_expr("(RE3 - R) / (RE1 / RE2)"),
    "MCARI" = rlang::parse_expr("((RE1 - R) - 0.2 * (RE1 - G)) * (RE1 / R)"),
    "MCARI1" = rlang::parse_expr("1.2 * (2.5 * (N - R) - 1.3 * (N - G))"),
    "MCARI2" = rlang::parse_expr("(1.5 * (2.5 * (N - R) - 1.3 * (N - G))) / ((((2.0 * N + 1)^2) - (6.0 * N - 5 * (R^0.5)) - 0.5)^0.5)"),
    "MCARI705" = rlang::parse_expr("((RE2 - RE1) - 0.2 * (RE2 - G)) * (RE2 / RE1)"),
    "MCARIOSAVI" = rlang::parse_expr("(((RE1 - R) - 0.2 * (RE1 - G)) * (RE1 / R)) / (1.16 * (N - R) / (N + R + 0.16))"),
    "MCARIOSAVI705" = rlang::parse_expr("(((RE2 - RE1) - 0.2 * (RE2 - G)) * (RE2 / RE1)) / (1.16 * (RE2 - RE1) / (RE2 + RE1 + 0.16))"),
    "MGRVI" = rlang::parse_expr("(G^2.0 - R^2.0) / (G^2.0 + R^2.0)"),
    "MNDVI" = rlang::parse_expr("(N - S2) / (N + S2)"),
    "MNLI" = rlang::parse_expr("(1 + L_savi) * ((N^2) - R) / ((N^2) + R + L_savi)"),
    "MRBVI" = rlang::parse_expr("(R^2.0 - B^2.0) / (R^2.0 + B^2.0)"),
    "MSAVI" = rlang::parse_expr("0.5 * (2.0 * N + 1 - (((2 * N + 1)^2) - 8 * (N - R))^0.5)"),
    "MSI" = rlang::parse_expr("S1 / N"),
    "MSR" = rlang::parse_expr("(N / R - 1) / ((N / R + 1)^0.5)"),
    "MSR705" = rlang::parse_expr("(RE2 / RE1 - 1) / ((RE2 / RE1 + 1)^0.5)"),
    "MTCI" = rlang::parse_expr("(RE2 - RE1) / (RE1 - R)"),
    "MTVI1" = rlang::parse_expr("1.2 * (1.2 * (N - G) - 2.5 * (R - G))"),
    "MTVI2" = rlang::parse_expr("(1.5 * (1.2 * (N - G) - 2.5 * (R - G))) / ((((2.0 * N + 1)^2) - (6.0 * N - 5 * (R^0.5)) - 0.5)^0.5)"),
    "ND705" = rlang::parse_expr("(RE2 - RE1) / (RE2 + RE1)"),
    "NDDI" = rlang::parse_expr("(((N - R) / (N + R)) - ((G - N) / (G + N))) / (((N - R) / (N + R)) + ((G - N) / (G + N)))"),
    "NDGI" = rlang::parse_expr("(((lambda_N - lambda_R) / (lambda_N - lambda_G)) * G + (1.0 - ((lambda_N - lambda_R) / (lambda_N - lambda_G))) * N - R) / (((lambda_N - lambda_R) / (lambda_N - lambda_G)) * G + (1.0 - ((lambda_N - lambda_R) / (lambda_N - lambda_G))) * N + R)"),
    "NDII" = rlang::parse_expr("(N - S1) / (N + S1)"),
    "NDMI" = rlang::parse_expr("(N - S1) / (N + S1)"),
    "NDPI" = rlang::parse_expr("(N - (alpha * R + (1.0 - alpha) * S1)) / (N + (alpha * R + (1.0 - alpha) * S1))"),
    "NDREI" = rlang::parse_expr("(N - RE1) / (N + RE1)"),
    "NDVI" = rlang::parse_expr("(N - R) / (N + R)"),
    "NDVI705" = rlang::parse_expr("(RE2 - RE1) / (RE2 + RE1)"),
    "NDYI" = rlang::parse_expr("(G - B) / (G + B)"),
    "NGRDI" = rlang::parse_expr("(G - R) / (G + R)"),
    "NIRv" = rlang::parse_expr("((N - R) / (N + R)) * N"),
    "NIRvH2" = rlang::parse_expr("N - R - k * (lambda_N - lambda_R)"),
    "NIRvP" = rlang::parse_expr("((N - R) / (N + R)) * N * PAR"),
    "NLI" = rlang::parse_expr("((N^2) - R) / ((N^2) + R)"),
    "NMDI" = rlang::parse_expr("(N - (S1 - S2)) / (N + (S1 - S2))"),
    "NRFIg" = rlang::parse_expr("(G - S2) / (G + S2)"),
    "NRFIr" = rlang::parse_expr("(R - S2) / (R + S2)"),
    "NormG" = rlang::parse_expr("G / (N + G + R)"),
    "NormNIR" = rlang::parse_expr("N / (N + G + R)"),
    "NormR" = rlang::parse_expr("R / (N + G + R)"),
    "OCVI" = rlang::parse_expr("(N / G) * (R / G)^c"),
    "OSAVI" = rlang::parse_expr("(N - R) / (N + R + 0.16)"),
    "PSRI" = rlang::parse_expr("(R - B) / RE2"),
    "RCC" = rlang::parse_expr("R / (R + G + B)"),
    "RDVI" = rlang::parse_expr("(N - R) / ((N + R)^0.5)"),
    "REDSI" = rlang::parse_expr("((705.0 - 665.0) * (RE3 - R) - (783.0 - 665.0) * (RE1 - R)) / (2.0 * R)"),
    "RENDVI" = rlang::parse_expr("(RE2 - RE1) / (RE2 + RE1)"),
    "RGBVI" = rlang::parse_expr("(G^2.0 - B * R) / (G^2.0 + B * R)"),
    "RGRI" = rlang::parse_expr("R / G"),
    "RI" = rlang::parse_expr("(R - G) / (R + G)"),
    "RVI" = rlang::parse_expr("RE2 / R"),
    "S2REP" = rlang::parse_expr("705.0 + 35.0 * ((((RE3 + R) / 2.0) - RE1) / (RE2 - RE1))"),
    "SARVI" = rlang::parse_expr("(1 + L_savi) * (N - (R - (R - B))) / (N + (R - (R - B)) + L_savi)"),
    "SAVI" = rlang::parse_expr("(1.0 + L_savi) * (N - R) / (N + R + L_savi)"),
    "SAVI2" = rlang::parse_expr("N / (R + (slb / sla))"),
    "SEVI" = rlang::parse_expr("(N / R) + fdelta * (1.0 / R)"),
    "SI" = rlang::parse_expr("sign((1.0-B)*(1.0-G)*(1.0-R)) * abs((1.0-B)*(1.0-G)*(1.0-R))^(1/3)"),
    "SIPI" = rlang::parse_expr("(N - A) / (N - R)"),
    "SLAVI" = rlang::parse_expr("N / (R + S2)"),
    "SR" = rlang::parse_expr("N / R"),
    "SR2" = rlang::parse_expr("N / G"),
    "SR3" = rlang::parse_expr("N2 / (G * RE1)"),
    "SR555" = rlang::parse_expr("RE2 / G"),
    "SR705" = rlang::parse_expr("RE2 / RE1"),
    "SeLI" = rlang::parse_expr("(N2 - RE1) / (N2 + RE1)"),
    "TCARI" = rlang::parse_expr("3 * ((RE1 - R) - 0.2 * (RE1 - G) * (RE1 / R))"),
    "TCARIOSAVI" = rlang::parse_expr("(3 * ((RE1 - R) - 0.2 * (RE1 - G) * (RE1 / R))) / (1.16 * (N - R) / (N + R + 0.16))"),
    "TCARIOSAVI705" = rlang::parse_expr("(3 * ((RE2 - RE1) - 0.2 * (RE2 - G) * (RE2 / RE1))) / (1.16 * (RE2 - RE1) / (RE2 + RE1 + 0.16))"),
    "TCI" = rlang::parse_expr("1.2 * (RE1 - G) - 1.5 * (R - G) * (RE1 / R)^0.5"),
    "TDVI" = rlang::parse_expr("1.5 * ((N - R) / ((N^2.0 + R + 0.5)^0.5))"),
    "TGI" = rlang::parse_expr("-0.5 * (190 * (R - G) - 120 * (R - B))"),
    "TRRVI" = rlang::parse_expr("((RE2 - R) / (RE2 + R)) / (((N - R) / (N + R)) + 1.0)"),
    "TSAVI" = rlang::parse_expr("sla * (N - sla * R - slb) / (sla * N + R - sla * slb)"),
    "TTVI" = rlang::parse_expr("0.5 * ((865.0 - 740.0) * (RE3 - RE2) - (N2 - RE2) * (783.0 - 740))"),
    "TVI" = rlang::parse_expr("(((N - R) / (N + R)) + 0.5)^0.5"),
    "TriVI" = rlang::parse_expr("0.5 * (120 * (N - G) - 200 * (R - G))"),
    "VARI" = rlang::parse_expr("(G - R) / (G + R - B)"),
    "VARI700" = rlang::parse_expr("(RE1 - 1.7 * R + 0.7 * B) / (RE1 + 1.3 * R - 1.3 * B)"),
    "VI700" = rlang::parse_expr("(RE1 - R) / (RE1 + R)"),
    "VIG" = rlang::parse_expr("(G - R) / (G + R)"),
    "WDRVI" = rlang::parse_expr("(alpha * N - R) / (alpha * N + R)"),
    "WDVI" = rlang::parse_expr("N - sla * R"),
    "mND705" = rlang::parse_expr("(RE2 - RE1) / (RE2 + RE1 - A)"),
    "mSR705" = rlang::parse_expr("(RE2 - A) / (RE2 + A)")
  )

  # --- 3. Filter Calculations Based on Available Bands ---

  available_bands <- names(spectra_df)
  calculations_to_run <- list()
  calculated_indices_names <- c()

  # Loop through all possible indices
  for (index_name in names(.index_band_requirements)) {
    required <- .index_band_requirements[[index_name]]

    # Check if all required bands are in the input data frame
    if (all(required %in% available_bands)) {

      # If yes, add the formula to our list of calculations
      calculations_to_run[[index_name]] <- .index_formulas[[index_name]]
      calculated_indices_names <- c(calculated_indices_names, index_name)
    }
  }

  # --- 4. Execute Calculations ---

  # If no indices can be calculated, warn the user and return the original data
  if (length(calculations_to_run) == 0) {
    warning("No vegetation indices could be calculated. Check column names. ",
            "Expected names like N, R, G, B, RE1, etc.")
    return(spectra_df)
  }

  # Print a helpful message
  message(paste("Calculating", length(calculated_indices_names), "indices:",
                paste(sort(calculated_indices_names), collapse = ", ")))

  # Use dplyr::mutate with the unquote-splice operator (!!!)
  # This dynamically adds all valid calculations as new columns
  indices_df <- dplyr::mutate(
    spectra_df,
    !!!calculations_to_run
  )

  return(indices_df)
}
