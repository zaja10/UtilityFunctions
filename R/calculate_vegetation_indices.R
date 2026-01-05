#' Calculate Vegetation Spectral Indices
#'
#' Calculates vegetation indices from the 'Awesome Spectral Indices' list based
#' on columns present in the input data frame.
#'
#' @param spectra_df A data frame where rows are observations and columns are spectral bands.
#'   Expected column names include: N, R, G, B, RE1, RE2, RE3, S1, S2, PAR, A, G1, N2.
#' @param g,C1,C2,L,L_savi,gamma,alpha,sla,slb,n,lambda_N,lambda_R,lambda_G,k,c,epsilon,fdelta
#'   Parameters for specific indices (see defaults).
#'
#' @return A data frame with the original columns plus new columns for calculated indices.
#'
#' @importFrom dplyr mutate
#' @importFrom rlang parse_expr !!!
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(N = 0.8, R = 0.1, G = 0.2, B = 0.1)
#' res <- calculate_vegetation_indices(df)
#' names(res)
#' }
calculate_vegetation_indices <- function(spectra_df,
                                         g = 2.5, C1 = 6.0, C2 = 7.5, L = 1.0,
                                         L_savi = 0.5, gamma = 1.0, alpha = 0.1,
                                         sla = 1.0, slb = 0.0, n = 2.0,
                                         lambda_N = 800.0, lambda_R = 670.0,
                                         lambda_G = 550.0, k = 0.0, c = 1.0,
                                         epsilon = 0.1, fdelta = 0.0) {
    # 1. Input Validation
    if (!is.data.frame(spectra_df)) stop("spectra_df must be a data frame.")

    # 2. Define Requirements (Mapping Index -> Bands)
    reqs <- list(
        AFRI1600 = c("N", "S1"), AFRI2100 = c("N", "S2"), ARI = c("G", "RE1"),
        ARI2 = c("G", "N", "RE1"), ARVI = c("B", "N", "R"), ATSAVI = c("N", "R"),
        AVI = c("N", "R"), BCC = c("B", "G", "R"), BNDVI = c("B", "N"),
        BWDRVI = c("B", "N"), CCI = c("G1", "R"), CIG = c("G", "N"),
        CIRE = c("N", "RE1"), CVI = c("G", "N", "R"), DSI = c("N", "S1"),
        DSWI1 = c("N", "S1"), DSWI2 = c("G", "S1"), DSWI3 = c("R", "S1"),
        DSWI4 = c("G", "R"), DSWI5 = c("G", "N", "R", "S1"), DVI = c("N", "R"),
        DVIplus = c("G", "N", "R"), EBI = c("B", "G", "R"), EVI = c("B", "N", "R"),
        EVI2 = c("N", "R"), ExG = c("B", "G", "R"), ExGR = c("B", "G", "R"),
        ExR = c("G", "R"), FCVI = c("B", "G", "N", "R"), GARI = c("B", "G", "N", "R"),
        GBNDVI = c("B", "G", "N"), GCC = c("B", "G", "R"), GDVI = c("N", "R"),
        GEMI = c("N", "R"), GLI = c("B", "G", "R"), GM1 = c("G", "RE2"),
        GM2 = c("RE1", "RE2"), GNDVI = c("G", "N"), GOSAVI = c("G", "N"),
        GRNDVI = c("G", "N", "R"), GRVI = c("G", "N"), GSAVI = c("G", "N"),
        GVMI = c("N", "S2"), IAVI = c("B", "N", "R"), IKAW = c("B", "R"),
        IPVI = c("N", "R"), IRECI = c("R", "RE1", "RE2", "RE3"),
        MCARI = c("G", "R", "RE1"), MCARI1 = c("G", "N", "R"),
        MCARI2 = c("G", "N", "R"), MCARI705 = c("G", "RE1", "RE2"),
        MCARIOSAVI = c("G", "N", "R", "RE1"), MCARIOSAVI705 = c("G", "RE1", "RE2"),
        MGRVI = c("G", "R"), MNDVI = c("N", "S2"), MNLI = c("N", "R"),
        MRBVI = c("B", "R"), MSAVI = c("N", "R"), MSI = c("N", "S1"),
        MSR = c("N", "R"), MSR705 = c("RE1", "RE2"), MTCI = c("R", "RE1", "RE2"),
        MTVI1 = c("G", "N", "R"), MTVI2 = c("G", "N", "R"), ND705 = c("RE1", "RE2"),
        NDDI = c("G", "N", "R"), NDGI = c("G", "N", "R"), NDII = c("N", "S1"),
        NDMI = c("N", "S1"), NDPI = c("N", "R", "S1"), NDREI = c("N", "RE1"),
        NDVI = c("N", "R"), NDVI705 = c("RE1", "RE2"), NDYI = c("B", "G"),
        NGRDI = c("G", "R"), NIRv = c("N", "R"), NIRvH2 = c("N", "R"),
        NIRvP = c("N", "PAR", "R"), NLI = c("N", "R"), NMDI = c("N", "S1", "S2"),
        NRFIg = c("G", "S2"), NRFIr = c("R", "S2"), NormG = c("G", "N", "R"),
        NormNIR = c("G", "N", "R"), NormR = c("G", "N", "R"), OCVI = c("G", "N", "R"),
        OSAVI = c("N", "R"), PSRI = c("B", "R", "RE2"), RCC = c("B", "G", "R"),
        RDVI = c("N", "R"), REDSI = c("R", "RE1", "RE3"), RENDVI = c("RE1", "RE2"),
        RGBVI = c("B", "G", "R"), RGRI = c("G", "R"), RI = c("G", "R"),
        RVI = c("R", "RE2"), S2REP = c("R", "RE1", "RE2", "RE3"),
        SARVI = c("B", "N", "R"), SAVI = c("N", "R"), SAVI2 = c("N", "R"),
        SEVI = c("N", "R"), SI = c("B", "G", "R"), SIPI = c("A", "N", "R"),
        SLAVI = c("N", "R", "S2"), SR = c("N", "R"), SR2 = c("G", "N"),
        SR3 = c("G", "N2", "RE1"), SR555 = c("G", "RE2"), SR705 = c("RE1", "RE2"),
        SeLI = c("N2", "RE1"), TCARI = c("G", "R", "RE1"),
        TCARIOSAVI = c("G", "N", "R", "RE1"), TCARIOSAVI705 = c("G", "RE1", "RE2"),
        TCI = c("G", "R", "RE1"), TDVI = c("N", "R"), TGI = c("B", "G", "R"),
        TRRVI = c("N", "R", "RE2"), TSAVI = c("N", "R"),
        TTVI = c("N2", "RE2", "RE3"), TVI = c("N", "R"), TriVI = c("G", "N", "R"),
        VARI = c("B", "G", "R"), VARI700 = c("B", "R", "RE1"),
        VI700 = c("R", "RE1"), VIG = c("G", "R"), WDRVI = c("N", "R"),
        WDVI = c("N", "R"), mND705 = c("A", "RE1", "RE2"), mSR705 = c("A", "RE2")
    )

    # 3. Define Formulas (R Expressions)
    forms <- list(
        AFRI1600 = "(N - 0.66 * S1) / (N + 0.66 * S1)",
        AFRI2100 = "(N - 0.5 * S2) / (N + 0.5 * S2)",
        ARI = "(1 / G) - (1 / RE1)",
        ARI2 = "N * ((1 / G) - (1 / RE1))",
        ARVI = "(N - (R - gamma * (R - B))) / (N + (R - gamma * (R - B)))",
        ATSAVI = "sla * (N - sla * R - slb) / (sla * N + R - sla * slb + 0.08 * (1 + sla^2.0))",
        AVI = "sign(N*(1.0-R)*(N-R)) * abs(N*(1.0-R)*(N-R))^(1/3)",
        BCC = "B / (R + G + B)",
        BNDVI = "(N - B) / (N + B)",
        BWDRVI = "(alpha * N - B) / (alpha * N + B)",
        CCI = "(G1 - R) / (G1 + R)",
        CIG = "(N / G) - 1.0",
        CIRE = "(N / RE1) - 1",
        CVI = "(N * R) / (G^2.0)",
        DSI = "S1 / N",
        DSWI1 = "N / S1",
        DSWI2 = "S1 / G",
        DSWI3 = "S1 / R",
        DSWI4 = "G / R",
        DSWI5 = "(N + G) / (S1 + R)",
        DVI = "N - R",
        DVIplus = "((lambda_N - lambda_R) / (lambda_N - lambda_G)) * G + (1.0 - ((lambda_N - lambda_R) / (lambda_N - lambda_G))) * N - R",
        EBI = "(R + G + B) / ((G / B) * (R - B + epsilon))",
        EVI = "g * (N - R) / (N + C1 * R - C2 * B + L)",
        EVI2 = "g * (N - R) / (N + 2.4 * R + L)",
        ExG = "2 * G - R - B",
        ExGR = "(2.0 * G - R - B) - (1.3 * R - G)",
        ExR = "1.3 * R - G",
        FCVI = "N - ((R + G + B) / 3.0)",
        GARI = "(N - (G - (B - R))) / (N - (G + (B - R)))",
        GBNDVI = "(N - (G + B)) / (N + (G + B))",
        GCC = "G / (R + G + B)",
        GDVI = "((N^n) - (R^n)) / ((N^n) + (R^n))",
        GEMI = "((2.0 * ((N^2.0) - (R^2.0)) + 1.5 * N + 0.5 * R) / (N + R + 0.5)) * (1.0 - 0.25 * ((2.0 * ((N^2.0) - (R^2)) + 1.5 * N + 0.5 * R) / (N + R + 0.5))) - ((R - 0.125) / (1 - R))",
        GLI = "(2.0 * G - R - B) / (2.0 * G + R + B)",
        GM1 = "RE2 / G",
        GM2 = "RE2 / RE1",
        GNDVI = "(N - G) / (N + G)",
        GOSAVI = "(N - G) / (N + G + 0.16)",
        GRNDVI = "(N - (G + R)) / (N + (G + R))",
        GRVI = "N / G",
        GSAVI = "(1.0 + L_savi) * (N - G) / (N + G + L_savi)",
        GVMI = "((N + 0.1) - (S2 + 0.02)) / ((N + 0.1) + (S2 + 0.02))",
        IAVI = "(N - (R - gamma * (B - R))) / (N + (R - gamma * (B - R)))",
        IKAW = "(R - B) / (R + B)",
        IPVI = "N / (N + R)",
        IRECI = "(RE3 - R) / (RE1 / RE2)",
        MCARI = "((RE1 - R) - 0.2 * (RE1 - G)) * (RE1 / R)",
        MCARI1 = "1.2 * (2.5 * (N - R) - 1.3 * (N - G))",
        MCARI2 = "(1.5 * (2.5 * (N - R) - 1.3 * (N - G))) / ((((2.0 * N + 1)^2) - (6.0 * N - 5 * (R^0.5)) - 0.5)^0.5)",
        MCARI705 = "((RE2 - RE1) - 0.2 * (RE2 - G)) * (RE2 / RE1)",
        MCARIOSAVI = "(((RE1 - R) - 0.2 * (RE1 - G)) * (RE1 / R)) / (1.16 * (N - R) / (N + R + 0.16))",
        MCARIOSAVI705 = "(((RE2 - RE1) - 0.2 * (RE2 - G)) * (RE2 / RE1)) / (1.16 * (RE2 - RE1) / (RE2 + RE1 + 0.16))",
        MGRVI = "(G^2.0 - R^2.0) / (G^2.0 + R^2.0)",
        MNDVI = "(N - S2) / (N + S2)",
        MNLI = "(1 + L_savi) * ((N^2) - R) / ((N^2) + R + L_savi)",
        MRBVI = "(R^2.0 - B^2.0) / (R^2.0 + B^2.0)",
        MSAVI = "0.5 * (2.0 * N + 1 - (((2 * N + 1)^2) - 8 * (N - R))^0.5)",
        MSI = "S1 / N",
        MSR = "(N / R - 1) / ((N / R + 1)^0.5)",
        MSR705 = "(RE2 / RE1 - 1) / ((RE2 / RE1 + 1)^0.5)",
        MTCI = "(RE2 - RE1) / (RE1 - R)",
        MTVI1 = "1.2 * (1.2 * (N - G) - 2.5 * (R - G))",
        MTVI2 = "(1.5 * (1.2 * (N - G) - 2.5 * (R - G))) / ((((2.0 * N + 1)^2) - (6.0 * N - 5 * (R^0.5)) - 0.5)^0.5)",
        ND705 = "(RE2 - RE1) / (RE2 + RE1)",
        NDDI = "(((N - R) / (N + R)) - ((G - N) / (G + N))) / (((N - R) / (N + R)) + ((G - N) / (G + N)))",
        NDGI = "(((lambda_N - lambda_R) / (lambda_N - lambda_G)) * G + (1.0 - ((lambda_N - lambda_R) / (lambda_N - lambda_G))) * N - R) / (((lambda_N - lambda_R) / (lambda_N - lambda_G)) * G + (1.0 - ((lambda_N - lambda_R) / (lambda_N - lambda_G))) * N + R)",
        NDII = "(N - S1) / (N + S1)",
        NDMI = "(N - S1) / (N + S1)",
        NDPI = "(N - (alpha * R + (1.0 - alpha) * S1)) / (N + (alpha * R + (1.0 - alpha) * S1))",
        NDREI = "(N - RE1) / (N + RE1)",
        NDVI = "(N - R) / (N + R)",
        NDVI705 = "(RE2 - RE1) / (RE2 + RE1)",
        NDYI = "(G - B) / (G + B)",
        NGRDI = "(G - R) / (G + R)",
        NIRv = "((N - R) / (N + R)) * N",
        NIRvH2 = "N - R - k * (lambda_N - lambda_R)",
        NIRvP = "((N - R) / (N + R)) * N * PAR",
        NLI = "((N^2) - R) / ((N^2) + R)",
        NMDI = "(N - (S1 - S2)) / (N + (S1 - S2))",
        NRFIg = "(G - S2) / (G + S2)",
        NRFIr = "(R - S2) / (R + S2)",
        NormG = "G / (N + G + R)",
        NormNIR = "N / (N + G + R)",
        NormR = "R / (N + G + R)",
        OCVI = "(N / G) * (R / G)^c",
        OSAVI = "(N - R) / (N + R + 0.16)",
        PSRI = "(R - B) / RE2",
        RCC = "R / (R + G + B)",
        RDVI = "(N - R) / ((N + R)^0.5)",
        REDSI = "((705.0 - 665.0) * (RE3 - R) - (783.0 - 665.0) * (RE1 - R)) / (2.0 * R)",
        RENDVI = "(RE2 - RE1) / (RE2 + RE1)",
        RGBVI = "(G^2.0 - B * R) / (G^2.0 + B * R)",
        RGRI = "R / G",
        RI = "(R - G) / (R + G)",
        RVI = "RE2 / R",
        S2REP = "705.0 + 35.0 * ((((RE3 + R) / 2.0) - RE1) / (RE2 - RE1))",
        SARVI = "(1 + L_savi) * (N - (R - (R - B))) / (N + (R - (R - B)) + L_savi)",
        SAVI = "(1.0 + L_savi) * (N - R) / (N + R + L_savi)",
        SAVI2 = "N / (R + (slb / sla))",
        SEVI = "(N / R) + fdelta * (1.0 / R)",
        SI = "sign((1.0-B)*(1.0-G)*(1.0-R)) * abs((1.0-B)*(1.0-G)*(1.0-R))^(1/3)",
        SIPI = "(N - A) / (N - R)",
        SLAVI = "N / (R + S2)",
        SR = "N / R",
        SR2 = "N / G",
        SR3 = "N2 / (G * RE1)",
        SR555 = "RE2 / G",
        SR705 = "RE2 / RE1",
        SeLI = "(N2 - RE1) / (N2 + RE1)",
        TCARI = "3 * ((RE1 - R) - 0.2 * (RE1 - G) * (RE1 / R))",
        TCARIOSAVI = "(3 * ((RE1 - R) - 0.2 * (RE1 - G) * (RE1 / R))) / (1.16 * (N - R) / (N + R + 0.16))",
        TCARIOSAVI705 = "(3 * ((RE2 - RE1) - 0.2 * (RE2 - G) * (RE2 / RE1))) / (1.16 * (RE2 - RE1) / (RE2 + RE1 + 0.16))",
        TCI = "1.2 * (RE1 - G) - 1.5 * (R - G) * (RE1 / R)^0.5",
        TDVI = "1.5 * ((N - R) / ((N^2.0 + R + 0.5)^0.5))",
        TGI = "-0.5 * (190 * (R - G) - 120 * (R - B))",
        TRRVI = "((RE2 - R) / (RE2 + R)) / (((N - R) / (N + R)) + 1.0)",
        TSAVI = "sla * (N - sla * R - slb) / (sla * N + R - sla * slb)",
        TTVI = "0.5 * ((865.0 - 740.0) * (RE3 - RE2) - (N2 - RE2) * (783.0 - 740))",
        TVI = "(((N - R) / (N + R)) + 0.5)^0.5",
        TriVI = "0.5 * (120 * (N - G) - 200 * (R - G))",
        VARI = "(G - R) / (G + R - B)",
        VARI700 = "(RE1 - 1.7 * R + 0.7 * B) / (RE1 + 1.3 * R - 1.3 * B)",
        VI700 = "(RE1 - R) / (RE1 + R)",
        VIG = "(G - R) / (G + R)",
        WDRVI = "(alpha * N - R) / (alpha * N + R)",
        WDVI = "N - sla * R",
        mND705 = "(RE2 - RE1) / (RE2 + RE1 - A)",
        mSR705 = "(RE2 - A) / (RE2 + A)"
    )

    # 4. Filter Calculations
    avail_bands <- names(spectra_df)
    todo <- list()

    for (idx in names(reqs)) {
        if (all(reqs[[idx]] %in% avail_bands)) {
            # Parse the string into an expression
            todo[[idx]] <- rlang::parse_expr(forms[[idx]])
        }
    }

    if (length(todo) == 0) {
        warning("No indices could be calculated. Missing required bands?")
        return(spectra_df)
    }

    cli::cli_alert_info("Calculating {length(todo)} indices: {paste(names(todo), collapse=', ')}")

    # 5. Apply
    # Splice the expressions into mutate
    res <- dplyr::mutate(spectra_df, !!!todo)

    return(res)
}
