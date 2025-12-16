test_that("FA model parsing works with Base R regex", {
    # 1. Create Data
    df <- data.frame(Site = c("A", "B"), Genotype = c("G1", "G2"), Yield = c(1, 2))

    # 2. Create Mock Object
    # fit_fa_model pulls data name from call. We assume it can find 'df' in the calling environment.
    # We construct a call that references 'df'.

    mock_model <- list(
        call = call("asreml", data = quote(df)),
        varcomp = matrix(0, nrow = 6, ncol = 1, dimnames = list(c(
            "fa(Site, 2)!SiteA!var", "fa(Site, 2)!SiteA!fa_1", "fa(Site, 2)!SiteA!fa_2",
            "fa(Site, 2)!SiteB!var", "fa(Site, 2)!SiteB!fa_1", "fa(Site, 2)!SiteB!fa_2"
        ), "component")),
        coefficients = list(
            fixed = matrix(0, nrow = 1, dimnames = list("Genotype_A", NULL)),
            random = matrix(0, nrow = 2, dimnames = list(c("Comp_1_Genotype_A", "Comp_2_Genotype_A"), NULL))
        ),
        vcoeff = list(fixed = c(1))
    )
    class(mock_model) <- "asreml"

    # Run the function. fit_fa_model attempts to get(data_name).
    # If data_name is "df", get("df") will look in parent frames.
    # Since test_that body is a function, and we define 'df' here,
    # 'get' *should* find it if 'fit_fa_model' is not strictly locked to package namespace
    # in a way that skips the call stack.
    # However, if it fails, we know we need to mock 'get' or pass data explicitly.
    # Based on standard R scoping, this usually works for testing.

    expect_error(
        {
            fit_fa_model(mock_model, classify = "fa(Site, 2):Genotype")
        },
        NA
    )
})
