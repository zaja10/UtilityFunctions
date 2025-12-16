test_that("Universal parser handles RR and interactions", {
    # 1. Mock RR Model Structure (Reduced Rank)
    # RR has loadings but NO internal !var. Specific variances are usually external (e.g. diag)

    mock_rr <- list(
        varcomp = matrix(c(
            # Loadings for Factor 1 (Site A, B)
            0.5, 0.6,
            # Specific variances (diag structure)
            0.1, 0.2
        ), ncol = 1, dimnames = list(c(
            "rr(Site, 1)!SiteA!rr_1", "rr(Site, 1)!SiteB!rr_1",
            "diag(Site)!SiteA", "diag(Site)!SiteB"
        ), "component")),
        coefficients = list(random = matrix(0, nrow = 0, ncol = 1)) # No scores for simplicity
    )
    class(mock_rr) <- "asreml"

    # Test RR parsing
    res_rr <- fit_fa_model(mock_rr, classify = "rr(Site, 1):Genotype", psi_term = "diag(Site)")

    expect_equal(res_rr$meta$type, "Reduced Rank (RR)")
    expect_equal(nrow(res_rr$loadings$raw), 2)
    expect_equal(res_rr$var_comp$psi$Psi, c(0.1, 0.2))

    # 2. Mock Multi-Trait Model
    # fa(Trait, 1):Genotype

    mock_mt <- list(
        varcomp = matrix(c(
            0.8, 0.9, # Loadings Trait1, Trait2
            0.05, 0.05 # Variances
        ), ncol = 1, dimnames = list(c(
            "fa(Trait, 1)!Trait1!var", "fa(Trait, 1)!Trait2!var",
            "fa(Trait, 1)!Trait1!fa_1", "fa(Trait, 1)!Trait2!fa_1"
        ), "component")),
        coefficients = list(random = matrix(0, nrow = 0, ncol = 1))
    )
    class(mock_mt) <- "asreml"

    res_mt <- fit_fa_model(mock_mt, classify = "fa(Trait, 1):Genotype")

    expect_equal(res_mt$meta$group, "Trait")
    expect_equal(rownames(res_mt$loadings$raw), c("Trait1", "Trait2"))

    # 3. Mock Interaction Model (vm)
    # fa(Site,1):vm(Genotype, G)

    mock_vm <- list(
        varcomp = matrix(c(0.5, 0.1, 0.5), ncol = 1, dimnames = list(c(
            "fa(Site, 1)!SiteA!fa_1", "fa(Site, 1)!SiteA!var", "fa(Site, 1)!SiteB!fa_1" # Incomplete for brevity
        ), "component")),
        coefficients = list(random = matrix(0, nrow = 0, ncol = 1))
    )
    class(mock_vm) <- "asreml"

    # Should run without error and extract correct genotype label
    # (Though we mock empty coefs, so scores will be NULL, checking metadata)
    # We test the regex logic by passing a mock coef with a VM name

    mock_vm$coefficients$random <- matrix(1, nrow = 1, dimnames = list("fa(Site, 1)_Comp_1:vm(Genotype, G)_G1", NULL))

    res_vm <- fit_fa_model(mock_vm, classify = "fa(Site, 1):vm(Genotype, G)")
    expect_equal(res_vm$meta$genotype, "Genotype")
    # Scores should have Genotype = G1
    if (!is.null(res_vm$scores)) {
        expect_equal(rownames(res_vm$scores$raw), "G1")
    }
})
