test_that("calculate_h2_from_list validation", {
    expect_error(calculate_h2_from_list("not a list"), "Input 'model_list' must be a list")

    # Valid list but invalid objects
    mock_list <- list(mdl1 = "string", mdl2 = NULL)

    # Should message skipping and return list of NAs
    suppressMessages({
        res <- calculate_h2_from_list(mock_list)
    })

    expect_true(is.vector(res))
    expect_equal(length(res), 2)
    expect_true(all(is.na(res)))
    expect_named(res, c("mdl1", "mdl2"))
})
