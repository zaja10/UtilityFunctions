test_that("makeNm converts factors/characters to numeric", {
    vec_char <- c("1", "2.5", "3")
    vec_fac <- factor(c("10", "20"))

    expect_equal(makeNm(vec_char), c(1, 2.5, 3))
    expect_equal(makeNm(vec_fac), c(10, 20))
    expect_warning(makeNm(c("A", "1")), "NAs introduced by coercion")
})
