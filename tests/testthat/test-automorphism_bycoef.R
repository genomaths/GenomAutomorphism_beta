test_that("automorphism_bycoef function test", {
    data(autm, package = "GenomAutomorphism")
    test1 <- automorphism_bycoef(x = autm[1:10])
    test1 <- test1$autm[9] == 49
    expect_true(test1)
})
