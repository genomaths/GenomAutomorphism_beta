test_that("automorphism_bycoef function test", {
    data(autm, package = "GenomAutomorphism")
    test1 <- automorphism_bycoef(x = autm[1:5])
    test1 <- test1$autm[5] == 27
    expect_true(test1)
})
