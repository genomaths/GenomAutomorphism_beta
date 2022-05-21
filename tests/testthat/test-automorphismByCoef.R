test_that("automorphismByCoef function test", {
    data(autm, package = "GenomAutomorphism")
    test1 <- automorphismByCoef(x = autm[1:5])
    test1 <- test1$autm[5] == 27
    expect_true(test1)
})
