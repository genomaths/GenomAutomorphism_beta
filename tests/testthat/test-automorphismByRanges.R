test_that("automorphismByRanges function test", {
    data(autm, package = "GenomAutomorphism")
    test1 <- automorphismByRanges(x = autm[c(24:35)])
    test1 <- test1$cube[2] == "TGCA"
    expect_true(test1)
})
