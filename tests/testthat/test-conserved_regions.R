test_that("conserved_regions function test", {
    data(autm, package = "GenomAutomorphism")
    test1 <- conserved_regions(autm[24:30])
    test1 <- length(test1) == 6
    expect_true(test1)
})
