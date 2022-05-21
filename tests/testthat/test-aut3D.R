test_that("aut3D function test", {
    data(aln)
    autms <- aut3D(
        seq = aln,
        start = 25,
        end = 30,
        verbose = FALSE
    )
    test1 <- autms$autm[2] == "2,1,1"
    expect_true(test1)
})
