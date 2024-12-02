test_that("automorphism_prob works", {
        data("autby_coef", package = "GenomAutomorphism")
        post_prob <- automorphism_prob(autby_coef[1:10])
        expect_true(round(sum(post_prob$Posterior_Probability[1:10]), 3) == 1)
    }
)
