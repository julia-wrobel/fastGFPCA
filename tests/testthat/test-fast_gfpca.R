test_that("Y dataframe has correct structure", {
  expect_error(fast_gfpca(mtcars))
})


test_that("Function works for variety of binwidths", {
  df_gfpca <- sim_gfpca(N = 50, J = 200, case = 2)$df_gfpca
  gfpca_mod <- fast_gfpca(df_gfpca, overlap = FALSE, binwidth = 4, family = binomial)
  gfpca_mod <- fast_gfpca(df_gfpca, overlap = TRUE, binwidth = 4, family = binomial)

  expect_message(fast_gfpca(df_gfpca, overlap = FALSE, binwidth = 30, family = binomial),
                 "binwidth should be no more than J/10. Converting to a new binwidth of 20.")
  expect_message(fast_gfpca(df_gfpca, overlap = FALSE, binwidth = 5, family = binomial))

})

test_that("Function works for different links and EF families",{
  df_gfpca <- sim_gfpca(N = 50, J = 100, case = 1, mu = TRUE, family = "poisson")$df_gfpca
  gfpca_mod <- fast_gfpca(df_gfpca, overlap = TRUE, binwidth = 10,
                          npc = 6,
                          family = poisson(link = "log"),
                          periodicity = TRUE)

  expect_equal(gfpca_mod$family$link, "log")
  expect_equal(gfpca_mod$family$family, "poisson")

  gfpca_mod <- fast_gfpca(df_gfpca, overlap = FALSE, binwidth = 10,
                          npc = 3,
                          family = poisson(link = "log"))

  expect_equal(gfpca_mod$family$link, "log")
  expect_equal(gfpca_mod$family$family, "poisson")

})

