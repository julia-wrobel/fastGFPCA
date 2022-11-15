test_that("binomial returns binary outcome", {
  N = 200
  J = 100
  dat_sim <- sim_gfpca(N = N, J = J, case = 1, family = "binomial")
  df = dat_sim$df_gfpca
  expect_equal(sort(unique(df$value)), c(0,1))

  expect_equal(dim(dat_sim$scores), c(N, 4))
  expect_equal(dim(dat_sim$phi), c(J, 4))
})

test_that("Poisson returns count outcome", {
  N = 100
  J = 200
  dat_sim <- sim_gfpca(N = N, J = J, case = 2, family = "poisson")
  df = dat_sim$df_gfpca

  expect_true(all(unique(df$value) %% 1 == 0))
  expect_equal(dim(dat_sim$scores), c(N, 4))
  expect_equal(dim(dat_sim$phi), c(J, 4))
})



test_that("Gaussian works", {
  N = 300
  J = 50
  dat_sim <- sim_gfpca(N = N, J = J, case = 1, mu = TRUE, family = "gaussian")
  df = dat_sim$df_gfpca

  expect_equal(dim(dat_sim$scores), c(N, 4))
  expect_equal(dim(dat_sim$phi), c(J, 4))
})



test_that("error is thrown if family argument is incorrect", {
  expect_error(sim_gfpca(N = N, J = J, case = 2, family = "whatever"))
})
