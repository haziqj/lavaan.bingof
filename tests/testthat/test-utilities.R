test_that("T2 matrix works", {
  fit <- fit_facmod_pml(1)$fit
  Sigma2 <- lavaan.bingof:::create_Sigma2_matrix(fit)

  pi <-
    create_resp_pattern(5) %>%
    mutate(across(starts_with("y"), as.numeric)) %>%
    unite("pattern", starts_with("y"), sep = "", remove = FALSE) %>%
    left_join(lavaan::lavTables(fit, dimension = 0), by = "pattern") %>%
    pull(est.prop) %>%
    suppressWarnings()
  pi[is.na(pi)] <- 0
  Sigma <- diag(pi) - tcrossprod(pi)
  T2 <- create_T2_mat(p = 5)
  Sigma2_created <- T2 %*% Sigma %*% t(T2)

  expect_equal(Sigma2_created[upper.tri(Sigma2_created)],
               Sigma2[upper.tri(Sigma2)], tolerance = 0.05)
})
