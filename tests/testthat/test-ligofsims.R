test_that("Sim wrapper works", {
  sink("/dev/null")
  res <- run_ligof_sims(
    model_no = 1,
    nsim = 1,
    samp_size = 500,
    samp = "srs",
    simtype = "type1",
    ncores = 1,
    starting_seed = 197
  )
  sink()
  expect_snapshot(res)
})
