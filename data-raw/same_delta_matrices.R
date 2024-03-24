# Are the Delta matrices same or no?

dat <- gen_data_bin_wt()
fit0 <- lavaan::sem(
  model = txt_mod(1),
  data = dat,
  std.lv = TRUE,
  estimator = "PML"
  # start = get_true_values(1),
  # do.fit = FALSE
)
fit1 <- lavaan::sem(
  model = txt_mod(1),
  data = dat,
  std.lv = TRUE,
  estimator = "PML",
  # start = get_true_values(1),
  # do.fit = FALSE,
  sampling.weights = "wt"
)

D0 <- get_Delta_mats(fit0)
D1 <- get_Delta_mats(fit1)

D0$Delta_til |> round(3)
D1$Delta_til |> round(3)

# Check the H matrix
get_H_inv_mat <- function(.lavobject) {
  list2env(extract_lavaan_info(.lavobject), environment())

  lavaan:::lav_model_information(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavcache = lavcache,
    lavimplied = lavimplied,
    lavh1 = NULL,
    lavoptions = lavoptions,
    extra = FALSE,
    augmented = TRUE,
    inverted = TRUE,
    use.ginv = FALSE
  )
}
get_H_inv_mat(fit0)
get_H_inv_mat(fit1)
