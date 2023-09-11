# create_Sigma2_matrix_complex <- function(.lavobject, .svy_design,
#                                          bootstrap = FALSE, nboot = 100) {
#   list2env(extract_lavaan_info(.lavobject), environment())
#   list2env(get_uni_bi_moments(.lavobject), environment())
#
#   y_form <- paste0("~ ", paste0("y", 1:p, collapse = " + "), sep = "")
#   v <-
#     srvyr::as_survey(.svy_design) %>%
#     mutate(across(starts_with("y"), \(y) as.numeric(y) - 1))
#   ystart <- min(grep("^y[0-9]", names(v$variables))) - 1
#   idx <- combn(p, 2)
#   for (k in seq_len(ncol(idx))) {
#     i <- idx[1, k]
#     j <- idx[2, k]
#     varname <- paste0("y", i, ".", j, collapse = "")
#     yi <- v$variables[, i + ystart, drop = TRUE]
#     yj <- v$variables[, j + ystart, drop = TRUE]
#     yij <- (yi == 1) * (yj == 1)  # both positive
#     v$variables[[varname]] <- yij
#     y_form <- paste0(y_form, paste0(" + ", varname))
#   }
#
#   x <- v$variables[, -(1:ystart)]
#   xbar <- c(pidot1, pidot2)  # pi2 (model probs)
#   # xbar <- c(pdot1, pdot2)  #p2 (proportions)
#   # x <- t(t(x) - xbar)
#
#   if (isTRUE(bootstrap)) {
#     bootsvy <- svrep::as_bootstrap_design(.svy_design,
#                                           type = "Rao-Wu-Yue-Beaumont",
#                                           replicates = nboot)
#     wt <- bootsvy$repweights
#     res <- NULL
#     for (b in seq_len(ncol(wt))) {
#       res[[b]] <- cov.wt(x, wt[, b], center = FALSE)$cov
#     }
#     res <- apply(simplify2array(res), 1:2, mean)
#   } else {
#     wt <- .svy_design$allprob$wt
#     res <- cov.wt(x, wt, center = xbar)$cov
#   }
#
#   res
# }
#
# create_Sigma2_matrix_complex_old <- function(.lavobject, .svy_design) {
#   list2env(extract_lavaan_info(.lavobject), environment())
#   list2env(get_uni_bi_moments(.lavobject), environment())
#
#   # Get strata information
#   strata  <- .svy_design$strata[, 1, drop = FALSE]
#   colnames(strata) <- "strata"
#   a <- unique(strata[, 1])
#
#   # Get cluster information
#   cluster <- .svy_design$cluster
#   colnames(cluster)[1] <- "cluster"
#
#   dat <-
#     dat %>%
#     bind_cols(strata, cluster)
#
#   # Create the data set of univariate and bivariate responses (1/0)
#   tmp <- dat %>%
#     select(starts_with("y")) %>%
#     mutate(across(starts_with("y"), function(x) as.numeric(x == 2)))
#
#   id <- combn(p, 2)
#   the_names <- paste0("y", 1:p)
#   the_names <- c(the_names, paste0(the_names[id[1, ]], the_names[id[2, ]]))
#   tmp2 <- as_tibble(((tmp[, id[1, ]] == 1) + (tmp[, id[2, ]] == 1)) == 2) %>%
#     mutate(across(everything(), as.numeric)) %>% suppressWarnings()
#   colnames(tmp2) <- the_names[-(1:p)]
#   tmp <- bind_cols(dat %>% select(-starts_with("y")), tmp, tmp2)
#
#   pi2 <- c(pidot1, pidot2)
#   Sigma <- list()
#   for (i in a) {
#     # Filter the strata a
#     current_dat <- tmp %>%
#       filter(strata == i)
#
#     # Obtain the responses in this strata (i.e. according to S_ab)
#     y <- current_dat %>%
#       select(starts_with("y")) %>%
#       mutate(across(starts_with("y"), function(x) as.numeric(x == 1))) %>%
#       as.matrix()
#
#     uab <- {sweep(y, 2, pi2, "-") * current_dat$wt} %>%
#       as_tibble() %>%
#       bind_cols(current_dat %>% select(-starts_with("y")), .) %>%
#       group_by(cluster) %>%
#       summarise(across(starts_with("y"), sum), .groups = "drop") %>%
#       select(starts_with("y")) %>%
#       as.matrix()
#     na <- nrow(uab)
#
#     ubar <- apply(uab, 2, mean)
#
#     Sigma[[i]] <- crossprod(sweep(uab, 2, ubar, "-")) * na / (na - 1)
#   }
#   Reduce("+", Sigma) / N
# }
#
# create_Sigma2_matrix_complex2 <- function(.lavobject, .svy_design) {
#   list2env(extract_lavaan_info(.lavobject), environment())
#   list2env(get_uni_bi_moments(.lavobject), environment())
#
#   y_form <- paste0("~ ", paste0("y", 1:p, collapse = " + "), sep = "")
#   v <-
#     srvyr::as_survey(.svy_design) %>%
#     mutate(across(starts_with("y"), \(y) as.numeric(y) - 1))
#   ystart <- min(grep("^y[0-9]", names(v$variables))) - 1
#   idx <- combn(p, 2)
#   for (k in seq_len(ncol(idx))) {
#     i <- idx[1, k]
#     j <- idx[2, k]
#     varname <- paste0("y", i, ".", j, collapse = "")
#     yi <- v$variables[, i + ystart, drop = TRUE]
#     yj <- v$variables[, j + ystart, drop = TRUE]
#     yij <- (yi == 1) * (yj == 1)  # both positive
#     v$variables[[varname]] <- yij
#     y_form <- paste0(y_form, paste0(" + ", varname))
#   }
#
#   x <- v$variables[, -(1:ystart)]
#   survey::svyvar(x, .svy_design) %>%
#     as.matrix()  # NOT GOOD BECAUSE svyvar() not as efficient as cov.wt()
#   # wt <- 1 / .svy_design$prob
#   # res <- cov.wt(x, wt, center = c(pidot1, pidot2))$cov
#   # res
# }
#
# create_Sigma_univariate_matrix_complex <- function(.lavobject, .svy_design) {
#   list2env(extract_lavaan_info(.lavobject), environment())
#   list2env(get_uni_bi_moments(.lavobject), environment())
#
#   # Get strata information
#   strata  <- .svy_design$strata[, 1, drop = FALSE]
#   colnames(strata) <- "strata"
#   a <- unique(strata[, 1])
#
#   # Get cluster information
#   cluster <- .svy_design$cluster
#   colnames(cluster)[1] <- "cluster"
#
#   dat <-
#     dat %>%
#     bind_cols(strata, cluster)
#
#   # Create the data set of univariate and bivariate responses (1/0)
#   tmp <- dat %>%
#     select(starts_with("y")) %>%
#     mutate(across(starts_with("y"), function(x) as.numeric(x == 2)))
#   tmp <- bind_cols(dat %>% select(-starts_with("y")), tmp)
#
#   Sigma <- list()
#   for (i in a) {
#     # Filter the strata a
#     current_dat <- tmp %>%
#       filter(strata == i)
#
#     # Obtain the responses in this strata (i.e. according to S_ab)
#     y <- current_dat %>%
#       select(starts_with("y")) %>%
#       mutate(across(starts_with("y"), function(x) as.numeric(x == 1))) %>%
#       as.matrix()
#
#     uab <- {sweep(y, 2, pidot1, "-") * current_dat$wt} %>%
#       as_tibble() %>%
#       bind_cols(current_dat %>% select(-starts_with("y")), .) %>%
#       group_by(cluster) %>%
#       summarise(across(starts_with("y"), sum), .groups = "drop") %>%
#       select(starts_with("y")) %>%
#       as.matrix()
#     na <- nrow(uab)
#
#     ubar <- apply(uab, 2, mean)
#
#     Sigma[[i]] <- crossprod(sweep(uab, 2, ubar, "-")) * na / (na - 1)
#   }
#   Reduce("+", Sigma) / N
# }
#
# get_w_for_delta <- function(.lavobject, .svy_design) {
#   if (is.null(.svy_design)) return(NA)
#
#   list2env(extract_lavaan_info(.lavobject), environment())
#   list2env(get_uni_bi_moments(.lavobject), environment())
#   y_form <- paste0("~ ", paste0("y", 1:p, collapse = " + "), sep = "")
#   v <-
#     srvyr::as_survey(.svy_design) %>%
#     mutate(across(starts_with("y"), \(y) as.numeric(y) - 1))
#   ystart <- min(grep("^y[0-9]", names(v$variables))) - 1
#   idx <- combn(p, 2)
#   for (k in seq_len(ncol(idx))) {
#     i <- idx[1, k]
#     j <- idx[2, k]
#     varname <- paste0("y", i, ".", j, collapse = "")
#     yi <- v$variables[, i + ystart, drop = TRUE]
#     yj <- v$variables[, j + ystart, drop = TRUE]
#     yij <- (yi == 1) * (yj == 1)  # both positive
#     v$variables[[varname]] <- yij
#     y_form <- paste0(y_form, paste0(" + ", varname))
#   }
#   x <- v$variables[, -(1:ystart)]
#   wt <- 1 / .svy_design$prob
#   apply(x, 2, \(x) sum(x * wt ^ 2)) / apply(x, 2, \(x) sum(x * wt))
# }
#
# Pearson_test_v3 <- function(object, approx_Omega2 = FALSE, svy_design = NULL,
#                             .order = 2) {
#   .order <- match.arg(as.character(.order), c("1", "2"))
#   list2env(test_begin(object, approx_Omega2, svy_design), environment())
#
#   Xi <- diag(1 / pi2_hat)
#   X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
#
#   # Rao-Scott adjustment
#   tmp <- eigen(Omega2)
#   U <- diag(tmp$val)
#   V <- tmp$vec
#   Omegahalf <- V %*% sqrt(U) %*% t(V)
#   # same as below
#   Omegahalf <- t(chol(Omega2))
#   mat <- t(Omegahalf) %*% Xi %*% (Omegahalf)
#   delta <- eigen(mat)$values
#
#   X2 <- X2 / mean(delta)
#   df <- S
#   if (.order == "2") {
#     a_sq <- mean((delta - mean(delta)) ^ 2) / mean(delta) ^ 2
#     # a_sq <- (sd(delta) / mean(delta)) ^ 2
#     X2 <- X2 / (1 + a_sq)
#     df <- S / (1 + a_sq)
#   }
#
#   data.frame(X2 = X2, df = df, name = "PearsonRS") %>%
#     after_test(., Xi, S)
# }
