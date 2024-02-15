get_samp <- function(pop) {
  schools_selected <- tibble(
    type = c("A", "B", "C"),
    sch_id = list(
      sort(sample(400, size = 50, replace = FALSE)),
      sort(sample(1000, size = 50, replace = FALSE)),
      sort(sample(600, size = 50, replace = FALSE))
    )
  ) |>
    unnest(sch_id)

  out <-
    left_join(schools_selected, pop, by = join_by(type, sch_id)) |>
    slice_sample(n = 10, by = c(type, sch_id)) |>
    arrange(type, sch_id, tch_id)

  svy <- svydesign(id = ~sch_id + tch_id, data = out, nest = TRUE,
                   strata = ~type + NULL, fpc = ~fpc1 + fpc2, deff = TRUE)

  out$wt <- weights(svy)
  out$wt.norm <- out$wt / sum(out$wt) * nrow(out)

  attr(out, "svy") <- svy
  out
}
