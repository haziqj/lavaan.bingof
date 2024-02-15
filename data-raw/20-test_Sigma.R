library(tidyverse)
library(survey)
library(SDaA)
letters <- c(letters, sapply(letters, function(x) paste0(x, letters)))

teacher.sample <- merge(teachers, teachmi, by = c("dist", "school"))
teacher.sample$teacher <- 1:nrow(teacher.sample)
teacher.sample$N <- NA
teacher.sample$N[teacher.sample$dist == "sm/me"] <- 66
teacher.sample$N[teacher.sample$dist == "large"] <- 245

# Two-stage stratified cluster sampling specification:
#
# - Population is teachers
# - [PSU] Strata are school districts (school = sm/me or large) SRS
# - [SSU] Teachers within schools SRS
teacher.design <-
  svydesign(id = ~school + teacher, data = teacher.sample,
            strata = ~dist + NULL, fpc = ~N + popteach)
summary(teacher.design)

# Try new population ----------------------------------------------------------
pop <-
  tibble(
    type = c("A", "B", "C"),
    fpc1 = c(400, 1000, 600),
    fpc2 = map(fpc1, \(x) sample(30, size = x, replace = TRUE)),
    mod = c(10, 0, -10)
  ) |>
  unnest(fpc2) |>
  group_by(type) |>
  mutate(
    sch_id = row_number(),
    tch_id = map(fpc2, \(x) letters[1:x])
  ) |>
  unnest(tch_id) |>
  ungroup() |>
  mutate(mod2 = map_dbl(mod, \(x) {
    if (x == 10) return(rnorm(1, mean = 10, sd = 1))
    if (x == -10) return(rnorm(1, mean = -10, sd = 0.5))
    if (x == 0) return(0)
  })) |>
  mutate(y = 0 + rnorm(n(), mean = 50, sd = 10)) |>
  select(-mod)

# Generate sample
res <- list()
for (i in 1:100) {

  schools_selected <- tibble(
    type = c("A", "B", "C"),
    sch_id = list(
      sort(sample(400, size = 50, replace = FALSE)),
      sort(sample(1000, size = 50, replace = FALSE)),
      sort(sample(600, size = 50, replace = FALSE))
    )
  ) |>
    unnest(sch_id)

  samp <-
    left_join(schools_selected, pop, by = join_by(type, sch_id)) |>
    slice_sample(n = 10, by = c(type, sch_id)) |>
    arrange(type, sch_id, tch_id)

  svy <- svydesign(id = ~sch_id + tch_id, data = samp, nest = TRUE,
                   strata = ~type + NULL, fpc = ~fpc1 + fpc2)
  res[[i]] <- c(
  svy = svymean(~ y, svy),
  # mean(pop$y)
  srs = mean(samp$y)
  )
}
res <- bind_rows(res)
pivot_longer(res, everything()) |>
  group_by(name) |>
  summarise(mean = mean(value),
            sd = sd(value))

