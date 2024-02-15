library(tidyverse)
letters <- c(letters, sapply(letters, function(x) paste0(x, letters)))

# Create population ------------------------------------------------------------
set.seed(123)
pop <-
  tibble(
    type = c("A", "B", "C"),
    fpc1 = c(400, 1000, 600),
    fpc2 = map(fpc1, \(x) { 5 + sample(45, size = x, replace = TRUE)})
  ) |>
  unnest(fpc2) |>
  group_by(type) |>
  mutate(
    sch_id = row_number(),
    tch_id = map(fpc2, \(x) letters[1:x])
  ) |>
  unnest(tch_id) |>
  ungroup()

pop <- bind_cols(pop, gen_data_bin(model_no = 5, n = nrow(pop)))

dat <- get_samp(pop)
svy <- attr(dat, "svy")
fit <- lavaan::sem(model = txt_mod(5), data = dat, estimator = "PML",
                   std.lv = TRUE)

# Fit the model ----------------------------------------------------------------
res1 <- run_ligof_sims_alt(1, nsim = 50)
res3 <- run_ligof_sims_alt(3, nsim = 50)
res5 <- run_ligof_sims_alt(5, nsim = 50)

# Try 2 ------------------------------------------------------------------------
set.seed(123)
n_individuals <- 10000
n_items <- 5

# Generate latent variable Z from N(0, 1)
latent_variable <- rnorm(n_individuals)
loadings <- c(0.8, 1.0, 1.2, 0.9, 1.1)
latent_influenced <- sweep(matrix(latent_variable, ncol = n_items, nrow = n_individuals), 2, loadings, `*`)


# Generate binary items based on thresholds
thresholds <- qnorm(c(0.05, 0.1, 0.4, 0.3, 0.2))
binary_items <- sweep(latent_influenced, 2, thresholds, FUN = ">")
binary_items <- +binary_items # Convert logical to binary

# Assign individuals to strata (e.g., based on quartiles of some correlated variable)
# Here, we simply use the latent variable for illustration
strata <- cut(latent_variable, breaks = quantile(latent_variable, probs = seq(0, 1, length.out = 5)), labels = FALSE, include.lowest = TRUE)

# Combine all data
population_data <- data.frame(ID = 1:n_individuals, Stratum = strata,  binary_items)

clusters_per_stratum <- 50
# Initialize a vector for cluster assignments
cluster_assignment <- numeric(n_individuals)

tmp <- NULL

# For each stratum, assign individuals to clusters
for (stratum in 1:4) {
  # Identify individuals in this stratum
  ind_in_stratum <- which(strata == stratum)

  # Determine the number of individuals in this stratum
  n_ind_stratum <- length(ind_in_stratum)

  # Generate random cluster sizes for this stratum that sum to the number of individuals in the stratum
  # This is a simplification and might need adjustment to ensure every individual is assigned
  sizes <- round(runif(clusters_per_stratum, 40, 60))
  sizes <- round(sizes / sum(sizes) * n_ind_stratum)

  # Ensure sum of sizes matches exactly the number of individuals in the stratum
  diff <- n_ind_stratum - sum(sizes)
  sizes[1] <- sizes[1] + diff

  # Assign individuals to clusters
  cluster_ids <- rep(1:clusters_per_stratum, times = sizes)

  # If sizes do not exactly match, adjust this process to ensure all individuals are assigned
  cluster_assignment[ind_in_stratum] <- cluster_ids[1:n_ind_stratum]

  tmp <- rbind(tmp, data.frame(
    Stratum = stratum,
    Cluster = 1:50,
    fpc2 = sizes
  ))
}

# Update the data frame with the correct cluster assignments
population_data$Cluster <- cluster_assignment
population_data$fpc1 <- 50
population_data <- merge(population_data, tmp, by = c("Stratum", "Cluster"))



# Get sample
sampled_clusters_per_stratum <- 10
individuals_per_cluster_sample <- 25

# Sample clusters from each stratum
sampled_data <- lapply(split(population_data, population_data$Stratum), function(stratum_data) {
  unique_clusters <- unique(stratum_data$Cluster)
  sampled_clusters <- sample(unique_clusters, sampled_clusters_per_stratum) # Sample 10 clusters from each stratum
  sampled_stratum_data <- stratum_data[stratum_data$Cluster %in% sampled_clusters, ]

  # Sample individuals within sampled clusters
  sampled_individuals <- do.call(rbind, lapply(split(sampled_stratum_data, sampled_stratum_data$Cluster), function(cluster_data) {
    if(nrow(cluster_data) > individuals_per_cluster_sample) {
      sample_n(cluster_data, individuals_per_cluster_sample)
    } else {
      cluster_data
    }
  }))

  return(sampled_individuals)
})

# Combine sampled data back into a single data frame
final_sample_data <- do.call(rbind, sampled_data)
svy <- svydesign(id = ~Cluster + ID, strata = ~Stratum, data = final_sample_data, fpc = ~fpc1 + fpc2, nest = TRUE)

dat <-
  as_tibble(final_sample_data) |>
  mutate(across(starts_with("X"), ordered)) |>
  rename(y1 = X1, y2 = X2, y3 = X3, y4 = X4, y5 = X5) |>
  mutate(wt = weights(svy))

fit <- lavaan::sem(model = txt_mod(1), data = dat, estimator = "PML",
                   std.lv = TRUE, sampling.weights = "wt")
