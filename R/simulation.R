#' Generate R value combinations
#'
#' @param start_value Starting value for R_w (R_z will be the negative)
#' @param step_size Step size for R value changes
#' @param num_steps Number of steps to take
#'
#' @return List of R_w and R_z combinations
#' @export
#'
#' @examples
#' generate_r_combinations(-0.5, 0.1, 10)
generate_r_combinations <- function(start_value, step_size, num_steps) {
  R_combinations <- list()

  for(i in 0:num_steps) {
    R_w_val <- start_value + i * step_size
    R_z_val <- -start_value - i * step_size
    R_combinations[[i+1]] <- c(R_w = R_w_val, R_z = R_z_val)
  }

  return(R_combinations)
}

#' Compare algorithm with different R values
#'
#' @param data Dataset to use for simulation
#' @param x_var Name of X variable
#' @param y_var Name of Y variable
#' @param z_var Name of Z variable
#' @param w_var Name of W variable
#' @param covariates Optional vector of covariate variable names
#' @param start_value Starting value for R_w (R_z will be the negative)
#' @param step_size Step size for R value changes
#' @param num_steps Number of steps to take
#' @param num_iterations Number of bootstrap iterations per R combination
#' @param parallel Whether to use parallel processing (default: TRUE)
#' @param num_cores Number of cores to use for parallel processing (default: NULL, uses detectCores())
#'
#' @return List of results for each R combination
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read.csv("your_data.csv")
#' results <- compare_algorithm_with_r(data, "X_var", "Y_var", "Z_var", "W_var",
#'                                    start_value = -0.5,
#'                                    step_size = 0.1, num_steps = 10)
#' }
compare_algorithm_with_r <- function(data, x_var, y_var, z_var, w_var,
                                     covariates = NULL,
                                     start_value = -0.5, step_size = 0.1,
                                     num_steps = 10, num_iterations = 500,
                                     parallel = TRUE, num_cores = NULL) {
  # Prepare data
  prepared_data <- prepare_data(data, x_var, y_var, z_var, w_var, covariates)

  # Generate R combinations
  R_combinations <- generate_r_combinations(start_value, step_size, num_steps)

  # Set up parallel processing if requested
  if (parallel) {
    if (is.null(num_cores)) {
      num_cores <- parallel::detectCores()
    }

    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)

    # Export necessary objects to the cluster
    parallel::clusterExport(
      cl,
      c("bi_tsls_estimation", "R_combinations", "prepared_data"),
      envir = environment()
    )
  }

  # Prepare list for storing results
  all_results <- list()

  # For each combination of R_w and R_z
  for (r_idx in seq_along(R_combinations)) {
    R_w <- R_combinations[[r_idx]]["R_w"]
    R_z <- R_combinations[[r_idx]]["R_z"]

    if (parallel) {
      # Initialize results for this combination using foreach
      results <- foreach::foreach(
        i = 1:num_iterations,
        .combine = 'cbind',
        .packages = c("MASS", "stats")
      ) %dopar% {
        # Bootstrap resampling
        set.seed(123 + i)
        bootstrap_indices <- sample(1:nrow(prepared_data), size = nrow(prepared_data), replace = TRUE)
        bootstrap_data <- prepared_data[bootstrap_indices, ]

        # Apply algorithm to the bootstrap sample
        bi_tsls_estimation(bootstrap_data, R_w, R_z)
      }
    } else {
      # Sequential processing
      results <- matrix(0, nrow = 2, ncol = num_iterations)

      for (i in 1:num_iterations) {
        set.seed(123 + i)
        bootstrap_indices <- sample(1:nrow(prepared_data), size = nrow(prepared_data), replace = TRUE)
        bootstrap_data <- prepared_data[bootstrap_indices, ]

        # Apply algorithm to the bootstrap sample
        results[, i] <- bi_tsls_estimation(bootstrap_data, R_w, R_z)
      }
    }

    all_results[[r_idx]] <- results
  }

  # Stop parallel cluster if used
  if (parallel) {
    parallel::stopCluster(cl)
  }

  return(all_results)
}

#' Process simulation results into a data frame
#'
#' @param results List of results from compare_algorithm_with_r
#' @param R_combinations List of R combinations used
#'
#' @return Data frame with processed results
#' @export
process_simulation_results <- function(results, R_combinations) {
  df <- data.frame()

  for (r_idx in seq_along(R_combinations)) {
    current_data <- results[[r_idx]]

    temp_df <- data.frame(
      R_w = R_combinations[[r_idx]]["R_w"],
      R_z = R_combinations[[r_idx]]["R_z"],
      R_idx = r_idx - 1,  # Index 0 to n
      beta_XY = current_data[1, ],
      beta_YX = current_data[2, ]
    )

    df <- rbind(df, temp_df)
  }

  # Clean data by removing extreme outliers
  df_cleaned <- df %>%
    dplyr::group_by(R_w, R_z) %>%
    dplyr::mutate(
      beta_XY = remove_extreme_outliers(beta_XY),
      beta_YX = remove_extreme_outliers(beta_YX)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::ungroup()

  return(df_cleaned)
}

#' Remove extreme outliers from a vector
#'
#' @param x Numeric vector
#' @param probs Probability quantiles to use as cutoffs (default: c(0.05, 0.95))
#'
#' @return Vector with extreme values replaced by NA
#' @keywords internal
remove_extreme_outliers <- function(x, probs = c(0.05, 0.95)) {
  qnt <- stats::quantile(x, probs = probs, na.rm = TRUE)
  y <- x
  y[x < qnt[1]] <- NA
  y[x > qnt[2]] <- NA
  return(y)
}

#' Run a complete Bi-TSLS simulation analysis
#'
#' @param data Dataset to use
#' @param x_var Name of X variable
#' @param y_var Name of Y variable
#' @param z_var Name of Z variable
#' @param w_var Name of W variable
#' @param covariates Optional vector of covariate variable names
#' @param start_value Starting value for R_w (R_z will be the negative)
#' @param step_size Step size for R value changes
#' @param num_steps Number of steps to take
#' @param num_iterations Number of bootstrap iterations per R combination
#' @param parallel Whether to use parallel processing
#' @param plot Whether to generate plots
#'
#' @return List containing results and optionally plots
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read.csv("your_data.csv")
#' analysis <- run_bi_tsls_analysis(
#'   data, "X_var", "Y_var", "Z_var", "W_var",
#'   start_value = -0.5, step_size = 0.1, num_steps = 10
#' )
#' }
run_bi_tsls_analysis <- function(data, x_var, y_var, z_var, w_var,
                                 covariates = NULL, start_value = -0.5,
                                 step_size = 0.1, num_steps = 10,
                                 num_iterations = 500, parallel = TRUE,
                                 plot = TRUE) {
  # Generate R combinations
  R_combinations <- generate_r_combinations(start_value, step_size, num_steps)

  # Run the simulation
  results <- compare_algorithm_with_r(
    data, x_var, y_var, z_var, w_var,
    covariates = covariates,
    start_value = start_value,
    step_size = step_size,
    num_steps = num_steps,
    num_iterations = num_iterations,
    parallel = parallel
  )

  # Process the results
  df_cleaned <- process_simulation_results(results, R_combinations)

  # Create output list
  output <- list(
    results = results,
    processed_data = df_cleaned,
    R_combinations = R_combinations
  )

  # Generate plots if requested
  if (plot) {
    xy_plot <- plot_bi_tsls_strength_XY(df_cleaned)
    yx_plot <- plot_bi_tsls_strength_YX(df_cleaned)

    output$plots <- list(
      xy_plot = xy_plot,
      yx_plot = yx_plot
    )
  }

  return(output)
}
