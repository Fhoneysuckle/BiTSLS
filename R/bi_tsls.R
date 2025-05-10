#' @name BiTSLS-package
#' @title Bi-directional Two-Stage Least Squares Estimation
#' @description A package for bi-directional Two-Stage Least Squares (Bi-TSLS) estimation.
#'
#' @importFrom stats as.formula coef lm na.omit quantile sd
#' @importFrom utils head
#' @importFrom dplyr group_by summarize ungroup mutate
#' @importFrom tidyr drop_na
#' @importFrom scales label_number
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores clusterExport
#' @importFrom magrittr %>%
#' @import ggplot2
NULL

#' Bi-directional Two-Stage Least Squares Estimation Function
#'
#' @description
#' The main function for bi-directional Two-Stage Least Squares (Bi-TSLS) estimation.
#' It allows for the estimation of causal effects in both directions between variables.
#'
#' @param data A data frame containing the necessary variables.
#' @param x_var Name of the X variable (primary variable in one direction).
#' @param y_var Name of the Y variable (primary variable in the other direction).
#' @param z_var Name of the Z (negative control exposure variable).
#' @param w_var Name of the W (negative control outcome variable).
#' @param covariates Optional vector of covariate variable names. If NULL, all other variables in the data will be used as covariates.
#' @param standardize_vars Logical, whether to standardize variables before analysis. Default is TRUE.
#' @param R_w Sensitivity parameter R_w. Default is 0.
#' @param R_z Sensitivity parameter R_z. Default is 0.
#'
#' @return A named vector containing the estimated causal effects in both directions (beta_xy and beta_yx).
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read.csv("your_data.csv")
#' result <- bi_tsls(data, "X_variable", "Y_variable", "Z_variable", "W_variable")
#' print(result)
#' }
bi_tsls <- function(data, x_var, y_var, z_var, w_var, covariates = NULL,
                    standardize_vars = TRUE, R_w = 0, R_z = 0) {
  # Validate inputs
  if(!all(c(x_var, y_var, z_var, w_var) %in% names(data))) {
    stop("One or more specified variables not found in the data.")
  }

  # Create a copy of the data to avoid modifying the original
  df <- data.frame(data)

  # Rename the key variables for internal processing
  names(df)[names(df) == x_var] <- "X"
  names(df)[names(df) == y_var] <- "Y"
  names(df)[names(df) == z_var] <- "Z"
  names(df)[names(df) == w_var] <- "W"

  # Handle covariate variables
  if(is.null(covariates)) {
    # Use all other variables as covariates
    all_vars <- names(df)
    main_vars <- c("X", "Y", "Z", "W")
    covariate_vars <- setdiff(all_vars, main_vars)

    # Remove any unwanted variables (like ID variables)
    if(length(covariate_vars) > 0) {
      # Keep the covariate variables
      df <- df[, c(main_vars, covariate_vars)]
    }
  } else {
    # Use specified covariate variables
    if(!all(covariates %in% names(data))) {
      stop("One or more covariate variables not found in the data.")
    }
    covariate_vars <- covariates

    # Rename covariate variables to V1, V2, ...
    for(i in seq_along(covariates)) {
      names(df)[names(df) == covariates[i]] <- paste0("V", i)
    }

    # Keep only the relevant variables
    df <- df[, c("X", "Y", "Z", "W", paste0("V", seq_along(covariates)))]
  }

  # Check for missing values
  if(any(is.na(df))) {
    warning("Missing values detected in the data. Rows with missing values will be removed.")
    df <- na.omit(df)
  }

  # Standardize variables if requested
  if(standardize_vars) {
    df$Z <- (df$Z - mean(df$Z, na.rm = TRUE)) / sd(df$Z, na.rm = TRUE)

    # Standardize covariate variables
    for(i in seq_along(covariate_vars)) {
      v_name <- paste0("V", i)
      if(v_name %in% names(df)) {
        df[[v_name]] <- (df[[v_name]] - mean(df[[v_name]], na.rm = TRUE)) /
          sd(df[[v_name]], na.rm = TRUE)
      }
    }
  }

  # Run the Bi-TSLS estimation
  result <- bi_tsls_estimation(df, R_w, R_z)

  # Return named result
  names(result) <- c("beta_xy", "beta_yx")
  return(result)
}

#' Bi-TSLS Estimation Algorithm
#'
#' @param data Prepared data frame with standardized variables
#' @param R_w Sensitivity parameter R_w
#' @param R_z Sensitivity parameter R_z
#'
#' @return Vector of estimated beta_xy and beta_yx
#' @keywords internal
bi_tsls_estimation <- function(data, R_w, R_z) {
  # Get covariate variable names
  V_names <- names(data)[grep("^V", names(data))]

  # Create V terms and interactions
  V_terms <- paste(V_names, collapse = " + ")

  # Create interaction terms
  ZV_interactions <- paste("Z:", V_names, collapse = " + ")
  WV_interactions <- paste("W:", V_names, collapse = " + ")

  # First stage regressions with interactions
  formula_W <- as.formula(paste("W ~ Z +", V_terms, "+", ZV_interactions))
  formula_Z <- as.formula(paste("Z ~ W +", V_terms, "+", WV_interactions))

  fit_W <- lm(formula_W, data = data)
  fit_Z <- lm(formula_Z, data = data)

  # Second stage regressions
  formula_X_Z <- as.formula(paste("X ~ Z + fit_W$fitted.values +", V_terms))
  formula_Y_Z <- as.formula(paste("Y ~ Z + fit_W$fitted.values +", V_terms))
  formula_X_W <- as.formula(paste("X ~ fit_Z$fitted.values + W +", V_terms))
  formula_Y_W <- as.formula(paste("Y ~ fit_Z$fitted.values + W +", V_terms))

  lm_X_Z <- lm(formula_X_Z, data = data)
  lm_Y_Z <- lm(formula_Y_Z, data = data)
  lm_X_W <- lm(formula_X_W, data = data)
  lm_Y_W <- lm(formula_Y_W, data = data)

  k1 <- coef(lm_X_W)["W"] / coef(lm_Y_W)["W"]
  k2 <- coef(lm_Y_Z)["Z"] / coef(lm_X_Z)["Z"]

  beta_xy_Bi_TSLS <- ((k2 * (1 + k1 * R_z - R_w * R_z)) - R_z) / (1 - k1 * k2 * R_w * R_z)
  beta_yx_Bi_TSLS <- ((k1 * (1 + k2 * R_w - R_w * R_z)) - R_w) / (1 - k1 * k2 * R_w * R_z)

  return(c(beta_xy_Bi_TSLS, beta_yx_Bi_TSLS))
}

#' Prepare data for Bi-TSLS analysis
#'
#' @param data Original data frame
#' @param x_var Name of X variable
#' @param y_var Name of Y variable
#' @param z_var Name of Z variable
#' @param w_var Name of W variable
#' @param covariates Optional vector of covariate variable names
#' @param standardize_vars Logical, whether to standardize Z and covariate variables
#'
#' @return Prepared data frame
#' @keywords internal
prepare_data <- function(data, x_var, y_var, z_var, w_var,
                         covariates = NULL, standardize_vars = TRUE) {
  # Create a copy of the data
  df <- data.frame(data)

  # Rename the key variables
  names(df)[names(df) == x_var] <- "X"
  names(df)[names(df) == y_var] <- "Y"
  names(df)[names(df) == z_var] <- "Z"
  names(df)[names(df) == w_var] <- "W"

  # Handle covariate variables
  if(is.null(covariates)) {
    # Use all other variables as covariates
    all_vars <- names(df)
    main_vars <- c("X", "Y", "Z", "W")
    covariates <- setdiff(all_vars, main_vars)
  }

  # Rename covariate variables
  for(i in seq_along(covariates)) {
    names(df)[names(df) == covariates[i]] <- paste0("V", i)
  }

  # Keep only relevant variables
  df <- df[, c("X", "Y", "Z", "W", paste0("V", seq_along(covariates)))]

  # Remove rows with missing values
  df <- na.omit(df)

  # Standardize variables if requested
  if(standardize_vars) {
    standardize <- function(x) {
      (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    }

    df$Z <- standardize(df$Z)

    # Standardize covariate variables
    for(i in seq_along(covariates)) {
      v_name <- paste0("V", i)
      df[[v_name]] <- standardize(df[[v_name]])
    }
  }

  return(df)
}
