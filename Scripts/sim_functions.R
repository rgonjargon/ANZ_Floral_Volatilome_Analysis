###########################################################################################
# Simulation function
###########################################################################################

#' Generate Spatial Autocorrelation Matrices
#'
#' Creates spatially autocorrelated matrices for two parent nodes (a and b) using
#' multivariate normal distributions with distance-based covariance.
#'
#' @param n Integer. Grid dimension (creates n x n matrix).
#' @param s_autocorr_a Numeric. Spatial autocorrelation parameter for variable a. Default is 0.
#' @param s_autocorr_b Numeric. Spatial autocorrelation parameter for variable b. Default is 0.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return A list containing two matrices:
#'   \item{a}{n x n matrix of spatially autocorrelated values for variable a}
#'   \item{b}{n x n matrix of spatially autocorrelated values for variable b}
#'
#' @examples
#' spatial <- generate_spatial_autocorr(n = 9, s_autocorr_a = 1, s_autocorr_b = 1)
#'
#' @export
generate_spatial_autocorr <- function(n, s_autocorr_a = 0, s_autocorr_b = 0, seed = 123) {
  set.seed(seed)
  
  rows <- rep(1:n, each = n)
  cols <- rep(1:n, times = n)
  distances <- as.matrix(dist(cbind(rows, cols)))
  
  a <- MASS::mvrnorm(1, mu = rep(0, n^2), Sigma = exp(-distances / s_autocorr_a))
  b <- MASS::mvrnorm(1, mu = 2 * (rows / n) - 1, Sigma = exp(-distances / s_autocorr_b))
  
  list(
    a = matrix(a, nrow = n, ncol = n),
    b = matrix(b, nrow = n, ncol = n)
  )
}

#' Generate Correlated Error Terms
#'
#' Creates correlated error matrices for three traits with customizable means,
#' standard deviations, and correlation structure.
#'
#' @param n Integer. Grid dimension (creates n x n matrices).
#' @param error_sd Numeric. Standard deviation of base error term. Default is 1.0.
#' @param error_t1_mean Numeric. Mean for trait 1 error. Default is 0.0.
#' @param error_t1_sd Numeric. Standard deviation for trait 1 error. Default is 1.0.
#' @param error_t2_mean Numeric. Mean for trait 2 error. Default is 0.5.
#' @param error_t2_sd Numeric. Standard deviation for trait 2 error. Default is 1.0.
#' @param error_t3_mean Numeric. Mean for trait 3 error. Default is -0.5.
#' @param error_t3_sd Numeric. Standard deviation for trait 3 error. Default is 1.0.
#' @param error_t3_sign Numeric. Sign multiplier for trait 3 base error. Default is -1.0.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return A list containing three n x n error matrices:
#'   \item{error_t1}{Error matrix for trait 1}
#'   \item{error_t2}{Error matrix for trait 2}
#'   \item{error_t3}{Error matrix for trait 3}
#'
#' @examples
#' errors <- generate_errors(n = 9, error_sd = 1.5)
#'
#' @export
generate_errors <- function(n, 
                            error_sd = 1.0,
                            error_t1_mean = 0.0, error_t1_sd = 1.0,
                            error_t2_mean = 0.5, error_t2_sd = 1.0,
                            error_t3_mean = -0.5, error_t3_sd = 1.0,
                            error_t3_sign = -1.0,
                            seed = 123) {
  set.seed(seed)
  
  error <- rnorm(n^2, mean = 0, sd = error_sd)
  error_t1 <- error + rnorm(n^2, mean = error_t1_mean, sd = error_t1_sd)
  error_t2 <- error + rnorm(n^2, mean = error_t2_mean, sd = error_t2_sd)
  error_t3 <- error_t3_sign * error + rnorm(n^2, mean = error_t3_mean, sd = error_t3_sd)
  
  list(
    error_t1 = matrix(error_t1, nrow = n, ncol = n),
    error_t2 = matrix(error_t2, nrow = n, ncol = n),
    error_t3 = matrix(error_t3, nrow = n, ncol = n)
  )
}

#' Generate Site Assignment Matrix
#'
#' Creates a matrix assigning spatial positions to experimental sites.
#' Supports 1, 2, or 4 site configurations.
#'
#' @param n Integer. Grid dimension (creates n x n matrix).
#' @param n_sites Integer. Number of sites (1, 2, or 4). Default is 4.
#'
#' @return An n x n matrix with integer site assignments.
#'
#' @details
#' Site configurations:
#' \itemize{
#'   \item 1 site: All positions assigned to site 1
#'   \item 2 sites: Left half = site 1, right half = site 2
#'   \item 4 sites: Four quadrants (NW=1, NE=2, SW=3, SE=4)
#' }
#'
#' @examples
#' site_matrix <- generate_site_matrix(n = 10, n_sites = 4)
#'
#' @export
generate_site_matrix <- function(n, n_sites = 4) {
  site_matrix <- matrix(0, nrow = n, ncol = n)
  
  if (n_sites == 1) {
    site_matrix[] <- 1
  } else if (n_sites == 2) {
    mid_col <- n / 2
    site_matrix[, 1:mid_col] <- 1
    site_matrix[, (mid_col + 1):n] <- 2
  } else if (n_sites == 4) {
    mid_row <- n / 2
    mid_col <- n / 2
    site_matrix[1:mid_row, 1:mid_col] <- 1
    site_matrix[1:mid_row, (mid_col + 1):n] <- 2
    site_matrix[(mid_row + 1):n, 1:mid_col] <- 3
    site_matrix[(mid_row + 1):n, (mid_col + 1):n] <- 4
  } else {
    stop("n_sites must be 1, 2, or 4")
  }
  
  site_matrix
}

#' Generate PID (Genotype) Assignment Matrix
#'
#' Creates a randomized matrix assigning PIDs (genotype identifiers) and checks
#' to spatial positions in a field trial.
#'
#' @param n Integer. Grid dimension (creates n x n matrix).
#' @param n_pid Integer. Number of unique PIDs. Default is 11.
#' @param n_chks Integer. Number of check plots. Default is 4.
#' @param custom_pid_matrix Matrix or NULL. Custom PID matrix to use instead of generating.
#'   Must be n x n. Default is NULL.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return An n x n matrix with randomized PID assignments (1 to n_pid for regular plots,
#'   n_pid + 1 for checks).
#'
#' @examples
#' pid_matrix <- generate_pid_matrix(n = 9, n_pid = 11, n_chks = 4)
#'
#' @export
generate_pid_matrix <- function(n, n_pid = 11, n_chks = 4, 
                                custom_pid_matrix = NULL, seed = 123) {
  if (!is.null(custom_pid_matrix)) {
    if (!is.matrix(custom_pid_matrix)) {
      stop("custom_pid_matrix must be a matrix")
    }
    if (nrow(custom_pid_matrix) != n || ncol(custom_pid_matrix) != n) {
      stop(paste0("custom_pid_matrix dimensions must be ", n, "x", n))
    }
    return(custom_pid_matrix)
  }
  
  set.seed(seed)
  
  n_regular <- n^2 - n_chks
  pid_matrix <- matrix(0, nrow = n, ncol = n)
  
  if (n_regular > 0) {
    pid_matrix[1:n_regular] <- rep(1:n_pid, length.out = n_regular)
  }
  
  if (n_chks > 0) {
    pid_matrix[(n_regular + 1):n^2] <- n_pid + 1
  }
  
  matrix(sample(pid_matrix), nrow = n, ncol = n)
}

#' Generate Breeding Value Matrix
#'
#' Creates a matrix of breeding values for each spatial position based on PID,
#' site effects, and random variation.
#'
#' @param n Integer. Grid dimension (creates n x n matrix).
#' @param pid_matrix Matrix. n x n matrix of PID assignments.
#' @param site_matrix Matrix. n x n matrix of site assignments.
#' @param n_pid Integer. Number of unique PIDs. Default is 11.
#' @param eff_pid Numeric. Effect size for PID differences. Default is 1.0.
#' @param eff_site Numeric. Effect size for site differences. Default is 0.3.
#' @param per_chks Numeric. Percentile (0-1) for check mean relative to PIDs. Default is 0.9.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return An n x n matrix of breeding values with site-specific variation.
#'
#' @examples
#' bv <- generate_bv_matrix(n = 9, pid_matrix, site_matrix)
#'
#' @export
generate_bv_matrix <- function(n, pid_matrix, site_matrix, 
                               n_pid = 11, eff_pid = 1.0, 
                               eff_site = 0.3, per_chks = 0.9,
                               seed = 123) {
  set.seed(seed)
  
  # Calculate PID means
  pid_means <- seq(0, 0 + eff_pid * n_pid, length.out = n_pid)
  
  # Calculate check mean at specified percentile
  chks_mean <- quantile(pid_means, probs = per_chks)
  
  # Extended means vector including check
  pid_means_extended <- c(pid_means, chks_mean)
  
  # Generate site-level means
  n_sites <- max(site_matrix)
  site_means <- seq(-eff_site/2, eff_site/2, length.out = n_sites)
  
  # Generate BV values with site-level variation
  bv <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n_pid + 1)) {
    mask <- pid_matrix == i
    bv[mask] <- rnorm(sum(mask), 
                      mean = pid_means_extended[i] + site_means[site_matrix[mask]], 
                      sd = 1)
  }
  
  bv
}

#' Simulate Field Trial Data
#'
#' Main function to simulate complete field trial data with spatial structure,
#' multiple traits, and a specified causal DAG.
#'
#' @param n Integer. Grid dimension (creates n x n spatial layout). Default is 9.
#' @param n_pid Integer. Number of unique PIDs (genotypes). Default is 11.
#' @param eff_pid Numeric. Effect size for PID differences. Default is 1.0.
#' @param n_chks Integer. Number of check plots. Default is 4.
#' @param per_chks Numeric. Percentile for check mean (0-1). Default is 0.9.
#' @param n_sites Integer. Number of sites (1, 2, or 4). Default is 4.
#' @param eff_site Numeric. Effect size for site differences. Default is 0.3.
#' @param custom_pid_matrix Matrix or NULL. Custom PID matrix. Default is NULL.
#' @param s_autocorr_a Numeric. Spatial autocorrelation for variable a. Default is 0.
#' @param s_autocorr_b Numeric. Spatial autocorrelation for variable b. Default is 0.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param error_sd Numeric. Standard deviation of base error. Default is 1.0.
#' @param error_t1_mean Numeric. Mean for trait 1 error. Default is 0.0.
#' @param error_t1_sd Numeric. SD for trait 1 error. Default is 1.0.
#' @param error_t2_mean Numeric. Mean for trait 2 error. Default is 0.5.
#' @param error_t2_sd Numeric. SD for trait 2 error. Default is 1.0.
#' @param error_t3_mean Numeric. Mean for trait 3 error. Default is -0.5.
#' @param error_t3_sd Numeric. SD for trait 3 error. Default is 1.0.
#' @param error_t3_sign Numeric. Sign multiplier for trait 3. Default is -1.0.
#' @param x_intercept Numeric. Intercept for variable x. Default is 0.0.
#' @param x_beta_a Numeric. Effect of a on x. Default is 0.0.
#' @param x_beta_b Numeric. Effect of b on x. Default is 0.0.
#' @param x_sd Numeric. SD for x noise. Default is 1.0.
#' @param t1_intercept Numeric. Intercept for trait 1. Default is 0.0.
#' @param t1_beta_x Numeric. Effect of x on trait 1. Default is 0.7.
#' @param t1_beta_a Numeric. Effect of a on trait 1. Default is 0.2.
#' @param t1_beta_b Numeric. Effect of b on trait 1. Default is 0.4.
#' @param t2_intercept Numeric. Intercept for trait 2. Default is 0.0.
#' @param t2_beta_x Numeric. Effect of x on trait 2. Default is 0.2.
#' @param t2_beta_a Numeric. Effect of a on trait 2. Default is 0.1.
#' @param t2_beta_b Numeric. Effect of b on trait 2. Default is 0.8.
#' @param t3_intercept Numeric. Intercept for trait 3. Default is 0.0.
#' @param t3_beta_x Numeric. Effect of x on trait 3. Default is -0.7.
#' @param t3_beta_a Numeric. Effect of a on trait 3. Default is -0.2.
#' @param t3_beta_b Numeric. Effect of b on trait 3. Default is -0.7.
#'
#' @return A list containing:
#'   \item{array}{3D array (n x n x 9) with all variables}
#'   \item{df_long}{Long-format tibble of the data}
#'   \item{df_wide}{Wide-format tibble with factors for modeling}
#'
#' @examples
#' sim <- sim_trial(n = 9, n_pid = 11, seed = 123)
#' print_matrices(sim$array)
#'
#' @export
sim_trial <- function(n = 9, 
                      n_pid = 11,
                      eff_pid = 1.0,
                      n_chks = 4,     
                      per_chks = 0.9,
                      n_sites = 4,
                      eff_site = 0.3,
                      custom_pid_matrix = NULL,
                      s_autocorr_a = 0, 
                      s_autocorr_b = 0,
                      seed = 123,
                      error_sd = 1.0,
                      error_t1_mean = 0.0, error_t1_sd = 1.0,
                      error_t2_mean = 0.5, error_t2_sd = 1.0,
                      error_t3_mean = -0.5, error_t3_sd = 1.0,
                      error_t3_sign = -1.0,
                      x_intercept = 0.0, x_beta_a = 0.0, x_beta_b = 0.0, x_sd = 1.0,
                      t1_intercept = 0.0, t1_beta_x = 0.7, t1_beta_a = 0.2, t1_beta_b = 0.4,
                      t2_intercept = 0.0, t2_beta_x = 0.2, t2_beta_a = 0.1, t2_beta_b = 0.8,
                      t3_intercept = 0.0, t3_beta_x = -0.7, t3_beta_a = -0.2, t3_beta_b = -0.7) {
  
  set.seed(seed)
  
  # Generate all component matrices
  spatial <- generate_spatial_autocorr(n, s_autocorr_a, s_autocorr_b, seed)
  errors <- generate_errors(n, error_sd, error_t1_mean, error_t1_sd, 
                            error_t2_mean, error_t2_sd, error_t3_mean, 
                            error_t3_sd, error_t3_sign, seed)
  site_matrix <- generate_site_matrix(n, n_sites)
  pid_matrix <- generate_pid_matrix(n, n_pid, n_chks, custom_pid_matrix, seed)
  bv <- generate_bv_matrix(n, pid_matrix, site_matrix, n_pid, 
                           eff_pid, eff_site, per_chks, seed)
  
  # Create array
  data_array <- array(
    dim = c(n, n, 9),
    dimnames = list(
      row = 1:n, col = 1:n,
      variable = c("pid", "site", "bv", "t1", "t2", "t3", "x", "a", "b")
    )
  )
  
  # Populate array with matrices
  data_array[, , "pid"] <- pid_matrix
  data_array[, , "site"] <- site_matrix
  data_array[, , "bv"] <- bv
  data_array[, , "a"] <- spatial$a
  data_array[, , "b"] <- spatial$b
  
  # Apply DAG
  set.seed(seed)
  data_array[, , "x"] <- x_intercept + bv + 
    x_beta_a * spatial$a + x_beta_b * spatial$b + 
    rnorm(n^2, mean = 0, sd = x_sd)
  
  data_array[, , "t1"] <- t1_intercept + 
    t1_beta_x * data_array[, , "x"] + 
    t1_beta_a * spatial$a + t1_beta_b * spatial$b + 
    errors$error_t1
  
  data_array[, , "t2"] <- t2_intercept + 
    t2_beta_x * data_array[, , "x"] + 
    t2_beta_a * spatial$a + t2_beta_b * spatial$b + 
    errors$error_t2
  
  data_array[, , "t3"] <- t3_intercept + 
    t3_beta_x * data_array[, , "x"] + 
    t3_beta_a * spatial$a + t3_beta_b * spatial$b - 
    errors$error_t3
  
  # Round values
  for (var in c("bv", "t1", "t2", "t3", "x", "a", "b")) {
    data_array[, , var] <- round(data_array[, , var], 2)
  }
  
  # Convert to dataframes
  df <- data_array |> 
    as.data.frame.table() |> 
    dplyr::as_tibble() |>
    dplyr::rename(value = "Freq") |>
    dplyr::mutate(
      row = as.numeric(row),
      col = as.numeric(col)
    )
  
  df_mod <- df |> 
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::mutate(
      pid = ifelse(pid == n_pid + 1, "chks", paste0("pid_", pid)),
      pid = as.factor(pid),
      site = as.factor(site)
    )
  
  list(
    array = data_array,
    df_long = df,
    df_wide = df_mod
  )
}

###########################################################################################
# Print matrices function
###########################################################################################

#' Print Matrices from Simulation Array
#'
#' Displays selected matrices from a simulation result array with optional
#' truncation for large grids.
#'
#' @param data_array Array. 3D array from sim_trial output.
#' @param matrices Character vector. Names of matrices to print. Default includes all.
#' @param max_rows Integer. Maximum rows to display. Default is 10.
#' @param max_cols Integer. Maximum columns to display. Default is 10.
#'
#' @return NULL (prints to console).
#'
#' @examples
#' sim <- sim_trial(n = 9)
#' print_matrices(sim$array, matrices = c("pid", "t1", "t2"))
#'
#' @export
print_matrices <- function(data_array, 
                           matrices = c("pid", "site", "bv", "x", "t1", "t2", "t3", "a", "b"),
                           max_rows = 10,
                           max_cols = 10) {
  n <- dim(data_array)[1]
  available_vars <- dimnames(data_array)$variable
  
  # Check that requested matrices exist
  invalid_vars <- setdiff(matrices, available_vars)
  if (length(invalid_vars) > 0) {
    warning("The following matrices are not available: ", paste(invalid_vars, collapse = ", "))
  }
  
  valid_vars <- intersect(matrices, available_vars)
  
  # Print each requested matrix
  for (var in valid_vars) {
    cat("\n", var, "\n")
    print(data_array[1:min(max_rows, n), 1:min(max_cols, n), var])
  }
}

###########################################################################################
# Plotting functions
###########################################################################################

#' Plot Simulated Trial Data and DAG
#'
#' Creates visualization of simulated spatial data and/or the causal DAG structure.
#'
#' @param sim_result List. Output from sim_trial function.
#' @param plot_data Logical. Whether to plot the spatial data. Default is TRUE.
#' @param plot_dag Logical. Whether to plot the causal DAG. Default is TRUE.
#' @param standardize Logical. Whether to standardize values for plotting. Default is TRUE.
#'
#' @return NULL (displays plots).
#'
#' @details
#' Data plot shows spatial heatmaps of all variables (except PID) in a faceted layout.
#' DAG plot visualizes the causal structure with color-coded nodes.
#'
#' @examples
#' sim <- sim_trial(n = 9)
#' plot_sim_trial(sim, plot_data = TRUE, plot_dag = TRUE)
#'
#' @export
plot_sim_trial <- function(sim_result, 
                           plot_data = TRUE, 
                           plot_dag = TRUE,
                           standardize = TRUE) {
  # Extract data from simulation result
  df <- sim_result$df_long
  
  # Create data plot if requested
  if (plot_data) {
    # Prepare data based on standardization option
    plot_df <- df |> 
      filter(!(variable %in% c("pid"))) |>
      group_by(variable)
    
    if (standardize) {
      plot_df <- plot_df |>
        mutate(plot_value = (value - mean(value)) / sd(value))
      fill_label <- "Z-score"
      title_suffix <- "(Standardized)"
    } else {
      plot_df <- plot_df |>
        mutate(plot_value = value)
      fill_label <- "Value"
      title_suffix <- "(Raw Values)"
    }
    
    data_plot <- ggplot(plot_df,
                        aes(x = col, y = row, fill = plot_value)) +
      geom_tile() +
      facet_wrap(~variable, nrow = 2, ncol = 4) +
      scale_fill_viridis_c() +
      theme_minimal() +
      labs(title = paste("Simulated Spatial Data", title_suffix),
           x = "Column",
           y = "Row",
           fill = fill_label)
    print(data_plot)
  }
  
  # Create DAG plot if requested
  if (plot_dag) {
    basic_dag <- dagify(
      t1 ~ x + a + b,
      t2 ~ x + a + b,
      t3 ~ x + a + b,
      x ~ a + b,
      exposure = "x",
      outcome = c("t1", "t2", "t3"),
      coords = list(
        x = c(x = 1, t1 = 3, t2 = 3, t3 = 3, a = 2, b = 2),
        y = c(x = 1, t1 = 2, t2 = 1, t3 = 0, a = 2, b = 0)
      )
    )
    
    basic_dag_plot <- basic_dag |> 
      tidy_dagitty() |> 
      mutate(var_type = case_when(
        name == "x" ~ "Exposure",
        name %in% c("t1", "t2", "t3") ~ "Outcome",
        name == "a" ~ "Climate",
        name == "b" ~ "Soil"
      )) |> 
      mutate(var_label = case_match(name,
                                    "x" ~ "Pedigree",
                                    "t1" ~ "Trait 1",
                                    "t2" ~ "Trait 2",
                                    "t3" ~ "Trait 3",
                                    "a" ~ "Climate",
                                    "b" ~ "Soil"
      ))
    
    dag_plot <- ggplot(basic_dag_plot, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_dag_edges() +
      geom_dag_point(aes(color = var_type), size = 12) +
      geom_richtext(
        aes(label = name), fontface = "bold", color = "white",
        fill = NA, label.color = NA,
        label.padding = grid::unit(c(3, 0, 0, 0), "pt")
      ) + 
      geom_dag_text(
        data = filter(basic_dag_plot, var_type %in% c("Exposure", "Outcome")), 
        aes(label = var_label), 
        nudge_y = -0.25, color = "black", size = 11, size.unit = "pt"
      ) +
      geom_dag_text(
        data = filter(basic_dag_plot, var_type %in% c("Climate", "Soil")), 
        aes(label = var_label), 
        nudge_y = 0.26, color = "black", size = 11, size.unit = "pt"
      ) +
      scale_color_manual(
        values = c("Climate" = "#44AA99", "Exposure" = "#6699CC", 
                   "Outcome" = "#882255", "Soil" = "#999933"), 
        guide = "none"
      ) +
      theme_dag() +
      labs(title = "Causal DAG: Simulated Data Structure")
    print(dag_plot)
  }
  
  # Return invisibly so nothing prints to console
  invisible(NULL)
}

###########################################################################################
# Bayesian Model Fitting Function
###########################################################################################

#' Fit Bayesian Mixed Model with brms
#'
#' Fits a multivariate Bayesian mixed model to simulated trial data using brms.
#'
#' @param formula Formula or NULL. Custom brms formula. If NULL, formula is constructed
#'   based on covariates and site_effects arguments. Default is NULL.
#' @param model_name Character. Name for saved model file.
#' @param covariates Logical. Include spatial covariates (a, b) in model. Default is TRUE.
#' @param site_effects Logical. Include random site effects. Default is TRUE.
#'
#' @return A brmsfit object.
#'
#' @details
#' Fits multivariate model for traits t1, t2, t3 with:
#' \itemize{
#'   \item Random effects for PID with correlation structure
#'   \item Optional spatial covariates a and b
#'   \item Optional random site effects
#'   \item No residual correlation (set_rescor = FALSE)
#' }
#'
#' @examples
#' \dontrun{
#' sim <- sim_trial(n = 9)
#' model <- fit_brm(model_name = "model1", covariates = TRUE, site_effects = TRUE)
#' }
#'
#' @export
fit_brm <- function(formula = NULL, 
                    model_name,
                    covariates = TRUE,
                    site_effects = TRUE) {
  if (is.null(formula)) {
    if (covariates & site_effects) {
      mod_formula <- bf(mvbind(t1, t2, t3) ~ a + b + (1+site|q|pid), 
                        family = gaussian)
    } else if (covariates & !site_effects) {
      mod_formula <- bf(mvbind(t1, t2, t3) ~ a + b + (1|q|pid), 
                        family = gaussian)
    } else if (!covariates & site_effects) {
      mod_formula <- bf(mvbind(t1, t2, t3) ~ (1+site|q|pid), 
                        family = gaussian)
    } else {
      mod_formula <- bf(mvbind(t1, t2, t3) ~ (1|q|pid), 
                        family = gaussian)
    }
    formula <- mod_formula + set_rescor(FALSE)
  }
  
  brm(formula,
      family = gaussian(link = "identity"),
      data = sim$df_wide,
      # file = paste0(here("Analysis/Models/", model_name)),
      backend = "cmdstanr",
      # threads = threading(3),
      silent = 2,
      refresh = 500,
      chains = 4,
      cores = 4,
      warmup = 1000,
      iter = 2000)
}

###########################################################################################
# Conditional Effects Plotting Function
###########################################################################################

#' Plot PID Conditional Effects from brms Model
#'
#' Creates a coefficient plot showing the estimated random effects of PIDs
#' (genotypes) on a specified trait from a fitted brms model.
#'
#' @param model A brmsfit object from fit_brm function.
#' @param resp Character. Response variable name (e.g., "t1", "t2", "t3"). Default is "t1".
#'
#' @return A ggplot object showing PID effects with error bars.
#'
#' @details
#' The plot displays:
#' \itemize{
#'   \item Point estimates of PID effects
#'   \item Credible intervals (error bars)
#'   \item PIDs ordered by estimated effect size
#' }
#'
#' @examples
#' \dontrun{
#' sim <- sim_trial(n = 9)
#' model <- fit_brm(model_name = "model1")
#' plot_pid_effects(model, resp = "t1")
#' }
#'
#' @export
plot_pid_effects <- function(model, resp = "t1") {
  conditional_effects(model, 
                      re_formula = NULL, 
                      effects = "pid",
                      resp = resp)[[1]] |>
    mutate(pid = fct_reorder(pid, estimate__)) |>
    ggplot(aes(x = pid, y = estimate__)) +
    geom_errorbar(aes(ymin = lower__, ymax = upper__), alpha = 1) +
    geom_point(size = 2) +
    coord_flip() +
    labs(
      title = paste("Random effect of Pedigree on", resp),
      x = "Pedigree",
      y = "Estimated Trait Value"
    )
}