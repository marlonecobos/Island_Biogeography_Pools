# function to check if matrix is positive definite
is_pos_def <- function(x, tol = 1e-8) {
  eigenvalues <- eigen(x, only.values = TRUE)$values
  eigenvalues <- ifelse(eigenvalues < tol, 0, eigenvalues)
  return(all(eigenvalues > 0))
}

# function to check if matrix is non-singular
non_sing <- function(x, tol = 1e-8) {
  return(det(x) > tol)
}


# function to assess sensitivity of ellipsoids
ellipsoid_sensitivity <- function(data, variable_columns, level = 0.95, 
                                  iterations = 4, train_proportion = 0.7, 
                                  tol = 1e-60) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data.frame.")
  }
  if (!is.numeric(variable_columns) && !is.character(variable_columns)) {
    stop("Argument 'variable_columns' must be numeric indices or character names.")
  }
  if (is.numeric(variable_columns) && any(variable_columns > ncol(data) | variable_columns < 1)) {
    stop("Numeric 'variable_columns' are out of bounds.")
  }
  if (is.character(variable_columns) && !all(variable_columns %in% names(data))) {
    stop("One or more 'variable_columns' not found in data.")
  }
  if (!is.numeric(level) || level <= 0 || level >= 1) {
    stop("Argument 'level' must be a number between 0 and 1.")
  }
  if (!is.numeric(iterations) || iterations < 1) {
    stop("Argument 'iterations' must be an integer greater than 0.")
  }
  if (!is.numeric(train_proportion) || train_proportion <= 0 || train_proportion >= 1) {
    stop("Argument 'train_proportion' must be a number between 0 and 1.")
  }
  
  iteration_results <- list()
  
  for (i in 1:iterations) {
    # 1. Split data into training and testing sets
    n_rows <- nrow(data)
    train_size <- floor(train_proportion * n_rows)
    train_indices <- sample(1:n_rows, size = train_size)
    
    train_data <- data[train_indices, variable_columns, drop = FALSE]
    test_data <- data[-train_indices, variable_columns, drop = FALSE]
    
    # Handle case with only one variable column
    if (is.vector(train_data)) {
      train_data <- as.data.frame(train_data)
      test_data <- as.data.frame(test_data)
      names(train_data) <- names(test_data) <- if(is.character(variable_columns)) variable_columns else names(data)[variable_columns]
    }
    
    # 2. In a loop use training data to get centroid and covar_matrix
    centroid <- colMeans(train_data)
    covar_matrix <- cov(train_data)
    
    # 3. Id which remaining data is inside or outside the limit.
    maha_dist <- stats::mahalanobis(x = test_data,
                                    center = centroid,
                                    cov = covar_matrix,
                                    tol = tol)
    
    cutoff <- stats::qchisq(level, df = length(variable_columns))
    
    inside_count <- sum(maha_dist <= cutoff)
    outside_count <- sum(maha_dist > cutoff)
    
    total_test <- nrow(test_data)
    sensitivity <- inside_count / total_test
    omission <- outside_count / total_test
    
    # 4. For each iteration return how many records not used in training are inside and outside the limit.
    iteration_results[[i]] <- data.frame(iteration = i,
                                         inside = inside_count,
                                         outside = outside_count,
                                         sensitivity = sensitivity,
                                         omission = omission)
  }
  
  all_iterations <- do.call(rbind, iteration_results)
  
  # 5. return mean and SD of results from each iteration.
  summary_stats <- data.frame(
    mean_inside = mean(all_iterations$inside),
    sd_inside = sd(all_iterations$inside),
    mean_outside = mean(all_iterations$outside),
    sd_outside = sd(all_iterations$outside),
    mean_sensitivity = mean(all_iterations$sensitivity),
    sd_sensitivity = sd(all_iterations$sensitivity),
    mean_omission = mean(all_iterations$omission),
    sd_omission = sd(all_iterations$omission)
  )
  
  return(list(iteration_results = all_iterations,
              summary = summary_stats))
}
