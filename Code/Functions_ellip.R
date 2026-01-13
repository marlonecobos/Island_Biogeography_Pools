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
