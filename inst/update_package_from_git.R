# Check if 'devtools' is installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  message("devtools is not installed. Installing now.")
  install.packages("devtools")
} else {
  message("devtools is already installed.")
}

# Update package
devtools::install_github("vloubiere/vlfunction")