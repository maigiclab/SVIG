sanitize_filepath <- function(x) {
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)  # Replace illegal characters with "_"
  x <- trimws(x)  # Trim leading/trailing whitespace
  x <- gsub("\\s+", "_", x)  # Replace spaces with "_"
  return(x)
}