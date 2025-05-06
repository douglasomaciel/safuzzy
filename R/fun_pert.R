# --- Funções de pertinência ---
zmf <- function(X, mfParams) {
  a <- mfParams[1]; b <- mfParams[2]
  res <- sapply(X, function(x) {
    ifelse(x <= a, 1,
           ifelse(x <= (b + a) / 2, 1 - 2 * ((x - a) / (b - a))^2,
                  ifelse(x <= b, 2 * ((x - b) / (b - a))^2, 0))
    )
  })
  list(x = X, y = res)
}

smf <- function(X, mfParams) {
  a <- mfParams[1]; b <- mfParams[2]
  res <- sapply(X, function(x) {
    ifelse(x <= a, 0,
           ifelse(x <= (b + a) / 2, 2 * ((x - a) / (b - a))^2,
                  ifelse(x <= b, 1 - 2 * ((x - b) / (b - a))^2, 1))
    )
  })
  list(x = X, y = res)
}

pimf <- function(X, mfParams) {
  a <- mfParams[1]; b <- mfParams[2]; c <- mfParams[3]; d <- mfParams[4]
  res <- sapply(X, function(x) {
    ifelse(x <= a, 0,
           ifelse(x <= (a + b)/2, 2*((x - a)/(b - a))^2,
                  ifelse(x <= b, 1 - 2*((x - b)/(b - a))^2,
                         ifelse(x <= c, 1,
                                ifelse(x <= (c + d)/2, 1 - 2*((x - c)/(d - c))^2,
                                       ifelse(x <= d, 2*((x - d)/(d - c))^2, 0)))))
    )
  })
  list(x = X, y = res)
}
