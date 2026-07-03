#' Rice Lines Adaptability and Stability Data Using Fuzzy Logic
#'
#' A real dataset containing the grain yield and plant heigth performance of upland rice lines
#' evaluated across multiple environments. This dataset is used to demonstrate
#' the application of the fuzzy logic methodology for adaptability and stability
#' analysis implemented in the \code{safuzzy} package.
#'
#' The data consists of phenotypic evaluations of elite lines and commercial
#' cultivars of upland rice. The analysis provides membership degrees that assist
#' breeders in selecting stable and high-yielding genotypes for target environments.
#'
#' @format A data frame (or tibble) with columns representing the experimental factors:
#' \describe{
#'   \item{\code{genotype}}{Factor representing the evaluated upland rice lines/genotypes.}
#'   \item{\code{environment}}{Factor representing the test environments (combinations of locations and crop years).}
#'   \item{\code{block}}{Factor representing the local control.}
#'   \item{\code{gy}}{Numeric variable containing the grain yield performance (e.g., kg/ha).}
#'   \item{\code{ph}}{Numeric variable containing the plant heigth performance (e.g., cm).}
#' }
#'
#' @source Data obtained from the breeding trials conducted and published by Maciel et al. (2025).
#'
#' @references Maciel, D. D. O., Guimarães, P. H. R., & Melo, P. G. S. (2025).
#' Harnessing fuzzy logic for adaptive and stable selection of upland rice lines.
#' \emph{Crop Breeding and Applied Biotechnology}, 25(2), e527425213.
#' \doi{10.1590/1984-70332025v25n2a28}
#'
#' @examples
#' data(ge_data)
#' head(ge_data)
"ge_data"
