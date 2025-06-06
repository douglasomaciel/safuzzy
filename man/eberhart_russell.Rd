% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eberhart_russel.R
\name{eberhart_russell}
\alias{eberhart_russell}
\title{Eberhart & Russel (1966) method}
\usage{
eberhart_russell(data, env, gen, rep, var)
}
\arguments{
\item{data}{Data file (data.frame) with variables.}

\item{env}{Column with environment information.}

\item{gen}{Coluna contendo informações de genótipo.}

\item{rep}{Column containing genotype information.}

\item{var}{Variable to be analyzed.}
}
\value{
A data frame containing the following estimates:
  \itemize{
    \item{\code{Gen}}{: Genotype.}
    \item{\code{B_0}}{: Mean of the variable for each genotype.}
    \item{\code{B_1}}{: Regression coefficient (B_1) for each genotype.}
    \item{\code{R2}}{: Coefficient of determination (R²) for each genotype.}
    \item{\code{GE}}{: Membership (\%) to the general stability genotypes group.}
    \item{\code{PA}}{: Membership (\%) to the  poorly adapted geotypes group.}
    \item{\code{FAV}}{: Membership (\%) to the favorable adaptade genotypes group.}
    \item{\code{UNF}}{: Membership (\%) to the unffavorable adaptade genotypes group.}
  }
}
\description{
Adaptability and Stability Analysis based on the interpretation of the Eberhart & Russel 1966 methodology, developed by Carneiro et al. 2018.
}
\examples{
# Caminho do arquivo de dados
caminho <- system.file("extdata", "Dados.csv", package = "safuzzy")

# Leitura do arquivo
dados = read.csv2(caminho)

# Uso da função com a variável de interesse
eberhart_russell(data = dados, env = environment, gen = treatment, rep = block, var = gy)

}
\references{
Carneiro, V. Q., Prado, A. L. D., Cruz, C. D., Carneiro, P. C. S., Nascimento, M., & Carneiro, J. E. D. S. (2018). Fuzzy control systems for decision-making in cultivars recommendation. Acta Scientiarum. Agronomy, 40, e39314.
}
\seealso{
\code{\link{hybrid}}

\code{\link{lin_binns}}
}
\author{
Douglas de Oliveira Maciel \email{douglasmaciel@discente.ufg.br}
}
