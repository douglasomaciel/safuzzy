% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lins_bins_mod.R
\name{lin_binns}
\alias{lin_binns}
\title{Modified method of Lin & Binns (1988)}
\usage{
lin_binns(data, env, gen, rep, var)
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
    \item{\code{PIF}}{: Performance index in favorable environments.}
    \item{\code{PID}}{: Performance index in unfavorable environments.}
    \item{\code{GE}}{: Membership (\%) to the general stability genotypes group.}
    \item{\code{PA}}{: Membership (\%) to the  poorly adapted geotypes group.}
    \item{\code{FAV}}{: Membership (\%) to the favorable adaptade genotypes group.}
    \item{\code{UNF}}{: Membership (\%) to the unffavorable adaptade genotypes group.}
  }
}
\description{
Stability and Adaptability Analysis based on the interpretation of the modified Lins and Bins (1988) methodology.
}
\examples{
# Caminho do arquivo de dados
caminho <- system.file("extdata", "Dados.csv", package = "safuzzy")

# Leitura do arquivo
dados = read.csv2(caminho)

# Uso da função com a variável de interesse
lin_binns(data = dados, env = environment, gen = treatment, rep = block, var = gy)

}
\references{
Carneiro, P. C. S. (1998). Novas metodologias de análise da adaptabilidade e estabilidade de comportamento (Doctoral dissertation, Universidade Federal de Viçosa.).
}
\seealso{
\code{\link{eberhart_russell}}

\code{\link{hybrid}}
}
\author{
Douglas de Oliveira Maciel \email{douglasmaciel@discente.ufg.br}
}
