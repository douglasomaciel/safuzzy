#' Eberhart & Russel método
#'
#' @description
#' Análise de Estabilidade e Adaptabilidade baseado na interpretação da metodologia de Eberhart & Russel 1966, desenvolvido por Carneiro et al. 2018.
#'
#' @references Carneiro, V. Q., Prado, A. L. D., Cruz, C. D., Carneiro, P. C. S., Nascimento, M., & Carneiro, J. E. D. S. (2018). Fuzzy control systems for decision-making in cultivars recommendation. Acta Scientiarum. Agronomy, 40, e39314.
#'
#' @param Dados Arquivo de dados com as variáveis. Deve ser organizado na seguinte ordem: Ambiente, Genótipo, Bloco.
#'
#' @param variavel Variável a ser análisada.
#'
#' @return Um data frame contendo as seguintes estimativas:
#'    \item{\code{Media(\beta_{0})}}{Para cada genótipo.}
#'
#' @examples
#' eberhart_russel(Dados = dados, varivavel = prod)

#' @export

eberhart_russel <- function(Dados, variavel) {
  # Captura o nome da variável passado sem aspas
  var_name <- deparse(substitute(variavel))

  # Renomear colunas iniciais
  colnames(Dados)[1:3] <- c("Amb", "Gen", "Rep")
  Dados$Amb <- as.factor(Dados$Amb)
  Dados$Gen <- as.factor(Dados$Gen)
  Dados$Rep <- as.factor(Dados$Rep)
  Dados <- Dados[order(Dados$Rep), ]
  Dados <- Dados[order(Dados$Gen), ]
  Dados <- Dados[order(Dados$Amb), ]

  if (ncol(Dados) == 4) {
    Dados <- cbind(Dados, Dados[[var_name]])
    colnames(Dados)[5] <- paste0(var_name, "_copy")
  }

  ngen <- length(unique(Dados$Gen))
  namb <- length(unique(Dados$Amb))
  nrep <- length(unique(Dados$Rep))

  vars <- colnames(Dados)[-(1:3)]
  MediaGen <- apply(Dados[, vars], 2, function(x) tapply(x, Dados$Gen, mean))
  MediaAmb <- apply(Dados[, vars], 2, function(x) tapply(x, Dados$Amb, mean))
  IJ <- apply(MediaAmb, 2, function(x) x - mean(x))

  DadosGA <- apply(Dados[, vars], 2, function(x) tapply(x, Dados$Gen:Dados$Amb, mean))
  var_index <- which(vars == var_name)
  Y <- matrix(DadosGA[, var_index], ncol = namb, byrow = TRUE)
  X <- IJ[, var_index]

  D <- t(apply(Y, 1, function(y) {
    a <- summary(lm(y ~ X))
    c(coefficients(a)[,1], R2 = a$r.squared)
  }))

  Y_obs <- Dados[[var_name]]
  modelo <- summary(aov(Y_obs ~ Dados$Gen + Dados$Rep / Dados$Amb + Dados$Amb + Dados$Gen:Dados$Amb))
  QMR <- modelo[[1]]$`Mean Sq`[6]
  GLR <- modelo[[1]]$`Df`[6]
  se <- QMR / nrep
  VarBo <- se * sum(X^2) / (namb * sum(X^2))
  VarB1 <- se / sum(X^2)


  Med <- D[,1]
  BMax <- mean(Med) + 3 * sd(Med)
  BMin <- mean(Med) - 3 * sd(Med)
  MedPad <- sapply(Med, function(x) 100 * (x - BMin) / (BMax - BMin))

  B1 <- D[,2]
  BMax <- 1 + qt(0.975, GLR) * sqrt(VarB1)
  BMin <- 1 + qt(0.025, GLR) * sqrt(VarB1)
  B1Pad <- sapply(B1, function(x) 6 * (x - BMin) / (BMax - BMin) - 2)

  R2 <- 100 * D[,3]
  D2 <- cbind(MedPad, B1Pad, R2)

  Pert_B0 <- cbind(MFBaixa = zmf(D2[,1], c(0, 100)),
                   MFAlta = smf(D2[,1], c(0, 100)))

  Pert_B1 <- cbind(MFMenor = zmf(D2[,2], c(-5, 1)),
                   MFIgual = pimf(D2[,2], c(-5, 1, 1, 7)),
                   MFMaior = smf(D2[,2], c(1, 7)))

  Pert_R2 <- cbind(MFBaixo = zmf(D2[,3], c(60, 100)),
                   MFAlto = smf(D2[,3], c(60, 100)))

  Regras <- matrix(c(
    1,1,1,3, 1,1,2,3, 1,2,1,3, 1,2,2,3,
    1,3,1,3, 1,3,2,3, 2,1,1,3, 2,1,2,2,
    2,2,1,3, 2,2,2,1, 2,2,2,5, 2,3,1,3,
    2,3,2,4
  ), nrow = 13, byrow = TRUE)

  PertSaida <- t(sapply(1:ngen, function(i)
    sapply(1:nrow(Regras), function(j)
      min(c(Pert_B0[i, Regras[j,1]], Pert_B1[i, Regras[j,2]], Pert_R2[i, Regras[j,3]])))
  ))

  GE <- apply(cbind(PertSaida[,10]), 1, max)
  Des <- apply(cbind(PertSaida[,8]), 1, max)
  PA <- apply(cbind(PertSaida[,c(1:7,9,12)]), 1, max)
  FAV <- apply(cbind(PertSaida[,13]), 1, max)

  Pertinencias <- round(100 * cbind(GE, PA, FAV, Des), 0)
  ResultadoER <- cbind(D, Pertinencias)
  colnames(ResultadoER)[4:7] <- c("GE", "PA", "FAV", "Des")

  return(ResultadoER)
}
