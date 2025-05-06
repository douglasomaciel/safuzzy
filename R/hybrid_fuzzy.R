#' Híbrido método
#'
#' @description
#' Análise de Estabilidade e Adaptabilidade baseado na interpretação da metodologia de Eberhart & Russel 1966, Associado ao método de Lins & Bins modificado, desenvolvido por Carneiro et al. (2020).
#'
#' @references Carneiro, A. R. T., Sanglard, D. A., Azevedo, A. M., Souza, T. L. P. O. D., Pereira, H. S., Melo, L. C., & Carneiro, P. C. S. (2020). Fuzzy logic applied to different adaptability and stability methods in common bean. Pesquisa agropecuária brasileira, 55, e01609.
#'
#' @param Dados Arquivo de dados com as variáveis. Deve ser organizado na seguinte ordem: Ambiente, Genótipo, Bloco.
#'
#' @param variavel Variável a ser análisada.

#' @export
hibrido <- function(Dados, variavel) {
  # --- Pré-processamento ---
  dados$Var <- dados[[deparse(substitute(variavel))]]
  colnames(dados)[1:3] <- c("Amb", "Gen", "Rep")
  dados$Amb <- as.factor(dados$Amb)
  dados$Gen <- as.factor(dados$Gen)
  dados$Rep <- as.factor(dados$Rep)
  dados <- dados[order(dados$Amb, dados$Gen, dados$Rep), ]

  ngen <- length(unique(dados$Gen))
  namb <- length(unique(dados$Amb))
  nrep <- length(unique(dados$Rep))

  # --- Médias ---
  MediaGen <- tapply(dados$Var, dados$Gen, mean)
  MediaAmb <- tapply(dados$Var, dados$Amb, mean)
  IJ <- MediaAmb - mean(MediaAmb)

  DadosGA <- tapply(dados$Var, list(dados$Gen, dados$Amb), mean)
  Y <- DadosGA
  X <- IJ

  D <- t(apply(Y, 1, function(y) {
    a <- summary(lm(y ~ X))
    c(coefficients(a)[, 1], R2 = a$r.squared)
  }))
  colnames(D) <- c("Bo", "B1", "R2")

  # --- Erro da regressão ---
  modelo <- summary(aov(Var ~ Gen + Rep/Amb + Amb + Gen:Amb, data = dados))
  QMR <- modelo[[1]]$`Mean Sq`[6]
  GLR <- modelo[[1]]$`Df`[6]
  se <- QMR / nrep
  VarB1 <- se / sum(X^2)

  B1Pad <- padronizar(D[, "B1"], -2, 4)
  R2Pad <- D[, "R2"] * 100

  Pert_B1 <- cbind(
    MFMenor = zmf(B1Pad, c(-5, 1))$y,
    MFIgual = pimf(B1Pad, c(-5, 1, 1, 7))$y,
    MFMaior = smf(B1Pad, c(1, 7))$y
  )

  Pert_R2 <- cbind(
    MFBaixo = zmf(R2Pad, c(60, 100))$y,
    MFAlto = smf(R2Pad, c(60, 100))$y
  )

  # --- PIF e PID ---
  AmbFav <- IJ >= 0
  AmbDesf <- IJ < 0
  YIJ <- matrix(DadosGA, ncol = namb, byrow = FALSE)

  PIF <- sapply(1:ngen, function(g) sum((YIJ[g, AmbFav] - apply(YIJ[, AmbFav], 2, max))^2) / (2 * sum(AmbFav)))
  PID <- sapply(1:ngen, function(g) sum((YIJ[g, AmbDesf] - apply(YIJ[, AmbDesf], 2, max))^2) / (2 * sum(AmbDesf)))

  PIFpad <- padronizar(PIF, 0, 100)
  PIDpad <- padronizar(PID, 0, 100)

  Pert_PIF <- cbind(
    MFBaixa = zmf(PIFpad, c(0, 100))$y,
    MFAlta = smf(PIFpad, c(0, 100))$y
  )
  Pert_PID <- cbind(
    MFBaixa = zmf(PIDpad, c(0, 100))$y,
    MFAlta = smf(PIDpad, c(0, 100))$y
  )

  # --- Regras e Inferência ---
  Regras <- matrix(c(
    1, 1, 1, 1, 3,
    1, 1, 1, 2, 2,
    1, 2, 1, 1, 3,
    1, 2, 1, 2, 3,
    1, 1, 2, 1, 3,
    1, 1, 2, 2, 1,
    1, 2, 2, 1, 3,
    1, 2, 2, 2, 3,
    1, 1, 3, 1, 3,
    1, 1, 3, 2, 4,
    1, 2, 3, 1, 3,
    1, 2, 3, 2, 4,
    2, 1, 1, 1, 3,
    2, 1, 1, 2, 2,
    2, 2, 1, 1, 3,
    2, 2, 1, 2, 3,
    2, 1, 2, 1, 3,
    2, 1, 2, 2, 2,
    2, 2, 2, 1, 3,
    2, 2, 2, 2, 3,
    2, 1, 3, 1, 3,
    2, 1, 3, 2, 3,
    2, 2, 3, 1, 3,
    2, 2, 3, 2, 3
  ), ncol = 5, byrow = TRUE)

  PertSaida <- t(sapply(1:ngen, function(i) sapply(1:nrow(Regras), function(j) {
    min(c(
      Pert_PIF[i, Regras[j, 1]],
      Pert_PID[i, Regras[j, 2]],
      Pert_B1[i, Regras[j, 3]],
      Pert_R2[i, Regras[j, 4]]
    ))
  })))

  GE <- apply(cbind(PertSaida[, 6]), 1, max)
  Des <- apply(cbind(PertSaida[, c(2, 14, 18)]), 1, max)
  PA <- apply(cbind(PertSaida[, c(1, 3:5, 7:9, 11, 13, 15:17, 19:24)]), 1, max)
  FAV <- apply(cbind(PertSaida[, c(10, 12)]), 1, max)
  GF <- apply(cbind(PertSaida[, 6]), 1, max)

  Pertinencias <- round(100 * cbind(GE, PA, FAV, Des), 0)
  Resultado <- cbind(PIF, PID, B1 = D[, "B1"], R2 = D[, "R2"], Pertinencias)
  return(Resultado)
}
