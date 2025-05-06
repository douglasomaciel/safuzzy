padronizar_b0 <- function(Dados, var_name, D, nrep) {
  # Extrai a variável de interesse
  Y_obs <- Dados[[var_name]]

  # Ajusta o modelo de ANOVA conforme Eberhart & Russell
  modelo <- summary(aov(Y_obs ~ Dados$Gen + Dados$Rep / Dados$Amb + Dados$Amb + Dados$Gen:Dados$Amb))

  # Cálculo dos limites para padronização da média
  Med <- D[,1]
  BMax <- mean(Med) + 3 * sd(Med)
  BMin <- mean(Med) - 3 * sd(Med)

  # Padronização da média
  MedPad <- sapply(Med, function(x) 100 * (x - BMin) / (BMax - BMin))

  return(MedPad)
}

