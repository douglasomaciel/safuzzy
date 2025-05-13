#' Eberhart & Russel método
#'
#' @description
#' Análise de Adaptabilidade e Estabilidade baseado na interpretação da metodologia de Eberhart & Russel 1966, desenvolvido por Carneiro et al. 2018.
#'
#' @references Carneiro, V. Q., Prado, A. L. D., Cruz, C. D., Carneiro, P. C. S., Nascimento, M., & Carneiro, J. E. D. S. (2018). Fuzzy control systems for decision-making in cultivars recommendation. Acta Scientiarum. Agronomy, 40, e39314.
#'
#' @param dados Arquivo de dados (data.frame) com as variáveis.
#' @param ambiente Coluna com as informações de embiente.
#' @param genotipo Coluna contendo informações de genótipo.
#' @param bloco coluna contendo informações de bloco.
#' @param variavel Variável a ser análisada.
#'
#' @return Um data frame contendo as seguintes estimativas:
#'   \itemize{
#'     \item{\code{Gen}}{: Genótipo.}
#'     \item{\code{B_0}}{: Média da variável para cada genótipo.}
#'     \item{\code{B_1}}{: Coeficiente de regressão (B_1i) para cada genótipo.}
#'     \item{\code{R2}}{: Coeficiente de determinação (R²) padronizado para cada genótipo.}
#'     \item{\code{GE}}{: Pertinência (\%) para o grupo de genótipos de estabilidade geral.}
#'     \item{\code{PA}}{: Pertinência (\%) para o grupo de genótipos pouco adaptados.}
#'     \item{\code{FAV}}{: Pertinência (\%) para o grupo de genótipos com adaptabilidade favorável.}
#'     \item{\code{Des}}{: Pertinência (\%) para o grupo de genótipos adaptados a ambientes desfavoráveis.}
#'   }
#'
#' @author Douglas de Oliveira Maciel \email{douglasmaciel@discente.ufg.br}
#'
#' @examples
#' # Caminho do arquivo de dados
#' caminho <- system.file("extdata", "Dados.csv", package = "safuzzy")
#'
#' # Leitura do arquivo
#' dados = read.csv2(caminho)
#'
#' # Uso da função com a variável de interesse
#' eberhart_russel(dados, ambiente = environment, genotipo = treatment, bloco = block, variavel = gy)
#'

#' @export

eberhart_russel = function(dados, ambiente, genotipo, bloco, variavel){
  library(dplyr)
  Dados <- dados %>%
    rename(Amb = {{ambiente}},
           Gen = {{genotipo}},
           Rep = {{bloco}},
           Yvar = {{variavel}})

  Dados$Amb <- as.factor(Dados$Amb)
  Dados$Gen <- as.factor(Dados$Gen)
  Dados$Rep <- as.factor(Dados$Rep)

  media_amb <- Dados %>%
    group_by(Amb) %>%
    summarise(Yvar = mean(Yvar, na.rm = TRUE)) # Média dos ambientes

  IJ <- media_amb %>%
    mutate(IJ = Yvar - mean(Yvar)) %>%  # Índice ambiental
    select(Amb, IJ)

  media_GxA <- Dados %>%
    group_by(Gen, Amb) %>%
    summarise(Yvar = mean(Yvar, na.rm = TRUE))

  # Juntar as médias GxA com o índice ambiental
  GxA_Env <- media_GxA %>%
    left_join(IJ, by = "Amb")

  # Realizar a regressão para cada genótipo
  reg <- GxA_Env %>%
    group_by(Gen) %>%
    summarise(
      med = mean(Yvar),
      b = coef(lm(Yvar ~ IJ))[2], # Coeficiente de regressão (b_i)
      R2 = summary(lm(Yvar ~ IJ))$r.squared # R-quadrado (R^2_i)
      # Para o desvio da regressão (S^2_di), seria mais complexo e envolveria os resíduos
    )

  # modelo
  mod = summary(aov(Yvar ~ Gen+Rep/Amb+Amb+Gen:Amb, data = Dados))
  nrep <- length(unique(Dados$Rep)) # Calcular nrep
  namb <- length(unique(Dados$Amb))
  QMR=mod[[1]]$`Mean Sq`[6]
  GLR=mod[[1]]$`Df`[6]
  se=QMR/nrep
  VarBo=se*sum(IJ$IJ^2)/(namb*sum(IJ$IJ^2))
  VarB1=se/sum(IJ$IJ^2)

  # padronizar a média
  LI=0
  LS=100
  BMax=mean(reg$med)+3*sd(reg$med)
  BMin=mean(reg$med)-3*sd(reg$med)

  reg <- reg %>%
    mutate(med_pad = LS + (LS - LI) * (med - BMax) / (BMax - BMin))


  # beta 1 padronizado
  LI=-2
  LS=4
  BMax=1+qt(0.975,GLR)*sqrt(VarB1)
  BMin=1+qt(0.025,GLR)*sqrt(VarB1)

  reg <- reg %>%
    mutate(b_pad = LS + (LS - LI) * (b - BMax) / (BMax - BMin))

  # R2 padronizado
  reg <- reg %>%
    mutate(r_pad = R2 * 100)

  # pertinencias
  pert_med <- reg %>%
    mutate(med_baixa = zmf(med_pad, 0, 100),
           med_alta = smf(med_pad, 0, 100)) %>%
    select(Gen, med_baixa, med_alta) # Incluir Gen

  pert_b <- reg %>%
    mutate(b_menor = zmf(b_pad, -5, 1),
           b_igual = pimf(b_pad, -5, 1, 1, 7),
           b_maior = smf(b_pad, 1, 7)) %>%
    select(Gen, b_menor, b_igual, b_maior) # Incluir Gen

  pert_r <- reg %>%
    mutate(r_baixa = zmf(r_pad, 60, 100),
           r_alta = smf(r_pad, 60, 100)) %>%
    select(Gen, r_baixa, r_alta) # Incluir Gen


  # regras
  Regras=matrix(c(1,1,1,3,
                  1,1,2,3,
                  1,2,1,3,
                  1,2,2,3,
                  1,3,1,3,
                  1,3,2,3,
                  2,1,1,3,
                  2,1,2,2,
                  2,2,1,3,
                  2,2,2,1,
                  2,2,2,5,
                  2,3,1,3,
                  2,3,2,4),nrow=13,byrow=T)

  # Aplicar as regras
  PertSaida <- sapply(1:nrow(reg), function(i) {
    sapply(1:nrow(Regras), function(j) {
      min(c(
        pert_med[[i, Regras[j, 1] + 1]], # +1 porque la 1ª columna es Gen
        pert_b[[i, Regras[j, 2] + 1]],   # +1 porque la 1ª columna es Gen
        pert_r[[i, Regras[j, 3] + 1]]    # +1 porque la 1ª columna es Gen
      ))
    })
  })

  PertSaida <- t(PertSaida) # Transpor para que as linhas sejam os genótipos

  GE <- apply(PertSaida[, 10, drop = FALSE], 1, max)
  Des <- apply(PertSaida[, 8, drop = FALSE], 1, max)
  PA <- apply(PertSaida[, c(1:7, 9, 12), drop = FALSE], 1, max)
  FAV <- apply(PertSaida[, 13, drop = FALSE], 1, max)
  GF <- apply(PertSaida[, 11, drop = FALSE], 1, max)

  Pertinencias <- data.frame(
    Gen = reg$Gen, # Agregar Gen aquí
    GE = GE,
    PA = PA,
    FAV = FAV,
    Des = Des
  )

  ResultadoER <- left_join(reg, Pertinencias, by = "Gen")

  saida = ResultadoER%>%
    mutate(B_0 = round(med,2))%>%
    mutate(B_1 = round(b,4))%>%
    mutate(R2 = r_pad,2)%>%
    mutate(GE = round(GE*100,0))%>%
    mutate(PA = round(PA*100,0))%>%
    mutate(FAV = round(FAV*100,0))%>%
    mutate(Des = round(Des*100,0))%>%
    select(Gen,B_0,B_1,R2,GE,PA,FAV,Des)

  return(saida)

}

