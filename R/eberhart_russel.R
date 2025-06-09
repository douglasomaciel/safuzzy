#' Eberhart & Russel (1966) method
#'
#' @description
#' Adaptability and Stability Analysis based on the interpretation of the Eberhart & Russel 1966 methodology, developed by Carneiro et al. 2018.
#'
#' @references Carneiro, V. Q., Prado, A. L. D., Cruz, C. D., Carneiro, P. C. S., Nascimento, M., & Carneiro, J. E. D. S. (2018). Fuzzy control systems for decision-making in cultivars recommendation. Acta Scientiarum. Agronomy, 40, e39314.
#'
#' @param data Data file (data.frame) with variables.
#' @param env Column with environment information.
#' @param gen Coluna contendo informações de genótipo.
#' @param rep Column containing genotype information.
#' @param var Variable to be analyzed.
#'
#' @import dplyr
#'
#' @return A data frame containing the following estimates:
#'   \itemize{
#'     \item{\code{Gen}}{: Genotype.}
#'     \item{\code{B_0}}{: Mean of the variable for each genotype.}
#'     \item{\code{B_1}}{: Regression coefficient (B_1) for each genotype.}
#'     \item{\code{R2}}{: Coefficient of determination (R²) for each genotype.}
#'     \item{\code{GE}}{: Membership (\%) to the general stability genotypes group.}
#'     \item{\code{PA}}{: Membership (\%) to the  poorly adapted geotypes group.}
#'     \item{\code{FAV}}{: Membership (\%) to the favorable adaptade genotypes group.}
#'     \item{\code{UNF}}{: Membership (\%) to the unffavorable adaptade genotypes group.}
#'   }
#'
#'
#' @seealso \code{\link{hybrid}}
#' @seealso \code{\link{lin_binns}}
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
#' eberhart_russell(data = dados, env = environment, gen = treatment, rep = block, var = gy)
#'

#' @export

eberhart_russell = function(data, env, gen, rep, var){
  Dados <- data %>%
    rename(Amb = {{env}},
           Gen = {{gen}},
           Rep = {{rep}},
           Yvar = {{var}})%>%
    mutate(across(c(Amb, Gen, Rep), as.factor))


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
      r2 = summary(lm(Yvar ~ IJ))$r.squared # R-quadrado (R^2_i)
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
    mutate(R2 = r2 * 100)

  # pertinencias
  pert_med <- reg %>%
    mutate(baixa = zmf(med_pad, 0, 100),
           alta = smf(med_pad, 0, 100)) %>%
    select(Gen, baixa, alta) # Incluir Gen

  pert_b <- reg %>%
    mutate(menor = zmf(b_pad, -5, 1),
           igual = pimf(b_pad, -5, 1, 1, 7),
           maior = smf(b_pad, 1, 7)) %>%
    select(Gen, menor, igual, maior) # Incluir Gen

  pert_r <- reg %>%
    mutate(baixo = zmf(R2, 60, 100),
           alto = smf(R2, 60, 100)) %>%
    select(Gen, baixo, alto) # Incluir Gen

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
                  2,3,1,3,
                  2,3,2,4),nrow=12,byrow=T)

  # Aplicar as regras
  PertSaida <- sapply(1:nrow(reg), function(i) {
    sapply(1:nrow(Regras), function(j) {
      min(c(
        pert_med[[i, Regras[j, 1] + 1]], # +1 porque a 1ª columna é Gen
        pert_b[[i, Regras[j, 2] + 1]],
        pert_r[[i, Regras[j, 3] + 1]]
      ))
    })
  }) %>% t()

  Pertinencias <- tibble(
    Gen = reg$Gen, # Agregar Gen aqui
    GE = apply(PertSaida[, 10, drop = FALSE], 1, max),
    PA = apply(PertSaida[, c(1:7, 9, 11), drop = FALSE], 1, max),
    FAV = apply(PertSaida[, 12, drop = FALSE], 1, max),
    UNF = apply(PertSaida[, 8, drop = FALSE], 1, max)
  )

  ResultadoER <- left_join(reg, Pertinencias, by = "Gen")

  saida = ResultadoER%>%
    rename(B_0 = med,
           B_1 = b)%>%
    mutate(B_1 = round(B_1,4),
           across(c(B_0,R2), ~ round(.x,2)),
      across(c(R2,GE,PA,FAV,UNF), ~round(.x *100,0))
    )%>%
    select(Gen,B_0,B_1,R2,GE,PA,FAV,UNF)

  return(saida)
}

