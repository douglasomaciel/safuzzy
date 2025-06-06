#' Hybrid method
#'
#' @description
#' Stability and Adaptability Analysis based on the interpretation of the Eberhart & Russel 1966 methodology, associated with the modified Lins & Bins (1988) methodology, developed by Carneiro et al. (2020).
#'
#' @references Carneiro, A. R. T., Sanglard, D. A., Azevedo, A. M., Souza, T. L. P. O. D., Pereira, H. S., Melo, L. C., & Carneiro, P. C. S. (2020). Fuzzy logic applied to different adaptability and stability methods in common bean. Pesquisa agropecuária brasileira, 55, e01609.
#'
#' @param data Data file (data.frame) with variables.
#' @param env Column with environment information.
#' @param gen Coluna contendo informações de genótipo.
#' @param rep Column containing genotype information.
#' @param var Variable to be analyzed.
#'
#' @import dplyr
#'
#' @return Um data frame contendo as seguintes estimativas:
#'   \itemize{
#'     \item{\code{Gen}}{: Genotype.}
#'     \item{\code{PIF}}{: Performance index in favorable environments.}
#'     \item{\code{PID}}{: Performance index in unfavorable environments.}
#'     \item{\code{B_1}}{: Regression coefficient (B_1) for each genotype.}
#'     \item{\code{R2}}{: Standardized coefficient of determination (R²) for each genotype.}
#'     \item{\code{GE}}{: Membership (\%) to the general stability genotypes group.}
#'     \item{\code{PA}}{: Membership (\%) to the  poorly adapted geotypes group.}
#'     \item{\code{FAV}}{: Membership (\%) to the favorable adaptade genotypes group.}
#'     \item{\code{UNF}}{: Membership (\%) to the unffavorable adaptade genotypes group.}
#'   }
#'
#' @author Douglas de Oliveira Maciel \email{douglasmaciel@discente.ufg.br}
#'
#' @seealso \code{\link{eberhart_russell}}
#' @seealso \code{\link{lin_binns}}
#' @seealso \code{\link{annicchiarico}}
#'
#' @examples
#' # Caminho do arquivo de dados
#' caminho <- system.file("extdata", "Dados.csv", package = "safuzzy")
#'
#' # Leitura do arquivo
#' dados = read.csv2(caminho)
#'
#' # Uso da função com a variável de interesse
#' hybrid(data = dados, env = environment, gen = treatment, rep = block, var = gy)
#'
#' @export

hybrid = function(data, env, gen, rep, var) {

  Dados <- data %>%
    rename(Amb = {{env}},
           Gen = {{gen}},
           Rep = {{rep}},
           Yvar = {{var}})%>%
    mutate(across(c(Amb, Gen, Rep), as.factor))


  media_amb <- Dados %>%
    group_by(Amb) %>%
    summarise(Media_Amb = mean(Yvar, na.rm = TRUE)) # Média dos ambientes

  media_geral <- mean(Dados$Yvar, na.rm = TRUE) # Média geral

  IJ <- media_amb %>%
    mutate(IJ = Media_Amb - media_geral) %>%  # Índice ambiental
    select(Amb, IJ)

  media_GxA <- Dados %>%
    group_by(Gen, Amb) %>%
    summarise(Yvar = mean(Yvar, na.rm = TRUE))

  # Encontrar a resposta máxima em cada ambiente (Mj)
  max_amb <- media_GxA %>%
    group_by(Amb) %>%
    summarise(Mj = max(Yvar, na.rm = TRUE))

  # Juntar as médias GxA com o índice ambiental e a resposta máxima
  GxA_Env <- media_GxA %>%
    left_join(IJ, by = "Amb") %>%
    left_join(max_amb, by = "Amb")

  # Calcular o número de ambientes favoráveis e desfavoráveis
  num_amb_fav <- sum(IJ$IJ >= 0)
  num_amb_desf <- sum(IJ$IJ < 0)

  # Calcular PIF i PID pera cada genotipo
  pif_pid <- GxA_Env %>%
    group_by(Gen) %>%
    summarise(
      PIF = sum((Yvar - Mj)^2 * (IJ >= 0)) / (2 * num_amb_fav),
      PID = sum((Yvar - Mj)^2 * (IJ < 0)) / (2 * num_amb_desf)
    )

  # Realizar a regressão para cada genótipo E JUNTAR PIF e PID
  reg <- GxA_Env %>%
    group_by(Gen) %>%
    summarise(
      med = mean(Yvar),
      b = coef(lm(Yvar ~ IJ))[2],
      r2 = summary(lm(Yvar ~ IJ))$r.squared
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

  # Padronizar pif e pid
  LI=0
  LS=100

  PIFMax=mean(pif_pid$PIF)+3*sd(pif_pid$PIF)
  PIFMin=mean(pif_pid$PIF)-3*sd(pif_pid$PIF)
  PIDMax=mean(pif_pid$PID)+3*sd(pif_pid$PID)
  PIDMin=mean(pif_pid$PID)-3*sd(pif_pid$PID)

  pif_pid <- pif_pid %>%
    mutate(pif_pad = LS + (LS - LI) * (PIF - PIFMax) / (PIFMax - PIFMin),
           pid_pad = LS + (LS - LI) * (PID - PIDMax) / (PIDMax - PIDMin))

  reg <- pif_pid %>%
    left_join(reg, by = "Gen")

  #Fuzificação
  pert_pif <- reg %>%
    mutate(baixo = zmf(pif_pad, 0, 100),
           alto = smf(pif_pad, 0, 100)) %>%
    select(Gen, baixo, alto)

  pert_pid <- reg %>%
    mutate(baixo = zmf(pid_pad, 0, 100),
           alto = smf(pid_pad, 0, 100)) %>%
    select(Gen, baixo, alto)

  pert_b <- reg %>%
    mutate(menor = zmf(b_pad, -5, 1),
           igual = pimf(b_pad, -5, 1, 1, 7),
           maior = smf(b_pad, 1, 7)) %>%
    select(Gen, menor, igual, maior)

  pert_r <- reg %>%
    mutate(baixo = zmf(R2, 60, 100),
           alto = smf(R2, 60, 100)) %>%
    select(Gen, baixo, alto)

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

  # Aplicar as regras
  PertSaida <- sapply(1:nrow(reg), function(i) {
    sapply(1:nrow(Regras), function(j) {
      min(c(
        pert_pif[[i, Regras[j, 1] + 1]], # +1 porque a 1 coluna é Gen
        pert_pid[[i, Regras[j, 2] + 1]],
        pert_b[[i, Regras[j, 3] + 1]],
        pert_r[[i, Regras[j, 4] + 1]]
      ))
    })
  }) %>% t()

  GE=apply(cbind(PertSaida[,6]),1,max)
  UNF=apply(cbind(PertSaida[,c(2,14,18)]),1,max)
  PA=apply(cbind(PertSaida[,c(1,3:5,7:9,11,13,15:17,19:24)]),1,max)
  FAV=apply(cbind(PertSaida[,c(10,12)]),1,max)

  Pertinencias <- data.frame(
    Gen = reg$Gen, # Agregar Gen aquí
    GE = GE,
    PA = PA,
    FAV = FAV,
    UNF = UNF
  )

  Resultado <- left_join(reg, Pertinencias, by = "Gen")

  saida = Resultado%>%
    rename(B_1 = b)%>%
    mutate(
      across(c(PIF,PID,R2),~round(.x,2)),
      B_1 = round(B_1,4),
      across(c(GE,PA,FAV,UNF),~round(.x*100,0))
    )%>%
    select(Gen,PIF,PID,B_1,R2,GE,PA,FAV,UNF)


  return(saida) # Retornar um data.frame com as estimativas de parâmetros de adpt. estab. e pertinências
}

