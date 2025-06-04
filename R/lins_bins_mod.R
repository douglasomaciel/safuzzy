#' Modified method of Lin & Binns (1988)
#'
#' @description Stability and Adaptability Analysis based on the interpretation of the modified Lins and Bins (1988) methodology.
#'
#' @references Carneiro, P. C. S. (1998). Novas metodologias de análise da adaptabilidade e estabilidade de comportamento (Doctoral dissertation, Universidade Federal de Viçosa.).
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
#'     \item{\code{PIF}}{: Performance index in favorable environments.}
#'     \item{\code{PID}}{: Performance index in unfavorable environments.}
#'     \item{\code{GE}}{: Membership (\%) to the general stability genotypes group.}
#'     \item{\code{PA}}{: Membership (\%) to the  poorly adapted geotypes group.}
#'     \item{\code{FAV}}{: Membership (\%) to the favorable adaptade genotypes group.}
#'     \item{\code{UNF}}{: Membership (\%) to the unffavorable adaptade genotypes group.}
#'   }
#'
#' @seealso \code{\link{eberhart_russell}}
#' @seealso \code{\link{hybrid}}
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
#' lin_binns(data = dados, env = environment, gen = treatment, rep = block, var = gy)
#'
#' @export

lin_binns = function(data, env, gen, rep, var) {

  Dados <- data %>%
    rename(Amb = {{env}},
           Gen = {{gen}},
           Rep = {{rep}},
           Yvar = {{var}})

  Dados$Amb <- as.factor(Dados$Amb)
  Dados$Gen <- as.factor(Dados$Gen)
  Dados$Rep <- as.factor(Dados$Rep)

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

  # Padronizar pif i pid

  LI=0
  LS=100

  PIFMax=mean(pif_pid$PIF)+3*sd(pif_pid$PIF)
  PIFMin=mean(pif_pid$PIF)-3*sd(pif_pid$PIF)
  PIDMax=mean(pif_pid$PID)+3*sd(pif_pid$PID)
  PIDMin=mean(pif_pid$PID)-3*sd(pif_pid$PID)

  pif_pid <- pif_pid %>%
    mutate(pif_pad = LS + (LS - LI) * (PIF - PIFMax) / (PIFMax - PIFMin))%>%
    mutate(pid_pad = LS + (LS - LI) * (PID - PIDMax) / (PIDMax - PIDMin))

  # pertinences
  pert_pif <- pif_pid %>%
    mutate(pif_baixa = zmf(pif_pad, 0, 60),
           pif_alta = smf(pif_pad, 0, 60)) %>%
    select(Gen, pif_baixa, pif_alta)

  pert_pid <- pif_pid %>%
    mutate(pid_baixa = zmf(pid_pad, 0, 100),
           pid_alta = smf(pid_pad, 0, 100)) %>%
    select(Gen, pid_baixa, pid_alta)


  # --- Regras e inferência
  Regras <- matrix(c(
    1, 1, 1,
    1, 2, 4,
    2, 1, 2,
    2, 2, 3
  ), ncol = 3, byrow = TRUE)

  # Aplicar les regles
  PertSaida <- sapply(1:nrow(pif_pid), function(i) {
    sapply(1:nrow(Regras), function(j) {
      min(c(
        pert_pif[[i, Regras[j, 1] + 1]], # +1 por a 1 coluna é Gen
        pert_pid[[i, Regras[j, 2] + 1]]
      ))
    })
  })

  PertSaida <- t(PertSaida) # Transpor para que as linhas sejam os genótipos

  GE=apply(cbind(PertSaida[,1]),1,max)
  UNF=apply(cbind(PertSaida[,3]),1,max)
  PA=apply(cbind(PertSaida[,4]),1,max)
  FAV=apply(cbind(PertSaida[,2]),1,max)
  GF=apply(cbind(PertSaida[,1]),1,max)

  Pertinencias <- data.frame(
    Gen = pif_pid$Gen, # Agregar Gen aquí
    GE = GE,
    PA = PA,
    FAV = FAV,
    UNF = UNF
  )

  Resultado <- left_join(pif_pid, Pertinencias, by = "Gen")

  saida = Resultado%>%
    mutate(PIF = round(PIF,2))%>%
    mutate(PID = round(PID,2)) %>%
    mutate(GE = round(GE*100,0))%>%
    mutate(PA = round(PA*100,0))%>%
    mutate(FAV = round(FAV*100,0))%>%
    mutate(UNF = round(UNF*100,0))%>%
    select(Gen,PIF,PID,GE,PA,FAV,UNF)


  return(saida) # Retornar um data.frame com as estimativas de parâmetros de adpt. estab. e pertinências
}

