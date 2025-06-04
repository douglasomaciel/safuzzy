#' Annicchiarico (1992) method
#'
#' @description
#' Adaptability and Stability Analysis based on the interpretation of the Annicchiarico 1992 methodology, developed by Carneiro et al. 2019.
#'
#' @references Carneiro, A. R. T., Sanglard, D. A., Azevedo, A. M., Souza, T. L. P. O. D., Pereira, H. S., & Melo, L. C. (2019). Fuzzy logic in automation for interpretation of adaptability and stability in plant breeding studies. Scientia Agricola, 76(2), 123-129.
#'
#' @param data Data file (data.frame) with variables.
#' @param env Column with environment information.
#' @param gen Coluna contendo informações de genótipo.
#' @param rep Column containing genotype information.
#' @param var Variable to be analyzed.
#'
#' @import dplyr
#' @import tidyr
#'
#' @return A data frame containing the following estimates:
#'   \itemize{
#'     \item{\code{Gen}}{: Genotype.}
#'     \item{\code{Wg}}{: General recommendation index for environments.}
#'     \item{\code{Wd}}{: Recommendation index for unfavorable environments.}
#'     \item{\code{Wf}}{: Recommendation index for favorable environments.}
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
#' annicchiarico(data = dados, env = environment, gen = treatment, rep = block, var = gy)
#'

#' @export

annicchiarico = function(data, env, gen, rep, var){
  Dados <- data %>%
    rename(Amb = {{env}},
           Gen = {{gen}},
           Rep = {{rep}},
           Yvar = {{var}})

  Dados <- Dados %>%
    mutate(Amb = as.factor(Amb),
           Gen = as.factor(Gen),
           Rep = as.factor(Rep))

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

  Y = media_GxA %>%
    left_join(media_amb, by = "Amb", suffix = c("", "_ambiente")) %>% # Junta a média do ambiente
    mutate(Yp = 100 * Yvar / Yvar_ambiente) %>% # Yvar é a média GxA, Yvar_ambiente é a média do ambiente
    select(-Yvar_ambiente) # Remove a coluna extra se não for mais necessária

  # Cálculo do Wg
  Wg <- Y %>%
    group_by(Gen) %>% # Agrupar por Genótipo
    summarise(
      # Média dos valores de Yp para cada genótipo
      mean_Yp_gen = mean(Yp, na.rm = TRUE),
      # Desvio padrão dos valores de Yp para cada genótipo
      sd_Yp_gen = sd(Yp, na.rm = TRUE)
    ) %>%
    # Aplicar a fórmula do Wg
    mutate(Wg = mean_Yp_gen - sd_Yp_gen * qnorm(0.75)) %>%
    select(Gen, Wg) # Selecionar apenas o Gen e o valor final de Wg

  # Passo 1: Identificar os ambientes desfavoráveis
  ambientes_desfavoraveis <- IJ %>%
    filter(IJ <= 0) %>%
    pull(Amb) # pull() extrai a coluna como um vetor simples

  # Passo 2: Filtrar o data frame Y para incluir APENAS os ambientes desfavoráveis
  Y_desfavoravel <- Y %>%
    filter(Amb %in% ambientes_desfavoraveis)

  # Passo 3: Calcular Wd sobre os Yp dos ambientes desfavoráveis
  Wd <- Y_desfavoravel %>%
    group_by(Gen) %>% # Agrupar por Genótipo
    summarise(
      # Média dos valores de Yp para cada genótipo, APENAS nos ambientes desfavoráveis
      mean_Yp_desfavoravel = mean(Yp, na.rm = TRUE),
      # Desvio padrão dos valores de Yp para cada genótipo, APENAS nos ambientes desfavoráveis
      sd_Yp_desfavoravel = sd(Yp, na.rm = TRUE)
    ) %>%
    # Aplicar a fórmula do Wd
    mutate(Wd = mean_Yp_desfavoravel - sd_Yp_desfavoravel * qnorm(0.75)) %>%
    select(Gen, Wd) # Selecionar apenas o Gen e o valor final de Wd

  # Passo 1: Identificar os ambientes favoráveis
  ambientes_favoraveis <- IJ %>%
    filter(IJ > 0) %>%
    pull(Amb) # pull() extrai a coluna como um vetor simples

  # Passo 2: Filtrar o data frame Y para incluir APENAS os ambientes favoráveis
  Y_favoravel <- Y %>%
    filter(Amb %in% ambientes_favoraveis)

  # Passo 3: Calcular Wd sobre os Yp dos ambientes favoráveis
  Wf <- Y_favoravel %>%
    group_by(Gen) %>% # Agrupar por Genótipo
    summarise(
      # Média dos valores de Yp para cada genótipo, APENAS nos ambientes favoráveis
      mean_Yp_favoravel = mean(Yp, na.rm = TRUE),
      # Desvio padrão dos valores de Yp para cada genótipo, APENAS nos ambientes favoráveis
      sd_Yp_favoravel = sd(Yp, na.rm = TRUE)
    ) %>%
    # Aplicar a fórmula do Wf
    mutate(Wf = mean_Yp_favoravel - sd_Yp_favoravel * qnorm(0.75)) %>%
    select(Gen, Wf) # Selecionar apenas o Gen e o valor final de Wf

  W = Wg %>%
    left_join(Wd, by="Gen")%>%
    left_join(Wf, by="Gen")

  pert_wd <- W %>%
    mutate(wd_baixa = zmf(Wd, 0, 200),
           wd_alta = smf(Wd, 0, 200)) %>%
    select(Gen, wd_baixa, wd_alta) # Incluir Gen

  pert_wf <- W %>%
    mutate(wf_baixa = zmf(Wf, 0, 200),
           wf_alta = smf(Wf, 0, 200)) %>%
    select(Gen, wf_baixa, wf_alta) # Incluir Gen

  Regras=matrix(c(1,2,1,
                  2,1,2,
                  2,2,3,
                  1,1,4),ncol=3,byrow=T)

  # Aplicar as regras
  PertSaida <- sapply(1:nrow(W), function(i) {
    sapply(1:nrow(Regras), function(j) {
      min(c(
        pert_wd[[i, Regras[j, 1] + 1]], # +1 porque la 1ª columna es Gen
        pert_wf[[i, Regras[j, 2] + 1]]
      ))
    })
  })

  PertSaida <- t(PertSaida) # Transpor para que as linhas sejam os genótipos

  GE <- apply(PertSaida[, 3, drop = FALSE], 1, max)
  UNF <- apply(PertSaida[, 2, drop = FALSE], 1, max)
  PA <- apply(PertSaida[, 4, drop = FALSE], 1, max)
  FAV <- apply(PertSaida[, 1, drop = FALSE], 1, max)

  Pertinencias <- data.frame(
    Gen = W$Gen, # Agregar Gen aquí
    GE = GE,
    PA = PA,
    FAV = FAV,
    UNF = UNF
  )

  Resultado <- left_join(W, Pertinencias, by = "Gen")

  saida = Resultado%>%
    mutate(Wg = round(Wg,2))%>%
    mutate(Wd = round(Wd,2))%>%
    mutate(Wf = round(Wf,2))%>%
    mutate(GE = round(GE*100,0))%>%
    mutate(PA = round(PA*100,0))%>%
    mutate(FAV = round(FAV*100,0))%>%
    mutate(UNF = round(UNF*100,0))%>%
    select(Gen,Wg,Wd,Wf,GE,PA,FAV,UNF)

  return(saida)

}
