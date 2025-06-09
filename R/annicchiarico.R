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
#' @import purrr
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
    rename(Amb = {{env}}, Gen = {{gen}}, Rep = {{rep}}, Yvar = {{var}}) %>%
    mutate(across(c(Amb, Gen, Rep), as.factor))

  media_amb <- Dados %>%
    group_by(Amb) %>%
    summarise(Yvar_amb = mean(Yvar, na.rm = TRUE), .groups = "drop")

  IJ <- media_amb %>%
    mutate(IJ = Yvar_amb - mean(Yvar_amb)) %>%
    select(Amb, IJ)

  Y <- Dados %>%
    group_by(Gen, Amb) %>%
    summarise(Yvar = mean(Yvar, na.rm = TRUE), .groups = "drop") %>%
    left_join(media_amb, by = "Amb") %>%
    mutate(Yp = 100 * Yvar / Yvar_amb) %>%
    select(Gen, Amb, Yp)

  calc_W <- function(df) {
    df %>%
      group_by(Gen) %>%
      summarise(
        mean_Yp = mean(Yp, na.rm = TRUE),
        sd_Yp = sd(Yp, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(W = mean_Yp - sd_Yp * qnorm(0.75)) %>%
      select(Gen, W)
  }

  Wg <- calc_W(Y) %>% rename(Wg = W)

  desfav <- IJ %>% filter(IJ <= 0) %>% pull(Amb)
  favorav <- IJ %>% filter(IJ > 0) %>% pull(Amb)

  Wd <- calc_W(Y %>% filter(Amb %in% desfav)) %>% rename(Wd = W)
  Wf <- calc_W(Y %>% filter(Amb %in% favorav)) %>% rename(Wf = W)

  W <- reduce(list(Wg, Wd, Wf), left_join, by = "Gen")

# Fuzificação
  pert_wd <- W%>%mutate(
    baixa = zmf(Wd,0,200),
    alta = smf(Wd,0,200)
  )%>%select(Gen, baixa, alta)

  pert_wf <- W%>%mutate(
    baixa = zmf(Wf,0,200),
    alta = smf(Wf,0,200)
  )%>%select(Gen, baixa, alta)

  Regras=matrix(c(1,2,1,
                  2,1,2,
                  2,2,3,
                  1,1,4),ncol=3,byrow=T)

  # Aplicar as regras
  PertSaida <- sapply(1:nrow(W), function(i) {
    sapply(1:nrow(Regras), function(j) {
      min(c(
        pert_wd[[i, Regras[j, 1]+1]],
        pert_wf[[i, Regras[j, 2]+1]]
      ))
    })
  }) %>% t()

  Pertinencias <- tibble(
    Gen = W$Gen, # Agregar Gen aqui
    GE = apply(PertSaida[, 3, drop = FALSE], 1, max),
    PA = apply(PertSaida[, 4, drop = FALSE], 1, max),
    FAV = apply(PertSaida[, 1, drop = FALSE], 1, max),
    UNF = apply(PertSaida[, 2, drop = FALSE], 1, max)
  )

  Resultado <- left_join(W, Pertinencias, by = "Gen")

  saida = Resultado%>%
    mutate(
      across(c(Wg,Wd,Wf),~round(.x,2)),
      across(c(GE,PA,FAV,UNF),~round(.x*100,0))
    )%>%
    select(Gen,Wg,Wd,Wf,GE,PA,FAV,UNF)

  return(saida)

}
