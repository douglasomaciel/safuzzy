#' Cruz, Torres & Vencovsky (1989) method
#'
#' @description
#' Adaptability and Stability Analysis based on the interpretation of the Cruz, Torres & Vencovsky (1989), developed by Carneiro et al. 2019.
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
#'
#' @return A data frame containing the following estimates:
#'   \itemize{
#'     \item{\code{Gen}}{: Genotype.}
#'     \item{\code{B_0}}{: Mean of the variable for each genotype.}
#'     \item{\code{B_1}}{: Regression coefficient (B_1) for each genotype.}
#'     \item{\code{B1_B2}}{: Regression coefficient (B1+B2) for each genotype.}
#'     \item{\code{R2}}{: Coefficient of determination (R²) for each genotype.}
#'     \item{\code{MdAF}}{: Membership (\%) to the average adaptability to favorable environments genotypes group.}
#'     \item{\code{Nda}}{: Membership (\%) to the  poorly adapted geotypes group.}
#'     \item{\code{MdAG}}{: Membership (\%) to the general adaptability to favorable environments genotypes group.}
#'     \item{\code{MaxGF}}{: Membership (\%) to the maximum adaptability to favorable environments genotypes group.}
#'     \item{\code{MaxDes}}{: Membership (\%) to the maximum adaptability to unfavorable environments genotypes group.}
#'     \item{\code{BE}}{: Membership (\%) to the low stability genotypes group.}
#'     \item{\code{BP}}{: Membership (\%) to the low yield genotypes group.}
#'   }
#'
#'
#' @seealso \code{\link{hybrid}}
#' @seealso \code{\link{lin_binns}}
#' @seealso \code{\link{eberhart_russell}}
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
#' cruz_torres_vencovsky(data = dados, env = environment, gen = treatment, rep = block, var = gy)

#' @export

cruz_torres_vencovsky <- function(data, env, gen, rep, var){
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
    mutate(TIJ = ifelse(IJ<=0,0,IJ-mean(IJ[IJ>0])))%>%
    select(Amb, IJ, TIJ)

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
      b = coef(lm(Yvar ~ IJ + TIJ))[2], # Coeficiente de regressão (b_i)
      b2 = coef(lm(Yvar~IJ + TIJ))[3],
      R2 = summary(lm(Yvar ~ IJ))$r.squared # R-quadrado (R^2_i)
      # Para o desvio da regressão (S^2_di), seria mais complexo e envolveria os resíduos
    )%>%
    mutate(b1_b2 = b+b2)

  # modelo
  mod = summary(aov(Yvar ~ Gen+Rep/Amb+Amb+Gen:Amb, data = Dados))
  nrep <- length(unique(Dados$Rep)) # Calcular nrep
  namb <- length(unique(Dados$Amb))
  QMR=mod[[1]]$`Mean Sq`[6]
  GLR=mod[[1]]$`Df`[6]
  se=QMR/nrep
  VarBo=se*sum(IJ$IJ^2)/(namb*sum(IJ$IJ^2))
  VarB1=se/sum(IJ$IJ^2)
  VarB2=QMR*sum(IJ$IJ^2)/(nrep*sum(IJ$TIJ^2)*(sum(IJ$IJ^2)-sum(IJ$TIJ^2)))
  VarB1B2=QMR/(nrep*sum(IJ$TIJ^2))

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

  # beta 1 + beta 2 padronizado
  LI=-2
  LS=4
  BMax=1+qt(0.975,GLR)*sqrt(VarB1B2)
  BMin=1+qt(0.025,GLR)*sqrt(VarB1B2)

  reg <- reg %>%
    mutate(b1_b2_pad = LS + (LS - LI) * (b1_b2 - BMax) / (BMax - BMin))

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

  pert_b1_b2 <- reg %>%
    mutate(b1_b2_menor = zmf(b1_b2_pad, -5, 1),
           b1_b2_igual = pimf(b1_b2_pad, -5, 1, 1, 7),
           b1_b2_maior = smf(b1_b2_pad, 1, 7)) %>%
    select(Gen, b1_b2_menor, b1_b2_igual, b1_b2_maior) # Incluir Gen

  pert_r <- reg %>%
    mutate(r_baixa = zmf(r_pad, 60, 100),
           r_alta = smf(r_pad, 60, 100)) %>%
    select(Gen, r_baixa, r_alta) # Incluir Gen

  Regras=matrix(c(2,2,2,2,1,
                  2,2,3,2,2,
                  2,2,1,2,3,
                  2,3,2,2,1,
                  2,3,3,2,2,
                  2,3,1,2,3,
                  2,1,2,2,4,
                  2,1,3,2,5,
                  2,1,1,2,6,
                  2,2,2,1,7,
                  2,2,3,1,7,
                  2,2,1,1,7,
                  2,3,2,1,7,
                  2,3,3,1,7,
                  2,3,1,1,7,
                  2,1,2,1,7,
                  2,1,3,1,7,
                  2,1,1,1,7,
                  1,2,2,1,8,
                  1,2,2,2,8,
                  1,2,3,1,8,
                  1,2,3,2,8,
                  1,2,1,1,8,
                  1,2,1,2,8,
                  1,3,2,1,8,
                  1,3,2,2,8,
                  1,3,3,1,8,
                  1,3,3,2,8,
                  1,3,1,1,8,
                  1,3,1,2,8,
                  1,1,2,1,8,
                  1,1,2,2,8,
                  1,1,3,1,8,
                  1,1,3,2,8,
                  1,1,1,1,8,
                  1,1,1,2,8),ncol=5,byrow=T)

  # Aplicar as regras
  PertSaida <- sapply(1:nrow(reg), function(i) {
    sapply(1:nrow(Regras), function(j) {
      min(c(
        pert_med[[i, Regras[j, 1] + 1]], # +1 porque la 1ª columna es Gen
        pert_b[[i, Regras[j, 2] + 1]],   # +1 porque la 1ª columna es Gen
        pert_b1_b2[[i, Regras[j, 3] + 1]],    # +1 porque la 1ª columna es Gen
        pert_r[[i, Regras[j, 4] + 1]]    # +1 porque la 1ª columna es Gen
      ))
    })
  })

  PertSaida <- t(PertSaida) # Transpor para que as linhas sejam os genótipos

  MdAF <- apply(PertSaida[, c(1,4), drop = FALSE], 1, max)
  MaxAF <- apply(PertSaida[, c(2,5), drop = FALSE], 1, max)
  Nad <- apply(PertSaida[, c(3,6), drop = FALSE], 1, max)
  MdAG <- apply(PertSaida[, 7, drop = FALSE], 1, max)
  MaxAG <- apply(PertSaida[, 8, drop = FALSE], 1, max)
  MaxDes <- apply(PertSaida[, 9, drop = FALSE], 1, max)
  BE <- apply(PertSaida[, 10:18, drop = FALSE], 1, max)
  BP <- apply(PertSaida[, 19:36, drop = FALSE], 1, max)

  Pertinencias <- data.frame(
    Gen = reg$Gen,
    MdAF = MdAF,
    Nad = Nad,
    MdAG = MdAG,
    MaxAG = MaxAG,
    MaxDes = MaxDes,
    BE = BE,
    BP = BP
  )

  Resultado <- left_join(reg, Pertinencias, by = "Gen")

  saida = Resultado%>%
    mutate(B_0 = round(med,2))%>%
    mutate(B_1 = round(b,4))%>%
    mutate(B1_B2 = round(b1_b2,4))%>%
    mutate(R2 = r_pad,2)%>%
    mutate(MdAF = round(MdAF*100,0))%>%
    mutate(Nad = round(Nad*100,0))%>%
    mutate(MdAG = round(MdAG*100,0))%>%
    mutate(MaxAG = round(MaxAG*100,0))%>%
    mutate(MaxDes = round(MaxDes*100,0))%>%
    mutate(BE = round(BE*100,0))%>%
    mutate(BP = round(BP*100,0))%>%
    select(Gen,B_0,B_1,B1_B2,R2,MdAF,Nad,MdAG,MaxAG,MaxDes,BE,BP)


  return(saida)

}
