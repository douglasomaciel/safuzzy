# --- Padronização ---
padronizar <- function(x, LI, LS) {
  BMax <- mean(x) + 3 * sd(x)
  BMin <- mean(x) - 3 * sd(x)
  sapply(x, function(v) LS + (LS - LI) * (v - BMax) / (BMax - BMin))
}
