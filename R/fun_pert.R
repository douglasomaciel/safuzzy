# --- Funções de pertinência ---
zmf = function(x, a, b){
  y = ifelse(x <= a, 1,
             ifelse(x <= (b + a) / 2, 1 - 2 * ((x - a) / (b - a))^2,
                    ifelse(x <= b, 2 * ((x - b) / (b - a))^2, 0)
                    )
             )
  return(y)
}


smf = function(x, a, b){
  y = ifelse(x <= a, 0,
             ifelse(x <= (b + a) / 2, 2 * ((x - a) / (b - a))^2,
                    ifelse(x <= b, 1 - 2 * ((x - b) / (b - a))^2, 1)
                    )
             )
  return(y)
}

pimf = function(x, a, b, c, d){
  y = ifelse(x <= a, 0,
             ifelse(x <= (a + b)/2, 2*((x - a)/(b - a))^2,
                    ifelse(x <= b, 1 - 2*((x - b)/(b - a))^2,
                           ifelse(x <= c, 1,
                                  ifelse(x <= (c + d)/2, 1 - 2*((x - c)/(d - c))^2,
                                         ifelse(x <= d, 2*((x - d)/(d - c))^2, 0)
                                         )
                                  )
                           )
                    )
             )
  return(y)
}
