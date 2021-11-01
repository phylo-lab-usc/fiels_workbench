#functions used

#Used to convert arbutus pvals into useful data frames or tibbles
arbutus_transform <- function ( pval , len , tib) {
  df <- t(data.frame(pval))
  if(missing(tib)) tib = FALSE
  ifelse(tib == TRUE, {
    fin <- as_tibble(df)
  }
  , fin <- data.frame(df))
  row.names(fin) <- c(1:len)
  fin
}
