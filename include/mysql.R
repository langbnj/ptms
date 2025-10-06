library(RMySQL)
library(tibble)

if (exists("superreservedcon"))
{
  try(dbDisconnect(superreservedcon), silent=T)
  rm(superreservedcon)
  # rm(superreserveddrv)
}

# superreserveddrv = dbDriver("MySQL")

superreservedcon = dbConnect(dbDriver("MySQL"), host="[sql server]", dbname="blang")

Query <- function(query)
{
  # a = dbGetQuery(superreservedcon, statement=query)
  # a = suppressWarnings(dbGetQuery(superreservedcon, statement=query))
  a = as_tibble(dbGetQuery(superreservedcon, statement=query))
}
