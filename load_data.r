load_sc_data <- function(data_type) {
  data_path <- switch(data_type,
                      "rna" = "/home/big/Public/astro/0_data/WTsscRNA.rds",
                      "m6A" = "/home/big/Public/astro/0_data/WTsscm6A.rds",
  )
  readRDS(file.path(data_path))
}
