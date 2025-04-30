
dir <- "G:/Shared drives/BM shared/1. Projects/EcoTest/Databases"
t1_file <- "t1nc-20211002.xlsx"

get_t1 <- function(sp) {
  readxl::read_excel(file.path(dir, t1_file), skip = 3)[, 1:22] %>% dplyr::filter(Species == sp)
}
