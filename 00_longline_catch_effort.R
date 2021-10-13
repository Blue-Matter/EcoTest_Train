
# Use on 32-bit R
get_access_data <- function(file_name, query, ...) {
  require(RODBC)
  ch <- odbcConnectAccess(file_name)
  on.exit(odbcClose(ch))
  out <- sqlQuery(ch, query, ...)
  return(out)
}

dir <- "G:/Shared drives/BM shared/1. Projects/EcoTest/Databases"
file_name <- file.path(dir, "t2ce_20201218web.mdb")
query <- "SELECT * FROM t2ce LEFT JOIN Flags ON t2ce.FleetID = Flags.FleetID WHERE GearGrpCode = 'LL'"

LLCE <- get_access_data(file_name, query)
write.csv(LLCE, file.path(dir, "t2ce_LL.csv"))
