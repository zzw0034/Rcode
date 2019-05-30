
library(magrittr)

ori_dir <- setwd("S:/hao_shi/gcp18/trendyv7/output1700/gridpool/")
sub_dirs <- list.dirs()[-1]

for(sub_dir in sub_dirs) {
  cat("Mask: ", sub_dir, "\n")
  files_1700_full <- list.files(path=sub_dir, full.names = TRUE)
  files_1700 <- list.files(path=sub_dir)
  for(scen in paste0("s", 0:4)) {
    dest_dir <- paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/",
                       sub("\\./", "", sub_dir), scen)
    dest_files <- sub("s0", scen, files_1700) %>% file.path(dest_dir, .)
    file.copy(files_1700_full, dest_files)
    cat("Scen: ", scen, "\n")
  }
}
