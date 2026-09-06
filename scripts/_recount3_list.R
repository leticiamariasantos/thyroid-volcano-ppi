suppressPackageStartupMessages(library(recount3))
ap <- available_projects()
cat("total projects:", nrow(ap), "\n")
# find thyroid-relevant projects without shell-hostile pipes
hits <- ap[grepl("THCA", ap$project) | grepl("GTEX", ap$project) | grepl("THYROID", ap$project, ignore.case = TRUE), ]
print(hits)
cat("\n--- columns ---\n")
print(colnames(ap))
