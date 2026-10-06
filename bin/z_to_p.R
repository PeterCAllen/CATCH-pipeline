#!/usr/bin/env Rscript

con_in  <- file("stdin", "r")
con_out <- stdout()

header <- readLines(con_in, n = 1)
if (length(header) == 0) quit(status = 0)

fields <- toupper(strsplit(trimws(header), "[ \t]+")[[1]])
z_idx <- which(fields %in% c("Z", "ZSCORE", "Z_SCORE", "ZSTAT"))
if (length(z_idx) == 0) stop("ERROR: no Z column found in header: ", header)
z_idx <- z_idx[1]

fields[z_idx] <- "P"
writeLines(paste(fields, collapse = "\t"), con_out)

while (length(line <- readLines(con_in, n = 1)) > 0) {
  parts <- strsplit(trimws(line), "[ \t]+")[[1]]
  z <- suppressWarnings(as.numeric(parts[z_idx]))
  parts[z_idx] <- if (is.na(z)) "NA" else format(2 * pnorm(-abs(z)), scientific = TRUE, digits = 6)
  writeLines(paste(parts, collapse = "\t"), con_out)
}

close(con_in)
