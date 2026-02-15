# Convert all EDF files in the current folder to ASC files

edf_files <- list.files(pattern = "\\.[Ee][Dd][Ff]$")

for (edf in edf_files) {
  system2("edf2asc", args = c("-y", edf))
}
