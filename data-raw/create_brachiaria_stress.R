# Regenerates data/brachiaria_stress.rda.
# Run from the package root with: source("data-raw/create_brachiaria_stress.R")

pkgload::load_all(".", quiet = TRUE)

brachiaria_stress <- simulate_brachiaria_stress(n = 120, seed = 123)

save(
  brachiaria_stress,
  file = file.path("data", "brachiaria_stress.rda"),
  compress = "bzip2",
  version = 3
)
