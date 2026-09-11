# A tiny, fully deterministic trial used by most tests. Trait `up` rises under
# stress, `down` falls, `flat` is constant everywhere, and genotype A is pushed
# twice as hard as genotype B, so every expected sign and ordering is known in
# advance.
toy_trial <- function(n = 8) {
  data.frame(
    genotype = rep(c("A", "B"), each = 2 * n),
    trt = rep(rep(c("ctrl", "stress"), each = n), 2),
    up = c(
      seq(10, 11, length.out = n), seq(14, 15, length.out = n),
      seq(10, 11, length.out = n), seq(12, 13, length.out = n)
    ),
    down = c(
      seq(5, 6, length.out = n), seq(2, 3, length.out = n),
      seq(5, 6, length.out = n), seq(4, 5, length.out = n)
    ),
    flat = rep(1, 4 * n),
    stringsAsFactors = FALSE
  )
}
