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

# The same trial laid out in complete blocks, with a large additive block shift
# on top of the treatment effect. Any routine that handles blocks correctly must
# recover the treatment difference unchanged while removing the block spread.
toy_block_trial <- function(n_per_block = 4, n_blocks = 4, block_shift = 5) {
  base <- expand.grid(
    rep = seq_len(n_per_block),
    block = paste0("B", seq_len(n_blocks)),
    trt = c("ctrl", "stress"),
    genotype = c("A", "B"),
    stringsAsFactors = FALSE
  )
  shift <- block_shift * (as.integer(factor(base$block)) - 1)
  trt_effect <- ifelse(base$trt == "stress", 1, 0)
  geno_scale <- ifelse(base$genotype == "A", 2, 1)
  jitter <- (base$rep - mean(seq_len(n_per_block))) / 10

  data.frame(
    genotype = base$genotype,
    block = base$block,
    trt = base$trt,
    up = 10 + 3 * trt_effect * geno_scale + shift + jitter,
    down = 10 - 3 * trt_effect * geno_scale + shift + jitter,
    stringsAsFactors = FALSE
  )
}
