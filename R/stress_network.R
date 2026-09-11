# Shrink a correlation matrix towards the identity until it is safely
# invertible. Returns the shrinkage intensity actually applied.
shrink_correlation <- function(r, lambda = "auto") {
  p <- ncol(r)
  identity <- diag(p)

  is_usable <- function(m) {
    ev <- try(eigen(m, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
    if (inherits(ev, "try-error")) {
      return(FALSE)
    }
    min(ev) > 1e-8 && (max(ev) / min(ev)) < 1e6
  }

  if (identical(lambda, "auto")) {
    grid <- c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
    for (l in grid) {
      m <- (1 - l) * r + l * identity
      if (is_usable(m)) {
        return(list(matrix = m, lambda = l))
      }
    }
    ps_warn("The correlation matrix stayed ill-conditioned; using lambda = 0.5.")
    return(list(matrix = 0.5 * r + 0.5 * identity, lambda = 0.5))
  }

  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) ||
    lambda < 0 || lambda >= 1) {
    ps_abort('`lambda` must be "auto" or a numeric scalar in [0, 1).')
  }
  m <- (1 - lambda) * r + lambda * identity
  if (!is_usable(m)) {
    ps_warn("The shrunk correlation matrix is ill-conditioned; partial correlations may be unstable.")
  }
  list(matrix = m, lambda = lambda)
}

# Partial correlations from the precision matrix.
partial_correlation <- function(r) {
  theta <- try(solve(r), silent = TRUE)
  if (inherits(theta, "try-error")) {
    ps_abort("The correlation matrix could not be inverted; increase `lambda` or drop collinear traits.")
  }
  d <- sqrt(diag(theta))
  pc <- -theta / outer(d, d)
  diag(pc) <- 1
  dimnames(pc) <- dimnames(r)
  pc
}

#' Partial-Correlation Network of Physiological Traits
#'
#' @description
#' Estimates the conditional dependence structure among physiological traits and
#' detects the modules of traits that respond to stress together. Unlike a
#' marginal correlation network, every edge here is the association between two
#' traits *after removing the effect of all the other traits in the panel*, which
#' removes the dense clusters of indirect correlations that make raw correlation
#' networks unreadable.
#'
#' @details
#' The correlation matrix is shrunk towards the identity,
#' \eqn{R_\lambda = (1 - \lambda) R + \lambda I}, before inversion. With
#' `lambda = "auto"` the smallest intensity on a fixed grid that yields a
#' well-conditioned matrix is used, which keeps the estimate usable when the
#' number of traits approaches the number of plots. Partial correlations are then
#' read off the precision matrix as
#' \eqn{-\theta_{ij} / \sqrt{\theta_{ii}\theta_{jj}}}.
#'
#' Significance uses \eqn{t = \rho \sqrt{(n - p) / (1 - \rho^2)}} on \eqn{n - p}
#' degrees of freedom, where `p` is the number of traits. When the design has
#' fewer observations than that requires, p-values are returned as `NA` and edge
#' selection falls back to `threshold` alone; the shrinkage also makes the test
#' approximate, so `threshold` should be treated as the primary filter.
#'
#' Modules are found on the graph of absolute edge weights. `"louvain"`
#' (multi-level modularity) is the default; `"walktrap"` and `"fast_greedy"` are
#' available for comparison.
#'
#' @param data A data frame with one row per experimental unit.
#' @param traits Character vector of trait columns. Defaults to every numeric
#'   column that is not a design column.
#' @param treatment Optional name of a treatment column used to restrict the
#'   network to particular levels.
#' @param level Optional character vector of `treatment` levels to keep. A
#'   network estimated only on stressed plants describes the stressed
#'   physiology; one estimated on the whole trial mixes both regimes.
#' @param method `"partial"` (default), `"pearson"` or `"spearman"`. The last two
#'   produce a marginal correlation network.
#' @param lambda Shrinkage intensity for the partial-correlation estimate:
#'   `"auto"` or a scalar in `[0, 1)`.
#' @param threshold Minimum absolute correlation retained as an edge.
#' @param alpha Adjusted p-value threshold for retaining an edge; ignored when
#'   p-values cannot be computed.
#' @param p_adjust Multiple-testing adjustment passed to [stats::p.adjust()].
#' @param cluster Module detection algorithm: `"louvain"`, `"walktrap"`,
#'   `"fast_greedy"` or `"none"`.
#' @param sri Optional `plantstress_sri` object. Its indices are averaged per
#'   trait and attached to the node table as `sri`, so that node color can show
#'   how each hub responded to stress.
#'
#' @return An object of class `"plantstress_network"`: a list with `nodes`,
#'   `edges`, `graph` (an [igraph::igraph] object), `modules`, `matrix` (the
#'   partial or marginal correlation matrix) and `parameters`.
#'
#' @references
#' Epskamp S., Fried E.I. (2018). A tutorial on regularized partial correlation
#' networks. \doi{10.1037/met0000167}
#'
#' Blondel V.D., Guillaume J.-L., Lambiotte R., Lefebvre E. (2008). Fast
#' unfolding of communities in large networks.
#' \doi{10.1088/1742-5468/2008/10/P10008}
#'
#' @seealso [plot_stress_signature()], [calculate_sri()]
#'
#' @examples
#' data(brachiaria_stress)
#'
#' net <- stress_network(
#'   brachiaria_stress,
#'   traits = c("Fv_Fm", "PIabs", "A", "gs", "E", "RWC", "SPAD", "shoot_biomass"),
#'   treatment = "drought_level",
#'   level = c("moderate", "severe")
#' )
#' net
#' head(net$edges)
#' @export
stress_network <- function(data,
                           traits = NULL,
                           treatment = NULL,
                           level = NULL,
                           method = c("partial", "pearson", "spearman"),
                           lambda = "auto",
                           threshold = 0.1,
                           alpha = 0.05,
                           p_adjust = "BH",
                           cluster = c("louvain", "walktrap", "fast_greedy", "none"),
                           sri = NULL) {
  method <- match.arg(method)
  cluster <- match.arg(cluster)
  check_data(data)
  check_prob(threshold, "threshold")
  check_prob(alpha, "alpha")
  check_p_adjust(p_adjust)

  if (!is.null(treatment)) {
    check_column(data, treatment, "treatment")
    if (!is.null(level)) {
      keep <- as.character(data[[treatment]]) %in% as.character(level)
      if (!any(keep)) {
        ps_abort(paste0(
          "No row of `", treatment, "` matches `level`: ",
          paste(level, collapse = ", "), "."
        ))
      }
      data <- data[keep, , drop = FALSE]
    }
  }

  traits <- resolve_traits(data, traits, exclude = c(treatment))
  if (length(traits) < 3L) {
    ps_abort("At least three traits are required to build a network.")
  }

  x <- as.matrix(data[traits])
  n <- sum(stats::complete.cases(x))
  if (n < 5L) {
    ps_abort("At least five complete observations are required to estimate a network.")
  }

  cor_method <- if (method == "spearman") "spearman" else "pearson"
  r <- stats::cor(x, use = "pairwise.complete.obs", method = cor_method)
  if (anyNA(r)) {
    ps_abort("The correlation matrix contains missing values; drop traits with too few paired observations.")
  }

  lambda_used <- NA_real_
  if (method == "partial") {
    shrunk <- shrink_correlation(r, lambda)
    lambda_used <- shrunk$lambda
    mat <- partial_correlation(shrunk$matrix)
    df_res <- n - length(traits)
  } else {
    mat <- r
    df_res <- n - 2L
  }

  pairs <- utils::combn(traits, 2L, simplify = FALSE)
  edges <- dplyr::bind_rows(lapply(pairs, function(v) {
    rho <- mat[v[1], v[2]]
    tibble::tibble(from = v[1], to = v[2], correlation = rho)
  }))

  if (df_res >= 1L) {
    tstat <- edges$correlation * sqrt(df_res / pmax(1 - edges$correlation^2, .Machine$double.eps))
    edges$p_value <- 2 * stats::pt(-abs(tstat), df = df_res)
    edges$p_adj <- stats::p.adjust(edges$p_value, method = p_adjust)
  } else {
    ps_warn(paste0(
      "Only ", n, " complete observations for ", length(traits),
      " traits: p-values are not identifiable, edges are selected by `threshold` alone."
    ))
    edges$p_value <- NA_real_
    edges$p_adj <- NA_real_
  }

  edges$sign <- ifelse(edges$correlation >= 0, "positive", "negative")
  keep_edge <- abs(edges$correlation) >= threshold &
    (is.na(edges$p_adj) | edges$p_adj <= alpha)
  kept <- edges[keep_edge, , drop = FALSE]
  kept <- dplyr::arrange(kept, dplyr::desc(abs(.data$correlation)))

  vertices <- tibble::tibble(name = traits)
  graph <- if (nrow(kept) == 0L) {
    ps_warn("No edge passed the selection rules; returning an empty network.")
    igraph::add_vertices(
      igraph::make_empty_graph(directed = FALSE),
      nv = length(traits), name = traits
    )
  } else {
    g <- igraph::graph_from_data_frame(
      d = as.data.frame(kept[c("from", "to", "correlation", "sign", "p_adj")]),
      directed = FALSE,
      vertices = as.data.frame(vertices)
    )
    igraph::E(g)$weight <- abs(kept$correlation)
    g
  }

  membership <- rep(NA_integer_, length(traits))
  names(membership) <- traits
  modularity_value <- NA_real_
  if (cluster != "none" && nrow(kept) > 0L) {
    comm <- switch(cluster,
      louvain = igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight),
      walktrap = igraph::cluster_walktrap(graph, weights = igraph::E(graph)$weight),
      fast_greedy = igraph::cluster_fast_greedy(graph, weights = igraph::E(graph)$weight)
    )
    membership[igraph::V(graph)$name] <- as.integer(igraph::membership(comm))
    modularity_value <- igraph::modularity(comm)
  }

  nodes <- tibble::tibble(
    trait = igraph::V(graph)$name,
    module = unname(membership[igraph::V(graph)$name]),
    degree = as.numeric(igraph::degree(graph)),
    strength = as.numeric(igraph::strength(graph, weights = igraph::E(graph)$weight)),
    betweenness = as.numeric(igraph::betweenness(graph, weights = NA, normalized = TRUE))
  )

  if (!is.null(sri)) {
    if (!is.data.frame(sri) || !all(c("trait", "sri") %in% names(sri))) {
      ps_abort("`sri` must be a `plantstress_sri` object or a data frame with `trait` and `sri`.")
    }
    mean_sri <- dplyr::summarise(
      dplyr::group_by(tibble::as_tibble(sri), .data$trait),
      sri = mean(.data$sri, na.rm = TRUE), .groups = "drop"
    )
    nodes <- dplyr::left_join(nodes, mean_sri, by = "trait")
  }

  nodes <- dplyr::arrange(nodes, dplyr::desc(.data$strength))

  modules <- dplyr::summarise(
    dplyr::group_by(nodes[!is.na(nodes$module), , drop = FALSE], .data$module),
    n_traits = dplyr::n(),
    traits = paste(.data$trait, collapse = ", "),
    hub = .data$trait[which.max(.data$strength)][1],
    .groups = "drop"
  )

  out <- list(
    nodes = nodes,
    edges = kept,
    all_edges = edges,
    graph = graph,
    modules = modules,
    matrix = mat,
    parameters = list(
      method = method,
      lambda = lambda_used,
      threshold = threshold,
      alpha = alpha,
      p_adjust = p_adjust,
      cluster = cluster,
      n_complete = n,
      n_traits = length(traits),
      modularity = modularity_value,
      treatment = treatment,
      level = level
    )
  )
  class(out) <- "plantstress_network"
  out
}

#' @param x A `plantstress_network` object.
#' @param ... Ignored.
#' @rdname stress_network
#' @export
print.plantstress_network <- function(x, ...) {
  p <- x$parameters
  cat("<plantstress_network>\n")
  cat("  Method:      ", p$method,
    if (p$method == "partial") paste0(" (lambda = ", signif(p$lambda, 3), ")") else "",
    "\n",
    sep = ""
  )
  cat("  Traits:      ", p$n_traits, " | complete observations: ", p$n_complete, "\n", sep = "")
  cat("  Edges kept:  ", nrow(x$edges), " of ", nrow(x$all_edges),
    " (|r| >= ", p$threshold, ")\n",
    sep = ""
  )
  cat("  Modules:     ", nrow(x$modules),
    if (is.finite(p$modularity)) paste0(" (modularity = ", signif(p$modularity, 3), ")") else "",
    "\n",
    sep = ""
  )
  if (nrow(x$modules) > 0L) {
    for (i in seq_len(nrow(x$modules))) {
      cat("   - module ", x$modules$module[i], " [hub: ", x$modules$hub[i], "] ",
        x$modules$traits[i], "\n",
        sep = ""
      )
    }
  }
  invisible(x)
}

#' Plot a Trait Stress Network
#'
#' @description
#' Renders a `plantstress_network` as a [ggplot2::ggplot()]. Edge color shows
#' the sign of the partial correlation and edge width its magnitude; node color
#' shows either the detected module or the mean stress response index, when one
#' was attached in [stress_network()].
#'
#' @param x A `plantstress_network` object.
#' @param color_by `"module"` (default) or `"sri"`.
#' @param layout Name of an `igraph` layout function, without the `layout_`
#'   prefix (for example `"fr"`, `"kk"`, `"circle"`).
#' @param seed Random seed, so that force-directed layouts are reproducible.
#' @param label_size Size of the trait labels.
#' @param ... Ignored.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @examples
#' data(brachiaria_stress)
#' net <- stress_network(
#'   brachiaria_stress,
#'   traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD"),
#'   treatment = "drought_level", level = "severe"
#' )
#' plot(net)
#' @export
plot.plantstress_network <- function(x,
                                     color_by = c("module", "sri"),
                                     layout = "fr",
                                     seed = 1L,
                                     label_size = 3.2,
                                     ...) {
  color_by <- match.arg(color_by)
  if (color_by == "sri" && !"sri" %in% names(x$nodes)) {
    ps_abort('Node indices are not available; call `stress_network(sri = )` to attach them.')
  }

  layout_fun <- paste0("layout_with_", layout)
  if (!layout_fun %in% getNamespaceExports("igraph")) {
    layout_fun <- paste0("layout_in_", layout)
    if (!layout_fun %in% getNamespaceExports("igraph")) {
      ps_abort(paste0("Unknown igraph layout: ", layout, "."))
    }
  }

  ps_local_seed(seed)
  coords <- getExportedValue("igraph", layout_fun)(x$graph)
  coords <- as.data.frame(coords)
  names(coords) <- c("x", "y")
  coords$trait <- igraph::V(x$graph)$name

  nodes <- dplyr::left_join(coords, x$nodes, by = "trait")

  edge_df <- x$edges
  if (nrow(edge_df) > 0L) {
    edge_df$x <- nodes$x[match(edge_df$from, nodes$trait)]
    edge_df$y <- nodes$y[match(edge_df$from, nodes$trait)]
    edge_df$xend <- nodes$x[match(edge_df$to, nodes$trait)]
    edge_df$yend <- nodes$y[match(edge_df$to, nodes$trait)]
  }

  p <- ggplot2::ggplot()
  if (nrow(edge_df) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(
        x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
        colour = .data$sign, linewidth = abs(.data$correlation)
      ),
      alpha = 0.65, lineend = "round"
    ) +
      ggplot2::scale_colour_manual(
        values = c(positive = "#1B7837", negative = "#B2182B"),
        name = "Edge sign"
      ) +
      ggplot2::scale_linewidth_continuous(range = c(0.3, 2.2), name = "|correlation|")
  }

  if (color_by == "module") {
    nodes$fill <- factor(nodes$module)
    p <- p +
      ggplot2::geom_point(
        data = nodes,
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$fill, size = .data$strength),
        shape = 21, colour = "grey20", stroke = 0.4, na.rm = TRUE
      ) +
      ggplot2::scale_fill_viridis_d(name = "Module", option = "D", na.value = "grey80")
  } else {
    p <- p +
      ggplot2::geom_point(
        data = nodes,
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$sri, size = .data$strength),
        shape = 21, colour = "grey20", stroke = 0.4, na.rm = TRUE
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#2C7BB6", mid = "grey95", high = "#D7191C",
        midpoint = 0, name = "Mean SRI"
      )
  }

  p +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$trait),
      size = label_size, vjust = -1.4, fontface = "italic"
    ) +
    ggplot2::scale_size_continuous(range = c(2.5, 8), name = "Strength") +
    ggplot2::labs(
      title = "Trait stress network",
      subtitle = paste0(
        x$parameters$method, " correlations, ", nrow(x$edges), " edges, ",
        nrow(x$modules), " modules"
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "right")
}
