#' Launch the plantstressR Dashboard
#'
#' @description
#' Opens a Shiny application that walks through the whole package workflow
#' without writing any code: load a table, declare the design, and read the
#' signature, the ranking, the trait network, the ordination and the classical
#' tolerance indices from tabs.
#'
#' The app is a front end, not a second implementation. Every number and every
#' figure it shows comes from the same exported functions documented in this
#' manual, so anything seen on screen can be reproduced from a script.
#'
#' @details
#' `shiny` is an optional dependency. Install it with `install.packages("shiny")`
#' if the app refuses to start. `bslib`, if present, is used for the theme.
#'
#' @param launch.browser Logical. Open the app in the default browser.
#' @param ... Further arguments passed to [shiny::runApp()], such as `port` or
#'   `host`.
#'
#' @return Invisibly `NULL`; called for the side effect of running the app.
#'
#' @seealso [calculate_sri()], [integrated_stress_index()], [stress_network()]
#'
#' @examples
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   run_plantstress_app()
#' }
#' @export
run_plantstress_app <- function(launch.browser = TRUE, ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    ps_abort(
      "Package `shiny` is required for the dashboard. Install it with `install.packages(\"shiny\")`."
    )
  }

  app_dir <- system.file("shiny", package = "plantstressR")
  if (!nzchar(app_dir) || !file.exists(file.path(app_dir, "app.R"))) {
    ps_abort("The bundled Shiny app could not be found; reinstall `plantstressR`.")
  }

  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
  invisible(NULL)
}
