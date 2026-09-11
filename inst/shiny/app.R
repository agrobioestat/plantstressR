# plantstressR interactive dashboard.
#
# Every analysis in this file goes through an exported plantstressR function and
# is written with its `plantstressR::` prefix. That is deliberate: the test in
# tests/testthat/test-shiny-app.R extracts those prefixed calls and fails if any
# of them stops being exported, so the app cannot silently rot behind the API
# the way the 0.1.0 dashboard did.

library(shiny)

if (!requireNamespace("plantstressR", quietly = TRUE)) {
  stop("The plantstressR package must be installed to run this app.")
}

NONE <- "— none —"

app_theme <- function() {
  if (requireNamespace("bslib", quietly = TRUE)) {
    bslib::bs_theme(version = 5, primary = "#2C7BB6", base_font = "system-ui")
  } else {
    NULL
  }
}

app_css <- "
  .ps-hint { color: #6b7280; font-size: 0.85rem; margin-top: -6px; }
  .ps-step { font-weight: 600; margin-top: 14px; }
  .ps-panel { padding-top: 14px; }
  .shiny-output-error-validation { color: #b45309; font-weight: 500; }
"

# UI ---------------------------------------------------------------------

ui <- fluidPage(
  theme = app_theme(),
  tags$head(tags$style(HTML(app_css))),
  titlePanel("plantstressR — physiological stress signatures"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "ps-step", "1. Data"),
      radioButtons(
        "source", NULL,
        choices = c("Example trial" = "example", "Upload a CSV" = "upload"),
        selected = "example"
      ),
      conditionalPanel(
        "input.source == 'upload'",
        fileInput("file", "CSV file", accept = c(".csv", ".txt")),
        checkboxInput("semicolon", "Semicolon separated / decimal comma", FALSE)
      ),
      div(class = "ps-step", "2. Design"),
      selectInput("treatment", "Treatment column", choices = NULL),
      selectInput("control", "Control level", choices = NULL),
      selectInput("by", "Group by (genotype)", choices = NULL),
      selectInput("block", "Block column", choices = NULL),
      div(
        class = "ps-hint",
        "A block is a complete replicate of the trial: a bench, a strip, a run."
      ),
      div(class = "ps-step", "3. Traits"),
      selectizeInput("traits", NULL,
        choices = NULL, multiple = TRUE,
        options = list(plugins = list("remove_button"))
      ),
      div(class = "ps-step", "4. Options"),
      selectInput("method", "Effect size",
        choices = c("glass", "cohen", "hedges", "relative")
      ),
      selectInput("weights", "Index weighting",
        choices = c("equal", "precision", "pca")
      ),
      br(),
      actionButton("run", "Run analysis", class = "btn-primary", width = "100%")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "Data",
          div(
            class = "ps-panel",
            h4("Design checks"),
            verbatimTextOutput("checks"),
            h4("Preview"),
            tableOutput("preview")
          )
        ),
        tabPanel(
          "Signature",
          div(
            class = "ps-panel",
            fluidRow(
              column(4, selectInput("sig_type", "Plot", c("heatmap", "radar"))),
              column(4, selectInput("sig_unit", "Unit", choices = NULL)),
              column(4, checkboxInput("sig_sig", "Significant traits only", FALSE))
            ),
            plotOutput("signature_plot", height = "460px"),
            h4("Stress Response Index"),
            downloadButton("dl_sri", "Download table"),
            br(), br(),
            tableOutput("sri_table")
          )
        ),
        tabPanel(
          "Ranking",
          div(
            class = "ps-panel",
            fluidRow(
              column(6, selectInput("isi_type", "Plot", c("ranking", "contribution")))
            ),
            plotOutput("isi_plot", height = "420px"),
            h4("Integrated Stress Index"),
            downloadButton("dl_isi", "Download table"),
            br(), br(),
            tableOutput("isi_table"),
            h4("Trait weights actually used"),
            tableOutput("weights_table")
          )
        ),
        tabPanel(
          "Network",
          div(
            class = "ps-panel",
            fluidRow(
              column(4, numericInput("net_threshold", "|r| threshold", 0.1, 0, 1, 0.05)),
              column(4, selectInput("net_cluster", "Modules",
                c("louvain", "walktrap", "fast_greedy", "none")
              )),
              column(4, selectInput("net_levels", "Levels", choices = NULL, multiple = TRUE))
            ),
            plotOutput("network_plot", height = "460px"),
            h4("Modules"),
            tableOutput("modules_table"),
            h4("Module scores"),
            tableOutput("module_scores_table")
          )
        ),
        tabPanel(
          "Ordination",
          div(
            class = "ps-panel",
            plotOutput("ordination_plot", height = "480px"),
            tableOutput("eigen_table")
          )
        ),
        tabPanel(
          "Tolerance indices",
          div(
            class = "ps-panel",
            fluidRow(
              column(6, selectInput("sti_trait", "Productivity trait", choices = NULL)),
              column(6, selectInput("sti_indices", "Indices",
                choices = c(
                  "TOL", "MP", "GMP", "HM", "SSI",
                  "STI", "YI", "YSI", "RDI", "SSPI"
                ),
                selected = c("STI", "GMP", "SSI", "TOL"), multiple = TRUE
              ))
            ),
            div(
              class = "ps-hint",
              "Needs a grouping column: these indices compare units within a trial."
            ),
            br(),
            downloadButton("dl_sti", "Download table"),
            br(), br(),
            tableOutput("sti_table")
          )
        )
      )
    )
  )
)

# Server -----------------------------------------------------------------

server <- function(input, output, session) {
  # Data ------------------------------------------------------------------

  dataset <- reactive({
    if (identical(input$source, "example")) {
      return(plantstressR::simulate_brachiaria_stress())
    }
    file <- input$file
    validate(need(file, "Choose a CSV file to begin."))
    read_fun <- if (isTRUE(input$semicolon)) utils::read.csv2 else utils::read.csv
    out <- try(read_fun(file$datapath, stringsAsFactors = FALSE), silent = TRUE)
    validate(need(!inherits(out, "try-error"), "That file could not be read as CSV."))
    validate(need(nrow(out) > 0, "The file has no rows."))
    out
  })

  numeric_cols <- reactive({
    dat <- dataset()
    names(dat)[vapply(dat, is.numeric, logical(1))]
  })

  # Keep the design pickers in step with whatever table is loaded.
  observeEvent(dataset(), {
    dat <- dataset()
    all_cols <- names(dat)
    nums <- numeric_cols()
    categorical <- setdiff(all_cols, nums)

    trt_guess <- if (length(categorical) > 0) categorical[length(categorical)] else all_cols[1]
    updateSelectInput(session, "treatment", choices = all_cols, selected = trt_guess)
    updateSelectInput(session, "by",
      choices = c(NONE, categorical),
      selected = if (length(categorical) > 1) categorical[1] else NONE
    )
    updateSelectInput(session, "block", choices = c(NONE, categorical), selected = NONE)
    updateSelectizeInput(session, "traits",
      choices = nums,
      selected = utils::head(nums, 8), server = FALSE
    )
    updateSelectInput(session, "sti_trait", choices = nums,
      selected = utils::tail(nums, 1)
    )
  })

  observeEvent(input$treatment, {
    dat <- dataset()
    req(input$treatment %in% names(dat))
    levels_present <- sort(unique(as.character(dat[[input$treatment]])))
    updateSelectInput(session, "control",
      choices = levels_present, selected = levels_present[1]
    )
    updateSelectInput(session, "net_levels",
      choices = levels_present, selected = levels_present
    )
  })

  opt <- function(x) if (is.null(x) || identical(x, NONE) || !nzchar(x)) NULL else x

  design <- reactive({
    list(
      treatment = input$treatment,
      control = input$control,
      by = opt(input$by),
      block = opt(input$block),
      traits = input$traits
    )
  })

  output$preview <- renderTable(
    {
      utils::head(dataset(), 8)
    },
    striped = TRUE, spacing = "xs", digits = 3
  )

  output$checks <- renderPrint({
    d <- design()
    req(d$treatment, d$control)
    validate(need(length(d$traits) > 0, "Select at least one trait."))
    res <- try(
      plantstressR::validate_stress_data(
        dataset(),
        treatment = d$treatment, control = d$control,
        traits = d$traits, by = d$by, block = d$block, verbose = FALSE
      ),
      silent = TRUE
    )
    if (inherits(res, "try-error")) {
      cat("The design is not usable yet:\n\n")
      cat(conditionMessage(attr(res, "condition")), "\n")
    } else {
      print(res)
    }
  })

  # Analyses --------------------------------------------------------------
  # Nothing recomputes until the button is pressed, so the app never thrashes
  # while the design pickers are still being filled in.

  guard <- function(expr, what) {
    out <- try(expr, silent = TRUE)
    if (inherits(out, "try-error")) {
      showNotification(
        paste0(what, ": ", conditionMessage(attr(out, "condition"))),
        type = "error", duration = 10
      )
      return(NULL)
    }
    out
  }

  sri <- eventReactive(input$run, {
    d <- design()
    validate(need(length(d$traits) > 0, "Select at least one trait."))
    guard(
      plantstressR::calculate_sri(
        dataset(),
        treatment = d$treatment, control = d$control, traits = d$traits,
        by = d$by, block = d$block, method = input$method, verbose = FALSE
      ),
      "Stress Response Index"
    )
  })

  observeEvent(sri(), {
    req(sri())
    units <- unique(sri()$unit)
    updateSelectInput(session, "sig_unit", choices = units, selected = units[1])
  })

  isi <- reactive({
    req(sri())
    guard(
      plantstressR::integrated_stress_index(sri(), weights = input$weights),
      "Integrated index"
    )
  })

  network <- reactive({
    req(sri())
    d <- design()
    validate(need(length(d$traits) >= 3, "A network needs at least three traits."))
    guard(
      plantstressR::stress_network(
        dataset(),
        traits = d$traits, treatment = d$treatment, level = input$net_levels,
        block = d$block, threshold = input$net_threshold,
        cluster = input$net_cluster, sri = sri()
      ),
      "Network"
    )
  })

  ordination <- reactive({
    req(sri())
    d <- design()
    validate(need(length(d$traits) >= 2, "An ordination needs at least two traits."))
    guard(
      plantstressR::stress_ordination(
        dataset(),
        traits = d$traits, treatment = d$treatment, block = d$block
      ),
      "Ordination"
    )
  })

  # Signature -------------------------------------------------------------

  output$signature_plot <- renderPlot({
    req(sri())
    p <- guard(
      plantstressR::plot_stress_signature(
        sri(),
        type = input$sig_type,
        units = input$sig_unit,
        significant_only = isTRUE(input$sig_sig)
      ),
      "Signature plot"
    )
    req(p)
    p
  })

  output$sri_table <- renderTable(
    {
      req(sri())
      out <- as.data.frame(sri())
      out[c("unit", "group", "trait", "sri", "conf_low", "conf_high", "p_adj")]
    },
    striped = TRUE, spacing = "xs", digits = 3
  )

  # Ranking ---------------------------------------------------------------

  output$isi_plot <- renderPlot({
    req(isi())
    p <- guard(plot(isi(), type = input$isi_type), "Ranking plot")
    req(p)
    p
  })

  output$isi_table <- renderTable(
    {
      req(isi())
      as.data.frame(isi())
    },
    striped = TRUE, spacing = "xs", digits = 3
  )

  output$weights_table <- renderTable(
    {
      req(isi())
      as.data.frame(plantstressR::stress_weights(isi()))
    },
    striped = TRUE, spacing = "xs", digits = 4
  )

  # Network ---------------------------------------------------------------

  output$network_plot <- renderPlot({
    req(network())
    p <- guard(plot(network(), color_by = "sri"), "Network plot")
    req(p)
    p
  })

  output$modules_table <- renderTable(
    {
      req(network())
      as.data.frame(network()$modules)
    },
    striped = TRUE, spacing = "xs"
  )

  output$module_scores_table <- renderTable(
    {
      req(network(), sri())
      out <- guard(
        plantstressR::stress_module_scores(network(), sri()),
        "Module scores"
      )
      req(out)
      as.data.frame(out)
    },
    striped = TRUE, spacing = "xs", digits = 3
  )

  # Ordination ------------------------------------------------------------

  output$ordination_plot <- renderPlot({
    req(ordination())
    p <- guard(plot(ordination()), "Ordination plot")
    req(p)
    p
  })

  output$eigen_table <- renderTable(
    {
      req(ordination())
      as.data.frame(ordination()$eigenvalues)
    },
    striped = TRUE, spacing = "xs", digits = 2
  )

  # Tolerance indices -----------------------------------------------------

  sti <- reactive({
    d <- design()
    validate(need(
      !is.null(d$by),
      "Choose a grouping column: tolerance indices compare units within a trial."
    ))
    req(input$sti_trait, length(input$sti_indices) > 0)
    guard(
      plantstressR::stress_tolerance_index(
        dataset(),
        trait = input$sti_trait, treatment = d$treatment, control = d$control,
        by = d$by, indices = input$sti_indices
      ),
      "Tolerance indices"
    )
  })

  output$sti_table <- renderTable(
    {
      req(sti())
      as.data.frame(sti())
    },
    striped = TRUE, spacing = "xs", digits = 3
  )

  # Downloads -------------------------------------------------------------

  csv_handler <- function(name, value) {
    downloadHandler(
      filename = function() paste0("plantstressR_", name, ".csv"),
      content = function(path) {
        out <- value()
        req(out)
        utils::write.csv(as.data.frame(out), path, row.names = FALSE)
      }
    )
  }

  output$dl_sri <- csv_handler("sri", sri)
  output$dl_isi <- csv_handler("isi", isi)
  output$dl_sti <- csv_handler("tolerance_indices", sti)
}

shinyApp(ui, server)
