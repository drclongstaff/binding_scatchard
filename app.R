# An app to analyse Michaelis Menten data

library(shiny)
library(tidyplots)

ThisApp <- "Michaelis Menten curves and linear transformations"
ThisVersion <- 0.5

# Load functions
source("./Functions/Functions_LoadandPlots.R")

# Define UI for the application
ui <- fluidPage(
  includeCSS("./www/styles3.css"), # make a few changes to the colours and fonts
  # supress unnecessary warnings
  # tags$style(type="text/css",
  #   ".shiny-output-error { visibility: hidden; }",
  # ".shiny-output-error:before { visibility: hidden; }"
  # ),

  titlePanel(h2(ThisApp, align = "center")),
  sidebarLayout(
    sidebarPanel(
      tags$h4("Load your data"),
      fluidRow(
        # Will call a function in the server to detect and load the user's data file
        column(6, fileInput("data", label = "Select data", accept = c(".csv", ".txt", ".xlsx"))),
        column(5, numericInput("sheet", "Excel sheet", value = 1, min = 1, step = 1))
      ),
      uiOutput("whatx"),
      uiOutput("whaty"),
      radioButtons(
        inputId = "raw", label = tags$h4("Select plot"),
        choices = c(
          "Non-linear",
          "Hanes",
          "Lineweaver-Burk",
          "Eadie-Hofstee",
          "Scatchard"
        ),
        selected = "Non-linear"
      ),
      tags$i("Please contact me with any comments on:"),
      helpText(h5(
        ThisApp, "version", ThisVersion, " last accessed", Sys.Date(), "at",
        tags$a(href = "mailto: drclongstaff@gmail.com", "drclongstaff@gmail.com")
      )),
      tags$br(),
      "Code can be found on my github site:",
      tags$a(href = "https://github.com/drclongstaff/", "here"),
      tags$br(),
      "Other apps and links for reproducible analysis in haemostasis assays are available",
      tags$a(href = "https://drclongstaff.github.io/shiny-clots/", "here")
    ),
    mainPanel(
      tabsetPanel(
        type = "tab",
        tabPanel("Plot",
          align = "center",
          plotOutput(outputId = "myplot"),
          h4(textOutput("text1")),
          h4("Results Table"),
          tableOutput("resultsTable"),
          helpText(h5("*Scatchard plots included for binding assays where V=Bound and S=Free ligand")),
          helpText(h5("In this case Vmax=Max bound and Km=Kd"))
        ),
        tabPanel("Raw data", DT::DTOutput("contents")),
        tabPanel("Transformed data", DT::DTOutput("resultsTable2")),
        tabPanel(
          "Help",
          tags$blockquote(h5(
            "►The app opens with an Excel data file from assays of plasminogen activation by streptokinase",
            tags$br(),
            "►Three independent sets of data are provided, A, B and C over a range of [plasminogen] substrate concentrations",
            tags$br(),
            "►There are 3 sheets: 1) All replicates; 2) Means; 3) 3 points [Pgn] <km, ~Km and approaching Vmax",
            tags$br(),
            "►You can explore nonlinear fits or linear transformations of these data sets",
            tags$br(),
            "►Load your own data for fitting as csv, txt or xlsx files (the app will detect the format)",
            tags$br(),
            "►The supplied data shows the expected data layout",
            tags$br(),
            "►An option is provided to show a Scatchard plot which is used in receptor-ligand binding analysis",
            tags$br(),
            "►However, the same nonlinear curve fitting of binding assays can be used to find Bmax and Kd (equivalent to Vmax and Km)",
            tags$br(),
            "►You can also see that the Scatchard plot is related to the Eadie-Hofstee plot, where the x and y axis have been switched"
          )),
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  readData <- reactive({
    inputFile <- input$data
    if (is.null(inputFile)) {
      readxl::read_excel("./Data/SK_MMdata.xlsx", sheet = input$sheet) |> as.data.frame()
    } else {
      req(inputFile)
      (
        load_file(input$data$name, input$data$datapath, input$sheet) |>
          janitor::remove_empty(
            which = c("rows", "cols"),
            cutoff = 1, quiet = TRUE
          ) |> # remove empty cols and rows
          sapply(\(x) replace(x, x %in% "", NA)) |> # replace empty cells with NA
          as.data.frame()
      )
    }
  })


  var <- reactive({
    mycols <- colnames(readData())
  })

  output$whatx <- renderUI({
    selectInput("colmnamesx",
      label = h4("Select x axis data"),
      choices = var(), selected = colnames(readData()[1])
    )
  })

  output$whaty <- renderUI({
    selectInput("colmnamesy",
      label = h4("Select y axis data"),
      choices = var(), selected = colnames(readData()[2])
    )
  })

  output$contents <- DT::renderDT({
    readData()
  })

  procDat <- reactive({
    if (is.null(input$colmnamesy)) {
      return()
    }
    if (is.null(input$colmnamesx)) {
      return()
    }
    req(input$colmnamesx, input$colmnamesy)
    req(readData())
    readData <- readData()

    S <- readData[[input$colmnamesx]]
    V <- readData[[input$colmnamesy]]

    allDat <- signif(data.frame(
      "S" = S, "V" = V,
      "Sh" = S, "Vh" = S / V,
      "Se" = V / S, "Ve" = V,
      "Sl" = 1 / S, "Vl" = 1 / V,
      "Ss" = V, "Vs" = V / S
    ), digits = 4)
    allDat
  })


  output$myplot <- renderPlot({
    req(input$colmnamesx, input$colmnamesy)
    req(procDat())
    req(tabData())
    if (is.null(input$colmnamesy)) {
      return()
    }
    if (is.null(input$colmnamesx)) {
      return()
    }
    procDat <- procDat()
    tabData <- tabData()
    switch(input$raw,
      "Lineweaver-Burk" = linPlot(procDat, Sl, Vl),
      "Hanes" = linPlot(procDat, Sh, Vh),
      "Eadie-Hofstee" = linPlot(procDat, Se, Ve),
      "Scatchard" = linPlot(procDat, Ss, Vs),
      "Non-linear" = mmPlot(procDat, S, V, as.numeric(tabData[1, 4]), as.numeric(tabData[1, 3]))
    )
  })

  S.dat <- reactive({
    readData <- readData()

    S <- readData()[[input$colmnamesx]]
    S
  })

  V.dat <- reactive({
    V <- readData()[[input$colmnamesy]]
    V
  })

  tabData <- reactive({
    # if(is.null(input$colmnamesx)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
    req(input$colmnamesx, input$colmnamesy)
    req(readData())
    readData <- readData()

    S <- readData[[input$colmnamesx]]
    V <- readData[[input$colmnamesy]]

    # Lineweaver-Burke
    X.l <- 1 / S
    Y.l <- 1 / V
    LModl <- lm(Y.l ~ X.l)
    slope.l <- coef(LModl)[2]
    int.l <- coef(LModl)[1]
    r.l <- signif(cor(X.l, Y.l), digits = 4)
    Kmlm.l <- signif(slope.l / int.l, digits = 4)
    Vmaxlm.l <- signif(1 / int.l, digits = 4)

    # Hanes
    X.h <- S
    Y.h <- S / V
    HModl <- lm(Y.h ~ X.h)
    slope.h <- coef(HModl)[2]
    int.h <- coef(HModl)[1]
    r.h <- signif(cor(X.h, Y.h), digits = 4)
    Kmlm.h <- signif(int.h / slope.h, digits = 4)
    Vmaxlm.h <- signif(1 / slope.h, digits = 4)

    # Eadie-Hofstee
    X.e <- V / S
    Y.e <- V
    EModl <- lm(Y.e ~ X.e)
    slope.e <- coef(EModl)[2]
    int.e <- coef(EModl)[1]
    r.e <- signif(cor(X.e, Y.e), digits = 4)
    Kmlm.e <- signif(-slope.e, digits = 4)
    Vmaxlm.e <- signif(int.e, digits = 4)

    # Scatchard
    X.s <- V
    Y.s <- V / S
    EModl <- lm(Y.s ~ X.s)
    slope.s <- coef(EModl)[2]
    int.s <- coef(EModl)[1]
    r.s <- signif(cor(X.s, Y.s), digits = 4)
    Kmlm.s <- signif(-1 / slope.s, digits = 4)
    Vmaxlm.s <- signif(int.s * Kmlm.s, digits = 4)

    fitMM <- nls(V ~ SSmicmen(S, Vm, K))
    fitted <- predict(fitMM)
    Vmax <- signif(coef(fitMM)[1], digits = 4)
    Km <- signif(coef(fitMM)[2], digits = 4)
    crcNls <- signif(cor(V, fitted), digits = 4)


    tabData <- matrix(c(
      "Non-linear fit", "[S] vs V", Vmax, Km, crcNls,
      "Linear fit Hanes", "[S] vs [S]/V", Vmaxlm.h, Kmlm.h, r.h,
      "Linear fit Lineweaver-Burk", "1/[S] vs 1/V", Vmaxlm.l, Kmlm.l, r.l,
      "Linear fit Eadie-Hofstee", "V/[S] vs V", Vmaxlm.e, Kmlm.e, r.e,
      "Linear fit Scatchard *", "Bound vs Bound/Free", Vmaxlm.s, Kmlm.s, r.s
    ), byrow = TRUE, nrow = 5)
    colnames(tabData) <- c("Fit", "plot x~y", "Vmax", "Km", "Correlation")

    # write.table(tabData, "clipboard", sep="\t", col.names=F, row.names=F)

    tabData
  })

  output$resultsTable <- renderTable({
    req(tabData())
    tabData()
  })

  output$resultsTable2 <- DT::renderDT({
    req(procDat())
    procDat()
  })
}
# Run the application
shinyApp(ui = ui, server = server)
