# An app to analyse Binding curves and Scatchard plots

library(shiny)
library(tidyplots)

ThisApp <- "Binding isotherms and linear Scatchard plots"
ThisVersion <- 0.8

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
  #helpText("Binding"),
  titlePanel(tags$h2(ThisApp)),# align = "center")),
  
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
          "Scatchard"
        ),
        selected = "Non-linear", inline = TRUE
      ),
      
      checkboxInput(inputId= "gen", label = h4("Make a new set of data"), value = FALSE),
      
      fluidRow(
        column(6,numericInput("concA",
                   label = h5("Conc of ligand A"), value = 0.1, min = 0, step = 0.05)),
        column(6,numericInput("Diss",
                              label = h5("Kd value"), value = 0.5, min=0, step = 0.01))
      ),
      
      fluidRow(
      column(6,numericInput("B1", label = h5("starting [B]"), value=0.1, min=0, step = 0.1)),
      column(6,numericInput("B2", label = h5("ending [B]"), value=10, min = 0, step = 0.5))),
      fluidRow(
      column(6,numericInput("npoints", label = h5("number of points"), value = 10, min=3, step = 1)),
      column(6, numericInput("jit", h5("add noise"), value = 0, min=0, max=1, step = 0.01))
      
      ),
      
      #sliderInput("jit", 
                 # label =h4( "Add noise"),
                 # min = 0, max = 0.1, value = 0),
      
      actionButton("goCalc", "Click to update the new data"),
      tags$br(),
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
          tableOutput("resultsTable")
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
    read.csv("./Data/Scatchard perfect.csv") |> as.data.frame()
      #newRes()
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

  selectData <- reactive({
    if(input$gen) selectData <- newRes()
    else selectData <- readData()
  })
  
  var <- reactive({
    mycols <- colnames(selectData()) #(readData())
  })

  output$whatx <- renderUI({
    selectInput("colmnamesx",
      label = h4("Select x axis data"),
      #choices = var(), selected = colnames(readData()[1])
      choices = var(), selected = colnames(selectData()[1])
    )
  })

  output$whaty <- renderUI({
    selectInput("colmnamesy",
      label = h4("Select y axis data"),
      #choices = var(), selected = colnames(readData()[2])
      choices = var(), selected = colnames(selectData()[2])
    )
  })

  output$contents <- DT::renderDT({
    #readData()
    #procDat()
    selectData()
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
    #if(input$gen) readData <- newRes()
    #else readData <- readData()
    
    #readData <-  readData()
    readData <- selectData()

    S <- readData[[input$colmnamesx]]
    V <- readData[[input$colmnamesy]]

    allDat <- signif(data.frame(
      "X" = S, "Y" = V,
      "Xs" = V, "Ys" = V / S
    ), digits = 4)
    allDat
  })

  newRes<- reactive({
    
    input$goCalc
    
    A<-isolate(input$concA)
    K<-isolate(input$Diss)
    G <- isolate(input$jit*input$concA)
    B1 <- isolate(input$B1)
    B2 <- isolate(input$B2)
    L <- isolate(input$npoints)
  
    B <- exp(seq(log(B1), log(B2), length.out = L))
    #if(input$gen) G <- G
    #else G <- 0
    #if(input$gen) X <- resGen(A, B, K, G)
    #else NULL
    X <- resGen(A, B, K, G)
    Res1 <- data.frame("Free"=round(B-X,8), "Bound"=round(X,8), "Added"=round(B,8))
    #Res1 <- data.frame("Free"=B-X, "Bound"=X, "Added"=B)
    
  })
  

  output$myplot <- renderPlot({
    req(input$colmnamesx, input$colmnamesy)
    req(procDat())
    #req(tabData())
    if (is.null(input$colmnamesy)) {
      return()
    }
    if (is.null(input$colmnamesx)) {
      return()
    }
    plotDat <- procDat()
    #plotDat <- readData()
    tabData <- tabData()
    switch(input$raw,
      #"Scatchard" = linPlot(plotDat, Ss, Vs, input$colmnamesx),
      #"Non-linear" = mmPlot(plotDat, S, V, as.numeric(tabData[1, 4]), as.numeric(tabData[1, 3]), input$colmnamesx)
      "Scatchard" = linPlot(plotDat, Xs, Ys, input$colmnamesx),
      "Non-linear" = mmPlot(plotDat, X, Y, as.numeric(tabData[1, 4]), as.numeric(tabData[1, 3]), input$colmnamesx)
    )
  })

  tabData <- reactive({
    if(is.null(input$colmnamesx)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
    req(input$colmnamesx, input$colmnamesy)
    req(readData())
    req(selectData())
    #readData <- procDat()
    #readData <- readData()
    readData <- selectData()

    S <- readData[[input$colmnamesx]]
    V <- readData[[input$colmnamesy]]

   

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
      "Non-linear fit", paste0(input$colmnamesx," vs Bound"), Vmax, Km, crcNls,
      "Linear fit Scatchard ", paste0("Bound vs Bound/", input$colmnamesx), Vmaxlm.s, Kmlm.s, r.s
    ), byrow = TRUE, nrow = 2)
    colnames(tabData) <- c("Fit", "plot x~y", "Bmax", "Kd", "Correlation")

    # write.table(tabData, "clipboard", sep="\t", col.names=F, row.names=F)

    tabData
  })

  output$resultsTable <- renderTable({
    #req(tabData())
    tabData()
  })

  output$resultsTable2 <- DT::renderDT({
    req(procDat())
    procDat()
    #selectData()
  })
}
# Run the application
shinyApp(ui = ui, server = server)
