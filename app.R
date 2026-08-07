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
        selected = "Non-linear"
      ),
      
      checkboxInput(inputId= "gen", label = h4("Make a new set of data"), value = FALSE),
      
      fluidRow(
        column(6,numericInput("concA",
                   label = h5("Choose a concentration of ligand A"), value = 1)),
        column(6,numericInput("Diss",
                              label = h5("Choose a Kd value"), value = 0.5))
      ),
      #sliderInput("concB", 
                 # label =h5( "Range of ligand B (log scale)"),
                 # min = -1.1, max = 2.2, value = c(-1, 1)),
      
      fluidRow(
      column(6,numericInput("B1", label = h5("starting [B]"), value=0.2)),
      column(6,numericInput("B2", label = h5("ending [B]"), value=10))),
      fluidRow(
      column(6,numericInput("npoints", label = h5("number of points"), value = 20)),
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
    #read.csv("./Data/Scatchard perfect.csv") |> as.data.frame()
      newRes()
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
    
    readData <- newRes() #adjData() #newRes() #readData()

    S <- readData[[input$colmnamesx]]
    V <- readData[[input$colmnamesy]]

    allDat <- signif(data.frame(
      "S" = S, "V" = V,
      "Ss" = V, "Vs" = V / S
    ), digits = 4)
    allDat
  })

  newRes<- reactive({
    
    input$goCalc
    
    A<-isolate(input$concA)
    K<-isolate(input$Diss)
    G <- isolate(input$jit)
    #Be<-seq(input$concB[1], input$concB[2], length.out = 20)
    #B<-10^(Be)
    B <- seq(input$B1, input$B2, length.out=input$npoints)
    xv<-1:length(B)
    for (i in 1:length(B)){
      b<-A+B[i]+K
      c<-A*B[i]
      x<-(b-(b^2-4*c)^0.5)/2
      xv[i]<-x
      xvn<-xv+rnorm(length(xv), 0, G)
      
      #X<-round(xv, 5)
      if(input$gen) X<-round(xvn, 5)
      else X <- round(xv, 5)
    }
    Res<-cbind((B-X), X, B)
    colnames(Res)<-c("Free", "Bound", "Added")
    #write.table(Res, "clipboard", sep="\t", col.names=TRUE, row.names=F)
    #write.table(Res, "Binding perfect.txt", sep="\t", col.names=TRUE, row.names=F)
    
    
    
    Res1<-data.frame(Res)
  })
  
  adjData<-reactive({
    if(is.null(input$colmnamesx)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
    #if(input$gen) selectedData<-newRes()
    if(input$gen) selectedData<-procDat()[,input$colmnamesy]+rnorm(length(selectedData[,1]), 0, input$jit)
    else selectedData<-readData()
    #selectedData[,input$colmnamesy]<-selectedData[,input$colmnamesy]+rnorm(length(selectedData[,1]), 0, input$jit)
    data.frame(selectedData)
    
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
    plotDat <- procDat()
    tabData <- tabData()
    switch(input$raw,
      "Scatchard" = linPlot(plotDat, Ss, Vs),
      "Non-linear" = mmPlot(plotDat, S, V, as.numeric(tabData[1, 4]), as.numeric(tabData[1, 3]))
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
    readData <- newRes()#readData()

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
      "Non-linear fit", "Free vs Bound", Vmax, Km, crcNls,
      "Linear fit Scatchard ", "Bound vs Bound/Free", Vmaxlm.s, Kmlm.s, r.s
    ), byrow = TRUE, nrow = 2)
    colnames(tabData) <- c("Fit", "plot x~y", "Bmax", "Kd", "Correlation")

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
