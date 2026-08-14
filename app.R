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
      #Select the columns of data from the file to use
      uiOutput("whatx"),
      uiOutput("whaty"),
      
      #Select non-linear plot or Scatchard transformation
      radioButtons(
        inputId = "raw", label = tags$h4("Select plot"),
        choices = c(
          "Non-linear",
          "Scatchard"
        ),
        selected = "Non-linear", inline = TRUE
      ),
      
      #This section calculates a new set of data from the numeric inputs for A+B=AB
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
      
      #New data is calculated only on pressing this button
      actionButton("goCalc", "Click to update the new data"),
      tags$br(),
      
      #Contact info and links to code
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
    
    #Main panel with plot and results summary plus extra tabs
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
        
        #Tab 2 is raw data
        tabPanel("Raw data", DT::DTOutput("contents")),
        
        #Tab 3 is transformed data
        tabPanel("Transformed data", DT::DTOutput("resultsTable2")),
        
        #Tab 4 is brief help
        tabPanel(
          "Help",
          tags$blockquote(h5(
            "►The app opens with a data file from from a binding assay",
            tags$br(),
            "►Alternatively check the box to generate a new set of data and adjust the binding parameters",
            tags$br(),
            "►Then explore the effects of changing Kd or [A] or [B] and errors on results and plots",
            tags$br(),
            "►Click the update data button to recalculate and replot results",
            tags$br(),
            "►It is also possible to look at the effect of using 'Added' rather than 'Free' reactants",
            tags$br(),
            "►'Added' is often used to approximate free, but is only valid when the conc of ligand [A] is << [B]",
            tags$br(),
            "►You can explore how selections affect the results and deviations in non-linear and linear plots",
            tags$br(),
            "►Load your own data for fitting as csv, txt or xlsx files (the app will detect the format)",
            tags$br(),
            "►The supplied data shows the expected data layout",
            tags$br(),
            "►There are tabs to access tables of raw and transformed data"
            
          )),
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  #Get and clean startup data or user data
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

  #Choose between startup or user data and generated data
  selectData <- reactive({
    if(input$gen) selectData <- newRes()
    else selectData <- readData()
  })
  
  #Get column names from selected data
  var <- reactive({
    mycols <- colnames(selectData())
  })

  #Identify the x column-usually free or added [B]
  output$whatx <- renderUI({
    selectInput("colmnamesx",
      label = h4("Select x axis data"),
      choices = var(), selected = colnames(selectData()[1])
    )
  })

  #Intentify the y data-usually the complex [AB]
  output$whaty <- renderUI({
    selectInput("colmnamesy",
      label = h4("Select y axis data"),
      choices = var(), selected = colnames(selectData()[2])
    )
  })

  #This is the raw data for the second tab
  output$contents <- DT::renderDT({
    selectData()
  })

  #Process the selected data to add the Scatchard transform
  procDat <- reactive({
    
    req(selectData())

    theData <- selectData()

    S <- theData[[input$colmnamesx]]
    V <- theData[[input$colmnamesy]]

    allDat <- signif(data.frame(
      "X" = S, "Y" = V,
      "Xs" = V, "Ys" = V / S
    ), digits = 4)
    allDat
  })

  #A new set of results is calculated from inputs and resGen function on goCalc activation
  newRes<- reactive({
    
    input$goCalc
    
    A<-isolate(input$concA)
    K<-isolate(input$Diss)
    G <- isolate(input$jit*input$concA)
    B1 <- isolate(input$B1)
    B2 <- isolate(input$B2)
    L <- isolate(input$npoints)
    #Make a log sequence of concentrations of B
    B <- exp(seq(log(B1), log(B2), length.out = L))
    #Call the external function resGen
    X <- resGen(A, B, K, G)
    #Set up the dataframe of results
    Res1 <- data.frame("Free"=round(B-X,8), "Bound"=round(X,8), "Added"=round(B,8))
    
  })
  
  #Plots are generated using external functions for non-linear or linear plots
  output$myplot <- renderPlot({
    
    req(tabData())
    req(procDat())
    
    #procDat for plotting and tabData for added lines
    plotDat <- procDat()
    tabData <- tabData()
    switch(input$raw,
      "Scatchard" = linPlot(plotDat, Xs, Ys, input$colmnamesx),
      "Non-linear" = mmPlot(plotDat, X, Y, as.numeric(tabData[1, 4]), as.numeric(tabData[1, 3]), input$colmnamesx)
    )
    
  })

  #This is the summary table below the plots
  tabData <- reactive({
    req(input$colmnamesx) #These definitely  seem to help prevent temporary error
    req(input$colmnamesy)
    req(selectData())
    
    theData <- selectData()

    S <- theData[[input$colmnamesx]]
    V <- theData[[input$colmnamesy]]

    # Linear Scatchard transformation
    X.s <- V
    Y.s <- V / S
    EModl <- lm(Y.s ~ X.s)
    slope.s <- coef(EModl)[2]
    int.s <- coef(EModl)[1]
    r.s <- signif(cor(X.s, Y.s), digits = 4)
    Kmlm.s <- signif(-1 / slope.s, digits = 4)
    Vmaxlm.s <- signif(int.s * Kmlm.s, digits = 4)

    #Non-linear fitting using SSmicmen
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

    tabData
    
  })

  #Output the results summary table
  output$resultsTable <- renderTable({
    req(tabData())
    tabData()
  })

  #Output the raw and transformed data in tab 3
  output$resultsTable2 <- DT::renderDT({
    req(procDat())
    procDat()
  })
}
# Run the application
shinyApp(ui = ui, server = server)
