library(shiny)


fluidPage(
  
  titlePanel (h2("Binding curves and Scatchard plots",  align = "right")),
 
 
                
  sidebarLayout(
    
    sidebarPanel(

      
      fileInput("data", label = "Select dataset:"),
     
      radioButtons(inputId = "sep", label = "File type", choices = c(csv = ",", txt ="\t"), selected = ","),
      
      checkboxInput(inputId= "header", label = "Columns have header text", value = TRUE),
      
      
      
      
      sliderInput("jit", 
                  label =h4( "Add noise"),
                  min = 0, max = 0.1, value = 0.01),
                 
      uiOutput("whatx"), 
      
      uiOutput("whaty"),
    
     
      
      radioButtons(inputId="raw", label = h5("Select plot"), choices = c("Raw data", "Scatchard", "Residual"
                                                                              ), selected = "Raw data"), 
      checkboxInput(inputId= "gen", label = h4("Make a new set of data"), value = FALSE),
      
      numericInput("concA",
                   label = h5("Choose a concentration of ligand A"), value = 1),
      
      sliderInput("concB", 
                  label =h5( "Range of ligand B (log scale)"),
                  min = -1.1, max = 2.2, value = c(-1, 1)),
      
      numericInput("Diss",
                   label = h5("Choose a Kd value"), value = 0.5),
      
      actionButton("goCalc", "Click to update the new data"),
      
      helpText("Please cite this page if you find it useful, Longstaff C, 2016, Shiny App for analysing binding data, version 0.7,
               URL address, last accessed", Sys.Date())
   
    
    
  ),
  mainPanel( 
    tabsetPanel(type="tab",
                tabPanel("Plot", 
                         
                         plotOutput(outputId = "myplot"),
                         
                         h5(textOutput("text3")),
                         
                         h4("Results Table"), tableOutput("resultsTable"), align = "center"),
                
                
                
                tabPanel("Raw data", dataTableOutput("contents")),
                
                tabPanel("Help",
                         
                  tags$blockquote(h5("►Load your own data file in csv or txt fomat (tab separator)- SET THE NOISE SLIDER TO ZERO",
                                  tags$br(),
                                  "►Avoid unusual characters such as % and ' in names",
                                  tags$br(),
                                  "►A set of data is supplied and noise can be added or removed",
                                  tags$br(),
                                  "►The supplied data has columns of [Free], [Bound] and [Added]-sometimes used to approximate [Free]",
                                  tags$br(), 
                                  "►Plots available are raw data and non-linear fit, Scatchard plot and fit and its residual",
                                  tags$br(), 
                                  "►You can generate a new set of data for A+B=AB",
                                  tags$br(),
                                  "►Specify [A], the range of [B] and Kd, and also the amount of noise",
                                  tags$br(),
                                  "►Results are presented in the table for non-linear fit and from the linearised Scatchard plot"
                                 
                                 )
  
  ),
  
                
                  
                  
                  
                  tags$img(src="screenCapBind.png", width=600, height=700)
                  
                )
                
               
                
    )
  )
  
)   
)
