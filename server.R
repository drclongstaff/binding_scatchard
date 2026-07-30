library(shiny)


function(input, output){


  readData <- reactive({

    inputFile <- input$data
    
   if  (input$gen) inputFile<-newRes()
   else{ if (is.null(inputFile)) 
      return(read.csv(file.path("./Data/Scatchard perfect.csv"))) 
      else(read.table(inputFile$datapath, sep = input$sep, header = input$header))
   }

  })

  var<- reactive({
    mycols<-colnames(readData())
 })
  
  

  output$whatx<-renderUI({
   
    selectInput("colmnamesx",
                label= h5("Select independent variable (free or added)"),
                choices = var(), selected = colnames(readData()[1]))
  })

  output$whaty<-renderUI({
   
    selectInput("colmnamesy",
                label= h5("Select dependent variable (bound)"),
                choices = var(), selected = colnames(readData()[2]))
                                 
  })

  newRes<- reactive({
    
    input$goCalc
    
    A<-isolate(input$concA)
    K<-isolate(input$Diss)
    Be<-seq(input$concB[1], input$concB[2], length.out = 20)
    B<-10^(Be)
    xv<-1:length(B)
    for (i in 1:length(B)){
      b<-A+B[i]+K
      c<-A*B[i]
      x<-(b-(b^2-4*c)^0.5)/2
      xv[i]<-x
      #xvn<-xv+rnorm(length(xv), 0, input$jit)
      
      X<-round(xv, 5)
    }
    Res<-cbind((B-X), X, B)
    colnames(Res)<-c("Free", "Bound", "Added")
    #write.table(Res, "clipboard", sep="\t", col.names=TRUE, row.names=F)
    #write.table(Res, "Binding perfect.txt", sep="\t", col.names=TRUE, row.names=F)
    
    
    
    Res1<-data.frame(Res)
  })
  
  adjData<-reactive({
    if(is.null(input$colmnamesx)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
    if(input$gen) selectedData<-newRes()
    else selectedData<-readData()
    selectedData[,input$colmnamesy]<-selectedData[,input$colmnamesy]+rnorm(length(selectedData[,1]), 0, input$jit)
    data.frame(selectedData)
    
  })
  
  
  output$myplot<-renderPlot({
    if(is.null(input$colmnamesx)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
    
    Fr<-adjData()[,input$colmnamesx]  
    Bo<-adjData()[,input$colmnamesy]
    X<-Bo
    Y<-Bo/Fr
    ScModl<-lm(Y~X)
    slope<- coef(ScModl)[2]
    int<- coef(ScModl)[1]
    r<- signif(cor(X,Y), digits=4)
    
    Kdlm<- signif(-1/slope,digits=4)
    Bmaxlm<-signif(Kdlm*int, digits=4)
    
    fitMM<-nls(Bo~Bmax*Fr/(Kd+Fr), start=list(Bmax=Bmaxlm, Kd=Kdlm))
    pred_x<-seq(0, max(Fr), length.out = 50)
    fitted<- predict(fitMM, list(Fr=pred_x))
    Bmax<- signif(coef(fitMM)[1], digits=4)
    Kd<- signif(coef(fitMM)[2], digits=4)
    crcNls<-signif(cor(Bo, predict(fitMM)), digits = 4)
    
    limy<-ifelse(int>max(Y), int, max(Y))
    limx<-ifelse((-int/slope)>max(X), -int/slope, max(X))
    
      switch(input$raw,
             "Raw data"=plot(Fr,Bo, main="Raw data", pch=21, cex=2, xlab=input$colmnamesx, ylab=input$colmnamesy, xlim= c(0, max(Fr)), ylim = c(0, Bmax*1.1), lines(pred_x, fitted,lwd=2)),
             #abline("h"=Bmax, col="red", lwd=2),
             #"Scatchard"=plot(X,Y, main="Scatchard plot", pch=17, xlab=input$colmnamesy, ylab=paste("Free/", input$colmnamesy), xlim= c(0, Bmaxlm), ylim = c(0, int), abline(lm(Y~X), col="red", lwd=2)),
             "Scatchard"=plot(X,Y, main="Scatchard plot", pch=17, cex=1.5, xlab=input$colmnamesy, ylab=paste(input$colmnamesy,"/",input$colmnamesx), xlim= c(0, limx), ylim = c(0, limy), abline(lm(Y~X), col="red", lwd=2)),
             "Residual"=plot(Fr, ScModl$resid, xlab=input$colmnamesx, ylab= "Residual", pch=4, col="blue", lwd=2)) 
      
  })
  
  output$resultsTable<-renderTable({
    if(is.null(input$colmnamesx)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
    
    Fr<-adjData()[,input$colmnamesx]  
    Bo<-adjData()[,input$colmnamesy]
    
    X<-Bo
    Y<-Bo/Fr
    ScModl<-lm(Y~X)
    slope<- coef(ScModl)[2]
    int<- coef(ScModl)[1]
    r<- signif(cor(X,Y), digits=4)
    
    Kdlm<- signif(-1/slope,digits=4)
    Bmaxlm<-signif(Kdlm*int, digits=4)
    
    fitMM<-nls(Bo~Bmax*Fr/(Kd+Fr), start=list(Bmax=Bmaxlm, Kd=Kdlm))
    fitted<- predict(fitMM)
    Bmax<- signif(coef(fitMM)[1], digits=4)
    Kd<- signif(coef(fitMM)[2], digits=4)
    crcNls<-signif(cor(Bo, predict(fitMM)), digits = 4)
    
    
    tabData<-matrix(c("Non-linear fit (binding isotherm)", Bmax, Kd, crcNls,
                      "Linear fit (Scatchard plot)", Bmaxlm, Kdlm, r
                      ), byrow=TRUE, nrow=2)
    colnames(tabData)<-c("Fit", "Max binding", "Kd", "Correlation")
    
    #write.table(tabData, "clipboard", sep="\t", col.names=F, row.names=F) 
    
    tabData
    
  })

  output$contents<-renderDataTable({
    myData<-adjData()
    myBound<-adjData()[,2]
    myFree<-adjData()[,1]
    myBound_Free<-myBound/myFree
    theseResults<-data.frame(myBound, myBound_Free)
    write.table(theseResults, "clipboard", sep="\t", col.names=F, row.names=F)
    adjData()
    

  })

  
}