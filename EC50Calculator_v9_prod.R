#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(reshape2)
library(dplyr)
library(drc)



#exampleInputFiles <- tempfile()
#file.copy("ExampleInputFiles.zip", exampleInputFiles)

#Define row/col names for validation of appropriate format
reqColNames <- c("1","2","3","4","5","6","7","8","9","10","11","12")
reqRowNames <- c("A","B","C","D","E","F","G","H")

#checks if row/col names match expected input
checkplatefile <- function(df){
        if(is.null(df))
            return(NULL)
        else if((ncol(df) < 13) == TRUE || (nrow(df) < 9) == TRUE){
            return(FALSE)
        }
        else if(identical(as.character(df[1,2:13]),reqColNames) == FALSE){
            return(FALSE)
        }
        else if(identical(df[2:9,1],reqRowNames) == FALSE){
            return(FALSE)
        }
        else{
            return(TRUE)
        }
}

#checks if plate name matches expected input
checkplatename <- function(df){
    name <- df[1,1]
    return(name)
}


#compare cell names for control (DMSO) list to cell names for samples with compound added to determine if correct control was chosen
checkcontrol <- function(compiled, control){
  Controlresults <- compiled[which(compiled$compound == control),]
  Control_concentrations <- length(unique(Controlresults$concentration))

  if((Control_concentrations > 1) == TRUE){
    return(FALSE)
  }
  else {
    return(TRUE)
  }

}




findLevels <- function(df){
  if(is.null(df))
    return(NULL)
  else {
    df_formatted <- formatplatefile(df)
    df_melt <- melt(df_formatted,id = "PlateRow")
    
    levels <- unique(df_melt$value)
    levels <- na.omit(levels)
    
    return(levels)
  }
}

#for finding and setting options for individual plots - pulls levels from input files
findLevels_noControl <- function(df, controlchoice){
  if(is.null(df))
    return(NULL)
  if(is.null(controlchoice))
    return(NULL)
  else {
    df_formatted <- formatplatefile(df)
    df_melt <- melt(df_formatted,id = "PlateRow")
    
    #removes DMSO from compound list
    if(controlchoice %in% df_melt$value){
      df_melt <- df_melt[-which(df_melt$value == controlchoice),]
    }

    
    levels <- unique(df_melt$value)
    levels <- na.omit(levels)
    
    return(levels)
  }
}



#format input csv to take in only expected rows/columns 
formatplatefile <- function(df){
    df_form <- df[1:9,1:13]
    colnames(df_form) <- df_form[1,]
    df_form <- df_form[-1,]
    rownames(df_form) <- df_form[,1]
    colnames(df_form)[1] <- "PlateRow"
    return(df_form)
}

#takes in unformatted files, runs through format plate file function then compiles as a melted table
generatecompileddata <- function(dose, cell, values, compound){
    #format table files
    dosedf <- formatplatefile(dose)
    celldf <- formatplatefile(cell)
    compounddf <- formatplatefile(compound)
    valuesdf <- formatplatefile(values)

    #melt and compile table files
    dosedf_melt <- melt(dosedf,id = "PlateRow")
    celldf_melt <- melt(celldf,id = "PlateRow")
    compounddf_melt <- melt(compounddf,id = "PlateRow")
    valuesdf_melt <- melt(valuesdf,id = "PlateRow")
    compiled <- dosedf_melt
    
    #for validations of correct order between plates
    if(identical(compiled[,1:2],celldf_melt[,1:2]) ==TRUE){
        compiled$cell_line <- celldf_melt$value
    }
    else {
        return(FALSE)
    }
    if(identical(compiled[,1:2],compounddf_melt[,1:2]) ==TRUE){
    compiled$compound <- compounddf_melt$value
    }
    else {
        return(FALSE)
    }
    if(identical(compiled[,1:2],valuesdf_melt[,1:2]) ==TRUE){
    compiled$cell_count <- valuesdf_melt$value
    }
    else {
        return(FALSE)
    }
    
    #more table formatting
    colnames(compiled)[2:3] <- c("PlateCol","concentration")
    compiled <- compiled[complete.cases(compiled[,c(4:6)]),]
    #compiled[is.na(compiled)] <- 0
    compiled <- compiled[,-c(1,2)]
    
    return(compiled)
}



generatenormalizeddata <- function(compiled, control){

  #extract DMSO control values and normalize compiled data to DMSO results by cell line
  Controlresults <- compiled[which(compiled$compound == control),]
  compiled <- compiled[-which(compiled$compound == control),]
  Controlresults <- aggregate(Controlresults$cell_count, list(Controlresults$cell_line), FUN=mean) 
  compiled$cell_count__ctrl <- Controlresults$x[match(compiled$cell_line, Controlresults$Group.1)]
  compiled$cell_count <- compiled$cell_count / compiled$cell_count__ctrl
  
  #make cell line and compound columns factors to determine levels for generating EC50
  compiled$cell_line <- factor(compiled$cell_line)
  compiled$compound <- factor(compiled$compound)
  
  
  return(compiled)
  #}
}



generateEC50results <- function(compiled){
  
  EC50_results <- as.data.frame(matrix(ncol = 3, nrow = (length(levels(compiled$cell_line))*length(levels(compiled$compound)))))
  colnames(EC50_results) <- c("Cell_line", "Compound","EC50")
  
  #for-loop to calculate EC50 based on cell line and compound
  x <- 1
  for (i in levels(compiled$cell_line)){
    for (j in levels(compiled$compound)){
      drm_df <- compiled[which(compiled$cell_line==as.character(i) & compiled$compound==as.character(j)),c(4,1)]
      model <- drm(drm_df, fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
      EC50 <- round(summary(model)$coefficients[4,1], digits = 2)
      EC50_results[x,1] <- i
      EC50_results[x,2] <- j
      EC50_results[x,3] <- EC50
      x <- x+1
    }
  }
  
  return(EC50_results)
  
}



plotGraph <- function(compiled, cellchoice, compoundchoice){
  
  #generate dataframe for specific cell and compound used
  drm_df <- compiled[which(compiled$cell_line==cellchoice & compiled$compound==compoundchoice),c(4,1)]
  model <- drm(drm_df, 
               fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
  EC50 <- round(summary(model)$coefficients[4,1], digits = 2)

  if((min(model$data$cell_count) < .5) == FALSE){
    warning <- "WARNING: EC50 is extrapolated"
  }
  else {
    warning <- ""
  }
  
  plot <- plot(model, xlab = "Concentration in uM", ylab = "Relative cell number", col="dark red", type = "all", ylim = c(0,1.5))
  plot_title <- title(main = paste("Cell: ", cellchoice, "    Compound: ", compoundchoice))
  EC50warning <- text(x=min(model$data$concentration), y=0.25, labels = c(paste0(warning)), adj = 0, col = "Firebrick", font = 2)
  EC50text <- text(x=min(model$data$concentration), y=0.15, labels = c(paste0("EC50 = ", EC50, " microM")), adj = 0)
  EC50plot <- list(plot, plot_title, EC50warning, EC50text)
  return(EC50plot)
}


# Define UI for application that draws a histogram
ui <- fillPage(
    
    # Application title
    titlePanel("EC50 Calculator"),
    
    # Sidebar to import 4 csv files -  must be in 96-well plate format and contain no data past 
    sidebarLayout(
        sidebarPanel(

            fileInput("cell_upload", "Upload Cell Line:",
                      accept = c(".csv")),
            
            fileInput("compound_upload", "Upload Compound:",
                      accept = c(".csv")),
            
            fileInput("dose_upload", "Upload Dose (in uM):",
                      accept = c(".csv")),
            
            fileInput("values_upload", "Upload Values:",
                      accept = c(".csv")),
            
            downloadButton("downloadData", "Download Compiled Results"),
            br(),
            br(),
            downloadButton("downloadPlot", "Download EC50 Plot"),
            br(),
            br(),
            br(),
            br(),
            br(),
            downloadButton("exampleData", "Download Example Input Files"),
            br(),
            br(),
            tags$button(
              id = 'close',
              type = "button",
              class = "btn action-button",
              onclick = "setTimeout(function(){window.close();},500);", "Close window"
            )
        ),


        # Show a plot of design rep. "table" can be used for troubleshooting
        mainPanel(
            #uiOutput("plotSelect"),
            sidebarPanel(
              selectInput(inputId = "controlchoice", 
                          label = strong("Select the control:"),
                          choices = NULL),
              selectInput(inputId = "cellchoice", 
                          label = strong("Select a cell line:"),
                          choices = NULL),
              selectInput(inputId = "compoundchoice", 
                          label = strong("Select a compound:"),
                          choices = NULL),
            ),
            br(),
            br(),
            tableOutput(outputId = "compiled"),
            br(),
            br(),
            plotOutput(outputId = "plot"),
 
        )
    )
  )





# Define server logic required to draw a histogram
server <- function(input, output) {
  
    #define and take in reactive variables
    cell <- reactive({
        inFile <- input$cell_upload
        if (is.null(inFile))
            return(NULL)
        cell <- read.csv(inFile$datapath, header = FALSE)
        return(cell)
    })
    
    observeEvent(compiledTable(), {
      freezeReactiveValue(input, "cellchoice")
      cellchoices <- findLevels(cell())
      updateSelectInput(inputId = "cellchoice", choices = cellchoices) 
    })
    
    compound <- reactive({
        inFile <- input$compound_upload
        if (is.null(inFile))
            return(NULL)
        compound <- read.csv(inFile$datapath, header = FALSE)
        return(compound)
    })
   
    observeEvent(compiledTable(), {
      freezeReactiveValue(input, "controlchoice")
      controlchoices <- findLevels(compound())
      updateSelectInput(inputId = "controlchoice", choices = controlchoices) 
    })
    

     
    observeEvent(input$controlchoice, {
      freezeReactiveValue(input, "compoundchoice")
      compoundchoices <- findLevels_noControl(compound(),input$controlchoice)
      updateSelectInput(inputId = "compoundchoice", choices = compoundchoices) 
    })
    
    dose <- reactive({
        inFile <- input$dose_upload
        if (is.null(inFile))
            return(NULL)
        dose <- read.csv(inFile$datapath, header = FALSE)
        return(dose)
    })
    
    values <- reactive({
      inFile <- input$values_upload
      if (is.null(inFile))
        return(NULL)
      values <- read.csv(inFile$datapath, header = FALSE)
      return(values)
    })
    
    compiledTable <- function(){
      req(input$dose_upload)
      req(input$cell_upload)
      req(input$values_upload)
      req(input$compound_upload)
      if(checkplatefile(dose()) == FALSE){
        return(NULL)
      }
      if(checkplatefile(cell()) == FALSE){
        return(NULL)
      }
      if(checkplatefile(values()) == FALSE){
        return(NULL)
      }
      if(checkplatefile(compound()) == FALSE){
        return(NULL)
      }
      generatecompileddata(dose(), cell(), values(), compound())
    }
    
    
    normalizedTable <- function(){
      req(compiledTable())
      req(input$controlchoice)
      generatenormalizeddata(compiledTable(), input$controlchoice)
    }
    
    
    output$compiled <- renderTable({
           validate(
              need(checkplatefile(values()) == TRUE, "Unexpected Values plate format found. Please upload a file in 96-well format."),                need(checkplatename(values()) == "Values", "Values plate is incorrect. Please fix label or upload correct file."),
              need(checkplatefile(cell()) == TRUE, "Unexpected Cell plate format found. Please upload a file in 96-well format."),
              need(checkplatename(cell()) == "Cell line", "Cell plate is incorrect. Please fix label or upload correct file."),
              need(checkplatefile(compound()) == TRUE, "Unexpected Compound plate format found. Please upload a file in 96-well format."),
              need(checkplatename(compound()) == "Compound", "Compound plate is incorrect. Please fix label or upload correct file."),
              need(checkplatefile(dose()) == TRUE, "Unexpected Dose plate format found. Please upload a file in 96-well format."),
              need(checkplatename(dose()) == "Dose (in uM)", "Dose plate is incorrect. Please fix label or upload correct file."),
              need(checkcontrol(compiledTable(),input$controlchoice) == TRUE, "Control is not correct. Control aligns with >1 concentration.")
            )
        
        #for troubleshooting
             #generatecompileddata(dose(), cell(), values(), compound())
      req(compiledTable())  
      req(normalizedTable())
      
      generateEC50results(normalizedTable())
    })
 
    
    output$plot <- renderPlot({
      req(normalizedTable())
      req(compiledTable())
      req(input$compoundchoice)
      req(input$cellchoice)
      validate(
        need(checkplatefile(values()) == TRUE, message = FALSE),
        need(checkplatename(values()) == "Values", message = FALSE),
        need(checkplatefile(cell()) == TRUE, message = FALSE),
        need(checkplatename(cell()) == "Cell line", message = FALSE),
        need(checkplatefile(compound()) == TRUE, message = FALSE),
        need(checkplatename(compound()) == "Compound", message = FALSE),
        need(checkplatefile(dose()) == TRUE, message = FALSE),
        need(checkplatename(dose()) == "Dose (in uM)", message = FALSE),
        need(checkcontrol(compiledTable(),input$controlchoice) == TRUE, message = FALSE)
      )

      #plotGraph(compiledTable(), input$cellchoice_plot, input$compoundchoice_plot)
      plotGraph(normalizedTable(), input$cellchoice, input$compoundchoice)
    })

    #output$cell <- renderText(print(input$cellchoice))
    #output$compound <- renderText(print(input$compoundchoice))
    
    
    
    
    
    
    output$downloadData <- downloadHandler(
      filename = function(){
        req(normalizedTable())
        paste("EC50_compiled_results", ".csv" , sep = "")
      },
      content = function(file){
        write.csv(generateEC50results(normalizedTable()), file, row.names = FALSE)
      }
    )
    
    output$downloadPlot <- downloadHandler(
      filename = function(){
        req(input$cellchoice)
        req(input$compoundchoice)
        paste("EC50_results_",input$cellchoice, "_", input$compoundchoice, ".png" , sep = "")
      },
      content = function(file){
        req(normalizedTable())
        req(input$cellchoice)
        req(input$compoundchoice)
        png(file)
        print(plotGraph(normalizedTable(), input$cellchoice, input$compoundchoice))
        dev.off()
      }
    )
    
    
    output$exampleData <- downloadHandler(
        filename <- function() {
          paste("ExampleInputFiles", "zip", sep=".")
        },
        
        content <- function(file) {
          file.copy("ExampleInputFiles.zip", file)
        },
        contentType = "application/zip"
      )
    

    observe({
      if (input$close > 0) stopApp() 
    })
    


}

# Run the application 
shinyApp(ui = ui, server = server)
