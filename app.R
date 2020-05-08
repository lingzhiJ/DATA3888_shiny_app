#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(ggplot2)
library(GEOquery)  ## go to https://github.com/seandavi/GEOquery for installation details
library(R.utils)
library(reshape2)
library(limma) 
library(tidyverse)
library(dplyr)
library(shiny)

data = read.csv('Gene_expression_matrix.csv') 
clinical_outcome <-getGEO("GSE120649")
clinical_outcome<- clinical_outcome$GSE120649_series_matrix.txt.gz
rejection_status  <- clinical_outcome$description
rejection_status <- unlist(lapply(strsplit(as.character(rejection_status), ": " ) , `[[` , 1)) 
rejection_status = rejection_status[rejection_status != 'TCMR']

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    titlePanel("How number of gene influence PCA"),
    sidebarLayout(
        sidebarPanel(
          numericInput("size", "number of genes:", 10, min = 100, max = 57819),
          verbatimTextOutput("value")
        ),
        mainPanel(
            tabsetPanel(
            tabPanel("PCA",shiny::plotOutput("PCA"))
            ),
            
        )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$PCA <- renderPlot({
      
        gse = data[1:input$size,2:length(data)]
        gse_pca <- prcomp(t(gse))
        df_toplot <- data.frame(rejection_status, 
                                pc1 = gse_pca$x[,1], pc2 = gse_pca$x[,2])
        
        g <- ggplot(df_toplot, aes(x = pc1, y = pc2, color = rejection_status)) + 
          geom_point() + 
          theme_minimal() 
        g
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
