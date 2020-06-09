library("shinyFiles")
library(shinyFiles)
library("shinyBS")
library("DT")
library(plyr)
library(dplyr)
#library(pryr)
library(data.table)   
library(stringr)
library(tidyr)
#library(xlsx)
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library(Biostrings)
library(plotly)
library(xtable)
library(plot3D)
library(gridExtra)
library(RColorBrewer)
library('stringdist')
library("parallel")
library(sjmisc)
library(magrittr)  
  shinyServer(function(input,output,session){
    
   
  shinyFileChoose(input, "files", root=c(root='.'), filetypes = c('', "txt", "csv"))
    
      original.data <- eventReactive(input$files, {
      
             file1 = parseFilePaths(roots = c(root = '.'), input$files)
      
             if(length(file1$datapath) == 0){ return() }
      
             data = read.csv(file = as.character(file1$datapath), sep = "\t")
             
             
             
        })
      
      output$table <- renderDataTable(original.data())
      


      filtered.Data <- eventReactive(input$FilterButton, { 
     
      a = as.numeric(input$select)
      
      b = as.character(input$text)
      
      c = colnames(original.data())
      
      d = as.character(c[[a]])
     
      TheData = as.data.frame(original.data())
    
      e = factor(TheData[,d]) 
      
      The.filtered = filter(TheData, e == b)
      
      } )
      
        
      output$filtered <- renderDataTable(filtered.Data())
        
          
        }
        
              )
      
      
      
      

    
  

  
