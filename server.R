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
library(shinyalert)
library(tidyverse)
 
 shinyServer(function(input,output,session) {
    
   
      shinyDirChoose(input,'directory', roots =  c(root = '.'), filetypes = c('', 'txt'))
   
   
      directory <- reactive(input$directory)
    
      path1 <- reactive( { parseDirPath(roots =  c(root = '.'), input$directory)  })
   
      output$uiInputFiles <- renderUI ( {
      
         List.of.files = list.files(path1())
         wellPanel(
         tags$div(class = "multicol", checkboxGroupInput(inputId = "inputFiles", label = "Select Files", inline=FALSE, choices = List.of.files, selected=c("1_Summary.txt","2_IMGT-gapped-nt-sequences.txt","4_IMGT-gapped-AA-sequences.txt","6_Junction.txt")))
                 )
                                     })
        
      output$uiLoadData <- renderUI({
         
            actionButton("LoadData", "Load Data", style="color: black; background-color: #bf6f62; border-color: #fff")
                            
                                    })  
     
      ThePath = file.path("C:","Users","tgala","Desktop","Biosciences","TOOLS - BIOINFORMATICS","clonal R evolution","Data","B1")
      
      The_File = list()
    
      original.data = eventReactive(input$LoadData, {
       
            List_of_Files = as.list(input$inputFiles)
       
            for (i in List_of_Files) { 
          
               File_Name_Without_Quotes = noquote(i)
          
               File_Path = paste(ThePath, File_Name_Without_Quotes, sep = "/")
          
               The_File[[i]] = read.csv(file = File_Path, header = TRUE, sep = "\t")
            
                                  }
        
            MergedFiles = merge(The_File[[1]], The_File[[2]], by="Sequence.ID")
        
                                                        } )
     
        
  
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
      
      
      
      
      dominant_Clone <- observeEvent(input$calculateDominantClone, {
        
        print("e")
        
        DataFrame_AfterFiltering = as.data.frame(filtered.Data())
        
        print(DataFrame_AfterFiltering)
        
        Selected_Columns = select(DataFrame_AfterFiltering, V.GENE.and.allele)
        
        print(Selected_Columns)
      
                                                       } )  
      
      
      output$clonoTable <- renderDataTable(dominant_Clone())
        
        
        
        
      
    
    } 
)
  
        
      
      
      

    
  

  
