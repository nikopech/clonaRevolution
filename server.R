library("shinyFiles")
library(shinyFiles)
library("DT")
library(plyr)
library(dplyr)
library(data.table)   
library(Biostrings)
library(xtable)
library(shinyalert)
library(tidyverse)
library(scales)
library(stringr)
library(e1071)
library(widyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(stringdist)
library(stringi)
library(reshape)
library(compare)
library(ggm)
library(visNetwork)
library(ggrepel)

shinyServer(function(input,output,session) {
    
   
      shinyDirChoose(input,'directory', roots =  c(root = '.'), filetypes = c('', 'txt'))
   
   
      directory <- reactive(input$directory)
    
      path1 <- reactive( { parseDirPath(roots =  c(root = '.'), input$directory)  })
   
      
      
      output$uiInputFiles <- renderUI ( {
      
         List.of.files = list.files(path1())
         
        
         
         
       wellPanel(
         tags$div(class = "multicol", checkboxGroupInput(inputId = "inputFiles", label = "Select Samples", inline=FALSE, choices = List.of.files, selected=c("3_Nt-sequences.txt",
                                                                                                                                                           "5_AA-sequences.txt")))
                 )
                                     })
        
      output$uiLoadData <- renderUI({
         
            actionButton("LoadData", "Load Data", style="color: black; background-color: #bf6f62; border-color: #fff")
                            
                                    })  
     
     
      
      
 ################################ LOADING DATA ######################################################     
    
      
      Dataframes_of_Samples <- function(z) { 
        
        
        The_Data_Frames = list()
        
        for (x in z) {
          
           Nt_File_path = paste0(x,"/3_Nt-sequences.txt") 
          
           AA_File_path = paste0(x,"/5_AA-sequences.txt")
          
           The_Nt = read.csv(file = Nt_File_path, header = TRUE, sep = "\t")
          
           The_AA = read.csv(file = AA_File_path, header = TRUE, sep = "\t")
          
           MergedFiles = merge(The_Nt, The_AA, by="Sequence.ID")
          
           Only_needed_columns = select(MergedFiles, Sequence.ID, Sequence.number.x, V.DOMAIN.Functionality.x, V.GENE.and.allele.x, D.GENE.and.allele.x, V.D.J.REGION.x, V.J.REGION.x, V.REGION.x, FR1.IMGT.x, CDR1.IMGT.x, FR2.IMGT.x, CDR2.IMGT.x, FR3.IMGT.x, CDR3.IMGT.x, JUNCTION.x, X3.V.REGION, X.N.D..J.REGION, X.N.D..REGION, D.REGION, N2.REGION, J.REGION.x, FR4.IMGT.x, J.GENE.and.allele.y, V.D.J.REGION.y, V.REGION.y, FR1.IMGT.y, CDR1.IMGT.y, FR2.IMGT.y, CDR2.IMGT.y, FR3.IMGT.y, CDR3.IMGT.y, JUNCTION.y, J.REGION.y, FR4.IMGT.y)
           
           id = as.character(substring(x, 17, 21))
           
           V_gene_column = as.character(Only_needed_columns$V.GENE.and.allele.x)
           
           D_gene_column = as.character(Only_needed_columns$D.GENE.and.allele.x)
           
           J_gene_column = as.character(Only_needed_columns$J.GENE.and.allele.y)
           
           V_Gene = unlist(strsplit(V_gene_column, "\\*"))[[1]]
           
           D_Gene = unlist(strsplit(D_gene_column, "\\*"))[[1]]
           
           J_Gene = unlist(strsplit(J_gene_column, "\\*"))[[1]]
           
           CDR3AA = as.character(Only_needed_columns$CDR3.IMGT.y)
           
           CDR3_length = as.numeric(nchar(CDR3AA))
           
           Only_needed_columns = cbind(Only_needed_columns, V_Gene, D_Gene, J_Gene, CDR3_length, id)
           
           The_Data_Frames[[x]]  = as.data.frame(Only_needed_columns)  
          
                      }
        
       return(The_Data_Frames)
        
    
                                        }
       
      
      
       original.data = eventReactive(input$LoadData, { 
         
               print(paste0(path1(),"/",input$inputFiles))
         
               List_of_samples = as.list(paste0(path1(),"/",input$inputFiles))  
               
               List_of_dfs_in_original.data = Dataframes_of_Samples(List_of_samples) 
       
               return(List_of_dfs_in_original.data)
               
       
                                               })
       
      
       observeEvent(input$number, {
         
         output$table = renderDataTable(original.data()[[input$number]])
       } )
             
       
         observeEvent(input$LoadData, { 
        
        shinyalert(title = "Loading Data. . .", type = "success") })
     
       
    ########################################################################################
         
         
    #################### FILTERING DATA ####################################################     
     
        filtered.Data <- eventReactive(input$FilterButton, { 
          
         
          dfs = original.data()
          
          List_of_filtered_dfs = list()
          
          
          for (y in 1:length(dfs)) {
          
              df_to_filter = as.data.frame(dfs[[y]]) 
        
              df_to_filter$V_Gene = substring(df_to_filter$V.GENE.and.allele.x, 1, 14)
           
              The_First_Filtering = df_to_filter
          
              if (input$checkfunctionalVgene == TRUE) {
            
                    
                    pattern = "P|ORF"
                    
                    The_First_Filtering = filter(The_First_Filtering, !grepl(pattern, df_to_filter$V.GENE.and.allele.x))
              }
              
             
          
          
              if (input$checkCDR3characters == TRUE) {
             
                        pattern = "\\*|X|#|\\."
             
                 The_First_Filtering = filter(The_First_Filtering, !grepl(pattern, The_First_Filtering$CDR3.IMGT.y))
             
                                                      } 
          
         
               if (input$checklandmarks == TRUE) {
            
                   
                      The_First_Filtering$firstLandmark = stri_sub(The_First_Filtering$JUNCTION.y, 1,1) 
              
                      The_First_Filtering$secondLandmark = stri_sub(The_First_Filtering$JUNCTION.y, -1)
            
                      The_First_Filtering = filter(The_First_Filtering, The_First_Filtering$firstLandmark == "C")
             
                      The_First_Filtering = filter(The_First_Filtering, The_First_Filtering$secondLandmark == "W")
                 
                                                  }
          
          
          
               if (input$checkproductive == TRUE) {
               
                       The_First_Filtering = filter(The_First_Filtering, V.DOMAIN.Functionality.x == "productive")
               
                                                    }
          
          
            
          
            
         #  
         # if (input$extrafiltering == TRUE) { 
         #  
         #  
         #   Selected_Column = as.numeric(input$select)
         #  
         #  Texted_Characteristic = as.character(input$text)
         #  
         #  The_Colnames = colnames(The_First_Filtering)
         #  
         #  Choose_from_the_colnames = as.character(The_Colnames[[Selected_Column]])
         #  
         #  TheData = as.data.frame(original.data())
         #  
         #  Column_to_choose = factor(TheData[,Choose_from_the_colnames]) 
         #  
         #  The.filtered = filter(TheData, Column_to_choose == Texted_Characteristic)
         #  
         #  The.filtered = as.data.frame(The.filtered)
         #   
         # }
          
          
       
           List_of_filtered_dfs[[y]] = as.data.frame(The_First_Filtering)
          
          } 
          
         return(List_of_filtered_dfs)
          
                                           
            } )  
        
           
           observeEvent(input$number2, {
           output$filtered <- renderDataTable(filtered.Data()[[input$number2]]) 
           } )
     
      
       observeEvent(input$FilterButton, { 
        
        shinyalert(title = "Filtering succeed", type = "success") })
      
      
       #################################################################################
       
       ##################### PIPELINE ################################################
       
       
      
      observeEvent(input$Executepipeline, {
        
        
      
        
        filtered_dfs = filtered.Data()
        
        List_of_Related = list()
        
         List_of_Clonotypes = list()
       
        List_of_Dominants = list()
        
        for (c in 1:length(filtered_dfs)) {
          
        df_to_execute = as.data.frame(filtered_dfs[[c]])
          
      
          if (input$select_clonotype == "V Nt sequence + CDR3 Amino Acids")

        { Selected_Columns = select(df_to_execute, V.REGION.x, CDR3.IMGT.y, V_Gene, D_Gene, J_Gene, CDR3_length, id)

        Number_of_reads = nrow(Selected_Columns)

        Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.REGION.x, CDR3.IMGT.y) %>% tally())

        Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)

        TheTable = cbind(Counts_of_unique_rows, Frequency)

        Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
        
        The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),]  }




         else if (input$select_clonotype == "V Gene and Allele + CDR3 Amino Acids")

        { Selected_Columns = select(df_to_execute, V.GENE.and.allele.x, CDR3.IMGT.y, V_Gene, D_Gene, J_Gene, CDR3_length, id)

      Number_of_reads = nrow(Selected_Columns)

      Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.GENE.and.allele.x, CDR3.IMGT.y) %>% tally())

      Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)

      TheTable = cbind(Counts_of_unique_rows, Frequency)
      
      Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
      
      The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }
      
         
     
         else if (input$select_clonotype == "V Nt sequence + CDR3 Nucleotide")

        { Selected_Columns = select(df_to_execute, V.REGION.x, CDR3.IMGT.x, V_Gene, D_Gene, J_Gene, CDR3_length, id)

          Number_of_reads = nrow(Selected_Columns)

          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.REGION.x, CDR3.IMGT.x) %>% tally())

          Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)

          TheTable = cbind(Counts_of_unique_rows, Frequency)

          Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
          
          The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }
      
          
          
          
          
      
          else if (input$select_clonotype == "V Gene and Allele + CDR3 Nucleotide")
      
         { Selected_Columns = select(df_to_execute, V.GENE.and.allele.x, CDR3.IMGT.x, V_Gene, D_Gene, J_Gene, CDR3_length, id)
      
           Number_of_reads = nrow(Selected_Columns)
    
           Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.GENE.and.allele.x, CDR3.IMGT.x) %>% tally())
      
           Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
      
           TheTable = cbind(Counts_of_unique_rows, Frequency)
      
           Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
           
           The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }
      

      
      
      

      
          else if (input$select_clonotype == "J Gene + CDR3 Amino Acids")
      
         { Selected_Columns = select(df_to_execute, J_Gene, CDR3.IMGT.y, V_Gene, D_Gene, CDR3_length, id)
      
           Number_of_reads = nrow(Selected_Columns)
      
           Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(J_Gene, CDR3.IMGT.y) %>% tally())
      
           Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
      
           TheTable = cbind(Counts_of_unique_rows, Frequency)
            
           Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
           
           The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }
      
     
           
           
           

        else if (input$select_clonotype == "J Gene + CDR3 Nucleotide")

        { Selected_Columns = select(df_to_execute, J_Gene, CDR3.IMGT.x, V_Gene, D_Gene, CDR3_length, id)

          Number_of_reads = nrow(Selected_Columns)

          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(J_Gene, CDR3.IMGT.x) %>% tally())

          Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)

          TheTable = cbind(Counts_of_unique_rows, Frequency)
          
          Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
          
          The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }







         else if (input$select_clonotype == "J Gene and Allele + CDR3 Amino Acids")

        { Selected_Columns = select(df_to_execute, J.GENE.and.allele.x, CDR3.IMGT.y, V_Gene, D_Gene, J_Gene, CDR3_length, id)

          Number_of_reads = nrow(Selected_Columns)

          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(J.GENE.and.allele.x, CDR3.IMGT.y) %>% tally())

          Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)

          TheTable = cbind(Counts_of_unique_rows, Frequency)

          Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
          
          The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }





         else if (input$select_clonotype == "J Gene and Allele + CDR3 Nucleotide")

        { Selected_Columns = select(df_to_execute, J.GENE.and.allele.x, CDR3.IMGT.x, V_Gene, D_Gene, J_Gene, CDR3_length, id)

          Number_of_reads = nrow(Selected_Columns)

          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(J.GENE.and.allele.x, CDR3.IMGT.x) %>% tally())

          Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)

          TheTable = cbind(Counts_of_unique_rows, Frequency)
          
          Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
          
          The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }






         else if (input$select_clonotype == "CDR3 Amino Acids")
      
         { Selected_Columns = select(df_to_execute, CDR3.IMGT.y, V_Gene, D_Gene, J_Gene, CDR3_length, id)
      
           Number_of_reads = nrow(Selected_Columns)
      
           Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(CDR3.IMGT.y) %>% tally())
      
           Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
      
           TheTable = cbind(Counts_of_unique_rows, Frequency)
           
           Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
           
           The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }
      
      


      
      
        else if (input$select_clonotype == "CDR3 Nucleotide")
      
         { Selected_Columns = select(df_to_execute, CDR3.IMGT.x, V_Gene, D_Gene, J_Gene, CDR3_length, id)
    
           Number_of_reads = nrow(Selected_Columns)
      
           Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(CDR3.IMGT.x) %>% tally())
    
           Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
      
           TheTable = cbind(Counts_of_unique_rows, Frequency)
            
           Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
           
           The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }
    
      
      
        
        
           if (input$select_clonotype == "V Gene + CDR3 Amino Acids") 
          
        { 
    
         Selected_Columns = select(df_to_execute, V_Gene, CDR3.IMGT.y, D_Gene, J_Gene, CDR3_length, id)
        
        Number_of_reads = nrow(Selected_Columns)
        
        Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V_Gene, CDR3.IMGT.y) %>% tally())
        
        Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
        
        TheTable = cbind(Counts_of_unique_rows, Frequency)

        Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
        
        The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }
        
        
         
        
        
       
        
        
          else if (input$select_clonotype == "V Gene + CDR3 Nucleotide")
        
         {
           Selected_Columns = select(df_to_execute, V_Gene, CDR3.IMGT.x, D_Gene, J_Gene, CDR3_length, id)
        
           Number_of_reads = nrow(Selected_Columns)
        
           Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V_Gene, CDR3.IMGT.x) %>% tally())
        
           Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
        
           TheTable = cbind(Counts_of_unique_rows, Frequency)
           
           Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
           
           The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }

        
        
        
          else if (input$select_clonotype == "V.D.J Nt sequence")
        
         { Selected_Columns = select(df_to_execute, V.D.J.REGION.x, V_Gene, D_Gene, J_Gene, CDR3_length, CDR3.IMGT.x, id)
        
         Number_of_reads = nrow(Selected_Columns)
        
         Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.D.J.REGION.x) %>% tally())
        
         Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
        
         TheTable = cbind(Counts_of_unique_rows, Frequency)
         
         Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
         
         The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }
        



        
        
        
        

          else if (input$select_clonotype == "V.D.J AA sequence")
        
         { Selected_Columns = select(df_to_execute, V.D.J.REGION.y, V_Gene, D_Gene, J_Gene, CDR3_length, CDR3.IMGT.y, id)
        
         Number_of_reads = nrow(Selected_Columns)
        
         Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.D.J.REGION.y) %>% tally())
        
         Frequency = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
        
         TheTable = cbind(Counts_of_unique_rows, Frequency)
         
         Clonotypes_Table = unique(left_join(TheTable, Selected_Columns))
         
         The_Dominant = Clonotypes_Table[which.max(Clonotypes_Table$Frequency),] }
        
        
        
        
        
        
        List_of_Clonotypes[[c]] = as.data.frame(Clonotypes_Table)
        
        List_of_Dominants[[c]] = as.data.frame(The_Dominant)
        
        ####################################################################################################
        
           
        #keep only the reads with the wanted frequency
         
        Clonotypes_with_desired_frequency = filter(Clonotypes_Table, Clonotypes_Table$Frequency >= input$numeric1 )
        
        
        
        Clonotypes_without_Dominant = Clonotypes_with_desired_frequency[which(Clonotypes_with_desired_frequency$Frequency != max(Clonotypes_with_desired_frequency$Frequency)),] 
       
     
          
        
          Clonotypes_with_same_CDR3_length_to_Dominant = filter(Clonotypes_without_Dominant, as.numeric(Clonotypes_without_Dominant$CDR3_length) == as.numeric(The_Dominant$CDR3_length))
          
      
        
        
           Column_of_V_Gene_others = as.character(Clonotypes_with_same_CDR3_length_to_Dominant$V_Gene)
        
          V_gene_family = unlist(strsplit(Column_of_V_Gene_others, "-"))[[1]]
          
       
          
          Clonotypes_with_same_CDR3_length_to_Dominant = cbind(Clonotypes_with_same_CDR3_length_to_Dominant, V_gene_family)
          
        
          
          Column_of_V_Gene_Dominant = as.character(The_Dominant$V_Gene)
          
          V_gene_family_Dominant = unlist(strsplit(Column_of_V_Gene_Dominant, "-"))[[1]]
          
          
          
          The_Dominant = cbind(The_Dominant, V_gene_family_Dominant)
          
          
        
          
          Clonotypes_with_same_CDR3_length_and_same_vgenefamily_to_Dom = Clonotypes_with_same_CDR3_length_to_Dominant[which(as.character(Clonotypes_with_same_CDR3_length_to_Dominant$V_gene_family) == as.character(The_Dominant$V_gene_family_Dominant)),]
          
  
          if ((input$select_clonotype == "V Gene + CDR3 Amino Acids") || (input$select_clonotype == "V Gene and Allele + CDR3 Amino Acids") || (input$select_clonotype == "J Gene + CDR3 Amino Acids") || (input$select_clonotype == "J Gene and Allele + CDR3 Amino Acids") || (input$select_clonotype == "CDR3 Amino Acids") || (input$select_clonotype =="V Nt sequence + CDR3 Amino Acids") || (input$select_clonotype == "V.D.J AA sequence"))

             
            {  CDR3_aminoacids_other = as.character(Clonotypes_with_same_CDR3_length_and_same_vgenefamily_to_Dom$CDR3.IMGT.y)

              CDR3_aminoacids_Dominant = as.character(The_Dominant$CDR3.IMGT.y)
              
              mismatches = stringdistmatrix(CDR3_aminoacids_other, CDR3_aminoacids_Dominant, method = "hamming")
              
              d = data.frame(CDR3_aminoacids_other, mismatches)
              
              Mismatches_filtering = filter(d, mismatches <= input$numeric2)
              
              The_related_clonos = Clonotypes_with_same_CDR3_length_and_same_vgenefamily_to_Dom[which(Clonotypes_with_same_CDR3_length_and_same_vgenefamily_to_Dom$CDR3.IMGT.y %in% Mismatches_filtering$CDR3_aminoacids_other), ]  
              
              Mismatches = as.character(Mismatches_filtering$mismatches)
              
              The_related_clonos = cbind(The_related_clonos, Mismatches)
              

              
              }
          
          else 
            
            { CDR3_nt_other = as.character(Clonotypes_with_same_CDR3_length_and_same_vgenefamily_to_Dom$CDR3.IMGT.x)
            
          
            
            CDR3_nt_Dominant = as.character(The_Dominant$CDR3.IMGT.x)
            
          
            mismatches = stringdistmatrix(CDR3_nt_other, CDR3_nt_Dominant, method = "hamming")
          
            d = data.frame(CDR3_nt_other, mismatches)
            
            Mismatches_filtering = filter(d, mismatches <= input$numeric2)
            
            The_related_clonos = Clonotypes_with_same_CDR3_length_and_same_vgenefamily_to_Dom[which(Clonotypes_with_same_CDR3_length_and_same_vgenefamily_to_Dom$CDR3.IMGT.x %in% Mismatches_filtering$CDR3_nt_other), ]  
            
            Mismatches = as.character(Mismatches_filtering$mismatches)
            
            The_related_clonos = cbind(The_related_clonos, Mismatches)
            
            
            
             }
              
              
              
              
              
              
              
              
              
             if (input$check == TRUE) {
              
             Mismatches_filtering = filter(Mismatches_filtering, str_detect(Mismatches_filtering$CDR3_aminoacids_other, input$pattern))
              
             }
             
             
               
              ################ CALCULATE  RELATED  READS #############################
              
              
              
                
                
              
            
               
                List_of_Related[[c]] = as.data.frame(The_related_clonos)
               
        }
        
        
              
                observeEvent(input$nameid3, {
        
                output$clonoTable = renderDataTable(List_of_Clonotypes[[input$nameid3]])
        
                output$thedominant = renderDataTable(List_of_Dominants[[input$nameid3]])    
                
                })
                
                
                output$related1 = renderDataTable(List_of_Related[[1]])
                
                Mean_Frequency_1 = as.numeric(mean(List_of_Related[[1]]$Frequency))
                
                output$meanF1 = renderText(paste("The Mean Frequency of Subclones is", Mean_Frequency_1))
                
                output$related2 = renderDataTable(List_of_Related[[2]])
                
                Mean_Frequency_2 =  as.numeric(mean(List_of_Related[[2]]$Frequency))
                
                output$meanF2 = renderText(paste("The Mean Frequency of Subclones is", Mean_Frequency_2))
                
                
                if ((input$select_clonotype == "V Gene + CDR3 Amino Acids") || (input$select_clonotype == "V Gene and Allele + CDR3 Amino Acids") || (input$select_clonotype == "J Gene + CDR3 Amino Acids") || (input$select_clonotype == "J Gene and Allele + CDR3 Amino Acids") || (input$select_clonotype == "CDR3 Amino Acids") || (input$select_clonotype =="V Nt sequence + CDR3 Amino Acids") || (input$select_clonotype == "V.D.J AA sequence"))
                
               
                  {
                  
                  
                 
                
        
                ALL_SAMPLES_RELATED_CLONOTYPES = merge_all(List_of_Related)
                
               
                
                output$relatedclonotypes = renderDataTable(ALL_SAMPLES_RELATED_CLONOTYPES)
      
                
                
                
                List_of_Frequencies = list()
                
                for (i in 1:length(List_of_Related)) {
                  
                  Number_of_related = as.numeric(nrow(List_of_Related[[i]]))
                  
                  Frequency = Number_of_related / as.numeric(nrow(ALL_SAMPLES_RELATED_CLONOTYPES))
                  id = unique(as.character(List_of_Related[[i]]$id))
                  Freq_df = data_frame(Frequency, id)
                  List_of_Frequencies[[i]] = Freq_df
                  
                }

                Frequencies = merge_all(List_of_Frequencies)
                
                
                 output$freq = renderDataTable(Frequencies)
                 
                
        
             
                  if (as.character(List_of_Dominants[[1]][,c(1,2)]) == as.character(List_of_Dominants[[2]][,c(1,2)]))
                     
                  { 
                    output$textdominant = renderText("SAME DOMINANT CLONE BETWEEN SAMPLES")  }
                 
                  else {
                    
                    Related_to_2nd_Dom_and_1st_Dominant = List_of_Related[[2]][which(as.character(List_of_Related[[2]][,c(1,2)]) == as.character(List_of_Dominants[[1]][,c(1,2)]))]
                    
                    Related_to_1st_Dom_and_2nd_Dominant = List_of_Related[[1]][which(as.character(List_of_Related[[1]][,c(1,2)]) == as.character(List_of_Dominants[[2]][,c(1,2)]))]
                  
                    output$ccc = renderDataTable(Related_to_2nd_Dom_and_1st_Dominant)
                    
                    output$ddd = renderDataTable(Related_to_1st_Dom_and_2nd_Dominant)
                    
                    }
            
               
             
               
               
               ############ 3 DATA FRAMES where (RELATED AND DOMINANT)
               
               
                List_for_graph_all = list(ALL_SAMPLES_RELATED_CLONOTYPES, List_of_Dominants[[1]], List_of_Dominants[[2]])
                
                ALL_RELATED_AND_DOMINANT = merge_all(List_for_graph_all)
               
               List_for_graph_1 = list(List_of_Related[[1]], List_of_Dominants[[1]])
               
               ONE_RELATED_AND_DOMINANT = merge_all(List_for_graph_1)
               
               List_for_graph_2 = list(List_of_Related[[2]], List_of_Dominants[[2]])
              
               TWO_RELATED_AND_DOMINANT = merge_all(List_for_graph_2)
               
              # Common_RELATED_AND_DOMINANT = inner_join(ONE_RELATED_AND_DOMINANT, TWO_RELATED_AND_DOMINANT, by)
               
              ######################################################### 
               
               
             
               ########### for 1 ##################################
               
                 
               Common_Subclones = inner_join(List_of_Related[[1]], List_of_Related[[2]], by = "CDR3.IMGT.y")
               
               
               
               # "%notin%" = Negate("%in%")
               # 
               #  No_common_one = List_of_Related[[1]][which(List_of_Related[[1]]$CDR3.IMGT.y %notin% Common_Subclones$CDR3.IMGT.y),]
               # 
               #  view(No_common_one)
               
               Frequency = as.character(ONE_RELATED_AND_DOMINANT$Frequency)
                 
               n = as.character(ONE_RELATED_AND_DOMINANT$n)
               
               V_Gene = as.character(ONE_RELATED_AND_DOMINANT$V_Gene)
               
               D_Gene = as.character(ONE_RELATED_AND_DOMINANT$D_Gene)
               
               J_Gene = as.character(ONE_RELATED_AND_DOMINANT$J_Gene)
               
               CDR3_length = as.character(ONE_RELATED_AND_DOMINANT$CDR3_length)
               
               id = as.character(ONE_RELATED_AND_DOMINANT$id)
               
               V_gene_family = as.character(ONE_RELATED_AND_DOMINANT$V_gene_family)
               
               Mismatches = as.character(ONE_RELATED_AND_DOMINANT$Mismatches)
             
               CDR3_AA_SUBCLONES = ONE_RELATED_AND_DOMINANT$CDR3.IMGT.y
               
               
               
               CDR3_AA_SUBCLONES_df = as.data.frame(CDR3_AA_SUBCLONES)
               
               CDR3_AA_SUBCLONES_df = cbind(CDR3_AA_SUBCLONES_df, Frequency, n,  V_Gene, D_Gene, J_Gene, CDR3_length, id, V_gene_family, Mismatches)
               
               
              
               
               mismatches_between_subclones = stringdistmatrix(CDR3_AA_SUBCLONES_df$CDR3_AA_SUBCLONES, CDR3_AA_SUBCLONES_df$CDR3_AA_SUBCLONES, method = "hamming")
               
               mismatches_between_subclones[mismatches_between_subclones > 1] = 0
               
               # write.table(mismatches_between_subclones, file = "Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE )
               
               mismatches_between_subclones = adjMatrix(mismatches_between_subclones)
               
               CDR3_AA_SUBCLONES_df$id = as.character(CDR3_AA_SUBCLONES_df$id)
               
               CDR3_AA_SUBCLONES_df[which(as.character(CDR3_AA_SUBCLONES_df$CDR3_AA_SUBCLONES) %in% as.character(Common_Subclones$CDR3.IMGT.y)),]$id  = "IN BOTH SAMPLES"
               
              
                if (as.character(List_of_Dominants[[1]][,c(1,2)]) == as.character(List_of_Dominants[[2]][,c(1,2)])) 
                
                 { CDR3_AA_SUBCLONES_df[which(as.character(CDR3_AA_SUBCLONES_df$CDR3_AA_SUBCLONES) %in% as.character(List_of_Dominants[[1]]$CDR3.IMGT.y)),]$id  = "DOMINANT" }
                
                
             
               output$graph1 = renderDataTable(CDR3_AA_SUBCLONES_df)
               
               output$igraph1 = renderPlot({  TheGraph = graph_from_adjacency_matrix(mismatches_between_subclones, mode = "directed")
               
             
               
               V(TheGraph)$color <- as.factor(CDR3_AA_SUBCLONES_df$id)
               
                V(TheGraph)$size <- CDR3_AA_SUBCLONES_df$Frequency
              
               plot(TheGraph, edge.arrow.size=.2, edge.color="orange",
                    # vertex.color="#e80000",
                    vertex.frame.color="#d18b71",
                    vertex.label.color="black", vertex.label=row_number(CDR3_AA_SUBCLONES)) 
               
               
               })  
              
               
                Frequency_ = as.character(TWO_RELATED_AND_DOMINANT$Frequency)

                n_ = as.character(TWO_RELATED_AND_DOMINANT$n)

                V_Gene_ = as.character(TWO_RELATED_AND_DOMINANT$V_Gene)

                D_Gene_ = as.character(TWO_RELATED_AND_DOMINANT$D_Gene)

                J_Gene_ = as.character(TWO_RELATED_AND_DOMINANT$J_Gene)

                CDR3_length_ = as.character(TWO_RELATED_AND_DOMINANT$CDR3_length)

                id_ = as.character(TWO_RELATED_AND_DOMINANT$id)

                V_gene_family_ = as.character(TWO_RELATED_AND_DOMINANT$V_gene_family)

                Mismatches_ = as.character(TWO_RELATED_AND_DOMINANT$Mismatches)

                CDR3_AA_SUBCLONES_2 = TWO_RELATED_AND_DOMINANT$CDR3.IMGT.y



                CDR3_AA_SUBCLONES_2_df = as.data.frame(CDR3_AA_SUBCLONES_2)

                CDR3_AA_SUBCLONES_2_df = cbind(CDR3_AA_SUBCLONES_2_df, Frequency_, n_,  V_Gene_, D_Gene_, J_Gene_, CDR3_length_, id_, V_gene_family_, Mismatches_)




                mismatches_between_subclones_ = stringdistmatrix(CDR3_AA_SUBCLONES_2_df$CDR3_AA_SUBCLONES_2, CDR3_AA_SUBCLONES_2_df$CDR3_AA_SUBCLONES_2, method = "hamming")

                mismatches_between_subclones_[mismatches_between_subclones_ > 1] = 0

               # # write.table(mismatches_between_subclones, file = "Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE )

                mismatches_between_subclones_ = adjMatrix(mismatches_between_subclones_)

                CDR3_AA_SUBCLONES_2_df$id = as.character(CDR3_AA_SUBCLONES_2_df$id)

                CDR3_AA_SUBCLONES_2_df[which(as.character(CDR3_AA_SUBCLONES_2_df$CDR3_AA_SUBCLONES_2) %in% as.character(Common_Subclones$CDR3.IMGT.y)),]$id  = "IN BOTH SAMPLES"


                if (as.character(List_of_Dominants[[1]][,c(1,2)]) == as.character(List_of_Dominants[[2]][,c(1,2)]))

                { CDR3_AA_SUBCLONES_2_df[which(as.character(CDR3_AA_SUBCLONES_2_df$CDR3_AA_SUBCLONES_2) %in% as.character(List_of_Dominants[[2]]$CDR3.IMGT.y)),]$id  = "DOMINANT" }



                output$graph2 = renderDataTable(CDR3_AA_SUBCLONES_2_df)

                output$igraph2 = renderPlot({  TheGraph2 = graph_from_adjacency_matrix(mismatches_between_subclones_, mode = "directed")

                V(TheGraph2)$color <- as.factor(CDR3_AA_SUBCLONES_2_df$id)

                V(TheGraph2)$size <- CDR3_AA_SUBCLONES_2_df$Frequency

                plot(TheGraph2, edge.arrow.size=.2, edge.color="orange",
                    # vertex.color="#e80000",
                     vertex.frame.color="#d18b71",
                     vertex.label.color="black", vertex.label=row_number(CDR3_AA_SUBCLONES_2))


                })
                
                
                DF_FOR_PLOT = full_join(List_of_Related[[1]], List_of_Related[[2]], by="CDR3.IMGT.y")
                
                  if (NA %in% DF_FOR_PLOT$Frequency.y)
                
               { DF_FOR_PLOT[which(is.na(DF_FOR_PLOT$Frequency.y)),]$Frequency.y = 0 } 
                
                if (NA %in% DF_FOR_PLOT$Frequency.x) 
                
                { DF_FOR_PLOT[which(is.na(DF_FOR_PLOT$Frequency.x)),]$Frequency.x = 0 }
                

               output$dotplot = renderPlot({ ggplot(DF_FOR_PLOT, aes(x = Frequency.x, y = Frequency.y) ) + 
                   geom_point(size = 3) + geom_label_repel(aes(label = rownames(DF_FOR_PLOT)), box.padding   = 0.35, point.padding = 0.5,
                 segment.color = 'red')  +  scale_x_continuous(limits = c(0, 0.005)) + scale_y_continuous(limits = c(0, 0.005)) + geom_smooth() 
                 })
               
               output$dottable = renderDataTable(DF_FOR_PLOT)
               
               DF_OF_DOMINANT = full_join(List_of_Dominants[[1]], List_of_Dominants[[2]], by = "CDR3.IMGT.y")
               
               output$domdotplot = renderPlot({ ggplot(DF_OF_DOMINANT, aes(x = Frequency.x, y = Frequency.y) ) + geom_point(size = 4)
                 })
               
               output$domdottable = renderDataTable(DF_OF_DOMINANT)
              
               
                 }
              
              
              
              
              
              
              
              
              
              
              
              
              
                
                
                else 
                  
                {
                  
                  
                  
                  
                  
                  
                  
                  ALL_SAMPLES_RELATED_CLONOTYPES = merge_all(List_of_Related)
                  
                  
                  
                  output$relatedclonotypes = renderDataTable(ALL_SAMPLES_RELATED_CLONOTYPES)
                  
                  
                  
                  
                  List_of_Frequencies = list()
                  
                  for (i in 1:length(List_of_Related)) {
                    
                    Number_of_related = as.numeric(nrow(List_of_Related[[i]]))
                    
                    Frequency = Number_of_related / as.numeric(nrow(ALL_SAMPLES_RELATED_CLONOTYPES))
                    id = unique(as.character(List_of_Related[[i]]$id))
                    Freq_df = data_frame(Frequency, id)
                    List_of_Frequencies[[i]] = Freq_df
                    
                  }
                  
                  Frequencies = merge_all(List_of_Frequencies)
                  
                  
                  output$freq = renderDataTable(Frequencies)
                  
                  
                  
                  
                  if (as.character(List_of_Dominants[[1]][,c(1,2)]) == as.character(List_of_Dominants[[2]][,c(1,2)]))
                    
                  { 
                    output$textdominant = renderText("SAME DOMINANT CLONE BETWEEN SAMPLES")  }
                  
                  else {
                    
                    Related_to_2nd_Dom_and_1st_Dominant = List_of_Related[[2]][which(as.character(List_of_Related[[2]][,c(1,2)]) == as.character(List_of_Dominants[[1]][,c(1,2)]))]
                    
                    Related_to_1st_Dom_and_2nd_Dominant = List_of_Related[[1]][which(as.character(List_of_Related[[1]][,c(1,2)]) == as.character(List_of_Dominants[[2]][,c(1,2)]))]
                    
                    output$ccc = renderDataTable(Related_to_2nd_Dom_and_1st_Dominant)
                    
                    output$ddd = renderDataTable(Related_to_1st_Dom_and_2nd_Dominant)
                    
                  }
                  
                  
                  List_for_graph_all = list(ALL_SAMPLES_RELATED_CLONOTYPES, List_of_Dominants[[1]], List_of_Dominants[[2]])
                  
                  ALL_RELATED_AND_DOMINANT = merge_all(List_for_graph_all)
                  
                  List_for_graph_1 = list(List_of_Related[[1]], List_of_Dominants[[1]])
                  
                  ONE_RELATED_AND_DOMINANT = merge_all(List_for_graph_1)
                  
                  List_for_graph_2 = list(List_of_Related[[2]], List_of_Dominants[[2]])
                  
                  TWO_RELATED_AND_DOMINANT = merge_all(List_for_graph_2)
                  
                  print(ONE_RELATED_AND_DOMINANT)
                  
                  print(TWO_RELATED_AND_DOMINANT)
                  
                  ######################################################### 
                  
                  
                  
                  ########### for 1 ##################################
                  
                  
                  Common_Subclones = inner_join(List_of_Related[[1]], List_of_Related[[2]], by = "CDR3.IMGT.x")
                  
                  
                  
                  Frequency = as.character(ONE_RELATED_AND_DOMINANT$Frequency)
                  
                  n = as.character(ONE_RELATED_AND_DOMINANT$n)
                  
                  V_Gene = as.character(ONE_RELATED_AND_DOMINANT$V_Gene)
                  
                  D_Gene = as.character(ONE_RELATED_AND_DOMINANT$D_Gene)
                  
                  J_Gene = as.character(ONE_RELATED_AND_DOMINANT$J_Gene)
                  
                  CDR3_length = as.character(ONE_RELATED_AND_DOMINANT$CDR3_length)
                  
                  id = as.character(ONE_RELATED_AND_DOMINANT$id)
                  
                  V_gene_family = as.character(ONE_RELATED_AND_DOMINANT$V_gene_family)
                  
                  Mismatches = as.character(ONE_RELATED_AND_DOMINANT$Mismatches)
                  
                  CDR3_nt_SUBCLONES = ONE_RELATED_AND_DOMINANT$CDR3.IMGT.x
                  
                  
                  
                  CDR3_nt_SUBCLONES_df = as.data.frame(CDR3_nt_SUBCLONES)
                  
                  CDR3_nt_SUBCLONES_df = cbind(CDR3_nt_SUBCLONES_df, Frequency, n,  V_Gene, D_Gene, J_Gene, CDR3_length, id, V_gene_family, Mismatches)
                  
                  
                  
                  
                  mismatches_between_subclones = stringdistmatrix(CDR3_nt_SUBCLONES_df$CDR3_nt_SUBCLONES, CDR3_nt_SUBCLONES_df$CDR3_nt_SUBCLONES, method = "hamming")
                  
                  mismatches_between_subclones[mismatches_between_subclones > 1] = 0
                  
                  # write.table(mismatches_between_subclones, file = "Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE )
                  
                  mismatches_between_subclones = adjMatrix(mismatches_between_subclones)
                  
                  CDR3_nt_SUBCLONES_df$id = as.character(CDR3_nt_SUBCLONES_df$id)
                  
                  CDR3_nt_SUBCLONES_df[which(as.character(CDR3_nt_SUBCLONES_df$CDR3_nt_SUBCLONES) %in% as.character(Common_Subclones$CDR3.IMGT.x)),]$id  = "IN BOTH SAMPLES"
                  
                  
                  if (as.character(List_of_Dominants[[1]][,c(1,2)]) == as.character(List_of_Dominants[[2]][,c(1,2)])) 
                    
                  { CDR3_nt_SUBCLONES_df[which(as.character(CDR3_nt_SUBCLONES_df$CDR3_nt_SUBCLONES) %in% as.character(List_of_Dominants[[1]]$CDR3.IMGT.x)),]$id  = "DOMINANT" }
                  
                  
                  
                  output$graph1 = renderDataTable(CDR3_nt_SUBCLONES_df)
                  
                  output$igraph1 = renderPlot({  TheGraph = graph_from_adjacency_matrix(mismatches_between_subclones, mode = "directed")
                  
                  V(TheGraph)$color <- as.factor(CDR3_nt_SUBCLONES_df$id)
                  
                  V(TheGraph)$size <- CDR3_nt_SUBCLONES_df$Frequency
                  
                  plot(TheGraph, edge.arrow.size=.2, edge.color="orange",
                       # vertex.color="#e80000",
                       vertex.frame.color="#d18b71",
                       vertex.label.color="black", vertex.label=row_number(CDR3_nt_SUBCLONES)) 
                  
                  
                  })  
                  
                  
                  Frequency_ = as.character(TWO_RELATED_AND_DOMINANT$Frequency)
                  
                  n_ = as.character(TWO_RELATED_AND_DOMINANT$n)
                  
                  V_Gene_ = as.character(TWO_RELATED_AND_DOMINANT$V_Gene)
                  
                  D_Gene_ = as.character(TWO_RELATED_AND_DOMINANT$D_Gene)
                  
                  J_Gene_ = as.character(TWO_RELATED_AND_DOMINANT$J_Gene)
                  
                  CDR3_length_ = as.character(TWO_RELATED_AND_DOMINANT$CDR3_length)
                  
                  id_ = as.character(TWO_RELATED_AND_DOMINANT$id)
                  
                  V_gene_family_ = as.character(TWO_RELATED_AND_DOMINANT$V_gene_family)
                  
                  Mismatches_ = as.character(TWO_RELATED_AND_DOMINANT$Mismatches)
                  
                  CDR3_nt_SUBCLONES_2 = TWO_RELATED_AND_DOMINANT$CDR3.IMGT.x
                  
                  
                  
                  CDR3_nt_SUBCLONES_2_df = as.data.frame(CDR3_nt_SUBCLONES_2)
                  
                  CDR3_nt_SUBCLONES_2_df = cbind(CDR3_nt_SUBCLONES_2_df, Frequency_, n_,  V_Gene_, D_Gene_, J_Gene_, CDR3_length_, id_, V_gene_family_, Mismatches_)
                  
                  
                  
                  
                  mismatches_between_subclones_ = stringdistmatrix(CDR3_nt_SUBCLONES_2_df$CDR3_nt_SUBCLONES_2, CDR3_nt_SUBCLONES_2_df$CDR3_nt_SUBCLONES_2, method = "hamming")
                  
                  mismatches_between_subclones_[mismatches_between_subclones_ > 1] = 0
                  
                  # # write.table(mismatches_between_subclones, file = "Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE )
                  
                  mismatches_between_subclones_ = adjMatrix(mismatches_between_subclones_)
                  
                  CDR3_nt_SUBCLONES_2_df$id = as.character(CDR3_nt_SUBCLONES_2_df$id)
                  
                  CDR3_nt_SUBCLONES_2_df[which(as.character(CDR3_nt_SUBCLONES_2_df$CDR3_nt_SUBCLONES_2) %in% as.character(Common_Subclones$CDR3.IMGT.x)),]$id  = "IN BOTH SAMPLES"
                  
                  
                  if (as.character(List_of_Dominants[[1]][,c(1,2)]) == as.character(List_of_Dominants[[2]][,c(1,2)]))
                    
                  { CDR3_nt_SUBCLONES_2_df[which(as.character(CDR3_nt_SUBCLONES_2_df$CDR3_nt_SUBCLONES_2) %in% as.character(List_of_Dominants[[2]]$CDR3.IMGT.x)),]$id  = "DOMINANT" }
                  
                  
                  
                  output$graph2 = renderDataTable(CDR3_nt_SUBCLONES_2_df)
                  
                  output$igraph2 = renderPlot({  TheGraph2 = graph_from_adjacency_matrix(mismatches_between_subclones_, mode = "directed")
                  
                  V(TheGraph2)$color <- as.factor(CDR3_nt_SUBCLONES_2_df$id)
                  
                  V(TheGraph2)$size <- CDR3_nt_SUBCLONES_2_df$Frequency
                  
                  plot(TheGraph2, edge.arrow.size=.2, edge.color="orange",
                       # vertex.color="#e80000",
                       vertex.frame.color="#d18b71",
                       vertex.label.color="black", vertex.label=row_number(CDR3_nt_SUBCLONES_2))
                  
                  
                  })
                  
                  DF_FOR_PLOT = full_join(List_of_Related[[1]], List_of_Related[[2]], by="CDR3.IMGT.x")
                  
                  if (NA %in% DF_FOR_PLOT$Frequency.y)
                    
                  { DF_FOR_PLOT[which(is.na(DF_FOR_PLOT$Frequency.y)),]$Frequency.y = 0 } 
                  
                  if (NA %in% DF_FOR_PLOT$Frequency.x) 
                    
                  { DF_FOR_PLOT[which(is.na(DF_FOR_PLOT$Frequency.x)),]$Frequency.x = 0 }
                  
                  
                  output$dotplot = renderPlot({ ggplot(DF_FOR_PLOT, aes(x = Frequency.x, y = Frequency.y) ) + 
                      geom_point(size = 3) + geom_label_repel(aes(label = rownames(DF_FOR_PLOT)), box.padding   = 0.35, point.padding = 0.5,
                                                              segment.color = 'red')  +  scale_x_continuous(limits = c(0, 0.005)) + scale_y_continuous(limits = c(0, 0.005)) + geom_smooth() 
                  })
                  
                  output$dottable = renderDataTable(DF_FOR_PLOT)
                  
                  DF_OF_DOMINANT = full_join(List_of_Dominants[[1]], List_of_Dominants[[2]], by = "CDR3.IMGT.x")
                  
                  output$domdotplot = renderPlot({ ggplot(DF_OF_DOMINANT, aes(x = Frequency.x, y = Frequency.y) ) + geom_point(size = 4)
                  })
                  
                  output$domdottable = renderDataTable(DF_OF_DOMINANT)
                  
                  
                  
                  
                  
                  
                  
                  
                }
              
              ####################
              
              ########## succed ##########
           
              observeEvent(input$Executepipeline, { 
                
                shinyalert(title = "Pipeline Executed!", type = "success") })
           
         
         
             observeEvent(input$downloadallclonotypes, {
                # write.table(TheTable, "All_Clonotypes.txt", row.names = FALSE, sep = "\t", quote = FALSE) 
               
               shinyalert(title = "All_Clonotypes Downloaded", type = "success")   } )
             
            
              observeEvent(input$downloadthedominantclone, {
              #  write.table(Dominant_Reads_all_info, "Dominant_Clone.txt", row.names = FALSE, sep = "\t", quote = FALSE)
               
               shinyalert(title = "Dominant_Clone Downloaded", type = "success")  } )
               
             
             observeEvent(input$downloadreadsrelatedtodominant, {
                write.table(ALL_SAMPLES_RELATED_READS, "Reads_related_to_the_Dominant.txt", row.names = FALSE, sep = "\t", quote = FALSE) 
               
               shinyalert(title = "Reads_related_to_the_Dominant Downloaded", type = "success")} )
             
             
            
              
              
              
              
              
              
             
           } )  
       
      
      
       
        
        
        
      
    
    } 
)
  
        
      
      
      

    
  

  
