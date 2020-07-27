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
     
      ThePath = file.path("C:","Users","tgala","Desktop","Biosciences","TOOLS - BIOINFORMATICS","clonal R evolution","patient-files","IG1P")
      
       The_File = list()
    
      original.data = eventReactive(input$LoadData, {
        
      
       
            List_of_Files = as.list(input$inputFiles)
       
             for (i in List_of_Files) { 
             
                File_Name_Without_Quotes = noquote(i)
             
                File_Path = paste(ThePath, File_Name_Without_Quotes, sep = "/")
             
                The_File[[i]] = read.csv(file = File_Path, header = TRUE, sep = "\t")
             
                                   }
        
            MergedFiles = merge(The_File[[1]], The_File[[2]], by="Sequence.ID")
            
            Only_needed_columns = select(MergedFiles, Sequence.ID, Sequence.number.x, V.DOMAIN.Functionality.x, V.GENE.and.allele.x, D.GENE.and.allele.x, V.D.J.REGION.x, V.J.REGION.x, V.REGION.x, FR1.IMGT.x, CDR1.IMGT.x, FR2.IMGT.x, CDR2.IMGT.x, FR3.IMGT.x, CDR3.IMGT.x, JUNCTION.x, X3.V.REGION, X.N.D..J.REGION, X.N.D..REGION, D.REGION, N2.REGION, J.REGION.x, FR4.IMGT.x, J.GENE.and.allele.y, V.D.J.REGION.y, V.REGION.y, FR1.IMGT.y, CDR1.IMGT.y, FR2.IMGT.y, CDR2.IMGT.y, FR3.IMGT.y, CDR3.IMGT.y, JUNCTION.y, J.REGION.y, FR4.IMGT.y)
        
                                                        } )
     
   
      output$table <- renderDataTable(original.data())
      
     
       observeEvent(input$LoadData, { 
        
        shinyalert(title = "Loading Data. . .", type = "success") })
     
     
        filtered.Data <- eventReactive(input$FilterButton, { 
        
         Selected_Column = as.numeric(input$select)
        
         Texted_Characteristic = as.character(input$text)
        
         The_Colnames = colnames(original.data())
        
         Choose_from_the_colnames = as.character(The_Colnames[[Selected_Column]])
        
         TheData = as.data.frame(original.data())
        
         Column_to_choose = factor(TheData[,Choose_from_the_colnames]) 
        
         The.filtered = filter(TheData, Column_to_choose == Texted_Characteristic)
         
          The.filtered$V_Gene = substring(The.filtered$V.GENE.and.allele.x, 1, 14)
          
          The.filtered = as.data.frame(The.filtered)
         
        } )  
      
     
      
      
      output$filtered <- renderDataTable(filtered.Data())
      
      observeEvent(input$FilterButton, { 
        
        shinyalert(title = "Filtering succeed", type = "success") })
      
      
      
      observeEvent(input$Executepipeline, {
        
     
      
          if (input$select_clonotype == "V region + CDR3 Amino Acids") 
          
        { Selected_Columns = select(filtered.Data(), V.REGION.x, CDR3.IMGT.y)
        
        Number_of_reads = nrow(Selected_Columns)
        
        Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.REGION.x, CDR3.IMGT.y) %>% tally())
        
        Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
        
        TheTable = cbind(Counts_of_unique_rows, Frequencyy)
        
        The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
        
        output$clonoTable = renderDataTable(TheTable)
        
        output$thedominant = renderDataTable(The_Dominant)
        
          } 
        
           
      else if (input$select_clonotype == "V Gene and Allele + CDR3 Amino Acids") 
        
        { Selected_Columns = select(filtered.Data(), V.GENE.and.allele.x, CDR3.IMGT.y)
      
      Number_of_reads = nrow(Selected_Columns)
      
      Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.GENE.and.allele.x, CDR3.IMGT.y) %>% tally())
      
      Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
      
      TheTable = cbind(Counts_of_unique_rows, Frequencyy)
      
      The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
      
      output$clonoTable = renderDataTable(TheTable)
      
      output$thedominant = renderDataTable(The_Dominant)
      
      }
        
        
       else if (input$select_clonotype == "V region + CDR3 Nucleotide")
        
        { Selected_Columns = select(filtered.Data(), V.REGION.x, CDR3.IMGT.x)
          
          Number_of_reads = nrow(Selected_Columns)
          
          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.REGION.x, CDR3.IMGT.x) %>% tally())
          
          Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
          
          TheTable = cbind(Counts_of_unique_rows, Frequencyy)
          
          The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
          
          output$clonoTable = renderDataTable(TheTable)
          
          output$thedominant = renderDataTable(The_Dominant) }
        
        
        
        else if (input$select_clonotype == "V Gene and Allele + CDR3 Nucleotide") 
          
        { Selected_Columns = select(filtered.Data(), V.GENE.and.allele.x, CDR3.IMGT.x)
          
          Number_of_reads = nrow(Selected_Columns)
          
          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V.GENE.and.allele.x, CDR3.IMGT.x) %>% tally())
          
          Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
          
          TheTable = cbind(Counts_of_unique_rows, Frequencyy)
          
          The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
          
          output$clonoTable = renderDataTable(TheTable)
          
          output$thedominant = renderDataTable(The_Dominant) }
        
        
      
        else if (input$select_clonotype == "J Gene + CDR3 Amino Acids")
          
        { Selected_Columns = select(filtered.Data(), J.REGION.x, CDR3.IMGT.y)
          
          Number_of_reads = nrow(Selected_Columns)
          
          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(J.REGION.x, CDR3.IMGT.y) %>% tally())
          
          Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
          
          TheTable = cbind(Counts_of_unique_rows, Frequencyy)
          
          The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
          
          output$clonoTable = renderDataTable(TheTable)
          
          output$thedominant = renderDataTable(The_Dominant) }
        
        
        
        else if (input$select_clonotype == "J Gene + CDR3 Nucleotide")
          
        { Selected_Columns = select(filtered.Data(), J.REGION.x, CDR3.IMGT.x)
          
          Number_of_reads = nrow(Selected_Columns)
          
          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(J.REGION.x, CDR3.IMGT.x) %>% tally())
          
          Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
          
          TheTable = cbind(Counts_of_unique_rows, Frequencyy)
          
          The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
          
          output$clonoTable = renderDataTable(TheTable)
          
          output$thedominant = renderDataTable(The_Dominant) }
        
        
        
        else if (input$select_clonotype == "J Gene and Allele + CDR3 Amino Acids") 
          
        { Selected_Columns = select(filtered.Data(), J.GENE.and.allele.x, CDR3.IMGT.y)
          
          Number_of_reads = nrow(Selected_Columns)
          
          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(J.GENE.and.allele.x, CDR3.IMGT.y) %>% tally())
          
          Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
          
          TheTable = cbind(Counts_of_unique_rows, Frequencyy)
          
          The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
          
          output$clonoTable = renderDataTable(TheTable)
          
          output$thedominant = renderDataTable(The_Dominant) }
        
        
        
        else if (input$select_clonotype == "J Gene and Allele + CDR3 Nucleotide") 
          
        { Selected_Columns = select(filtered.Data(), J.GENE.and.allele.x, CDR3.IMGT.x)
          
          Number_of_reads = nrow(Selected_Columns)
          
          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(J.GENE.and.allele.x, CDR3.IMGT.x) %>% tally())
          
          Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
          
          TheTable = cbind(Counts_of_unique_rows, Frequencyy)
          
          The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
          
          output$clonoTable = renderDataTable(TheTable)
          
          output$thedominant = renderDataTable(The_Dominant) }
        
        
        
        else if (input$select_clonotype == "CDR3 Amino Acids") 
          
        { Selected_Columns = select(filtered.Data(), CDR3.IMGT.y)
          
          Number_of_reads = nrow(Selected_Columns)
          
          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(CDR3.IMGT.y) %>% tally())
          
          Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
          
          TheTable = cbind(Counts_of_unique_rows, Frequencyy)
          
          The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
          
          output$clonoTable = renderDataTable(TheTable)
          
          output$thedominant = renderDataTable(The_Dominant) }
        
        
        
        else if (input$select_clonotype == "CDR3 Nucleotide") 
          
        { Selected_Columns = select(filtered.Data(), CDR3.IMGT.x)
          
          Number_of_reads = nrow(Selected_Columns)
          
          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(CDR3.IMGT.x) %>% tally())
          
          Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
          
          TheTable = cbind(Counts_of_unique_rows, Frequencyy)
          
          The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
          
          output$clonoTable = renderDataTable(TheTable)
          
          output$thedominant = renderDataTable(The_Dominant) }
        
        else if (input$select_clonotype == "V Gene + CDR3 Amino Acids") 
          
        { 
          
        Selected_Columns = select(filtered.Data(), V.GENE.and.allele.x, CDR3.IMGT.y)
        
        Selected_Columns$V_Gene = substring(Selected_Columns$V.GENE.and.allele.x, 1, 14)
        
        Selected_Columns = select(Selected_Columns, V_Gene, CDR3.IMGT.y)
        
        Number_of_reads = nrow(Selected_Columns)
        
        Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V_Gene, CDR3.IMGT.y) %>% tally())
        
        Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
        
        TheTable = cbind(Counts_of_unique_rows, Frequencyy)
        
        The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
        
        output$clonoTable = renderDataTable(TheTable)
        
        output$thedominant = renderDataTable(The_Dominant) }
        
        
        else if (input$select_clonotype == "V Gene + CDR3 Nucleotide") 
          
        { 
          
          Selected_Columns = select(filtered.Data(), V.GENE.and.allele.x, CDR3.IMGT.x)
          
          print("nikos")
          
          Selected_Columns$V_Gene = substring(Selected_Columns$V.GENE.and.allele.x, 1, 14)
          
          print("achilleas")
          
          Selected_Columns = select(Selected_Columns, V_Gene, CDR3.IMGT.x)
          
          Number_of_reads = nrow(Selected_Columns)
          
          Counts_of_unique_rows = as.data.frame(Selected_Columns %>% group_by(V_Gene, CDR3.IMGT.x) %>% tally())
          
          Frequencyy = as.numeric(Counts_of_unique_rows$n / Number_of_reads)
          
          TheTable = cbind(Counts_of_unique_rows, Frequencyy)
          
          The_Dominant = TheTable[which.max(TheTable$Frequencyy),]
          
          output$clonoTable = renderDataTable(TheTable)
          
          output$thedominant = renderDataTable(The_Dominant) }
        
        
        
        
        
        
        
        #keep only the reads with the wanted frequency
         
        Reads_with_desired_frequency = filter(TheTable, TheTable$Frequencyy >= input$numeric1 )
         
         
        #Dominant READS
         Dominant_Reads = Reads_with_desired_frequency[which.max(Reads_with_desired_frequency$Frequencyy),]
         
         write.table(Dominant_Reads, "Dominant_Reads.txt", row.names = FALSE, sep = "\t", quote = FALSE)
         
         Dominant_Reads_all_info = unique(inner_join(Dominant_Reads, filtered.Data()))
         
         write.table(Dominant_Reads_all_info, "Dominant_Reads_all_info.txt", sep = "\t", row.names = FALSE, quote = FALSE)
         
         
        
         #READS_exempt_Dominant
         
         Reads_without_Dominant = Reads_with_desired_frequency[which(Reads_with_desired_frequency$Frequencyy != max(Reads_with_desired_frequency$Frequencyy)),]
         
         write.table(Reads_without_Dominant, "Reads_without_Dominant.txt", sep = "\t", row.names = FALSE, quote = FALSE)
         
         Reads_without_Dominant_all_info = unique(inner_join(Reads_without_Dominant, filtered.Data()))
         
         write.table(Reads_without_Dominant_all_info, "Reads_without_Dominant_all_info.txt", sep = "\t", row.names = FALSE, quote = FALSE)
         
       
         
         ######### for CDR3 association to Dominant ########
         
          CDR3n_Dominant = select(Dominant_Reads_all_info, CDR3.IMGT.x)  
         
          Length_of_CDR3n_Dominant = as.numeric(unique(apply(CDR3n_Dominant, 2, nchar)))
          
          ########   #######
          
          
          Reads_with_same_CDR3n_length_to_Dominant = filter(Reads_without_Dominant_all_info, nchar(as.character(Reads_without_Dominant_all_info$CDR3.IMGT.x)) == Length_of_CDR3n_Dominant)
          
         
          
          ########## calculate V gene family of DOMINANT and  READS (without Dominant) #####
          
          V_Gene_and_allele_Dominant = as.data.frame(select(Dominant_Reads_all_info, V.GENE.and.allele.x))
          
          V_Gene_and_allele_Dominant$V_Gene_group_Dominant = substring(V_Gene_and_allele_Dominant$V.GENE.and.allele.x, 1, 12)
          
          
          
          
          
           V_Gene_and_allele_Reads_without_Dominant = as.data.frame(select(Reads_with_same_CDR3n_length_to_Dominant, V.GENE.and.allele.x))
         
           V_Gene_and_allele_Reads_without_Dominant$V_Gene_group_Reads_without_Dominant = substring(V_Gene_and_allele_Reads_without_Dominant$V.GENE.and.allele.x, 1, 12)
         
            
          
            
            ######## filtering ..to find wich reads have the same V gene group with Dominant Clone###########
          
            
           
           Related_READS = filter(V_Gene_and_allele_Reads_without_Dominant, V_Gene_group_Reads_without_Dominant == unique(V_Gene_and_allele_Dominant$V_Gene_group_Dominant))
           
           
           
          ############# give all the information of the Related to Dominant READS #############
             
           
              Related_READS_all_info = Reads_with_same_CDR3n_length_to_Dominant[which(Reads_with_same_CDR3n_length_to_Dominant$V.GENE.and.allele.x %in% unique(Related_READS$V.GENE.and.allele.x)), ]
               
              Related_READS_all_info = unique(Related_READS_all_info)
              
              write.table(Related_READS_all_info, "Related_READS_all_info.txt", sep = "\t", row.names = FALSE, quote = FALSE )
              
              
              
              CDR3_aminoacids_Related = as.character(Related_READS_all_info$CDR3.IMGT.y)
              
              CDR3_aminoacids_Dominant = unique(as.character(Dominant_Reads_all_info$CDR3.IMGT.y))
              
              
              
              n = stringdistmatrix(CDR3_aminoacids_Related, CDR3_aminoacids_Dominant, method = "hamming")
              
              
              
              d = data.frame(CDR3_aminoacids_Dominant, CDR3_aminoacids_Related, n)
              
              Mismatches_filtering = filter(d, n <= input$numeric2)
              
              
              
             if (input$check == TRUE) {
              
             Mismatches_filtering = filter(Mismatches_filtering, str_detect(Mismatches_filtering$CDR3_aminoacids_Related, input$pattern))
              
             }
             
              
              
               tel = Related_READS_all_info[which(Related_READS_all_info$CDR3.IMGT.y %in% Mismatches_filtering$CDR3_aminoacids_Related), ]
              
              
              
              
              
             #####tables#####
              
               #output$relatedtodominant = renderDataTable(Related_READS_all_info)
              
              #output$aaa = renderDataTable(d)
              
              #output$bbb = renderDataTable(Mismatches_filtering)
              
              
              
              output$ccc = renderDataTable(tel)
           
              ####################
              
              ########## succed ##########
           
              observeEvent(input$Executepipeline, { 
                
                shinyalert(title = "Pipeline Executed!", type = "success") })
           
         
         
             observeEvent(input$downloadallclonotypes, {
                write.table(TheTable, "All_Clonotypes.txt", row.names = FALSE, sep = "\t", quote = FALSE) 
               
               shinyalert(title = "All_Clonotypes Downloaded", type = "success")   } )
             
            
              observeEvent(input$downloadthedominantclone, {
               write.table(Dominant_Reads_all_info, "Dominant_Clone.txt", row.names = FALSE, sep = "\t", quote = FALSE)
               
               shinyalert(title = "Dominant_Clone Downloaded", type = "success")  } )
               
             
             observeEvent(input$downloadreadsrelatedtodominant, {
               write.table(Related_READS_all_info, "Reads_related_to_the_Dominant.txt", row.names = FALSE, sep = "\t", quote = FALSE) 
               
               shinyalert(title = "Reads_related_to_the_Dominant Downloaded", type = "success")} )
             
             
            
              
              
              
              
              
              
             
           } )  
       
      
      
      
        
        
        
      
    
    } 
)
  
        
      
      
      

    
  

  
