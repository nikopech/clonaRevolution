library(shiny)
library(shinyFiles)
library(shinyalert)
library(shinyWidgets)
library(shinyjs)
library(dplyr)
# Define UI
ui <- fluidPage(

  
        
  
          navbarPage(

              title = "clonalRevolution", # span( "clonal_R_evolution", style = "background-color: dark; color: #e80000"),
              
              id = "navbar",
 
              position = "fixed-top", 

              inverse=TRUE,

                  
                    tabPanel(
  
                       "Welcome", 

                        value = "welcome", 
        
                        tags$style(type="text/css", "body {padding-top: 70px;}"), #padding for navigation bar

                               
                              mainPanel(
                                  br(),
                                  br(),
                                  br(),
                                  h1("A tool for clonal evolution analysis" , align = "right"),
                                  br(),
                                  br(),
                                  h3("Approaching the dominant clone, subclones and their evolutionary interactions",style = "color:red", width = 20,  align = "right"),
                                  br(),br(),
                                  br(),br(),
                                  br(),br(),
                                  img(src = "dar.png", height = 500, width = 800, align = "right")
                                        )
                              ),

                   tabPanel(

                       "Data Analysis",

                        value = "data analysis",

                        tags$style(type="text/css", "body {padding-top: 70px;}"), #padding for navigation bar
 
                        
                       
                        fluidRow(
  
                            column(5,
                            br(),br(),br(),
                            
                            shinyDirButton('directory', 'Folder select', 'Please select a folder'),
                            
                            br(),
                            br(),
                            
                          
                            
                          
                
                             uiOutput("uiInputFiles"),
                            
                            br(),
                            
                            
                            br(),
                            useShinyjs(),
                            # tags$style(appCSS),
                            
                            
                            br(),
                            
                            uiOutput("uiLoadData"),
                            br(),
                            
                            
                           
                            
                            br(),
                            br(),
                            
                           
                           
                            
                       
                            
                            
                            ),
                            
                            br(),
                            
                             
                            
                            br(),br(),br(),br(), br(),br(),br(),br(), br()  
                                                                 
                            # h3("SubClones"), column(3, style = "background-color:#bf6f62", 
                            #        sliderInput("slider1", h3("IGHV similarity"),
                            #                    min = 0, max = 100, value = 50),
                            #        sliderInput("slider1", h3("CDR3 similarity"),
                            #                    min = 0, max = 100, value = 50),
                            #        
                            #        
                            #               
                            #              
                            #              
                            #              
                            #        ),
                            # 
                            # column(1, style = "background-color:#1f0a0a", br(),br(), actionButton("calculateDominantClone","calculate"), br(),br(),br() 
                            #        
                            # ),
                            # 
                            # column(3,actionButton("calculateDistances", "Calculate Distances"), 
                            #          br(), br(), 
                            #        actionButton("visgraph", "Visualize Graph")),
                            #        
                            # 
                            # 
                                  ),
                       
                     
 
                        h2("Filters") ,   
         
                        
                        fluidRow( 
                  
                          useShinyalert(),  
                          
                          
                            
                            
                          column(3,style = "background-color:#bf6f62",
                                 
                                 checkboxInput("checkfunctionalVgene", "Keep only reads with Functional V-Gene", value = FALSE, width = NULL),
                                 
                                 br(),
                                 
                                 checkboxInput("checkCDR3characters", "Keep only reads that CDR3 has no special characters (X,*,#,.) ", value = FALSE, width = NULL),
                                 
                                 br(),
                                 
                                 checkboxInput("checklandmarks", "Keep only reads with CDR3 valid start/end landmarks", value = FALSE, width = NULL),
                                 
                                 br(), 
                                 
                                 checkboxInput("checkproductive", "Keep only reads with productive functionality", value = FALSE, width = NULL), br()
                          ),   
                                  
         
                            column(1, style = "background-color:#1f0a0a", br(),br(), 
                                            
                                            
                           
                                            actionButton("FilterButton","Filter"),
                                   
                                            
                                   
                                            br(),br(),
                                  
                                            
                                tableOutput("filtertable"), br() 
                                  ),
                          
                          column(3,style = "background-color:#bf6f62",
                                 
                                 checkboxInput("extrafiltering", "extra user's filtering", value = FALSE, width = NULL),
                                 
                                 selectInput("select", h3("Select column"), 
                                             choices = list("Sequence.number" = 1, "Sequence.ID" = 2, "Functionality" = 3, "V.GENE.and.allele" = 4, "V.REGION.score" = 5, "V.REGION.identity.." = 6, "V.REGION.identity.nt" = 7, "V.REGION.identity....with.ins.del.events." = 8, "V.REGION.identity.nt..with.ins.del.events." = 9, "J.GENE.and.allele" = 10, "J.REGION.score" = 11, "J.REGION.identity.." = 12, "J.REGION.identity.nt" = 13, "D.GENE.and.allele" = 14, "D.REGION.reading.frame" = 15, "CDR1.IMGT.length" = 16, "CDR2.IMGT.length" = 17, "CDR3.IMGT.length" = 18, "CDR.IMGT.lengths" = 19, "FR.IMGT.lengths" = 20, "AA.JUNCTION",	"JUNCTION.frame" = 21, "Orientation" = 22, "Functionality.comment" = 23, "V.REGION.potential.ins.del" = 24, "J.GENE.and.allele.comment" = 25, "V.REGION.insertions" = 26, "V.REGION.deletions" = 27, "Sequence" = 28, "X" = 29), selected = 1),
                                 
                                 textInput("text", h3("Include String"), 
                                           value = ""
                                 ),
                                 br()) 
                          
                         
                          
                            
                            
                        
                                ),
                       
                           br(), br(),br(),
                       
                       
                           
                       
                       h2("Pipeline"),
                       
                       fluidRow(
                         
                         
                           
                            column(3, style = "background-color:#bf6f62", 
                                   selectInput("select_clonotype", "Calculate Clonotypes by :", c("V Gene + CDR3 Amino Acids","V Gene and Allele + CDR3 Amino Acids", 
                                                                                    "V Gene + CDR3 Nucleotide","V Gene and Allele + CDR3 Nucleotide",
                                                                                    "J Gene + CDR3 Amino Acids", "J Gene + CDR3 Nucleotide",
                                                                                    "J Gene and Allele + CDR3 Amino Acids", "J Gene and Allele + CDR3 Nucleotide",
                                                                                    "CDR3 Amino Acids", "CDR3 Nucleotide", "V Nt sequence + CDR3 Amino Acids", "V Nt sequence and CDR3 Nucleotide", "V.D.J Nt sequence", "V.D.J AA sequence") ), br(),
                                  
                                    
                                   ),
                            
                            column(2, style = "background-color:#bf6f62", 
                                   numericInput("numeric1", "keep clonotypes with frequency amount greater than...", value = 1, br())                                                                
                                   
                            ),
                            
                            column(2, style = "background-color:#bf6f62", 
                                   numericInput("numeric2", "mismatches allowed in Related to Dominant sequences", value = 1, br())                                                                
                                   
                            ),
                            
                             column(2, style = "background-color:#bf6f62", 
                                    checkboxInput("check", "specific pattern", value = FALSE, width = NULL),
                                    
                                    textInput("pattern", h3("motif"), 
                                              value = ""
                                   )
                                    
                             ),
                            
                            column(1, style = "background-color:#30100b", br(),br(), actionButton("Executepipeline","Execute"),br(),br() 
                  
                                   )
                    )
                ),
              
              
              tabPanel(
                
                "Loaded Data",
                
                value = "filtered sequences",
                
                tags$style(type="text/css", "body {padding-top: 70px;}"), #padding for navigation bar
                                 
                mainPanel(
                    
                         
                        h2("Raw Data"),br(),br(),
                        
                        dataTableOutput("table"), br(),
                        
                        numericInput("number", "print dataframe of sample..", value = 1, width = 100), br(), 
                        
                        h2("Filtered Data"),br(),br(),
                        
                        dataTableOutput("filtered"), br(), 
                        
                         
                               
                               numericInput("number2", "print dataframe of sample..", value = 1, width = 100)
                                
                  
                          )  
                
                       ),
              
              
              tabPanel(
                
                "Clonotypes",
                
                value = "dominant clone",
                
                tags$style(type="text/css", "body {padding-top: 70px;}"), #padding for navigation bar
                
                h2("All Clonotypes"),br(),br(),
                
                dataTableOutput("clonoTable"),br(),
                
                numericInput("nameid3", "print dataframe of sample..", value = 1, width = 100), br(),
                
                actionButton("downloadallclonotypes", "Download"), br(),br(),
                
                h2("Dominant Clone"),br(),br(),
              
               dataTableOutput("thedominant"), br(),
               
               textOutput("textdominant"),
               
               actionButton("downloadthedominantclone", "Download")
               
                ),
              
              
              
              
              tabPanel(
                
                "Related to Dominant",
                
                value = "reads_related",
                
                tags$style(type="text/css", "body {padding-top: 70px;}"), #padding for navigation bar
             
                h2("Related to Dominant Clone Clonotypes"),br(),br(),
                
                 h3("(from all samples)"), br(),br(),
                
              
                
                dataTableOutput("relatedclonotypes"), br(),
                
                
                
                
                
                
                h3("(from 1st sample)"), br(),br(),
                
                dataTableOutput("related1"),br(),
                
                textOutput("meanF1"),br(),
                
                
                
                
                h3("(from 2nd Sample)"), br(),br(),
                
                dataTableOutput("related2"),br(),
                
                textOutput("meanF2"),br(),
                
                h2("Relative Frequency of Subclones (between samples)"),
                
                dataTableOutput("freq"),br(),
                
               
                
                
                dataTableOutput("ccc"),
                dataTableOutput("ddd"),
                
                
                actionButton("downloadreadsrelatedtodominant", "Download")
                
               
                
                 ),
              
              tabPanel(
                
                "Visualization",
                
                value = "subclones",
                
                tags$style(type="text/css", "body {padding-top: 70px;}"), #padding for navigation bar
                
                
               h2("1st Sample Graph"),
                
                 plotOutput(
                  "igraph1",
                  width = "170%",
                  height = "750px",
                  click = NULL,
                  dblclick = NULL,
                  hover = NULL,
                  brush = NULL,
                  inline = FALSE
                ),
                
                h3("Graph Clarification / 1st Sample"),br(),
                
                
                
                dataTableOutput("graph1"), 
                
                br(), 
               
               h2("2nd Sample Graph"),
                
                plotOutput(
                  "igraph2",
                  width = "170%",
                  height = "750px",
                  click = NULL,
                  dblclick = NULL,
                  hover = NULL,
                  brush = NULL,
                  inline = FALSE
                ),
                
                h3("Graph Clarification / 2nd Sample"),br(),
                
                
                
                dataTableOutput("graph2"), br(), 
               
               h2("Subclones Frequency Dotplot"), br(),
                
                plotOutput("dotplot", width = "170%",
                           height = "750px"), br(), 
               
               h3("Dotplot Clarification"), br(),
               
               dataTableOutput("dottable"), br(), 
               
              
                h2("Dominant Clones Frequency Dotplot"), br(),
               
               plotOutput("domdotplot"), br(),
               
               h3("Dotplot Clarification"), br(),
               
               dataTableOutput("domdottable")
               
               
                
                
                
              )
              
              
              
              
              
          )
)




         
         
         
         
         
         
         
    



# Define server logic ----
# server <- function(input, output) {
#   
# }

# Run the app ----
# shinyApp(ui = ui, option = list(height = 1388))