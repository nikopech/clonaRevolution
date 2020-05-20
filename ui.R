library(shiny)

# Define UI
ui <- fluidPage(
          navbarPage(

              title = span( "clonal_R_evolution", style = "background-color: dark; color: red"),

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
                            br(),
                            fileInput("file", h3("File input")), 
                            br()
                                   ),
                            
                             
                            column(2, br(),br()),
                            column(2, br()),
                            
                            br(),
                            
                            
                           
                            
                            
                        
                            h3("SubClones"), column(3, style = "background-color:orange", 
                                   sliderInput("slider1", h3("IGHV similarity"),
                                               min = 0, max = 100, value = 50),
                                   sliderInput("slider1", h3("CDR3 similarity"),
                                               min = 0, max = 100, value = 50),
                                   
                                   
                                          
                                         
                                         
                                         
                                   ),
                            
                            column(1, style = "background-color:black", br(),br(), actionButton("calculateDominantClone","calculate"), br(),br(),br() 
                                   
                            ),
                            
                            column(3,actionButton("calculateDistances", "Calculate Distances"), 
                                     br(), br(), 
                                   actionButton("visgraph", "Visualize Graph")),
                                   
                            
                        
                                 ),
                       
                     
 
                        h3("Filters") ,   
         
                        
                        fluidRow( 
                  
                            column(3,style = "background-color:orange",
                                selectInput("select", h3("Select column"), 
                                        choices = list("Sequence.number" = 1, "Sequence.ID" = 2, "Functionality" = 3, "V.GENE.and.allele" = 4, "V.REGION.score" = 5, "V.REGION.identity.." = 6, "V.REGION.identity.nt" = 7, "V.REGION.identity....with.ins.del.events." = 8, "V.REGION.identity.nt..with.ins.del.events." = 9, "J.GENE.and.allele" = 10, "J.REGION.score" = 11, "J.REGION.identity.." = 12, "J.REGION.identity.nt" = 13, "D.GENE.and.allele" = 14, "D.REGION.reading.frame" = 15, "CDR1.IMGT.length" = 16, "CDR2.IMGT.length" = 17, "CDR3.IMGT.length" = 18, "CDR.IMGT.lengths" = 19, "FR.IMGT.lengths" = 20, "AA.JUNCTION",	"JUNCTION.frame" = 21, "Orientation" = 22, "Functionality.comment" = 23, "V.REGION.potential.ins.del" = 24, "J.GENE.and.allele.comment" = 25, "V.REGION.insertions" = 26, "V.REGION.deletions" = 27, "Sequence" = 28, "X" = 29), selected = 1),
                                 
                                textInput("text", h3("Exclude feature"), 
                                          value = " "
                                          ),
                                  br()),
                            
                            
                               
                                  
         
                            column(1, style = "background-color:black", br(),br(), actionButton("excludeButton","Exclude"), br(),br(),
                                  
                           
                                            actionButton("FilterButton","Filter"),br(),br(),
                                  
                                            actionButton("ResetButton","Reset"),br(),br(),br(),br()
                                  ),
                                   
                         
                            column(2,
                                h6("Filters Used"),
                                tableOutput("filtertable"), br() 
                                  )
                  
                        
                                ),
                       
                       br(), br(),br(),
                       
                       h3("Dominant Clone"),
                       
                       fluidRow(
                           
                            column(3, style = "background-color:orange", 
                                sliderInput("slider1", h3("IGHV identity"),
                                            min = 0, max = 100, value = 50),
                                sliderInput("slider1", h3("CDR3 identity"),
                                            min = 0, max = 100, value = 50)
                                   ),
                            
                            column(1, style = "background-color:black", br(),br(), actionButton("calculateDominantClone","calculate"), br(),br(),br() 
                  
                                   )
                    )
                ),
              
              
              tabPanel(
                
                "Filtered Sequences",
                
                value = "filtered sequences",
                
                tags$style(type="text/css", "body {padding-top: 70px;}"), #padding for navigation bar
                
                
              ),
              
              
              tabPanel(
                
                "Dominant Clone",
                
                value = "dominant clone",
                
                tags$style(type="text/css", "body {padding-top: 70px;}"), #padding for navigation bar
              ),
              
              
              tabPanel(
                
                "Subclones",
                
                value = "subclones",
                
                tags$style(type="text/css", "body {padding-top: 70px;}"), #padding for navigation bar
              )
              
              
              
          )
)




         
         
         
         
         
         
         
    



# Define server logic ----
server <- function(input, output) {

}

# Run the app ----
shinyApp(ui = ui, server = server)
