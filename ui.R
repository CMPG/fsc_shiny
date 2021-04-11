library(shiny)
library(shinythemes)
library(RColorBrewer)

shinyUI(
  navbarPage("fastsimcoal web app 1.1.1",
             
  ##### DOCUMENTATION TAB ------------------------------------------------------
  
  tabPanel("Home",
          fluidRow(
            column(width = 3),
            column(width = 6,
                   
            h4('Welcome to the fastsimcoal (FSC) Shiny app'),
            
            p('This application aims at helping fastsimcoal users writing
               the parameter file and interpreting the output of FSC.'),
            p('FSC takes as input a parameter file (.par) and the 
              observed SFS (.obs) and outputs various files,
              including one or several files containing the expected SFS.'),
             
            h4('Input of the application'),
             
            p('Different files can be analysed with this app:'),
              tags$ul(
                tags$li("the parameter file (.par);"), 
                tags$li("observed and expected SFS (depending on the number of populations,
                       one could consider marginal or joint SFS);"), 
                tags$li("likelihoods of different parameters values
                       (for single or multiple runs).")
              ),
  
             p('A single par file can be uploaded to plot the demographic model.
               All other files should be uploaded as a single ZIP archive.'),
             
             fileInput('mainfile', 'Upload a .par or .zip file:',
                       accept = c('.par', '.zip')
             ),
             
             h4('Files uploaded:'),
             textOutput('files.list')
            ),
            column(width = 3))
  ),
             
  ##### PAR FILE TAB -----------------------------------------------------------
  
    
  tabPanel("Model plot",
    fluidPage(theme = shinytheme("sandstone"),
      tags$head(tags$style(type="text/css", ".container-fluid {max-width: 1300px}")),
            
      # Sidebar
      sidebarPanel(
        h4("Graphical parameters"),
        sliderInput("widthMP",
                    "Plot width", 40, 100, 100,
                    ticks = FALSE, post = "%"),
        sliderInput("heightMP",
                    "Plot height", 300, 1000, 600,
                    ticks = FALSE, post = "px"),

        numericInput("gentimeMP", "Generation time", value = 1, min = 1, max = NA, step = 1),
        
        sliderInput("yL", "Zoom",
                    min = 0, max = 1, value = c(0, 1), step = 0.05),
        
        checkboxInput("logMP", "Log scale (time)", value = FALSE),
        
        checkboxInput("backMP", "Arrows go backwards", value = TRUE),
        checkboxInput("arrowsMP", "Represent splits with arrows", value = TRUE),
        checkboxInput("growthMP", "Indicate growth changes", value = TRUE),
        checkboxInput("migMP", "Plot migration events", value = TRUE),
        checkboxInput("ibMP", "Plot instantaneous bottlenecks", value = TRUE)
      ),

      # Par file plot
      
      mainPanel(
        h2("Demographic model"),
        
        p("Visualise the simulated model."),
        
        conditionalPanel(condition = "!output.parUp",
          p("Please upload a par file.")
        ),
        
        conditionalPanel(condition = "output.parUp",
          uiOutput("plotPar")
        )
      )
    )
  ),


  ##### SFS PLOT TAB -----------------------------------------------------------
  
  tabPanel("SFS",
           
     # Sidebar
     sidebarPanel(
       h4("Graphical parameters"),
       sliderInput("widthSFS",
                   "Plot width",40,100,100,
                   ticks = FALSE, post = "%"),
       sliderInput("heightSFS",
                   "Plot height",300,1000,600,
                   ticks = FALSE, post = "px"),
       sliderInput("cexAxisSFS",
                   "Font size",0.5, 2.5, 1.8,
                   ticks = FALSE),
       checkboxInput(inputId = "logSFS", label = "Log scale", value = TRUE),
       numericInput("numcols", "Number of columns", value = 2, min = 1, max = 3)
     ),

     # SFS plots
     mainPanel(
       h2("Marginal SFS"),
       
       conditionalPanel(condition = "!output.sfsUp",
                        p("Please upload a SFS.")
       ),
       conditionalPanel(condition = "output.sfsUp",
                        p("Select population"),
                        
                        h3("Observed vs. expected SFS"),
                        
                        uiOutput("plotSFSui"),
                        
                        h3("Scaled observed vs. expected SFS"),
                        
                        uiOutput("plotScaledSFSui"),
                        
                        h3("Entries contributing the most to differences in likelihood"),
                        
                        uiOutput("plotWorstEntriesui")
       ),

       h2("Joint SFS"),
       # conditional panel: if "joint" files, choose the joint SFS to plot
       conditionalPanel(condition = "!output.jsfsUp",
          conditionalPanel(condition = "output.sfsUp",
            uiOutput("jsfsPopus"),
            sliderInput("freqCutoff","Threshold for text display:",
                        min = -5, max = 0, step = 0.5, value = -3),
            selectInput("jsfscolor", "Pick a color palette:",
                        c("Red", rownames(brewer.pal.info[brewer.pal.info$category != "qual", ]))),
            uiOutput("plot2DSFSfrom1Dui")
          ),
          conditionalPanel(condition = "!output.sfsUp",
            p("Please upload a joint SFS.")
          )
       ),
       conditionalPanel(condition = "output.jsfsUp",
                        uiOutput("jsfsSelect"),
                        sliderInput("freqCutoff","Threshold for text display:",
                                    min = -5, max = 0, step = 0.5, value = -3),
                        selectInput("jsfscolor", "Pick a color palette:",
                                    c("Red", rownames(brewer.pal.info[brewer.pal.info$category != "qual", ]))),
                        uiOutput("plot2DSFSui")
       )
       
     )
  ),
  
  ##### OUTPUT TAB -------------------------------------------------------------
  
  tabPanel("Likelihoods",
           
           fluidPage(
             # Sidebar
             sidebarPanel(
               fileInput("multiRun", "Upload fastsimcoal bestLhoods",
                         accept = c(".txt")
               ),
               h4("Graphical parameters"),
               sliderInput("widthMR",
                           "Plot width",40,100,100,
                           ticks=FALSE,post="%"),
               sliderInput("heightMR",
                           "Plot height",500,1000,800,
                           ticks=FALSE,post="px"),
               sliderInput("cexAxisMR",
                           "Font size",1,3,2,step=0.1,
                           ticks=FALSE),
               sliderInput("yLlhood", "Zoom",
                           min = 0, max = 1, value = c(0, 1), step = 0.01),
               checkboxInput("logMR", "Log scale", value = FALSE),
               numericInput("selecMR", "Run number (0 = maximum likelihood):",
                            0, 1, 50, 1)
               
             ),
             
             # Par file plot
             mainPanel(
               
               h1("Likelihoods"),
               
               p("Upload the concatenated best likelihoods of different runs to
                 compare them."),
               
               uiOutput("varsMR"),
               uiOutput("plotbestLhoods")
               )
           )
  )
  
  )
)



