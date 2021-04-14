#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(shinyWidgets)
require(visNetwork, quietly = TRUE)
options(shiny.maxRequestSize = 30*1024^2)
R_BUILD_TAR=tar
# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("JUMPn: Comprehensive Proteomic Network Analysis"),
    theme = "style/style.css",
    fluid = TRUE,
    collapsible=TRUE,
    tags$style(HTML("
      #first {
          border: 4px double red;
      }
      #second {
          border: 4px double black;
      }
    ")),
    
    tabsetPanel(
        id='whichtab',
        tabPanel("Home",
                 value='Home',
                 includeHTML("JUMPn_Helpers/html_files/home.html"),
                 tags$script(src = "plugins/scripts.js"),
                 tags$head(
                     tags$link(rel = "stylesheet", 
                               type = "text/css", 
                               href = "plugins/font-awesome-4.7.0/css/font-awesome.min.css")
                 )
        ),
        # Sidebar with a slider input for number of bins 
        tabPanel(
            "Commence Analysis", 
            value='Commence',
            includeHTML("JUMPn_Helpers/html_files/network_analysis.html"),
            uiOutput('DynamicUI')
           
        ),
        tabPanel(
            title=uiOutput('WGCNATitlePage'),
            value='WGCNAPage',
            includeHTML("JUMPn_Helpers/html_files/net_output.html"),
            useShinyjs(),
            inlineCSS(list(.networkplotborder1 = "border: 2px solid black")),
            inlineCSS(list(.networkplotborder2 = "border: 2px solid black")),
            mainPanel(
                div(style='display:inline-block;vertical-align:top;',dropdownButton(icon=icon('download'),
                                                                                    switchInput(inputId='TotalData',
                                                                                                label = 'Download All Data',
                                                                                                labelWidth = '50px'),
                                                                                    switchInput(inputId='PublicationData',
                                                                                                label = 'Download Publication Tables',
                                                                                                labelWidth = '50px'),
                                                                                      switchInput(inputId='Figure1a',
                                                                                              label = 'Download Fig. 1',
                                                                                              labelWidth = '50px'),
                                                                                    switchInput(inputId='Figure2a',
                                                                                                label = 'Download Fig. 2',
                                                                                                labelWidth = '50px'),
                                                                                    switchInput(inputId='Figure3a',
                                                                                                label = 'Download Fig. 3',
                                                                                                labelWidth = '50px'),
                                                                                    switchInput(inputId='Figure4a',
                                                                                                label = 'Download Fig. 4',
                                                                                                labelWidth = '50px'),
                                                                                    switchInput(inputId='Figure5a',
                                                                                                label = 'Download Fig. 5',
                                                                                                labelWidth = '50px'),
                                                                                    downloadButton('downloadSelectedData1', 'Download Selected Figures',
                                                                                                 style='background-color:orange;')
                                                                                    )),
                div(style='display:inline-block;vertical-align:top;',
                    actionButton('SubmitNewSearch', 'Submit a New Search', style='position: absolute; right:10px;background-color:orange')),
                selectInput('ExpressionDisplay', h5("Select the Expression Format"),
                            choices=list("Boxplot" = 'Boxplot',
                                         "Trends" = "Trends"),
                            selected = "Boxplot"),
                
                column(6,
                       uiOutput('ExpressionDisplayObject')),
                column(6,
                       uiOutput('ClusterPathwayHeatmap')),
                
                imageOutput('buffeerex'),
                imageOutput('buffeerey'),
                
                uiOutput('ClusterDTSelectionBox'),
                uiOutput('DynamicClusterDTDisplay'),
                fluidRow(
                  column(6,
                         plotOutput('GeneExpression')),
                  column(6,
                         plotOutput('Cluster_Expression'))
                ),
                htmlOutput('ClusterPathwayTableTitle'),
                dataTableOutput('ClusterPathwayTable'),
                
                
            )
        ),
        tabPanel(
            title=uiOutput('PPITitlePage'),
            value='PPIPage',
            includeHTML("JUMPn_Helpers/html_files/path_output.html"),
            mainPanel(
              fluidRow(
                column(8,
                  uiOutput('DynamicClusterSelectionBox')),
                column(4,
                  div(style='display:inline-block;vertical-align:top;',selectInput('ClusterExpressionDisplay',
                              h5("Select the Expression Format"),
                              choices=list("Boxplot" = 'Boxplot',
                                           "Trends" = "Trends",
                                           "Pathway Barplot"="Pathway Barplot",
                                           "Pathway Circle Plot"="Pathway Circle Plot"),
                              selected = "Boxplot",
                              width = '50%')),
                  div(style='display:inline-block;vertical-align:top;',
                      actionButton('SubmitNewSearch', 'Submit a New Search', style='position: absolute; right:10px;background-color:orange'))
                )
              ),
              fluidRow(
                column(8,
                       div(style='display:inline-block;vertical-align:top;',uiOutput('DynamicClusterDisplay'))),
                                                                                         
                column(4,
                       uiOutput('ExpressionFormat'))
              ),
              #plotOutput('groupedclusterbars'),
              fluidRow(
                column(8,
                  uiOutput('DynamicModuleSelectionBox')),
                column(4,
                  uiOutput('ModulePathwaySelection'))
              ),
              fluidRow(
                column(8,
                       uiOutput('DynamicModularDisplay')),
                column(4,
                       uiOutput('EMAP'))
              )
            )
        ),
        #tabPanel(
         #   "Literature Search (LitFinder)",
          #  mainPanel(
           #     width = 9,
            #    #downloadButton("downloadData", "Download Data"),
             #   textOutput("inputFile1"),               #Returns the name of the input file the user has chosen to upload
              #  textOutput("Item_List"),                #Returns the list of gene names that have been searched for
               # textOutput("list"),                     
                #dataTableOutput('Table'),                      #Creates a data table which represents the most cited papers for each gene
                #plotOutput("bar",height = 500),         #Creates a plot which depicts the number of papers existing for each protein
                #plotOutput("newbar",height = 600),       #Creates a plot depicting the number of citations for each paper
                #visNetworkOutput("LitNetwork"),
                #dataTableOutput("top_papers"),
                #visNetworkOutput("LitNet1"),
                #dataTableOutput("LitTable1"),
                #visNetworkOutput("LitNet2"),
                #dataTableOutput("LitTable2"),
                #visNetworkOutput("LitNet3"),
                #dataTableOutput("LitTable3"),
                #visNetworkOutput("LitNet4"),
                #dataTableOutput("LitTable4"),
                #visNetworkOutput("LitNet5"),
                #dataTableOutput("LitTable5")
                
            #),
            #sidebarPanel(
             #   width=3,
              #  uiOutput("CondPanelSelection"),
                #selectInput('search_object', h4("Select the Cluster or Module You Would Like to Search")),
                #selectInput("select", h4("Select Organism"),                   #The selectInput() function creates a selection widget wherein the user selects the organism from which their data derives
                 #           choices = list("Human" = "Human", "Mouse" = "Mouse" ), selected = "Human"), #The selection widget defaults to 'Human' using the 'selected' attribute
               # textInput("Keyword", h4("Keyword Filter")),                    #The second widget utilizes the textInput function which allows the user to specify a keyword they intend to filter for
              #  textInput("Year", h4("Year Filter")),                          #The third widget utilizes the textInput function which allows the user to specify a year they intend to filter for
               # textInput("Journal", h4("Journal Filter")),                    #The fourth widget utilizes the textInput function which allows the user to specify a journal they intend to filter for
                #textInput("Citation", h4("The number of top most cited papers returned (A useful start is five or ten)"),                     #The fifth widget utilizes the textInput function which allows the user to specify the number of top cited papers to be returned
                 #         value = 5),
                #textInput("OutputFile", h4("Output File"),                     #The sixth widget utilizes the textInput function which allows the user to specify the name of the file to be exported after the program has finished
                 #         value = "myexamplefile"),        #Directions are provided for the outputFile function which is created using the 'value' attribute
                #selectInput('filestype',
                 #           h4("Choose File Upload Type"),
                  #          choices=list('csv Extension File' = 'a',
                   #                      'txt Extension File' = 'b',
                    #                     'Manual Text Upload' = 'c')),
                #uiOutput("CondPanel1"),
                #uiOutput("CondPanel2"),
                #actionButton('do2', 'Search Now'),
                #textOutput("Action")
            #)
        #),
        tabPanel("Help",
                 includeHTML("JUMPn_Helpers/html_files/help.html"),
                 tags$head(
                     tags$link(rel = "stylesheet", 
                               type = "text/css", 
                               href = "plugins/font-awesome-4.7.0/css/font-awesome.min.css")
                 )
        )
        
    )
)
