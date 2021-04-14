######################################
#### Define Modalsfor Dynamic UI #####
######################################

inputModal <- function(failed=FALSE){
  modalDialog(class = 'text-center',
    div(style='display: inline-block;vertical-align:top;',fileInput("background_file", label = strong(HTML("Upload a background file for pathway enrichment")))),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'BackgroundFileInfo', label=strong('i'),style='pill',color='default')),
    div(style='display: inline-block;vertical-align:top;',fileInput("ppi_file", label = strong(HTML("Upload a PPI database")))),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'PPIFileInfo', label=strong('i'),style='pill',color='default')),
    div(style='display: inline-block;vertical-align:top;',fileInput("annotation_file", label = strong(HTML("Upload an annotation file for pathway/ontological enrichment")))),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'AnnotationFileInfo', label=strong('i'),style='pill',color='default')),
    div(style='display: inline-block;vertical-align:top;',prettyCheckboxGroup('background', label=strong(HTML('Select the background to be used for pathway/ontological enrichment',)),
                choices=list('In-House Human/Mouse Background (Default)'='In-House Human/Mouse Background (Default)',
                             'User-Supplied Background'='User-Supplied Background',
                             'Merge In-House Background With User-Supplied Background'='Merge In-House Background With User-Supplied Background'),
                selected = c('In-House Human/Mouse Background (Default)'),
                status = 'primary')),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'BackgroundSelectInfo', label=strong('i'),style='pill',color='default')),
    div(style='display: inline-block;vertical-align:top;',prettyCheckboxGroup('PathwayEnrichment', label=strong(HTML('Select databases used for pathway/ontological enrichment')),
                choices=list('KEGG'='KEGG',
                             'GO'='GO',
                             'Reactome'='Reactome',
                             'HALLMARK'='HALLMARK',
                             'User-Supplied Database in .xlsx Format'),
                selected=c('KEGG', 'GO','Reactome','HALLMARK'),
                status='primary')),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'AnnotationSelectInfo', label=strong('i'),style='pill',color='default')),
  )
}

inputModal3 <- function(failed=FALSE){
  modalDialog(
    h4(strong('DISCLAIMER:',
    'Thank you for choosing JUMPn for your network and pathway analysis.',
    'JUMPn combines Co-Expression clustering with protein-protein interaction data',
    'to readily derive biological insights from complex data. JUMPn users have the choice of',
    'Executing integrated analysis (WGCNA+PPI Analysis; Default) or independent',
    'analysis (Either WGCNA or PPI Analysis separately). Before beginning, please choose between',
    'the available options for analysis.')),
    radioButtons('ModeOfAnalysis', label = h5(strong('Select a Mode of Analysis for JUMPn to Execute')),
                 choices=list('WGCNA+PPI'='WGCNA+PPI',
                              'Independent WGCNA'='Independent WGCNA',
                              'Independent PPI'='Independent PPI')),
    footer = actionButton('MODA_Begin', 'Begin')
  )
}
testSweetAlert <- function(session){
  sendSweetAlert(
    session=session,
    title='Information',
    type='info',
    btn_labels = NA,
    text = tags$div(class='text-center','Thank you for choosing JUMPn for your network and pathway analysis.JUMPn combines Co-Expression clustering with protein-protein interaction data
            to readily derive biological insights from complex data. JUMPn users have the choice of executing integrated analysis (WGCNA+PPI Analysis; Default) or independent
            analysis (Either WGCNA or PPI Analysis separately). Before beginning, please choose between the available options for analysis.',
            awesomeRadio('ModeOfAnalysis', label = h5(strong('Select a Mode of Analysis for JUMPn to Execute')),
                          choices=list('WGCNA+PPI'='WGCNA+PPI',
                                        'Independent WGCNA'='Independent WGCNA',
                                        'Independent PPI'='Independent PPI'),
                         width='100%'),
            div(style='display:inline-block;vertical-align:top;',actionButton('MODA_Begin', 'Begin'))
            
    )
    
  )
}

messageModal <- function(type,session){
  sendSweetAlert(
    session=session,
    type='info',
    title='Information',
    btn_labels = NA,
    text = tags$div(
      h4(strong(sprintf('Great, you can begin the %s mode of analysis. Check the proper file format with the "View Example File" button below', type))),
      actionButton('ViewExampleFile', 'View Example File'),
      actionButton('Continue', 'Continue')
      
    )
    
  )
  
}

imageModal <- function(output){
  modalDialog(
    h4(strong("DISCLAIMER: Gene List (Gene Symbols) must be in Column 2. Expression data must be in Column 4 onward; you can use as many expression columns as desired. 
    Columns 1 and 3 can include whatever information you'd like. Your file must have a header row, and you can use whatever labels you'd like.")),
    renderImage({
      list(src='www/images/WGCNAInput.png',
            contentType = 'image/png',
            width = 750,
            height = 600)
    }),
    easyClose = TRUE
  )
  
}

imageModal2 <- function(output){
  modalDialog(
    h4(strong("DISCLAIMER: PPI Analysis requires a different input than is needed for WGCNA. Column 1 must have the Cluster/Group to which the gene belongs. ",
              "Column 2 must have the accession ID (or some other misscelaneous categorical information), and Column 3 must contain the gene list (Gene Symbols),",
              'If you are conducting independent PPI analysis after running WGCNA, then this 3-column file is produced automatically, which you  can',
              'use directly as input for independent PPI after downloading the WGCNA data. Otherwise, you can make your own custom file which adheres to the described format.',
              'Additional columns may be provided beyond the standard 3-column format,  but these columns will not be evaluated.')),
    renderImage({
      list(src='www/images/PPIInput.png',
           contentType = 'image/png',
           width = 200,
           height = 450)
    }),
    easyClose = TRUE
  )
  
}
programFinishedModal <- function(session){
  sendSweetAlert(
    session = session,
    title='Success !!',
    type='success',
    btn_labels = NA,
    text = tags$div('Results have been  generated. Click the button below to be taken to the first results page!',
                    actionButton('TakeMeThere', 'Continue to Results'))
  )
}

ArrivedToResultsModal <- function(whereami,session){
  if (whereami=='WGCNA+PPI'){
    sendSweetAlert(
      session=session,
      type='info',
      title='Information',
      btn_labels = NA,
      text = tags$div(
        'Here are your WGCNA results. Below you will find co-expression trends (in the form of boxplot and linear figures), as well as pathway analysis (if selected). 
          Once you have finished investigating your results, you can navigate to the PPI results via the Page 2 tab above.',
        actionButton('ViewResults', 'View Results')
      )
    )
  } else if (whereami=='Independent WGCNA'){
    sendSweetAlert(
      session=session,
      type='info',
      title='Information',
      btn_labels = NA,
      text = tags$div(
        'Here are your WGCNA results. Below you will find co-expression trends (in the form of boxplot and linear figures), as well as pathway analysis (if selected). 
          Once you have finished investigating your results, you can submit another search  by navigating back to the  Commence Analysis Page. If you wish to autamatically submit
          a PPI search with these results generated from WGCNA, click the  button below',
        actionButton('ViewResults', 'View Results'),
        actionButton('SubmitQuickPPI', 'Start a PPI Search With WGCNA Results')
      )
    )
  } else if (whereami=='Independent PPI'){
    sendSweetAlert(
      session=session,
      type='info',
      title='Information',
      btn_labels = NA,
      text = tags$div(
        'Here are your PPI results. Below you will find co-expression trends (in the form of boxplot and linear figures), as well as pathway analysis (if selected). 
        Once you have finished investigating your results, you can submit another search  by navigating back to the  Commence Analysis Page.'
      )
    )
  }
}

QuickPPIModal <- function(){
  fluidRow(
    tags$style(HTML('#btn1 {
                        background-color:red
                    } #btn2 {
                        background-color:red
                    } #Advanced {
                        background-color:lightblue
                    } #PPISearch {
                        background-color:lightblue
                    } #WGCNASearch {
                        background-color:lightblue
                    } #do {
                        background-color:red
                    
                    
    ')),
    column(4),
    column(4, 
           h4(strong(HTML('PPI Analysis Parameters'))),
           box(status='primary', background = 'blue', width=40, height=10,
               id='second',
               fileInput("input_table", label = strong(HTML("Upload Input File"))),
               selectInput("skip_PPI_network_analysis", h5(strong("Should Co-Expression Clusters Be Modularized by Protein Interactions?")),
                           choices=list('Yes'='Yes',
                                        'No'='No')),
               sliderInput("TOM_triggered_module_size", label=h5(strong('Specify How Large a Module Must Be to Trigger TOM')),
                           min=1, max=200, value=100, step=1, ticks=TRUE),
               div(style='display:inline-block;vertical-align:top;',actionButton('Advanced', 'Advanced Parameters')),
               div(style='display:inline-block;vertical-align:top;',actionButton('do', label=h4(strong('Submit Independent PPI Analysis')))))
           
    ),
    column(4)
  )
}
ConfirmSearchModal <- function(session){
  sendSweetAlert(
    session=session,
    type = 'warning',
    title='Warning:',
    btn_labels = NA,
    text = tags$div('You are about to submit a search. Is your file of the correct format and are your parameters set. If so click "Confirm Search". Otherwise press "Back" and double check your parameters.',
                    actionButton("ConfirmSearch", "Confirm Search"),
                    actionButton('BackButton', 'Back'))
  )
}

ConfirmNewSearchModal <- function(session){
  sendSweetAlert(
    session=session,
    type = 'warning',
    title='Warning:',
    btn_labels = NA,
    text = tags$div('Submiting a new search will cause data from the previous run to be  overwritten. Wee highly recommend, if you have not already, to download the data generated with the
                    download handler below. Otherwise, you can confirm your new search with the Continue with New Search button',
                    actionButton("ConfirmNewSearch", "Confirm New Search"),
                    actionButton('DownloadData', 'Download Data'),
                    actionButton('CancelNewSearch', 'Cancel'))
  )
}

WidgetInfoModal <- function(session,messageinfo){
  sendSweetAlert(
    session=session,
    title='Information',
    text=messageinfo,
    type='info',
  )
}
  
successful_download <- function(session){
  sendSweetAlert(
    session=session,
    title = 'Success',
    text='Congrats, publciation data has been downloaded!',
    type='success'
  )
}


  