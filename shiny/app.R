library(plotly)
library(shiny)
library(DT)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(data.table)
library(stringr)

# Weird hack as R sometimes "forgets" its working directory
wd <- setwd(".")
setwd(wd)

# Make an empty Google analytics file (for dev version - not for production)
if (!file.exists('google-analytics.js'))
{
  file.create('google-analytics.js')
}

filename <- 'antigen.tsv'
filename <- normalizePath(filename)

fileInfo <- file.info(filename)
modifiedDate <- strsplit(as.character(fileInfo$mtime), ' ')[[1]][1]


sentences <- fread(filename,sep='\t',header=T,stringsAsFactors=TRUE)
sentences <- sentences[sentences$pmid!='None',]
sentences$pmid <- as.integer(as.character(sentences$pmid))


# Fill in some details for the citation table
sentences$pmidLink <- paste("<a target=\"_blank\" href='https://www.ncbi.nlm.nih.gov/pubmed/", sentences$pmid, "'>", sentences$pmid, "</a>", sep='')
sentences$journalShort <- strtrim(sentences$journal,51)
sentences[str_length(sentences$journalShort)==51,'journalShort'] <- paste(strtrim(sentences[str_length(sentences$journalShort)==51,]$journalShort,50),'...',sep='')
sentences$journalShort <- factor(sentences$journalShort)
sentences$probability <- round(sentences$probability,3)

sentences <- sentences[order(sentences$probability,decreasing=T),]

collated_normalized <- plyr::count(sentences[,c('normalizedTerm'),drop=F])
colnames(collated_normalized) <- c('normalizedTerm','sentenceCount')
collated_normalized <- collated_normalized[order(collated_normalized$normalizedTerm),]
collated_normalized <- collated_normalized[order(collated_normalized$sentenceCount,decreasing=T),]
#collated_normalized$normalizedTerm <- factor(as.character(collated_normalized$normalizedTerm),levels=unique(collated_normalized$normalizedTerm))
rownames(collated_normalized) <- 1:nrow(collated_normalized)

collated_intext <- plyr::count(sentences[,c('termInText'),drop=F])
colnames(collated_intext) <- c('termInText','sentenceCount')
collated_intext <- collated_intext[order(collated_intext$termInText),]
collated_intext <- collated_intext[order(collated_intext$sentenceCount,decreasing=T),]
#collated_intext$normalizedTerm <- factor(as.character(collated_intext$termInText),levels=unique(collated_intext$termInText))
rownames(collated_intext) <- 1:nrow(collated_intext)

#sentences$normalizedTerm <- factor(as.character(sentences$normalizedTerm),levels=unique(collated$normalizedTerm))


colors <- brewer.pal(3,'Set2')
color_Driver <- colors[1]
color_TumorSuppressor <- colors[2]
color_Oncogene <- colors[3]

citationTableExplanation <- "<br /><br /><br /><b>Citation Table:</b><br />Select row in table above to see associated citations and sentences<br /><br />"

ui <- function(req) {
  fluidPage(
    tags$head(
      includeHTML("google-analytics.js"),
      tags$style(".rightAlign{float:right; margin-left:5px; margin-bottom: 20px;}")
      ),
    titlePanel("",windowTitle="TumorAntigens"),
    helpText(includeHTML("header.html")),
    tabsetPanel(type = "tabs",
                tabPanel("By Normalized Term", 
                         mainPanel(
                           plotlyOutput("barchart_normalized"),
                           downloadButton("download_table_collated_normalized", label = "Download", class='rightAlign'),
                           DT::dataTableOutput("table_collated_normalized"),
                           HTML(citationTableExplanation),
                           downloadButton("download_table_sentences_normalized_all", label = "Download All", class='rightAlign'),
                           downloadButton("download_table_sentences_normalized_shown", label = "Download Shown", class='rightAlign'),
                           DT::dataTableOutput("table_sentences_normalized"),
                           helpText(paste("Last updated on",modifiedDate))
                         )
                )
                ,
                tabPanel("By Term in Text", 
                         mainPanel(
                           plotlyOutput("barchart_intext"),
                           downloadButton("download_table_collated_intext", label = "Download", class='rightAlign'),
                           DT::dataTableOutput("table_collated_intext"),
                           HTML(citationTableExplanation),
                           downloadButton("download_table_sentences_normalized_all2", label = "Download All", class='rightAlign'),
                           downloadButton("download_table_sentences_normalized_shown2", label = "Download Shown", class='rightAlign'),
                           DT::dataTableOutput("table_sentences_intext"),
                           helpText(paste("Last updated on",modifiedDate))
                         )
                )
    )
  )
  
}
row <- data.frame(normalizedTerm="ABO",stringsAsFactors=F)

server <- function(input, output, session) {
  output$table_collated_normalized <- DT::renderDataTable({
    DT::datatable(collated_normalized[,c('normalizedTerm','sentenceCount')],
                  selection = 'single',
                  rownames = FALSE,
                  colnames=c('Normalized Term','Sentence #'),
                  options = list(lengthMenu = c(5, 30, 50), pageLength = 20))
  })
  
  tableProxy_normalized = dataTableProxy('table_collated_normalized')
  
  
  output$table_sentences_normalized <- DT::renderDataTable({
    if(length(input$table_collated_normalized_rows_selected)>0) {
      row <- collated_normalized[input$table_collated_normalized_rows_selected,]
      entries <- sentences[sentences$normalizedTerm==row$normalizedTerm,]
    } else {
      entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
      colnames(entries) <- colnames(sentences)
    }
    DT::datatable(entries[,c('pmidLink','journalShort','year','HUGOorUniprotID','termInText','probability','sentence'),],
                  selection = 'none',
                  rownames = FALSE,
                  colnames=c('PMID','Journal','Year', 'HUGO/UniProt ID', 'Term in Text', 'Probability', 'Sentence'),
                  escape = FALSE,
                  options = list(pageLength = 20))
  })
  
  
  output$barchart_normalized <- renderPlotly({
    table <- collated_normalized
    if (nrow(table) > 0) {
      table$normalizedTerm <- factor(as.character(table$normalizedTerm),levels=unique(table$normalizedTerm))
      
      p <- plot_ly(table, x=~normalizedTerm, y=~sentenceCount, source='barchart_normalized', type = 'bar', name = 'Antigen', marker=list(color=color_TumorSuppressor)) %>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack', margin = list(b = 200), xaxis=list(title = "", tickangle = 45))%>% 
        config(displayModeBar = F)
      
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  observe({
    s <- event_data("plotly_click", source = "barchart_normalized")
    if (length(s) > 0) {
      table <- collated_normalized
      rowNumber <- rownames(table[table$normalizedTerm==s$x,])
      tableProxy_normalized %>% selectRows(as.numeric(rowNumber))
    }
  })
  
  output$download_table_collated_normalized <- downloadHandler(
    filename = function() {
      return("antigens_collatedByNormalized.tsv")
    },
    content = function(file) {
      write.table(collated_normalized, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_table_sentences_normalized_all <- downloadHandler(
    filename = function() {
      return("antigens_sentences.tsv")
    },
    content = function(file) {
      outdata <- sentences[,c("HUGOorUniprotID","normalizedTerm","termInText","pmid","journal","year","title","probability","sentence")]
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_table_sentences_normalized_shown <- downloadHandler(
    filename = function() {
      return("antigens_sentences_selected.tsv")
    },
    content = function(file) {
      if(length(input$table_collated_normalized_rows_selected)>0) {
        row <- collated_normalized[input$table_collated_normalized_rows_selected,]
        entries <- sentences[sentences$normalizedTerm==row$normalizedTerm,]
      } else {
        entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
        colnames(entries) <- colnames(sentences)
      }
      
      outdata <- entries[,c("HUGOorUniprotID","normalizedTerm","termInText","pmid","journal","year","title","probability","sentence")]
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  
  
  
  
  
  
  
  output$table_collated_intext <- DT::renderDataTable({
    DT::datatable(collated_intext[,c('termInText','sentenceCount')],
                  selection = 'single',
                  rownames = FALSE,
                  colnames=c('Term in Text','Sentence #'),
                  options = list(lengthMenu = c(5, 30, 50), pageLength = 20))
  })
  
  tableProxy_intext = dataTableProxy('table_collated_intext')
  
  
  output$table_sentences_intext <- DT::renderDataTable({
    if(length(input$table_collated_intext_rows_selected)>0) {
      row <- collated_intext[input$table_collated_intext_rows_selected,]
      entries <- sentences[sentences$termInText==row$termInText,]
      #entries <- sentences[sentences$matching_id==row$matching_id]
    } else {
      entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
      colnames(entries) <- colnames(sentences)
    }
    DT::datatable(entries[,c('pmidLink','journalShort','year','HUGOorUniprotID','normalizedTerm','probability','sentence'),],
                  selection = 'none',
                  rownames = FALSE,
                  colnames=c('PMID','Journal','Year', 'HUGO/UniProt ID', 'Normalized Term', 'Probability', 'Sentence'),
                  escape = FALSE,
                  options = list(pageLength = 20))
  })
  
  
  output$barchart_intext <- renderPlotly({
    table <- collated_intext
    if (nrow(table) > 0) {
      table$termInText <- factor(as.character(table$termInText),levels=unique(table$termInText))
      
      p <- plot_ly(table, x=~termInText, y=~sentenceCount, source='barchart_intext', type = 'bar', name = 'Antigen', marker=list(color=color_TumorSuppressor)) %>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack', margin = list(b = 200), xaxis=list(title = "", tickangle = 45))%>% 
        config(displayModeBar = F)
      
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  observe({
    s <- event_data("plotly_click", source = "barchart_intext")
    if (length(s) > 0) {
      table <- collated_intext
      rowNumber <- rownames(table[table$termInText==s$x,])
      tableProxy_intext %>% selectRows(as.numeric(rowNumber))
    }
  })
  
  output$download_table_collated_intext <- downloadHandler(
    filename = function() {
      return("antigens_collatedByTermInText.tsv")
    },
    content = function(file) {
      write.table(collated_intext, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_table_sentences_normalized_all2 <- downloadHandler(
    filename = function() {
      return("antigens_sentences.tsv")
    },
    content = function(file) {
      outdata <- sentences[,c("HUGOorUniprotID","normalizedTerm","termInText","pmid","journal","year","title","probability","sentence")]
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_table_sentences_normalized_shown2 <- downloadHandler(
    filename = function() {
      return("antigens_sentences_selected.tsv")
    },
    content = function(file) {
      if(length(input$table_collated_intext_rows_selected)>0) {
        row <- collated_intext[input$table_collated_intext_rows_selected,]
        entries <- sentences[sentences$termInText==row$termInText,]
      } else {
        entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
        colnames(entries) <- colnames(sentences)
      }
      
      outdata <- entries[,c("HUGOorUniprotID","normalizedTerm","termInText","pmid","journal","year","title","probability","sentence")]
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  
  
  
  
  
  
  
  
  
  
  
}

shinyApp(ui, server)
