# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/




library(shiny)
library(shinythemes)
library(bslib)
library(readr)
library(reshape2)
library(DT)
library(dplyr)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(biomaRt)
library(SummarizedExperiment)



# Set max upload size
options(shiny.maxRequestSize = 1000 * 1024^2)

ui <- fluidPage(
  theme = bs_theme(bootswatch = "sandstone"),
  br(),
  titlePanel("GENE EXPRESS"),
  
  sidebarLayout(
    sidebarPanel(
      br(),
      fileInput("gene_express_file", "Upload Gene Expression data (CSV UTF-8 Comma delimited)",
                accept = ".csv", buttonLabel = "Browse", placeholder = "Upload file"),
      
      fileInput("meta_data_file", "Upload Meta Data (CSV UTF-8 Comma delimited)",
                accept = ".csv", buttonLabel = "Browse", placeholder = "Upload file"),
      
      radioButtons("chooseorg", "Select Organism",
                   choices = c("Mouse" = "Mouse", "Human" = "Human")),
      
      selectInput("selectcol_metadata", "Select a column from Meta Data", choices = NULL),
      
      selectInput("selectgeneformat", "Select Gene Format:",
                  c("SYMBOL", "ENTREZID", "ENSEMBL")),
      
      textInput("input_genes", "Enter genes (comma-separated)", value = "Kras, Trp53"),
      
      actionButton("entergenes", "Enter")
    ),
    
    mainPanel(
      tabsetPanel(
        
        tabPanel("About",
                 br(), br(),
                 h4("About this Tool"),
                 p("This tool allows researchers to upload gene expression data and metadata to visualize selected genes.")
        ),
        
        tabPanel("Data Table",
                 br(), br(),
                 h4("Gene Expression Data"),
                 card(
                   card_header("Gene Expression Data"),
                   card_body(
                     DT::DTOutput("gene_express_data"),
                     br()
                   )
                 ),
                 br(),
                 
                 h4("Meta Data"),
                 card(
                   card_header("Meta Data"),
                   card_body(
                     DT::DTOutput("meta_data"),
                     br()
                   )
                 ),
                 br(),
                 
                 h4("Subset"),
                 card(
                   card_header("Subset"),
                   card_body(
                     DT::DTOutput("filtered_gene_table"),
                     actionButton("view_subset_modal", "View Full Table"),
                     downloadButton("download", "Download")
                   )
                 )
        ),
        
        tabPanel("Plots",
                 br(), br(),
                 card(
                   card_header("Log2 of Average Gene Expression"),
                   card_body(
                     plotOutput("barchart_plot"),
                     br(),
                     actionButton("view_log_modal", "Maximize"),
                     style = "padding: 10px;"
                   )
                 ),
                 
                 br(),
                 card_header("Gene Expression"),
                 card_body(
                   plotOutput("barchart_plot2"),
                   br(),
                   actionButton("view_plot2_modal", "Maximize"),
                   style = "padding: 10px;"
                 )
        ),
        
        tabPanel("Publications",
                 br(), p("Publications content goes here.")
        )
      )
    )
  )
)

















server <- function(input, output, session) {
  bs_themer()
  
  # Read uploaded gene expression file
  gene_data <- reactive({
    req(input$gene_express_file)
    input_gene_exp(input$gene_express_file$datapath)
  })
  
  # Read uploaded meta data file
  
  # Read uploaded meta data file and reorder based on gene_data column names
  meta_data <- reactive({
    req(input$meta_data_file)
    req(gene_data())  # Ensure gene data is available before using colnames
    raw_meta <- input_meta_data(input$meta_data_file$datapath)
    reordered_meta <- reorder_metadata(gene_data(), raw_meta)
    return(reordered_meta)
  })
  
  # Render gene expression table
  output$gene_express_data <- DT::renderDataTable({
    req(gene_data())
    datatable(gene_data(), options = list(scrollX = TRUE))
  })
  
  # Render meta data table
  output$meta_data <- DT::renderDataTable({
    req(meta_data())
    datatable(meta_data(), options = list(scrollX = TRUE))
  })
  
  # Update dropdown for selecting metadata column
  observeEvent(meta_data(), {
    updateSelectInput(
      session,
      "selectcol_metadata",
      choices = colnames(meta_data())
    )
  })
  
  
  
  # Create reactive SummarizedExperiment object
  se_object <- reactive({
    req(gene_data(), meta_data(), input$chooseorg)
    
    # Reorder metadata
    reordered_meta <- reorder_metadata(gene_data(), meta_data())
    
    # Create SE
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(Q3norm = gene_data()),
      colData = S4Vectors::DataFrame(reordered_meta)
    )
    
    # Add row annotations
    se <- rowData_se_obj(se, chooseorg = input$chooseorg)
    
    return(se)
  })
  
  
  
  # Reactive: Subset based on user genes + format
  filtered_gene_data <- reactive({
    req(se_object(), input$selectgeneformat, input$input_genes)
    
    subset_se_object(
      se_object = se_object(),
      selectgeneformat = input$selectgeneformat,
      input_genes = input$input_genes
    )
  })
  
  # Output: Table of subset expression
  output$filtered_gene_table <- DT::renderDT({
    req(filtered_gene_data())
    datatable(
      as.data.frame(assay(filtered_gene_data())),
      options = list(scrollX = TRUE)
    )
  })
  

  # Modal popup for full Subset data frame
  observeEvent(input$view_subset_modal, {
    showModal(modalDialog(
      title = "Subset Table",
      DTOutput("subset_data_full"),
      size = "l",
      easyClose = TRUE
    ))
  })
  
  output$subset_data_full <- DT::renderDT({
    datatable(as.data.frame(assay(filtered_gene_data())),
              options = list(scrollX = TRUE)
    )
  })
  
  
  
  #Using only filtered_gene_data
  
  #Reactive: Summary data frame with log2 average expression
  
  # Reactive: Summary data frame with log2 average expression per gene
  summary_log_expression <- reactive({
    req(filtered_gene_data())
    
    # Extract expression matrix and convert to data.frame
    expr_mat <- assay(filtered_gene_data())
    gene_ids <- rownames(expr_mat)
    
    # Pivot to long format: Gene, Sample, Expression
    long_df <- reshape2::melt(expr_mat)  # or use tidyr::pivot_longer if preferred
    colnames(long_df) <- c("Gene", "Sample", "Expression")
    
    # Summarise average and log2 expression
    summary_df <- long_df %>%
      group_by(Gene) %>%
      summarise(
        avg_expression = mean(Expression, na.rm = TRUE),
        log_of_average_expression = log2(avg_expression),
        .groups = "drop"
      ) %>%
      arrange(log_of_average_expression)
    
    # Ensure Gene is treated as an ordered factor
    summary_df$Gene <- factor(summary_df$Gene, levels = summary_df$Gene)
    return(summary_df)
  })
  
  # Output: Small in-page barplot
  output$barchart_plot <- renderPlot({
    summary_df <- summary_log_expression()
    gene_column <- colnames(summary_df)[1]
    
    ggplot(summary_df, aes(x = .data[[gene_column]], y = log_of_average_expression)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      geom_text(aes(label = round(log_of_average_expression, 2)),
                hjust = -0.1, size = 4) +
      labs(title = "Log2(Average Gene Expression)",
           x = "Gene", y = "Log2(Average Expression)") +
      theme_minimal() +
      coord_flip()
  })
  
  # Trigger modal popup on button click
  observeEvent(input$view_log_modal, {
    showModal(modalDialog(
      title = "Log2 Average Gene Expression - Full Plot",
      plotOutput("log_avg_plot_full", height = "600px"),
      size = "xl",
      easyClose = TRUE
    ))
  })
  
  # Output: Full-size barplot in modal
  output$log_avg_plot_full <- renderPlot({
    summary_df <- summary_log_expression()
    gene_column <- colnames(summary_df)[1]
    
    ggplot(summary_df, aes(x = .data[[gene_column]], y = log_of_average_expression)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      geom_text(aes(label = round(log_of_average_expression, 2)),
                hjust = -0.1, size = 4) +
      labs(title = "Log2(Average Gene Expression) – Full Plot",
           x = "Gene", y = "Log2(Average Expression)") +
      theme_minimal() +
      coord_flip()
  })
  
  
  #Create a merge dataset where: filtered_gene_data object - use the meta data and assay data 
  
  #Create a display of the contents of whats chosen in -  selectcol_metadata, 
  #it needs to pull out respective information from the subset 
  
  #generate a plot accordingly

  
  merged_object <- reactive({
    req(input$selectcol_metadata)
    req(filtered_gene_data())
  })
  
  
  
  
}


shinyApp(ui = ui, server = server)

