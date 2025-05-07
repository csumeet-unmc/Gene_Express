# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinythemes)
library(bslib)
library(readr)
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
                  c("Gene Name", "Entrez ID", "Ensembl ID")),
      
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
                     br(),
                     actionButton("view_subset_modal", "View Full Table"),
                     br(),
                     downloadButton("download", "Download")
                   )
                 )
        ),
        
        tabPanel("Plots",
                 br(), br(),
                 h4("Merged Gene Expression Table"),
                 card(
                   card_header("Merged Gene Expression Table"),
                   card_body(
                     DTOutput("barchart_table"),
                     style = "padding: 10px;"
                   )
                 ),
                 br(),
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
    read_csv(input$gene_express_file$datapath, show_col_types = FALSE)
  })
  
  # Read uploaded meta data file
  meta_data <- reactive({
    req(input$meta_data_file)
    read_csv(input$meta_data_file$datapath, show_col_types = FALSE)
  })
  
  
  
  # Output the gene expression table
  output$gene_express_data <- DT::renderDataTable({
    req(gene_data())
    datatable(gene_data(), options = list(scrollX = TRUE))
  })
  
  # Output the meta data table
  output$meta_data <- DT::renderDataTable({
    req(meta_data())
    datatable(meta_data(), options = list(scrollX = TRUE))
  })
  

  
  #Find column from meta_data
  observeEvent(meta_data(), {
    updateSelectInput(
      session,
      "selectcol_metadata",
      choices = colnames(meta_data())
    )
  })
  
  #Create filter based on user gene input 
  
  ## Subset the gene data based on user-typed genes
  
  filtered_gene_data <- reactive({
    req(gene_data())
    req(input$chooseorg)
    req(input$selectgeneformat)
    req(input$input_genes)
    
    genes_list <- trimws(unlist(strsplit(input$input_genes, ",")))
    
    data <- as.data.frame(gene_data())  # Force it to normal data.frame
    
    # Default empty table to prevent crash
    filtered_data <- data.frame()
    
    # Case 1: User selects (Mouse or Human) + (Gene Name)
    if ((input$chooseorg %in% c("Mouse", "Human")) && input$selectgeneformat == "Gene Name") {
      filtered_data <- data %>%
        dplyr::filter(.data[[names(data)[1]]] %in% genes_list)%>%
        dplyr::rename(!!"mgi_symbol" := !!sym(names(data)[1]))
      
      
      # Case 2: User selects (Mouse) + (Entrez ID)
    } else if (input$chooseorg == "Mouse" && input$selectgeneformat == "Entrez ID") {
      
      # Map Gene Symbols to ENTREZ IDs first
      data$entrez_id <- as.character(mapIds(
        x = org.Mm.eg.db,
        keys = data[[names(data)[1]]],   
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      ))
      
      # Bring entrez_id column to the front
      filtered_data <- data %>%
        dplyr::select(1 ,entrez_id, everything()) %>%
        dplyr::filter(entrez_id %in% genes_list)%>%
        dplyr::rename(!!"mgi_symbol" := !!sym(names(data)[1]))
      
      # Case 3: User selects (Mouse) + Ensembl ID 
    } else if (input$chooseorg == "Mouse" && input$selectgeneformat == "Ensembl ID") {
      colnames(data)[1] <- "mgi_symbol"  # Rename for mapping
      
      # Connect to Ensembl
      mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      
      # Get mapping from MGI to Ensembl
      gene_express_ann <- getBM(
        attributes = c("mgi_symbol", "ensembl_gene_id"),
        filters = "mgi_symbol",
        values = data$mgi_symbol,
        mart = mart
      )
      
      # Merge and filter
      filtered_data <- data %>%
        left_join(gene_express_ann, by = "mgi_symbol") %>%
        dplyr::relocate(ensembl_gene_id, .after = mgi_symbol) %>%
        dplyr::filter(ensembl_gene_id %in% genes_list)%>%
        dplyr::rename(!!"mgi_symbol" := !!sym(names(data)[1]))
    }
    
    filtered_data
  })
  
  
  # Output the filtered gene data
  output$filtered_gene_table <- renderDT({
    req(filtered_gene_data())
    filtered_gene_data()
  })
  
  
  
  # Modal popup for full Subset Data table
  observeEvent(input$view_subset_modal, {
    showModal(modalDialog(
      title = "Subset Table",
      DTOutput("subset_data_full"),
      size = "l",
      easyClose = TRUE
    ))
  })
  
  output$subset_data_full <- DT::renderDT({
    datatable(filtered_gene_data(), options = list(scrollX = TRUE))
  })
  
  
  
  #Wrangled Gene Expression Table output
  output$subset_data_full <- DT::renderDT({
    req(filtered_gene_data())
    datatable(filtered_gene_data(), options = list(scrollX = TRUE))
  })
  
  #Download the subset 
  output$downloadData <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0(input$filtered_gene_data, ".csv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      write.csv(data(), file)
      
    }
  )
  
  
  #Wrangle filtered_gene_data
  
  merge_df <- reactive({
    req(filtered_gene_data())
    req(meta_data())
    req(input$selectcol_metadata)
    req(input$selectgeneformat)
    
    md <- meta_data()
    names(md)[names(md) == "Scan name"] <- "Sample_id"
    req(input$selectcol_metadata %in% names(md))
    
    # Identify and exclude all ID/annotation columns from pivoting
    expression_data <- filtered_gene_data()
    id_cols <- c("mgi_symbol", "entrez_id", "ensembl_gene_id")
    sample_ids <- setdiff(colnames(expression_data), id_cols)
    
    # Pivot only the sample expression columns
    long_df <- expression_data %>%
      pivot_longer(
        cols = all_of(sample_ids),
        names_to = "Sample_id",
        values_to = "Expression"
      ) 
    
    # Clean up whitespace for accurate join
    long_df$Sample_id <- trimws(long_df$Sample_id)
    md$Sample_id <- trimws(md$Sample_id)
    
    # Perform the join
    joined_df <- left_join(
      long_df,
      md[, c("Sample_id", input$selectcol_metadata)],
      by = "Sample_id"
    )
    
    return(joined_df)
  })
      
    
    
  #Wrangled Gene Expression Table output
  output$barchart_table <- DT::renderDT({
    req(merge_df())
    datatable(merge_df(), options = list(scrollX = TRUE))
  })
  
  
  
  
  
  
  #Reactive: Summary data frame with log2 average expression
  
  summary_log_expression <- reactive({
    req(merge_df())
    
    data <- merge_df()
    gene_column <- if ("mgi_symbol" %in% colnames(data)) "mgi_symbol" 
    
    summary_df <- data %>%
      group_by(.data[[gene_column]]) %>%
      summarise(avg_expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
      mutate(log_of_average_expression = log2(avg_expression)) %>%
      arrange(log_of_average_expression)
    
    summary_df[[gene_column]] <- factor(summary_df[[gene_column]], levels = summary_df[[gene_column]])
    summary_df
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

  
}


shinyApp(ui = ui, server = server)

