# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
# Find out more about building applications with Shiny here:
#    https://shiny.posit.co/

library(shiny)
library(shinythemes)
library(bslib)
library(readr)
library(stringr)
library(reshape2)
library(DT)
library(dplyr)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(biomaRt)
library(SummarizedExperiment)
library(shinydashboard)
library(packer)
library(tibble)
library(sortable)
library(shinyalert)


source("func4.0.R", local = TRUE) 

# Set max upload size
options(shiny.maxRequestSize = 1000 * 1024^2)


ui <- page_sidebar(
  theme = bs_theme(bootswatch = "flatly"),
  title = "GENE EXPRESS",
  
  sidebar = sidebar(
    bg = "white",
    accordion(
      accordion_panel(
        #LAUNCHPAD
        
        "Launchpad",
        fileInput("gene_express_file", "Upload Gene Expression data (CSV UTF-8 Comma delimited)",
                  accept = ".csv", buttonLabel = "Browse", placeholder = "Upload file"),
        actionButton("view_exp_genedat", "View"),
        fileInput("meta_data_file", "Upload Meta Data (CSV UTF-8 Comma delimited)",
                  accept = ".csv", buttonLabel = "Browse", placeholder = "Upload file"),
        actionButton("view_meta_dat", "View")
      ),
      
      #CONTROL PANEL
      
      accordion_panel(
        "Control Panel",
        radioButtons("chooseorg", "Select Organism",
                     choices = c("Mouse" = "Mouse", "Human" = "Human")),
        
        selectInput("selectgeneformat", "Select Gene Format:",
                    c("SYMBOL", "ENTREZID", "ENSEMBL")),
        
        textInput("input_genes", "Enter genes (comma-separated)", value = "Kras, Trp53"),
        
        actionButton("viewsubsetdata", "View Subset Data")
      ), 
      
      #SWITCHBOARD
      
      accordion_panel("Switchboard",
      selectInput("selectcol_metadata", "Meta-Variable ", choices = NULL),  
      selectInput("x","X-axis", choices= NULL),
      selectInput("y", "Y-axis",choices = NULL),
      selectInput("facet_wrap", "Facet Wrap", choices = NULL),
      selectInput("shiny_sortable", "Sorting:",
                  choices = c("", "Rank List", "Bucket List"), selected = NULL),
      #verbatimTextOutput("ranked_result"),
      ) 
    )
  ),
  
  navset_card_underline(
    title = "VISUALISATION PLOTS",
    nav_panel(
      "Average Gene Expression",
      plotOutput("avg_gene_exp", height = "500px")  
    ),
    nav_panel("Heatmap",
              plotOutput("heatmap"), height= "500px")
  )
)






###SERVER

server <- function(input, output, session) {
  bs_themer()
  
  # Read and process gene expression data
  gene_data <- reactive({
    req(input$gene_express_file)
    input_gene_exp(input$gene_express_file$datapath)
  })
  
  # Read and process meta data, reordered to match gene_data columns
  meta_data <- reactive({
    req(input$meta_data_file, gene_data())
    raw_meta <- input_meta_data(input$meta_data_file$datapath)
  })
  
  # Modal: View uploaded gene expression data
  observeEvent(input$view_exp_genedat, {
    req(gene_data())
    showModal(modalDialog(
      title = "Gene Expression Data",
      DT::dataTableOutput("gene_express_data"),
      easyClose = TRUE,
      size = "xl"
    ))
  })
  
  output$gene_express_data <- DT::renderDataTable({
    req(gene_data())
    datatable(gene_data(), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Modal: View uploaded meta data
  observeEvent(input$view_meta_dat, {
    req(meta_data())
    showModal(modalDialog(
      title = "Meta Data",
      DT::dataTableOutput("meta_data_table"),
      easyClose = TRUE,
      size = "l"
    ))
  })
  
  output$meta_data_table <- DT::renderDataTable({
    req(meta_data())
    datatable(meta_data(), options = list(scrollX = TRUE, pageLength = 10))
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
  
  
  # Reactive: Subset based on user genes and gene format selected by user
  filtered_gene_data <- reactive({
    req(se_object(), input$selectgeneformat, input$input_genes)
    
    subset_se_object(
      se_object = se_object(),
      selectgeneformat = input$selectgeneformat,
      input_genes = input$input_genes
    )
  })
  
  
  
  
  # Output: Table of subset expression
  observeEvent(input$viewsubsetdata, {
    req(filtered_gene_data())
    showModal(modalDialog(
      title = "Subset",
      DT::dataTableOutput("subset_data"),
      easyClose = TRUE,
      size = "xl"
    ))
  })
  
  
  #Output: Subset Datatable
  output$subset_data <- DT::renderDataTable({
    req(filtered_gene_data())
    
    # Extract assay matrix and convert to data.frame
    subset_df <- as.data.frame(SummarizedExperiment::assay(filtered_gene_data()))
    datatable(subset_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  

  # Dropdown: Update based on meta data column names
  observeEvent(filtered_gene_data(), {
    req(filtered_gene_data())
    
    # Extract column names from the assay matrix (i.e., sample names)
    col_names <- colnames(SummarizedExperiment::colData(filtered_gene_data()))
    col_names <- stringr::str_replace_all(col_names, "\\.", " ")
    
    updateSelectInput(
      session,
      "selectcol_metadata",
      choices = col_names
    )
  })
  
  
  
  #Reactive summary data
  summary_df <- reactive({
    req(filtered_gene_data(), input$input_genes, input$selectcol_metadata)
    
    # Wrangle the data to bring expression and metadata into one data frame
    merge_df <- wrangle_data(
      se_object = filtered_gene_data(),
      input_genes = input$input_genes,
      selectcol_metadata = input$selectcol_metadata
    )
    
    # Filter by input$input_genes
    if (!is.null(input$input_genes)) {
      gene_list <- trimws(unlist(strsplit(input$input_genes, ",")))
      merge_df <- merge_df %>% filter(gene %in% gene_list)
    }
    
    # Ensure that merge_df is a data.frame and RETURN it
    return(merge_df)
    
    })
  
  

  #Compute Average gene expression
  compute_AGE <- reactive({
    req(summary_df(), input$input_genes, input$selectcol_metadata)
    
    df <- summary_df() %>%
      group_by(gene, .data[[input$selectcol_metadata]]) %>%
      summarise(
        avg_expression = mean(Expression, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        log2_avg_expression = log2(avg_expression + 1)
      )
    
    return(df)
  })
  
  
  #Sorting Shiny alert feature : Rank List
  observeEvent(input$shiny_sortable, {
    req(compute_AGE(), input$selectcol_metadata)
    
    if (input$shiny_sortable == "Rank List") {
      labels <- unique(compute_AGE()[[input$selectcol_metadata]])
      
      showModal(modalDialog(
        title = "Sort by Rank List",
        p("Drag to reorder"),
        rank_list(
          text = "Reorder the items below:",
          labels = as.character(labels),
          input_id = "ranked_items"
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("save_ranking", "Save Ranking")
        ),
        easyClose = TRUE
      ))
    }
  })
  
  # Display the result after ranking
  observeEvent(input$save_ranking, {
    removeModal()
    output$ranked_result <- renderPrint({
      input$ranked_items
    })
  })
  

  
  #X-axis and Y-axis Drop down for Average gene expression plot - switchboard
  observeEvent(compute_AGE(), {
    data <- compute_AGE()
    
    updateSelectInput(
      session,
      "x",
      choices = colnames(data)
    )
    
    updateSelectInput(
      session,
      "y",
      choices = colnames(data)
    )
    
    updateSelectInput(
      session,
      "facet_wrap",
      choices = colnames(data)
    )
  })
  
  
  
  
  ##Define: Plot Barchart
  plot_reactive <- reactive({
    df <- compute_AGE()
    req(df, input$x, input$y)
    
    # Safety: Check if x is numeric (since we're flipping it)
    if (!is.numeric(df[[input$x]])) {
      validate(need(FALSE, paste("Axis x must be numeric.")))
    }
    
    # Apply ranking order to metadata column
    if (!is.null(input$ranked_items)) {
      df[[input$selectcol_metadata]] <- factor(
        df[[input$selectcol_metadata]],
        levels = input$ranked_items
      )
    }
    
    # Plot
    ggplot(df, aes(x = .data[[input$x]], y = .data[[input$y]])) +
      geom_bar(stat = "identity", fill = "skyblue", position = position_dodge(width = 0.8)) +
      geom_text(
        aes(label = round(.data[[input$x]], 2)),
        position = position_dodge(width = 0.7),
        hjust = 1.00, size = 4, color = "black"
      ) +
      labs(
        title = "Average Gene Expression",
        x = input$x,
        y = input$y
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)
      ) +
      facet_wrap(as.formula(paste("~", input$facet_wrap)), nrow = 5) 
  })
  
   ## Output - Barchart
  output$avg_gene_exp <- renderPlot({
    plot_reactive() 
  })
  
}

shinyApp(ui = ui, server = server)



