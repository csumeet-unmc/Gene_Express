library(AnnotationDbi)
library(SummarizedExperiment)


#Function to read the gene expression data 
#1. check if the file has .csv 
#2. check if there are two columns minimum
#3. Convert first column into rownames (genes) 
#4. Return a Matrix

input_gene_exp <- function(input_gene_file) {
  # Check if file ends with .csv
  if (grepl("\\.csv$", input_gene_file)) {
    # Read the CSV as a data frame
    gene_express_file <- read.csv(input_gene_file, check.names = FALSE)
    
    # Check for at least two columns
    if (ncol(gene_express_file) < 2) {
      stop("The CSV must have at least two columns.")
    }
    
    # Make the 1st column the row names
    rownames(gene_express_file) <- gene_express_file[[1]]
    
    # Drop the first column
    gene_express_file <- gene_express_file[, -1, drop = FALSE]
    
    #Convert into Matrix
    gene_express_file <-as.matrix(gene_express_file)

    return(gene_express_file)
    
  } else {
    message("Please upload a valid Gene Expression CSV UTF-8 (Comma delimited) file.")
    return(NULL)
  }
}



#Function to read meta data  
#1. check if the file has .csv 
#2. Return a Data frame

input_meta_data <- function(input_meta_dat) {
  
  #Ensure that the user input file is a csv 
  if (grepl("\\.csv$", input_meta_dat)) {
    # Read the CSV as a data frame
    mdata <- read.csv(input_meta_dat, check.names = FALSE)
    
    # Check for at least two columns
    if (ncol(mdata) < 2) {
      stop("The CSV must have at least two columns.")
    }
  
  
  #Make sure its a data frame
  mdata = as.data.frame(mdata)
  
  # Make the 1st column the row names
  rownames(mdata) <- mdata[[1]]
  
  # Drop the first column
  mdata <- mdata[, -1, drop = FALSE]
  
  return(mdata)
  
  } else {
    print("Please upload a valid Meta Data CSV UTF-8 (Comma delimited) file ")
  } 
  
}



##Function to make sure the order of metadata is same as gene expression data

#input the gene expression data and the meta data 

reorder_metadata <- function(gene_exp, meta_dat) {
  
  # 1. Check if all rownames in meta_dat exist in colnames of gene_exp
  if (all(rownames(meta_dat) %in% colnames(gene_exp)))  {
    
    # 2. If order does not match exactly, then reorder meta_dat
    if (!all(colnames(gene_exp) == rownames(meta_dat))) {
      
      # 3. Reorder meta_dat to match the order of gene_exp columns
      meta_dat <- meta_dat[colnames(gene_exp), , drop = FALSE]
      
      # Final Check (optional)
      
      if (!all(colnames(gene_exp) == rownames(meta_dat))) {
        warning("Reordering failed. Samples are still misaligned.")
      } else {
        message("Metadata successfully reordered to match gene expression data.")
      }
    }
    
    return(meta_dat)
    
  } else {
    
    stop("Not all metadata rownames are present in gene expression column names.")
  }
}



#Function to ensure the rowData is mapped against Ensembl Id, Symbol and Entrez id

rowData_se_obj <- function(se_object, chooseorg) {
  # Extract gene symbols from the SE object
  gene_symbols <- rownames(assay(se_object))  # assumes rownames are SYMBOLs
  
  chooseorg <- "Mouse"  # or "Human"
  
  if (chooseorg == "Mouse") {
    # Use Mouse database
    gene_annotations <- AnnotationDbi::select(
      x = org.Mm.eg.db,
      keys = gene_symbols,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID", "ENSEMBL")
    )
    
  } else if (chooseorg == "Human") {
    # Use Human database
    gene_annotations <- AnnotationDbi::select(
      x = org.Hs.eg.db,
      keys = gene_symbols,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID", "ENSEMBL")
    )
  } else {
    stop ("Invalid organism selected: choose 'Mouse' or 'Human'")
  }
  
  # Match and reorder annotations to match rownames of SE object
  gene_annotations_ordered <- gene_annotations[match(gene_symbols, gene_annotations$SYMBOL), ]
  
  # Set rownames of annotation table to match SE object
  rownames(gene_annotations_ordered) <- gene_symbols
  
  # Attach annotations to SE object
  rowData(se_object) <- S4Vectors::DataFrame(gene_annotations_ordered)
  
  return(se_object)
}


##Function to Subset the object based on user input :
#1. Based on select gene format (selectgeneformat) : the user tells what to look at in rowData
#of the object and choose app symbol, entrez id and emsembl id 
#2. Now, it should map the genes accordingly and subset based on chose user input of genes (input_genes)
#3. Display accordingly 

subset_se_object <- function(se_object, selectgeneformat, input_genes) {
  # Safety checks
  if (is.null(se_object) || is.null(input_genes) || is.null(selectgeneformat)) {
    stop("Invalid input to subset_se_object()")
  }
  
  # Handle input_genes as either character vector or comma-separated string
  if (length(input_genes) == 1 && grepl(",", input_genes)) {
    gene_list <- trimws(unlist(strsplit(input_genes, ",")))
  } else {
    gene_list <- trimws(as.character(input_genes))  # Ensure it's character
  }
  
  # Get row annotations
  row_annots <- rowData(se_object)
  
  # Validate gene format
  if (!(selectgeneformat %in% colnames(row_annots))) {
    stop(paste0("Column '", selectgeneformat, "' not found in rowData."))
  }
  
  # Match input genes to annotation column
  matched_rows <- which(row_annots[[selectgeneformat]] %in% gene_list)
  
  # Return subset or NULL with warning
  if (length(matched_rows) > 0) {
    return(se_object[matched_rows, , drop = FALSE])
  } else {
    warning("No matching genes found for the selected gene format.")
    return(NULL)
  }
}



