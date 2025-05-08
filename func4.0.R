#Function to read the gene expression data 
#1. check if the file has .csv 
#2. check if there are two columns minimum
#3. Convert first column into rownames (genes) 
#4. Return a Matrix

input_gene_exp <- function(input_gene_file) {
  # Check if file ends with .csv
  if (grepl("\\.csv$", input_gene_file)) {
    # Read the CSV as a data frame
    gdata <- read.csv(input_gene_file, check.names = FALSE)
    
    # Check for at least two columns
    if (ncol(gdata) < 2) {
      stop("The CSV must have at least two columns.")
    }
    
    # Make the 1st column the row names
    rownames(gdata) <- gdata[[1]]
    
    # Drop the first column
    gdata <- gdata[, -1, drop = FALSE]
    
    #Convert into Matrix
    gdata <-as.matrix(gdata)

    return(gdata)
    
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



#Write a function to make sure the order of metdata is same as gene expression data
#input the gene expression data and the meta data 

reorder_metadata <- function(gene_exp, meta_dat) {
  
  # 1. Check if all rownames in meta_dat exist in colnames of gene_exp
  if (all(rownames(meta_dat) %in% colnames(gene_exp))) {
    
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