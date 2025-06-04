library(R6)
library(shiny)

# Read xlsx file
library(readxl)

library(SummarizedExperiment)

DataHandling <- R6Class("DataHandling",
    public = list(
        

    #' @return SummarizedExperiment for differential expression
    transform_xlsx = function(table_path) {
        
        samples_df <- readxl::read_excel(table_path, sheet = "Samples", col_types = "text")

        construct_metadata_df <- readxl::read_excel(table_path, sheet = "Construct_Metadata", col_types = "text")
        
        counts_df <- readxl::read_excel(table_path, sheet = "Construct_Counts")

        # Drops the first column, retaining only numeric/count data.
        counts_matrix <- as.matrix(counts_df[, -1])
        # Accesses the first column of counts_df as a vector to get barcodes
        rownames(counts_matrix) <- counts_df[[1]]

        dset <- SummarizedExperiment(
            assays = list(counts = counts_matrix),
            colData = samples_df,
            rowData = construct_metadata_df
        )
        return(dset)
    },
    
    server = function(input, output) {    }

    
  )
)

# Create an instance of the new class
dataHandling <- DataHandling$new()