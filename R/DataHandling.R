#' DataHandling R6 Class
#'
#' Class for handling data import and transformation from Excel files.
#'
#' @importFrom R6 R6Class
#' @importFrom readxl read_excel
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom R6 R6Class
#'
#' @return An R6 object of class DataHandling
#' @export
DataHandling <- R6Class("DataHandling",
    public = list(

    #' @return A new instance of DataHandling
    initialize = function() {
    },
    #' Transform Excel file into SummarizedExperiment
    #'
    #' Reads an Excel file with sheets "Samples", "Construct_Metadata", and "Construct_Counts" and
    #' converts it into a SummarizedExperiment object.
    #'
    #' @param table_path Character string. Path to the Excel file.
    #' @return A SummarizedExperiment object containing assays (counts), colData (samples), and rowData (construct metadata).
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
    #' Server logic placeholder
    #'
    #' @param input Shiny input object.
    #' @param output Shiny output object.
    server = function(input, output) {    }

    
  )
)

# Create an instance of the new class
dataHandling <- DataHandling$new()