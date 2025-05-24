library(R6)
library(shiny)

# Read xlsx file
library(readxl)

library(SummarizedExperiment)

# Define the Test class
DataHandling <- R6Class("DataHandling",
    public = list(
    
        writeTMPxlsx = function(report) {

            # Combine all data frames into a named list for writexl
            data_list <- lapply(names(report), function(sheet_name) {
                report[[sheet_name]]
            })
            names(data_list) <- names(report)

            # Save to a temporary file
            temp_file <- tempfile(fileext = ".xlsx")
            
            # Write the data to an Excel file using writexl
            writexl::write_xlsx(data_list, path = temp_file)

            return(temp_file)
        }, 

    #' @return Dataframe List
    transform_xlsx = function(table_path) {
        
        samples_df <- as.data.frame(
            suppressWarnings(readxl::read_excel(table_path, sheet = "Samples", col_types = "text"))
        )

        # Read "Construct_Metadata" and "Construct_Counts" with default guessing (numeric when possible)
        df1 <- suppressWarnings(
            suppressWarnings(readxl::read_excel(table_path, sheet = "Construct_Metadata"))
        )
        construct_metadata_df <- as.data.frame(df1)
        
        df2 <- suppressWarnings(
            suppressWarnings(readxl::read_excel(table_path, sheet = "Construct_Counts"))
        )
        counts_df <- as.data.frame(df2)


        # Set first column as rownames then delet first column
        rownames(counts_df) <- counts_df[, 1]
        counts_df <- counts_df[, -1]
        counts_matrix <- as.matrix(counts_df)

        print("NEW COUNTS HEAD DEBUG")
        print(head(counts_matrix))


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