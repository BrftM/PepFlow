library(R6)
library(shiny)

# Read xlsx file
library(readxl)

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
        required_sheets <- c("Samples", "Construct_Metadata", "Construct_Counts")
        sheet_names <- readxl::excel_sheets(table_path)

        tables_list <- lapply(sheet_names, function(sheet) {
            readxl::read_excel(table_path, sheet = sheet)
        })
        names(tables_list) <- sheet_names

        samples_df <- as.data.frame(tables_list[["Samples"]])
        construct_metadata_df <- as.data.frame(tables_list[["Construct_Metadata"]])
        counts_df <- as.data.frame(tables_list[["Construct_Counts"]])

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