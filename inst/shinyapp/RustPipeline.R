library(shiny)
library(R6)
library(DT)
library(shinyjs)
library(shinyFiles)
library(shinyalert)

# Plotting dependencys
library(dplyr)
library(ggplot2)
library(plotly)
library(patchwork)
library(crosstalk)

# Excel
library(readxl)  

# Pepitope
library(pepitope)



RustPipeline <- R6Class("RustPipeline",
  public = list(
    # Containts rv$all_construcs
    rv = reactiveValues(all_constructs = NULL),    

    prepare_peptide_table = function(peptide_table_paths) {
      print("Paths")
      # Print the paths to confirm
      print(peptide_table_paths)

      # Create a named list of data frames from each Excel file, using sheet names
      all_constructs <- lapply(peptide_table_paths, function(path) {
        sheet_names <- tryCatch({
          readxl::excel_sheets(path)
        }, error = function(e) {
          message("Error reading sheets from file: ", path)
          shinyalert(
            title = "Error reading sheets from peptide file", 
            text = paste("Error reading sheets from peptide file: ", paste(path, collapse = ", ")),
            type = "error"
          )
          return(NULL)
        })
        
        if (is.null(sheet_names)) return(NULL)

        # Read each sheet and store in a named list
        sheets_data <- lapply(sheet_names, function(sheet) {
          tryCatch({
            readxl::read_xlsx(path, sheet = sheet)
          }, error = function(e) {
            message("Error reading sheet: ", sheet, " from file: ", path)
            shinyalert(
              title = "Error reading sheets from peptide file", 
              text = paste("Error reading sheet: ", paste(sheet, collapse = ", ")),
              type = "error"
          )
            return(NULL)
          })
        })
        
        # Set the names of the list to the sheet names
        names(sheets_data) <- sheet_names
        return(sheets_data)
      })

      # Flatten the list structure if needed (combine lists from multiple files)
      all_constructs <- unlist(all_constructs, recursive = FALSE)

      print("Show peptide table head post mod")
      print(head(all_constructs))
      
      return(all_constructs)
    }, 




    ui = function() {
      tabPanel("Quality control",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            useShinyjs(),  # Include ShinyJS
            radioButtons("metadata_option", "Sample Metadata Source:", choices = c("Upload File", "Create Manually")),
            conditionalPanel(
              condition = "input.metadata_option == 'Upload File'",
              fileInput("samples_tsv", "Please select the samples.tsv file"),
            ),
            conditionalPanel(
              condition = "input.metadata_option == 'Create Manually'",
              textInput("sample_id", "Sample ID"),
              textInput("patient", "Patient"),
              textInput("rep", "Rep"),
              textInput("origin", "Origin"),
              textInput("barcode", "Barcode"),
              actionButton("add_row", "Add Row"),
              DTOutput("metadata_table"),
              textInput("export_filename", "Export Filename", value = "metadata.tsv"),
              downloadButton("export_table", "Export Table")
            ),
            
            #numericInput("max_mismatches", "Max Mismatches", value = 0, min = 0, max = 10),
            textInput("read_structures", "Read Structures", value = "7B+T"),
            
            fileInput("peptide_table", "Please select one or more peptide_table.xlsx files", multiple = TRUE, accept = c(".xlsx")),
            checkboxGroupInput("selected_tables", "Choose tables to keep:", choices = NULL),
            shinyFilesButton("fastq_file", "Select FASTQ File", "Please select a FASTQ file", multiple = FALSE),
            verbatimTextOutput("fastq_file_path"),
            actionButton("run_pipeline", "Run Pipeline"),
            div(id= "export_metrics", style = "display: none;",
              downloadButton("download_new_peptide_table", "Download Results: 2-all-metrics.xlsx"),
            ),
            verbatimTextOutput("status_2")
          ),
          mainPanel(
            width = 9,
              tabsetPanel(
                tabPanel("Barcode overlap", plotOutput("subplot_barcodes")),    
                tabPanel("Sample Metadata", tableOutput("sample_meta_data")),
                tabPanel("Construct Metadata", tableOutput("construct_meta_data")),
                tabPanel("Construct Counts", tableOutput("construct_counts")),
                tabPanel("Fqtk: Metrics", tableOutput("fqtk_metrics_df")),
                tabPanel("Read counts", plotlyOutput("subplot_plot_1")),
                tabPanel("Barcode Reads", plotlyOutput("subplot_plot_2")),
                tabPanel("Sample Correlation", plotlyOutput("subplot_plot_3"))
            ),
          )
        )
      )
    },

    server = function(input, output, session) {
      status_2 <- reactiveVal("Waiting for input...")
      output$status_2 <- renderText({status_2()})

      # 1. Step: Sample-handling
      # 1.1. get Path
      samples_tsv_path <- reactive({
        req(input$metadata_option == "Upload File")
        input$samples_tsv$datapath
      })
      
      # 1.2. Option-Upload: observe path if upload and start checks for required columns
      observeEvent(input$samples_tsv, {
        req(input$samples_tsv)  # Ensure file input is not NULL
        
        # Read the TSV file
        samples_tsv <- read.delim(input$samples_tsv$datapath, header = TRUE, sep = "\t")
        
        # Define required columns
        required_columns <- c("sample_id", "patient", "rep", "origin", "barcode")
        
        # Check for missing columns
        missing_columns <- setdiff(required_columns, colnames(samples_tsv))
    
        if (length(missing_columns) > 0) {
          # Log the error message and button update using runjs
        shinyalert(
                title = "Missing Columns!", 
                text = paste("Error: Missing columns in samples.tsv: ", paste(missing_columns, collapse = ", ")),
                type = "error"
          )
        } 
      })

      # 1.3. Option-Create: create sample via table 
      metadata_data <- reactiveVal(data.frame(sample_id = character(), patient = character(), rep = character(), origin = character(), barcode = character(), stringsAsFactors = FALSE))
      observeEvent(input$add_row, {
        new_row <- data.frame(sample_id = input$sample_id, patient = input$patient, rep = input$rep, origin = input$origin, barcode = input$barcode, stringsAsFactors = FALSE)
        metadata_data(rbind(metadata_data(), new_row))
        output$metadata_table <- renderDT({
          datatable(metadata_data(), editable = "cell", rownames = FALSE)
        }, server = FALSE)
      })
      output$metadata_table <- renderDT({
        datatable(metadata_data(), editable = "cell", rownames = FALSE)
      }, server = FALSE)

      output$export_table <- downloadHandler(
        filename = function() {
          input$export_filename
        },
        content = function(file) {
          write.table(metadata_data(), file, sep = "\t", row.names = FALSE, quote = FALSE)
        }
      )


      # 2. Step: Fastq-handling
      volumes = getVolumes()
      output$fastq_file_path <- renderText({ fastq_file_path() })
      shinyFiles::shinyFileChoose(input, "fastq_file",  roots = volumes, filetypes = c("gz"))
      fastq_file_path <- reactive({
        req(input$fastq_file)
        shinyFiles::parseFilePaths(volumes, input$fastq_file)$datapath 
      })



      # 3. Step: Peptide-handling
      # 3.1. get Path
      peptide_table_path <- reactive({
          req(input$peptide_table)
          input$peptide_table$datapath
      })

      # 3.2. read Constructs from tables (returns structured list of peptide tables)
      observeEvent(input$peptide_table, {
          
          # Errors resolved inside prepare_peptide_table
          #### Check here which tables shall be chosen?
          # Set all_construcs global via reative
         self$rv$all_constructs <- rust_pipeline$prepare_peptide_table(
            peptide_table = peptide_table_path()
          )

          table_names <- names(self$rv$all_constructs)
          print(table_names)
          # Show a shinyalert with a checkboxGroupInput
          showModal(modalDialog(
            title = "Select Tables",
            checkboxGroupInput("selected_tables_modal", 
                              "Choose tables to keep:", 
                              choices = table_names, 
                              selected = table_names),  # Preselect all tables
            footer = tagList(
              modalButton("Cancel"),
              actionButton("confirm_selection", "Confirm Selection")
            )
          ))
      })

      observeEvent(input$confirm_selection, {
            # Results of selection
            selected_tables <- input$selected_tables_modal

            if (length(selected_tables) > 0) {

              # Build subset with the selected
              self$rv$selected_constructs <- self$rv$all_constructs[selected_tables]

              # Close the dialog
              removeModal()
              # Display subset
              lib = "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
              self$rv$valid_barcodes = readr::read_tsv(lib, col_names=FALSE)$X1

              output$subplot_barcodes <- renderPlot({
                print(names(self$rv$selected_constructs))
                pepitope::plot_barcode_overlap(self$rv$selected_constructs, self$rv$valid_barcodes)
              })
              } else {
              # Show an error if no tables were selected
              shinyalert("Error", "Please select at least one table to proceed.", type = "error")
            }
          })

      # Step 4: Run pipeline
      observeEvent(input$run_pipeline, {
        req(fastq_file_path(), peptide_table_path())

        # Set tables from selection
        selected_tables <- self$rv$selected_constructs

        # Check if sample was uploaded or has to be created
        if (input$metadata_option == "Upload File") {
          req(samples_tsv_path())
          samples_tsv <- samples_tsv_path()
        } else {
          samples_tsv <- tempfile(fileext = ".tsv")
        }
      
         # Step 1: Run fqtk
        runjs("document.getElementById('status_2').innerText = 'Step 1/10 - Running fqtk demux...';")
        # Future optional: max_mismatches = input$max_mismatches,
        tmp_dir <- pepitope::demux_fq(fastq_file_path(), samples_tsv, input$read_structures)

        runjs("document.getElementById('status_2').innerText = 'Step 2/10 - Fqtk demux finished';")
        
        ### Check for results??? And alert if someting is missing?
        print(list.files(tmp_dir, pattern="\\.fq\\.gz$"))
        result_files = list.files(tmp_dir, pattern="\\.fq\\.gz$")
        print("All files")
        print(list.files(tmp_dir))
        shinyalert(
                title = "demux resulted in nice files", 
                text = paste("Counted files \\.fq\\.gz: ", paste(result_files, collapse = ", ")),
                type = "success"
          )
      


        #valid_barcodes = readr::read_tsv("test.txt", col_names=FALSE)
        runjs("document.getElementById('status_2').innerText = 'Step 4/10 - Run guide-counter count...';")

        dset <- tryCatch({
          ###  all_constructs = list(my_bc_type = my_data_frame)
          dset <- pepitope::count_bc(tmp_dir, selected_tables, self$rv$valid_barcodes)
                  shinyalert(
                  title = "demux resulted in nice files", 
                  text = paste("Counted files: ", paste(result_files, collapse = ", ")),
                  type = "success"
          )
        }, error = function(e) {
            # Handle the error gracefully
            shinyalert(
              title = "Error",
              text = paste("An error occurred: ", e$message),
              type = "error"
            )
            return(NULL)  # Return NULL in case of error
          })

        # colData(dset) – access the sample metadata as data.frame
        # rowData(dset) – access the construct metadata as data.frame
        # assay(dset) – access the construct counts as matrix
        # Step 2: Display fqtk metrics
        output$sample_meta_data <- renderTable({
          colData(dset)
        })
        runjs("document.getElementById('status_2').innerText = 'Step 3/10 - Display demux metrics...';")

   
        # Step 6: Output: {output}.counts.txt, {output}.-extended-counts.txt,  {output}.stats.txt display them in a table
        runjs("document.getElementById('status_2').innerText = 'Step 5/10 - Guide-counter successfully. Disyplay results';")
        
        output$construct_meta_data <- renderTable({
          rowData(dset)
        })

        # Add rownames = barcodes to count matrix
        construct_counts <- as.data.frame(assay(dset)) |> 
              mutate(barcode=rownames(dset)) |>
              dplyr::select(barcode, everything())

        output$construct_counts <- renderTable({
          construct_counts 
        })

        # "Files in tmp directory"
        # [1] "barcodes.counts.txt"
        # [2] "barcodes.extended-counts.txt"
        # [3] "barcodes.stats.txt"
        # [4] "demux-metrics.txt"
        # Read the fqtk metrics file if available
        fqtk_metrics_df <- if (file.exists(file.path(tmp_dir, "demux-metrics.txt"))) {
          read.table(file.path(tmp_dir, "demux-metrics.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        } else {
          data.frame(Message = "demux-metrics.txt not found")
        }

         # Step 2: Display metrics
        output$fqtk_metrics_df <- renderTable({
          fqtk_metrics_df
        })
   
        # Step 7: Perform plotting step with data of SummarizedExperiment 
        runjs("document.getElementById('status_2').innerText = 'Step 6/10 - Run plotting';")
        plot = pepitope::plot_reads(dset)

        output$subplot_plot_1 <- renderPlotly({
          # Convert ggplots to plotly
          subplot(ggplotly(plot[[1]], height=300), ggplotly(plot[[2]], height=300), nrows=1)
        })
        runjs("document.getElementById('status_2').innerText = 'Step 7/10 - Subplot_plot_1 done';")

        output$subplot_plot_2 <- renderPlotly({
          # Convert ggplots to plotly
          subplot(ggplotly(plot[[3]], height=300), ggplotly(plot[[4]], height=300), nrows=1)
        })


        runjs("document.getElementById('status_2').innerText = 'Step 9/10 - Subplot_plot_2 done';")

        output$subplot_plot_3 <- renderPlotly({
          # Convert ggplots to plotly
          ggplotly(pepitope::plot_distr(dset), height=500, tooltip="text")
        })
        runjs("document.getElementById('status_2').innerText = 'Step 10/10 - Subplot_plot_3 done';")

        # Show display tabs when all the plots are done not before
        shinyjs::show("export_metrics")
        #shinyjs::show("metrics_tables")
        shinyalert(
                title = "Count completed", 
                text = paste("Have fun checking the results! "),
                type = "success"
        )

        output$download_new_peptide_table <- downloadHandler(
          filename = function() {
            "2-all-metrics.xlsx"
          },
          content = function(file) {
            # 1. Export Sample Metadata (colData)
            sample_metadata <- as.data.frame(colData(dset))

            # 2. Export Construct Metadata (rowData)
            construct_metadata <- as.data.frame(rowData(dset))

            # 3. Export Construct Counts (assay)
            construct_counts <- as.data.frame(assay(dset))
            
            # Include row names in the Construct Counts data
            construct_counts <- tibble::rownames_to_column(construct_counts, var = "Barcode")

            # Combine data into a list of data frames, named as the sheet names
            data_list <- list(
              "Samples" = sample_metadata,
              "Construct_Metadata" = construct_metadata,
              "Construct_Counts" = construct_counts
            )

            # Write all sheets to an Excel file using writexl
            writexl::write_xlsx(data_list, file)
          }
        )
      })
    }
  )
)

rust_pipeline <- RustPipeline$new()

ui <- fluidPage(
  rust_pipeline$ui()
)

server <- function(input, output, session) {
  rust_pipeline$server(input, output, session)
}

shinyApp(ui, server)
