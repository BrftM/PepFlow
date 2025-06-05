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

library(Biostrings)
RustPipeline <- R6Class("RustPipeline",
  public = list(
    # Containts rv$all_construcs
    rv = reactiveValues(all_constructs = NULL),    
    selected_sheets = reactiveVal(NULL),
    
    prepare_peptide_table = function(peptide_table_paths) {

      # Create a named list of data frames from each Excel file, using sheet names
      all_constructs <- lapply(peptide_table_paths, function(path) {
        sheet_names <- tryCatch({
          suppressWarnings(readxl::excel_sheets(path))
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
            suppressWarnings(readxl::read_xlsx(path, sheet = sheet))
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
            actionButton("help_btn_2", "Upload info ℹ️", title = "Need help for what to upload?"),
            tags$h4("1. Metadata-section"),
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

            tags$h4("2. Peptide-table-section"),
            checkboxGroupInput("selected_tables", "", choices = NULL), 
            fileInput("peptide_table", "Please select one or more peptide_table.xlsx files", multiple = TRUE, accept = c(".xlsx")),
            div(id= "show_edit_sheets", style = "display: none;",
              actionButton("edit_sheets", "Edit sheet selection and reverse complement"),
            ),

            tags$h4("3. Fastq-section"),
            shinyFilesButton("fastq_file", "Select FASTQ File", "Please select a FASTQ file", multiple = FALSE),
            verbatimTextOutput("fastq_file_path"),
            textInput("read_structures", "Read Structures", value = "7B+T"),
            #numericInput("max_mismatches", "Max Mismatches", value = 0, min = 0, max = 10),
            tags$hr(),

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
                tabPanel("Construct Counts", h4("Top 50 Constructs (based on total counts)"), tableOutput("construct_counts")),
                tabPanel("Fqtk: Metrics", tableOutput("fqtk_metrics_df")),
                tabPanel("Read counts", plotlyOutput("subplot_plot_1")),
                tabPanel("Barcode Reads", plotlyOutput("subplot_plot_2")),
                tabPanel("Barcode distribution", plotlyOutput("subplot_plot_3"))
            ),
          )
        )
      )
    },

    server = function(input, output, session) {

      observeEvent(input$help_btn_2, {
            showModal(modalDialog(
              title = "📘 Help: Upload Instructions",
              HTML(
                "<div style='line-height: 1.5;'>

                  <h4>🧾 1. Sample Sheet</h4>
                  <p>Upload or \"Create Manually\" a sample sheet with the following format:</p>
                  <pre style='background:#f8f9fa; border:1px solid #dee2e6; padding:10px;'>
                      sample_id     patient       rep origin    barcode
                      lib1          pat2+common   1   Library   AAGACCA
                      lib2          pat3          1   Library   CCAGTGT
                      mock1         pat1+common   1   Mock      TGAGTCC
                      mock2         pat1+common   2   Mock      CAAGATG
                      screen1       pat1+common   1   Sample    AACCGAC
                      screen2       pat1+common   2   Sample    AGAATCG
                  </pre>
                  <p>The barcode in the table is not the Oligo barcode of the constructs but the sample barcode!</p>

                  <h4>🔬 2. Peptide Table</h4>
                  <p>Upload a peptide table (one sheet per patient), with columns:</p>
                  <code>var_id, mut_id, pep_id, pep_type, gene_name, gene_id, tx_id, tiled, n_tiles, nt, peptide, barcode</code>
                  <p>
                    Assure that the sheets of the peptide table don’t use overlapping oligo barcodes. 
                    You can check for overlaps in the <strong>“Barcode overlap”</strong> tab by selecting the relevant sheets 
                    and choosing whether barcodes should be used in reverse complement.
                  </p>
                  <ul>
                    <li><code>barcode</code> – if only one barcode per construct</li>
                    <li><code>barcode_1</code>, <code>barcode_2</code>, etc. – if multiple barcodes per construct</li>
                  </ul>
                  <p><em>Example — Sheetname: pat1</em></p>
                  <pre style='background:#f8f9fa; border:1px solid #dee2e6; padding:10px;'>
          var_id                mut_id     pep_id       pep_type   gene_name  ...  barcode_1    barcode_2
          chr1:114713908_T/A    NRAS_Q61L  NRAS_Q61     ref        NRAS       ...  AACAACAACACC AACA...
          chr1:114713908_T/A    NRAS_Q61L  NRAS_Q61L    alt        NRAS       ...  AACAACAACGGT AACA...
                  </pre>
                  <p><em>Example — Sheetname: pat2</em></p>
                  <pre style='background:#f8f9fa; border:1px solid #dee2e6; padding:10px;'>
          var_id                mut_id     pep_id       pep_type   gene_name  ...  barcode_1    barcode_2
          chr1:114713908_T/C    NRAS_Q61R  NRAS_Q61R    alt        NRAS       ...  AACAGTGGTCTT AACC...
          chr2:208248388_G/A    IDH1_R132L IDH1_R132    ref        IDH1       ...  AACATAACGAGG AACC...
                  </pre>

                  <h4>📂 3. FASTQ File</h4>
                  <p>Upload a FASTQ file containing sequencing results from co-culture experiments.</p>

                </div>"
              ),
              easyClose = TRUE,
              size = "l",
              footer = NULL
            ))
          })


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
      selected_sheets <- reactiveVal()
      reverse_flags <- reactiveVal(list())

      # 3.1. get Path
      peptide_table_path <- reactive({
          req(input$peptide_table)
          input$peptide_table$datapath
      })


      # 3.2. display plot
      observeEvent(input$peptide_table, {
       
          self$rv$all_constructs <- rust_pipeline$prepare_peptide_table(
            peptide_table = peptide_table_path()
          )

          # Define valid barcodes
          lib = "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
          self$rv$valid_barcodes <- readr::read_tsv(lib, col_names=FALSE)$X1

          output$subplot_barcodes <- renderPlot({
            pepitope::plot_barcode_overlap(self$rv$all_constructs, self$rv$valid_barcodes)
          })

          # Select all if the default is choosen
          self$rv$selected_constructs <- self$rv$all_constructs
          
          shinyjs::show("show_edit_sheets")
      })

      # 3.3. select sheets and if they shall be rev comped -> 3.2.
      observeEvent(input$edit_sheets, {
        req(self$rv$all_constructs)
        sheet_names <- names(self$rv$all_constructs)

        current_selection <- selected_sheets() %||% sheet_names

        showModal(modalDialog(
          title = "Edit Sheets and Reverse Complement",
          tagList(
            tags$h4("1. Select Sheets to Use"),
            checkboxGroupInput("modal_selected_sheets", NULL, 
                              choices = sheet_names, 
                              selected = current_selection),
            
            tags$hr(),
            tags$h4("2. Reverse Complement Options"),
            div(
              lapply(sheet_names, function(sheet) {
                checkboxInput(paste0("rev_", sheet), paste("Reverse Complement for", sheet), 
                              value = reverse_flags()[[sheet]] %||% FALSE)
              })
            )
          ),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_edit_sheets", "Apply Selections")
          ),
          size = "l",
          easyClose = FALSE
        ))
      })

      


      # 3.2. 
      observeEvent(input$confirm_edit_sheets, {
        req(input$modal_selected_sheets)
        removeModal() 
        # Delete pre-selection from default etc. to mitigate side effects
        self$rv$selected_constructs <- NULL
 
        sel <- input$modal_selected_sheets
        selected_sheets(sel) # ✅ update global value of selected sheets

        rf <- list()
        for (sheet in names(self$rv$all_constructs)) {
          rf[[sheet]] <- input[[paste0("rev_", sheet)]] %||% FALSE
        }
        reverse_flags(rf)

        # Re-apply rev-comp and subset constructs
        processed <- list()
        for (sheet in sel) {
          df <- self$rv$all_constructs[[sheet]]

          # Support for "barcode" or "barcode_1", "barcode_2", ...
          barcode_cols <- grep("^barcode(_\\d+)?$", names(df), value = TRUE)
          if (length(barcode_cols) == 0) {
            shinyalert(
              title = "Missing barcode column",
              text = paste0("Sheet '", sheet, "' must contain at least one column named 'barcode' or 'barcode_1', 'barcode_2', etc."),
              type = "error"
            )
            return(NULL)
          }

          # Apply reverse complement if flagged
          if (rf[[sheet]]) {
            for (bc_col in barcode_cols) {
              df[[bc_col]] <- as.character(reverseComplement(DNAStringSet(as.character(df[[bc_col]]))))
            }
          }

          processed[[sheet]] <- df
        }

        names(processed) <- sel
        self$rv$selected_constructs <- processed

        
          # Re-render the barcode overlap plot
          output$subplot_barcodes <- renderPlot({
            tryCatch({
              pepitope::plot_barcode_overlap(self$rv$selected_constructs, self$rv$valid_barcodes)

            }, error = function(e) {
              shinyalert(
                title = "Duplicate barcodes or processing error",
                text = paste0("An error occurred while processing one of the selected sheets: ", e$message,
                              "\nPlease review the sheet structure and ensure unique barcodes."),
                type = "error"
              )
            })
          })
      })

      # Step 4: Run pipeline
      observeEvent(input$run_pipeline, {
        req(fastq_file_path(), peptide_table_path())

        error_occurred <- FALSE  # flag to track error
        # Set tables from selection
        selected_tables <- self$rv$selected_constructs

        tryCatch({
          pepitope:::merge_constructs(selected_tables)
        }, error = function(e) {
          shinyalert(
            title = "Error",
            text = paste("An error occurred because of overlapping barcodes:", e$message),
            type = "error"
          )
          error_occurred <<- TRUE  # set flag from error handler
        })
        if (error_occurred) return()

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
        tmp_dir <- tempdir()  # Or your custom temp directory path
        # List all files and folders inside the temp directory
        files_to_delete <- list.files(tmp_dir, full.names = TRUE, recursive = TRUE)
        
        # Exclude samples_tsv
        files_to_delete <- files_to_delete[normalizePath(files_to_delete) != normalizePath(samples_tsv)]

        # Delete them
        unlink(files_to_delete, recursive = TRUE, force = TRUE)

        
        fastq_path <- shinyFiles::parseFilePaths(volumes, input$fastq_file)$datapath

        # Check for space(s) in path and replace if present
        if (grepl(" ", fastq_path)) {

          fastq_path <- shQuote(fastq_path)
        }

        tmp_dir <- pepitope::demux_fq(fastq_path, samples_tsv, input$read_structures)

        runjs("document.getElementById('status_2').innerText = 'Step 2/10 - Fqtk demux finished';")
        
        ### Check for results??? And alert if someting is missing?
        print(list.files(tmp_dir, pattern="\\.fq\\.gz$"))
        result_files = list.files(tmp_dir, pattern="\\.fq\\.gz$")
        print("All files")
        print(list.files(tmp_dir))
     

        runjs("document.getElementById('status_2').innerText = 'Step 4/10 - Run guide-counter count...';")

      
        dset <- tryCatch({
          pepitope::count_bc(tmp_dir, selected_tables, self$rv$valid_barcodes)
        }, error = function(e) {
          shinyalert(
            title = "Error",
            text = paste("An error occurred while counting barcodes:", e$message),
            type = "error"
          )
          error_occurred <<- TRUE  # set flag from error handler
        })
        if (error_occurred) return()
        # colData(dset) – access the sample metadata as data.frame
        # rowData(dset) – access the construct metadata as data.frame
        # assay(dset) – access the construct counts as matrix
        # Step 2: Display fqtk metrics
        print("head dset")
        print(head(dset))
        print("-------------")
        print(head(colData(dset)))
        print(head(rowData(dset)))
        print(head(assay(dset)))

        output$sample_meta_data <- renderTable({
          colData(dset)
        })
        runjs("document.getElementById('status_2').innerText = 'Step 3/10 - Display demux metrics...';")

   
        # Step 6: Output: {output}.counts.txt, {output}.-extended-counts.txt,  {output}.stats.txt display them in a table
        runjs("document.getElementById('status_2').innerText = 'Step 5/10 - Guide-counter successfully. Disyplay results';")

        output$construct_meta_data <- renderTable({
          construct_meta_data <- rowData(dset)
          head(construct_meta_data, 50) 
        })


        output$construct_counts <- renderTable({
          df <- as.data.frame(assay(dset))
          df$Total <- rowSums(df)  # assuming numeric values
          df <- df[order(-df$Total), ]  # descending sort
          df <- tibble::rownames_to_column(df, var = "barcode")
          head(df, 50)
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
            construct_counts <- tibble::rownames_to_column(construct_counts, var = "barcode")

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
