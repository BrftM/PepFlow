#' RustPipeline R6 Class
#'
#' This class handles the demultiplexing of fastq files and the guide counting of variants.
#'
#' @importFrom shiny fluidPage sidebarLayout sidebarPanel mainPanel tabsetPanel tabPanel
#' @importFrom shiny inputPanel textInput selectInput numericInput checkboxInput radioButtons
#' @importFrom shiny actionButton fileInput downloadButton uiOutput verbatimTextOutput textOutput
#' @importFrom shiny renderText renderUI renderPlot renderTable
#' @importFrom shiny reactive reactiveValues observe observeEvent eventReactive
#' @importFrom shiny runApp
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyalert useShinyalert
#' @importFrom shinyFiles shinyFilesButton shinyFileChoose
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom readxl read_excel excel_sheets read_xlsx
#' @importFrom writexl write_xlsx
#' @importFrom dplyr %>% filter select mutate
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom VariantAnnotation readVcf
#' @importFrom shinyFiles getVolumes
#' @importFrom Biostrings DNAStringSet reverseComplement
#' @importFrom pepitope example_peptides example_fastq demux_fq count_bc plot_distr plot_reads plot_barcode_overlap
#' @importFrom R6 R6Class
#'
#' @return An R6 object of class RustPipeline
#' @export
RustPipeline <- R6Class("RustPipeline",
  public = list(

    #' @field rv Reactive values container for internal state.
    rv = NULL,
    #' @field selected_sheets Reactive value storing user-selected Excel sheets.
    selected_sheets = NULL,
    #' @field test_mode_2 Logical flag to toggle test mode (default FALSE).
    test_mode_2 = FALSE,
    #' @field lib Character URL of barcode library file.
    lib = NULL,
    #' @field valid_barcodes Character vector of valid barcode strings.
    valid_barcodes = NULL,
    #' @field reverse_flags Reactive value storing flags for reverse complement processing.
    reverse_flags = NULL,

    #' Initialize RustPipeline object
    #'
    #' @param test_mode_2 Logical, optional. Whether to activate test mode. Default is FALSE.
    #' @return An initialized RustPipeline object.
    initialize = function(test_mode_2 = FALSE){
      self$rv <- shiny::reactiveValues(all_constructs = NULL)
      self$test_mode_2 <- test_mode_2
      self$lib <- "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
      self$valid_barcodes <- readr::read_tsv(self$lib, col_names=FALSE)$X1
      self$selected_sheets <- shiny::reactiveVal()
      self$reverse_flags <- shiny::reactiveVal()
    },
    
    #' Prepare peptide table from Excel files
    #'
    #' Reads Excel files and extracts all sheets as data frames.
    #'
    #' @param peptide_table_paths Character vector of file paths to Excel files.
    #' @return A named list of data frames with each sheet's content.
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
      
      return(all_constructs)
    }, 

    #' UI method placeholder
    #'
    #' @param ... Parameters passed to UI function (if any).
    ui = function() {
      tabPanel("Quality control",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            useShinyjs(),  
            shiny::actionButton("help_btn_2", "Upload info ℹ️", title = "Need help for what to upload?"),
            shiny::checkboxInput("use_test_data_2", "Use test data", value = FALSE),
            shiny::tags$h4("1. Peptide-table-section"),
            shiny::checkboxGroupInput("selected_tables", "", choices = NULL), 
            shiny::fileInput("peptide_table", "Please select one or more peptide_table.xlsx files", multiple = TRUE, accept = c(".xlsx")),
            shiny::div(id= "show_edit_sheets", style = "display: none;",
              actionButton("edit_sheets", "Edit sheet selection and reverse complement"),
            ),
            shiny::tags$h4("2. Metadata-section"),
            shiny::fileInput("samples_tsv", "Please select the samples.tsv file"),
            shiny::tags$h4("3. Fastq-section"),        
            shinyFilesButton("fastq_file", "Select FASTQ File", "Please select a FASTQ file", multiple = FALSE, accept = c(".fastq.gz", ".fastq")),
            shiny::verbatimTextOutput("fastq_file_path"),
            shiny::textInput("read_structures", "Read Structures", value = "7B+T"),
            shiny::tags$hr(),
            shiny::actionButton("run_pipeline", "Run Pipeline"),
           
            shiny::div(id= "export_metrics", style = "display: none;",
              shiny::downloadButton("download_new_peptide_table", "Download Results: all-metrics.xlsx"),
            ),
            shiny::verbatimTextOutput("status_2")
          ),
          mainPanel(
            width = 9,
              tabsetPanel(
                tabPanel("Barcode overlap", shiny::plotOutput("subplot_barcodes")),    
                tabPanel("Sample Metadata", shiny::tableOutput("sample_meta_data")),
                tabPanel("Construct Metadata", shiny::tableOutput("construct_meta_data")),
                tabPanel("Construct Counts", shiny::h4("Top 50 Constructs (based on total counts)"), shiny::tableOutput("construct_counts")),
                tabPanel("Fqtk: Metrics", shiny::tableOutput("fqtk_metrics_df")),
                tabPanel("Read counts", plotly::plotlyOutput("subplot_plot_1")),
                tabPanel("Barcode Reads", plotly::plotlyOutput("subplot_plot_2")),
                tabPanel("Barcode distribution", plotly::plotlyOutput("subplot_plot_3"))
            ),
          )
        )
      )
    },

    #' Server logic placeholder
    #'
    #' @param input Shiny input object.
    #' @param output Shiny output object.
    #' @param session Shiny session object.
    server = function(input, output, session) {

      observeEvent(input$help_btn_2, {
            showModal(modalDialog(
              title = "📘 Help: Upload Instructions or Use Test Data",
              HTML(
                "<div style='line-height: 1.5;'>
                  <h4>🚀 Workflow Trigger Options</h4>
                  <p>You can run the workflow in two ways:</p>
                  <ol>
                    <li><strong>Use your own data:</strong> Upload the peptide table, sample sheet, and FASTQ file.</li>
                    <li><strong>Use built-in test data:</strong> Enable the <code>Use test data</code> checkbox to run a demo pipeline.</li>
                  </ol>

                  <h4>🧪 Test Data from the following sources and description in the following section</h4>
                    <ul>
                      <li>🔬 1. Peptide Table: Generated via <code>pepitope::example_peptides()</code></li>
                      <li>🧾 2. Sample Sheet: <code>my_samples.tsv</code> included in <code>pepitope</code></li>
                      <li>📂 3. FASTQ File Simulated via <code>example_fastq()</code> function</li>
                    </ul>



                  <h4>🔬 1. Peptide Table</h4>
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

                  <h4>🧾 2. Sample Sheet</h4>
                  <p>Upload a sample sheet with the following format:</p>
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

                  

                  <h4>📂 3. FASTQ File</h4>
                  <p>Upload a FASTQ file containing sequencing results from co-culture experiments.  The read structure for the fastq sequences follows the demultiplexing format denoted as <strong>\"7B+T\"</strong>, indicating that the first 7 base pairs represent the sample barcode (<strong>\"B\"</strong>), followed by the target sequence (<strong>\"T\"</strong>)</p>
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
        req(input$samples_tsv)
        input$samples_tsv$datapath
      })
      
      # 1.2. Upload: observe path if upload and start checks for required columns
      observeEvent(input$samples_tsv, {
        req(input$samples_tsv)  
        
        # If test mode was enabled, reset it and show alert
        if (isTRUE(input$use_test_data_2)) {
          # Uncheck the test data checkbox
          updateCheckboxInput(inputId = "use_test_data_2", value = FALSE)

          self$test_mode_2 <- FALSE

          shinyalert(
              title = "Switched to real data",
              text = "A samples sheet was uploaded. Test mode is now disabled.",
              type = "info"
          )
        }

        # Read the TSV file
        samples_tsv <- tryCatch({
          read.delim(input$samples_tsv$datapath, header = TRUE, sep = "\t")
        }, error = function(e) {
          shinyalert(
            title = "File Read Error",
            text = paste("Could not read samples.tsv:", e$message),
            type = "error"
          )
          return(NULL)
        })

        required_columns <- c("sample_id", "patient", "rep", "origin", "barcode")
        
        missing_columns <- setdiff(required_columns, colnames(samples_tsv))
    
        if (length(missing_columns) > 0) {
          shinyalert(
                  title = "Missing Columns!", 
                  text = paste("Error: Missing columns in samples.tsv: ", paste(missing_columns, collapse = ", ")),
                  type = "error"
            )
        } 
      })

      # 2. Step: Fastq-handling
      ## This registers the input
      volumes = shinyFiles::getVolumes()
      shinyFiles::shinyFileChoose(input, "fastq_file",  roots = volumes, filetypes = c("gz"))

      fastq_file_path <- reactive({
        req(input$fastq_file)
        
        shinyFiles::parseFilePaths(volumes, input$fastq_file)$datapath           
      })

      observeEvent(input$fastq_file, {
        req(input$fastq_file) 
        
        output$fastq_file_path <- renderText({ fastq_file_path() })

        # If test mode was enabled, reset it and show alert
        if (isTRUE(input$use_test_data_2)) {

          updateCheckboxInput(inputId = "use_test_data_2", value = FALSE)

          self$test_mode_2 <- FALSE

          shinyalert(
              title = "Switched to real data",
              text = "A fastq file was uploaded. Test mode is now disabled.",
              type = "info"
          )
          }
      })

      # 3. Step: Peptide-handling


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

          # If test mode was enabled, reset it and show alert
          if (isTRUE(input$use_test_data_2)) {
            updateCheckboxInput(inputId = "use_test_data_2", value = FALSE)

            self$test_mode_2 <- FALSE

            shinyalert(
                title = "Switched to real data",
                text = "A peptide table was uploaded. Test mode is now disabled.",
                type = "info"
            )
          }

          output$subplot_barcodes <- renderPlot({
            pepitope::plot_barcode_overlap(self$rv$all_constructs, self$valid_barcodes)
          })
          
          shinyjs::show("show_edit_sheets")
      })

      # 3.3. select sheets and if they shall be rev comped -> 3.2.
      observeEvent(input$edit_sheets, {
        req(self$rv$all_constructs)
        sheet_names <- names(self$rv$all_constructs)

        current_selection <- self$selected_sheets() %||% sheet_names

        showModal(modalDialog(
          title = "Edit Sheets and Reverse Complement",
          tagList(
            shiny::tags$h4("1. Select Sheets to Use"),
            shiny::checkboxGroupInput("modal_selected_sheets", NULL, 
                              choices = sheet_names, 
                              selected = current_selection),
            
            shiny::tags$hr(),
            shiny::tags$h4("2. Reverse Complement Options"),
            shiny::div(
              lapply(sheet_names, function(sheet) {
                shiny::checkboxInput(paste0("rev_", sheet), paste("Reverse Complement for", sheet), 
                              value = self$reverse_flags()[[sheet]] %||% FALSE)
              })
            )
          ),
          footer = tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton("confirm_edit_sheets", "Apply Selections")
          ),
          size = "l",
          easyClose = FALSE
        ))
      })

      


      # 3.2. 
      observeEvent(input$confirm_edit_sheets, {
        req(input$modal_selected_sheets)
        shiny::removeModal() 
        # Delete pre-selection from default etc. to mitigate side effects
        self$rv$selected_constructs <- NULL
 
        sel <- input$modal_selected_sheets
        self$selected_sheets(sel) 

        rf <- list()
        for (sheet in names(self$rv$all_constructs)) {
          rf[[sheet]] <- input[[paste0("rev_", sheet)]] %||% FALSE
        }
        self$reverse_flags(rf)

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
              df[[bc_col]] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(as.character(df[[bc_col]]))))
            }
          }

          processed[[sheet]] <- df
        }

        names(processed) <- sel
        self$rv$selected_constructs <- processed

        output$subplot_barcodes <- renderPlot({
          tryCatch({
            pepitope::plot_barcode_overlap(self$rv$selected_constructs, self$valid_barcodes)

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
        use_test_2 <- isTRUE(input$use_test_data_2)

          if (use_test_2) {
            self$test_mode_2 <- TRUE
            runjs("document.getElementById('status_2').innerText = 'Step 0/10 - Loading test data...';")

            shinyalert(
              title = "Switched to test data",
              text = "Test mode is activated.",
              type = "info"
            )

            tryCatch({

              # Load peptide constructs
              self$rv$all_constructs_test <- pepitope::example_peptides(self$valid_barcodes)
              

              # Visualize barcode overlap
              output$subplot_barcodes <- renderPlot({
                pepitope::plot_barcode_overlap(self$rv$all_constructs_test, self$valid_barcodes)
              })

              # Load built-in sample sheet and fastq file
              sample_sheet <- system.file("my_samples.tsv", package = "pepitope")
              filtered_samples <- readr::read_tsv(sample_sheet, show_col_types = FALSE)

              fastq <- pepitope::example_fastq(sample_sheet, self$rv$all_constructs_test)
             
              # Store filtered samples TSV temporarily
              filtered_tsv <- tempfile(fileext = ".tsv")
              readr::write_tsv(filtered_samples, filtered_tsv)

              fastq_path <- fastq

              runjs("document.getElementById('status_2').innerText = 'Step 0/10 - Loading test data completed';")

            }, error = function(e) {
              shinyalert(
                title = "Failed to load test data",
                text = paste("An error occurred during test data setup:", e$message),
                type = "error"
              )
              runjs("document.getElementById('status_2').innerText = 'Error: Test data loading failed.';")

              return() 
            })
          } else {
            self$test_mode_2 <- FALSE

            if (!use_test_2 && is.null(samples_tsv_path()) || is.null(fastq_file_path()) || is.null(peptide_table_path()) ) {
              shinyalert(
                title = "Missing Input",
                text = "Please upload the missing file or enable 'Use test data'.",
                type = "warning"
              )
              runjs("document.getElementById('status_2').innerText = 'Waiting for input...';")
              return()
            }

            # Step 1: Read sample sheet
            samples_tsv <- samples_tsv_path()
            samples <- readr::read_tsv(samples_tsv)

            # Step 2: Update peptide table 
            self$rv$all_constructs <- rust_pipeline$prepare_peptide_table(
              peptide_table = peptide_table_path()
            )  

            if (is.null(self$rv$selected_constructs)) {
                  selected_tables <- self$rv$all_constructs
                  self$rv$selected_constructs <-  self$rv$all_constructs
            } else {
                  selected_tables <- self$rv$selected_constructs
            }

            fastq_path <- fastq_file_path()

            # Step 3: Filter only selected patients from sample sheet
            filtered_samples <- samples[samples$patient %in% names(selected_tables), ]

            # Step 4: Write filtered data to a temp file
            filtered_tsv <- tempfile(fileext = ".tsv")
            readr::write_tsv(filtered_samples, filtered_tsv)



            output$subplot_barcodes <- renderPlot({
              tryCatch({

                if (is.null(self$rv$selected_constructs)) {
                  pepitope::plot_barcode_overlap(self$rv$all_constructs, self$valid_barcodes)
                } else {
                  pepitope::plot_barcode_overlap(self$rv$selected_constructs, self$valid_barcodes)
                }

              }, error = function(e) {
                shinyalert(
                  title = "Duplicate barcodes or processing error",
                  text = paste0("An error occurred while processing one of the selected sheets: ", e$message,
                                "\nPlease review the sheet structure and ensure unique barcodes."),
                  type = "error"
                )
              })
            })
            
          }

         # Step 4: Run fqtk
        runjs("document.getElementById('status_2').innerText = 'Step 1/10 - Running fqtk demux...';")
       
        # Optional: Warn if no patients matched
        if (nrow(filtered_samples) == 0) {
          warning("No matching patients found in samples_tsv for selected_tables.")
          shinyalert(
            title = "Error",
            text = paste("No matching patients found in samples_tsv for selected_tables"),
            type = "error"
          )
        }

        # Step 5: Run demultiplexing with the filtered TSV
        tmp_dir <- tryCatch({
          pepitope::demux_fq(fastq_path, filtered_tsv, input$read_structures)
        }, error = function(e) {
          shinyalert(
            title = "Error during demultiplexing",
            text = paste("An error occurred in 'demux_fq()':", e$message),
            type = "error"
          )
          runjs("document.getElementById('status_2').innerText = 'Error during demultiplexing.';")
          return(NULL)
        })

        # Stop if demuxing failed
        if (is.null(tmp_dir)) return()

        runjs("document.getElementById('status_2').innerText = 'Step 2/10 - Fqtk demux finished';")

        runjs("document.getElementById('status_2').innerText = 'Step 4/10 - Run guide-counter count...';")

        # Step 6: Count barcodes
        dset <- tryCatch({
          # Choose the correct constructs table based on test mode
          selected_tables <- if (self$test_mode_2) {
            self$rv$all_constructs_test
          } else {
            self$rv$selected_constructs
          }

          pepitope::count_bc(tmp_dir, selected_tables, self$valid_barcodes)

        }, error = function(e) {
          shinyalert(
            title = "Error during barcode counting",
            text = paste("An error occurred in 'count_bc()':", e$message),
            type = "error"
          )
          runjs("document.getElementById('status_2').innerText = 'Error during barcode counting.';")
          return(NULL)
        })

        # Stop if counting failed
        if (is.null(dset)) return()


        # colData(dset) – access the sample metadata as data.frame
        # rowData(dset) – access the construct metadata as data.frame
        # assay(dset) – access the construct counts as matrix
        # Step 2: Display fqtk metrics
        output$sample_meta_data <- renderTable({
          SummarizedExperiment::colData(dset)
        })
        runjs("document.getElementById('status_2').innerText = 'Step 3/10 - Display demux metrics...';")

   
        # Step 6: Output: {output}.counts.txt, {output}.-extended-counts.txt,  {output}.stats.txt display them in a table
        runjs("document.getElementById('status_2').innerText = 'Step 5/10 - Guide-counter successfully. Disyplay results';")


        output$construct_meta_data <- renderTable({
          construct_meta_data <- SummarizedExperiment::rowData(dset)
          head(construct_meta_data, 50)
        })


        output$construct_counts <- renderTable({
          df <- as.data.frame(SummarizedExperiment::assay(dset))
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
          plotly::subplot(plotly::ggplotly(plot[[1]], height=300), plotly::ggplotly(plot[[2]], height=300), nrows=1)
        })
        runjs("document.getElementById('status_2').innerText = 'Step 7/10 - Subplot_plot_1 done';")

        output$subplot_plot_2 <- renderPlotly({
          plotly::subplot(plotly::ggplotly(plot[[3]], height=300), plotly::ggplotly(plot[[4]], height=300), nrows=1)
        })


        runjs("document.getElementById('status_2').innerText = 'Step 9/10 - Subplot_plot_2 done';")

        output$subplot_plot_3 <- renderPlotly({
          plotly::ggplotly(pepitope::plot_distr(dset), height=500, tooltip="text")
        })
        runjs("document.getElementById('status_2').innerText = 'Step 10/10 - Subplot_plot_3 done';")

        shinyjs::show("export_metrics")
        shinyalert(
                title = "Count completed", 
                text = paste("Have fun checking the results! "),
                type = "success"
        )
        output$download_new_peptide_table <- downloadHandler(
          filename = function() {
            paste0(Sys.Date(), "_all-metrics", ".xlsx")
          },
          content = function(file) {
            # 1. Export Sample Metadata (colData)
            sample_metadata <- as.data.frame(SummarizedExperiment::colData(dset))

            # 2. Export Construct Metadata (rowData)
            construct_metadata <- as.data.frame(SummarizedExperiment::rowData(dset))

            # 3. Export Construct Counts (assay)
            construct_counts <- as.data.frame(SummarizedExperiment::assay(dset))
      
            # Include row names in the Construct Counts data
            construct_counts <- tibble::rownames_to_column(construct_counts, var = "barcode")

            data_list <- list(
              "Samples" = sample_metadata,
              "Construct_Metadata" = construct_metadata,
              "Construct_Counts" = construct_counts
            )
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

shiny::shinyApp(ui, server)
