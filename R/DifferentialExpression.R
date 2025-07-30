#' DifferentialExpression Class
#'
#' This class handles differential expression analysis.
#'
#' @importFrom R6 R6Class
#' @importFrom shiny tabPanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny actionButton checkboxInput fileInput selectInput downloadButton
#' @importFrom shiny verbatimTextOutput textOutput plotOutput renderText renderPlot
#' @importFrom shiny observeEvent reactive reactiveVal renderUI updateCheckboxInput updateSelectInput
#' @importFrom shiny modalDialog showModal HTML req
#' @importFrom shinyalert shinyalert
#' @importFrom ggplot2 ggplot aes geom_point theme
#' @importFrom shinyjs useShinyjs runjs
#' 

#' @field test_mode_3 Logical flag for test mode.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#'
#' @return An R6 object of class DifferentialExpression.
#' @export
DifferentialExpression <- R6Class("DifferentialExpression",
    public = list(

    test_mode_3 = FALSE,
    
    #' Initialize DifferentialExpression
    #'
    #' @param test_mode_3 Logical, whether to use test data. Default FALSE.
    initialize = function(test_mode_3 = NULL) {
        self$test_mode_3 <- test_mode_3
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

    
    #' Run differential expression analysis
    #'
    #' @param dset Output from transform_xlsx()
    #' @param ref_group String name of the reference group (e.g., "CDK4")
    #' @param comp_group String name of the comparison group (e.g., "Bcell.only")
    #' @return List containing a SummarizedExperiment and the comparison config
    run_differential_expression = function(dset, ref_group, comp_group) {
        message("Step 2: Performing differential expression analysis")
        
        cfg <- list(c(comp_group, ref_group))

        runjs("document.getElementById('status_4').innerText = 'Step 3/9 - Perform diff_expr';")
        
        res <- tryCatch({
            pepitope::screen_calc(dset, cfg)
        }, error = function(e) {
            shinyalert(
                title = "Error",
                text = paste("An error occurred:", e$message),
                type = "error"
            )
        })
        
        runjs("document.getElementById('status_4').innerText = 'Step 4/9 - Export 3-screen.pdf...';")
                
        return(list(res = res, cfg=cfg))
      },
        
      #' Generate UI components
      #' @return shiny UI elements for tab panel
      ui = function() {
          tabPanel("Co-culture screen",
              sidebarLayout(
                  sidebarPanel(
                      width = 3,
                      actionButton("help_btn_3", "Upload info ℹ️", title = "Need help for what to upload?"),
                      checkboxInput("use_test_data_3", "Use test data", value = FALSE),
                      fileInput("final_peptide_table_1", "Please select all-metrics.xlsx file"),
                      selectInput("ref_group", "Reference Group", choices = NULL),
                      selectInput("comp_group", "Comparison Group", choices = NULL),
                      actionButton("differential", "Perform Analysis"),
                      verbatimTextOutput("status_4"),
                      div(id= "export_plot_data", style = "display: none;",
                        downloadButton("download_pdf", "Download whole graphics as .PDF"),
                        downloadButton("download_plot_data", "Download whole data as .xlsx"),
                      ),
                  ),
                  mainPanel(
                        width = 9,
                        plotlyOutput("de_plot") 
                  )
              )
          )
      },
      #' Server logic for differential expression module
      #'
      #' @param input Shiny input object
      #' @param output Shiny output object
      server = function(input, output) {
        wd <- normalizePath(".")

        status_4 <- reactiveVal("Waiting for input...")
        output$status_4 <- renderText({ status_4() })

        final_peptide_table_path_1 <- reactive({
            req(input$final_peptide_table_1$datapath)
            input$final_peptide_table_1$datapath
        }) 
        
        observeEvent(input$final_peptide_table_1, {
            req(final_peptide_table_path_1())       

            # If test mode was enabled, reset it and show alert
            if (isTRUE(input$use_test_data_3)) {
            # Uncheck the test data checkbox
            updateCheckboxInput(inputId = "use_test_data_3", value = FALSE)

            self$test_mode_3 <- FALSE

            # Show user notification
            shinyalert(
                title = "Switched to real data",
                text = "A Metrics file was uploaded. Test mode is now disabled, proceed with \"Perform Analysis\".",
                type = "info"
            )
            }     
        
            sheet_names <- suppressWarnings(readxl::excel_sheets(final_peptide_table_path_1()))
           
            required_sheets <- c("Samples", "Construct_Metadata", "Construct_Counts")
            missing_sheets <- setdiff(required_sheets, sheet_names)

            if (length(missing_sheets) > 0) {

                runjs(paste("document.getElementById('status_4').innerText = 'Error missing sheet: ", paste(missing_sheets, collapse = ", "), "';", sep = ""))

                shinyalert(
                        title = "Missing Sheet(s)!", 
                        text = paste("Error missing sheet in all-metrics file: ", paste(missing_sheets, collapse = ", ")),
                        type = "error"
                )
                # Clear the UI selections (in case there are any selections made)
                updateSelectInput(inputId = "ref_group", choices = NULL)
                updateSelectInput(inputId = "comp_group", choices = NULL)
            }

            samples_df <- as.data.frame(
                suppressWarnings(readxl::read_excel(final_peptide_table_path_1(), sheet = "Samples", col_types = "text"))
            )


            # Select only unique origins for display of comparison choices
            origins <- unique(samples_df$origin)
            
            if (length(origins) <2) {
                shinyalert(
                        title = "Not enough origins!", 
                        text = paste("Error missing distinct groups in 'origin' only one available: ", paste(origins, collapse = ", ")),
                        type = "error"
                )
                runjs("document.getElementById('status_4').innerText = 'Need at least two distinct groups in origin or short to compare.';")            
            } else {
                updateSelectInput(inputId = "ref_group", choices = origins, selected = origins[1])
                updateSelectInput(inputId = "comp_group", choices = origins, selected = origins[2])
                runjs("document.getElementById('status_4').innerText = 'Two distinct origins provided.';")
            }
        })

        observeEvent(input$help_btn_3, {
            showModal(modalDialog(
                title = "📘 Help: Upload Instructions",
                HTML(
                "<div style='line-height: 1.5;'>
                <h4>🚀 Test Mode</h4>
                <p>You can check the <code>Use test data</code> box to skip file uploads and run an example differential expression analysis based on internal data (pat1) otherwise supply the following described file.</p>


                    <h4>🗂️ Excel File Format (all-metrics.xlsx)</h4>
                    <p>Upload an Excel file, resulting from the previous step Quality control, containing the following three sheets:</p>
                    <ul>
                    <li><strong>Samples</strong></li>
                    <li><strong>Construct_Metadata</strong></li>
                    <li><strong>Construct_Counts</strong></li>
                    </ul>
                    <p>⚠️ All three sheets must contain only data from a single patient.</p>

                    <h4>🧾 Samples</h4>
                    <p>Must include one unique patient only and at least two different origins. Format example:</p>
                    <pre style='background:#f8f9fa; border:1px solid #dee2e6; padding:10px;'>
            sample_id  patient       rep  origin  barcode  total_reads  mapped_reads  smp        short                  label
            mock1      pat1+common   1    Mock    TGAGTCC  224687       224687        Mock-1     pat1+common Mock-1     pat1+common Mock-1 (mock1)
            screen1    pat1+common   1    Sample  AACCGAC  454355       454355        Sample-1   pat1+common Sample-1   pat1+common Sample-1 (screen1)
                    </pre>

                    <h4>🧬 Construct_Metadata</h4>
                    <p>This sheet must also contain only entries for one patient.</p>
                    <pre style='background:#f8f9fa; border:1px solid #dee2e6; padding:10px;'>
            barcode        bc_type  var_id             gene_name  mut_id         pep_id        pep_type  gene_id        tx_id           n_tiles  BbsI_replaced  tiled  nt  peptide
            AACAACCATCCA   pat1     chr1:46458643_T/C  BKLA1      KLN2A1_EdP     KLN2A1_Ed     ref       KNLG...        MDST...         1        0              ...    93  LEDDAA...
            AACAACCGCATT   pat1     chr1:56458644_C/T  BKLA1      KLN2A1_EdL     KLN2A1_EdL    alt       KNLG...        MDST...         1        0              ...    93  LEDDAA...
                    </pre>

                    <h4>📊 Construct_Counts</h4>
                    <p>Must contain count data for the same patient as in the Samples sheet.</p>
                    <ul>
                    <li>The number of columns must match the number of sample rows of the used sample sheet.</li>
                    <li>All samples must belong to the same patient.</li>
                    </ul>
                    <pre style='background:#f8f9fa; border:1px solid #dee2e6; padding:10px;'>
                           mock1  screen1   ...
            AACAACCATCCA   17     25    
            AACAACCGCATT   17     2     
            AACAACACAAGC   0      0     
            ...
                    </pre>

                </div>"
                ),
                easyClose = TRUE,
                size = "l",
                footer = NULL
            ))
            })

        observeEvent(input$differential , {

            use_test_3 <- isTRUE(input$use_test_data_3)
            self$test_mode_3 <- use_test_3

            if (!use_test_3 && is.null(final_peptide_table_path_1())) {
                shinyalert(
                title = "Missing Input",
                text = "Please upload a metrics file or enable 'Use test data'.",
                type = "warning"
                )
                runjs("document.getElementById('status_4').innerText = 'Waiting for input...';")
            }
            runjs("document.getElementById('status_4').innerText = 'Step 0/9 - Setup Required data...';")

            dset <- tryCatch({
                if (use_test_3) {
                    runjs("document.getElementById('status_4').innerText = 'Step 0/10 - Loading test data...';")
                    shinyalert(
                        title = "Switched to test data",
                        text = "Test mode is activated.",
                        type = "info"
                    )
                    lib <- "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
                    valid_barcodes <- readr::read_tsv(lib, col_names = FALSE)$X1
                    constructs <- pepitope::example_peptides(valid_barcodes)
                    sample_sheet <- system.file("my_samples.tsv", package = "pepitope")
                    fastq <- pepitope::example_fastq(sample_sheet, constructs)
                    tmp_dir <- pepitope::demux_fq(fastq, sample_sheet, read_structures = "7B+T")
                    dset <- pepitope::count_bc(tmp_dir, constructs, valid_barcodes)

                    # Keep only 'pat1' samples
                    dset <- dset[, grepl("pat1", dset$patient)]

                    runjs("document.getElementById('status_4').innerText = 'Step 0/10 - Loading test data completed';")
                } else {
                    req(final_peptide_table_path_1()) 
                    dset <- self$transform_xlsx(final_peptide_table_path_1())
                }
                dset 
                }, error = function(e) {
                shinyalert(
                    title = "Data Load Failed",
                    text = paste("Could not prepare input data:", e$message),
                    type = "error"
                )
                runjs("document.getElementById('status_4').innerText = 'Error: Failed to load input data.';")
                })

                if (is.null(dset)) return(runjs("document.getElementById('status_4').innerText = 'Error: Dset is NUll.';"))

            runjs("document.getElementById('status_4').innerText = 'Step 1/9 - Starting analysis...';")
           
            res_list <- tryCatch({        
                if (use_test_3) {  
                    res_list <- self$run_differential_expression(dset, "Mock", "Sample")
                    comparison_name <- "Sample vs Mock"
                    res_list
                } else {  
                    res_list <- self$run_differential_expression(dset, input$ref_group , input$comp_group )
                    comparison_name <- paste(input$comp_group, "vs", input$ref_group)
                    res_list
                }
            }, error = function(e) {
                shinyalert(
                title = "Differential Expression Error",
                text = paste("Error in differential expression analysis:", e$message),
                type = "error"
                )
                runjs("document.getElementById('status_4').innerText = 'Error: Failed to run differential expression.';")
            })

            if (is.null(res_list$res)) return()
                    
            runjs("document.getElementById('status_4').innerText = 'Step 7/9 - Analysis completed successfully';")
                       
            plt <- tryCatch({
                pepitope::plot_screen(res_list$res[[comparison_name]])
            }, error = function(e) {
                shinyalert(
                    title = "Plotting Error",
                    text = paste("Could not generate plot:", e$message),
                    type = "error"
                )
            })
            if (is.null(plt)) return()

            output$de_plot <- renderPlotly({
                plotly::ggplotly(plt, tooltip = "text")
            })
    
            shinyjs::show("export_plot_data")

            output$download_pdf <- downloadHandler(
            filename = function() {
                paste0("Differential_Expression_Result_", Sys.Date(), "_patient_", unique(dset$patient), "_", comparison_name, ".pdf")
            },
            content = function(file) {
                plot_with_title <- plt + ggplot2::ggtitle(paste("Differential Expression: Patient", unique(dset$patient), "-", comparison_name))
                pdf(file, width = 8, height = 6)
                print(plot_with_title)  
                dev.off()  
            }
            )

            output$download_plot_data <- downloadHandler(
                filename = function() {
                    paste0("diff_expression_summary_", Sys.Date(), "_patient_", unique(dset$patient), "_", comparison_name, ".xlsx")
                },
                content = function(file) {
                    sheet_list <- list()
                    sheet_list[[comparison_name]] <- res_list$res[[comparison_name]]
                    writexl::write_xlsx(sheet_list, path = file)
                }
            )
            runjs("document.getElementById('status_4').innerText = 'Step 9/9 - Plotting successful';")
        })
    }
  )
)
