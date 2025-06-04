## Class system (R6Class)
library(R6)
## UI/server logic
library(shiny)
## Data wrangling (dplyr, tibble, purrr)
library(dplyr)
## Static plotting
library(ggplot2)
library(ggpp)

## Interactive plots
library(plotly)

#Define the DifferentialExpression class
DifferentialExpression <- R6Class("DifferentialExpression",
    public = list(

    dataHandling = NULL,
    

    initialize = function(dataHandling) {
    self$dataHandling <- dataHandling
    },
    
    #' Wrapper method that builds config, triggers function calls to pepitope.
    #' @param table_list: Output from transform_xlsx()
    #' @param ref_group: String name of the reference group (e.g., "CDK4")
    #' @param comp_group: String name of the comparison group (e.g., "Bcell.only")
    #' Returns: List containing a SummarizedExperiment and the comparison
    run_differential_expression = function(dset, ref_group, comp_group) {
        message("Step 2: Performing differential expression analysis")
        
        # Apply make.names to match factor levels of DESeq2
        print(paste("Compare", ref_group, "vs.", comp_group))
        cfg <- list(c(comp_group, ref_group))

        runjs("document.getElementById('status_4').innerText = 'Step 3/9 - Perform diff_expr';")
        
        res <- tryCatch({
            # Attempt to run the function
            print(length(unique(dset$patient)))

            pepitope::screen_calc(dset, cfg)
        }, error = function(e) {
            # Show the error to the user via shinyalert
            shinyalert(
                title = "Error",
                text = paste("An error occurred:", e$message),
                type = "error"
            )
        # Return NULL (or another fallback value)
        return(NULL)
        })
        
        runjs("document.getElementById('status_4').innerText = 'Step 4/9 - Export 3-screen.pdf...';")
                
        return(list(res = res, cfg=cfg))
      },

      ui = function() {
          tabPanel("Co-culture screen",
              sidebarLayout(
                  sidebarPanel(
                      width = 3,
                      actionButton("help_btn_3", "Upload info ℹ️", title = "Need help for what to upload?"),
                      fileInput("final_peptide_table_1", "Please select 2-all-metrics.xlsx file"),
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
        
            sheet_names <- suppressWarnings(readxl::excel_sheets(final_peptide_table_path_1()))
           
            required_sheets <- c("Samples", "Construct_Metadata", "Construct_Counts")
            missing_sheets <- setdiff(required_sheets, sheet_names)

            if (length(missing_sheets) > 0) {
                # Display an error message indicating missing sheets
                runjs(paste("document.getElementById('status_4').innerText = 'Error missing sheet: ", paste(missing_sheets, collapse = ", "), "';", sep = ""))

                # Optionally, disable or hide further controls to prevent the user from continuing with the analysis
                shinyalert(
                        title = "Missing Sheet(s)!", 
                        text = paste("Error missing sheet in all-metrics file: ", paste(missing_sheets, collapse = ", ")),
                        type = "error"
                )
                # Clear the UI selections (in case there are any selections made)
                updateSelectInput(inputId = "ref_group", choices = NULL)
                updateSelectInput(inputId = "comp_group", choices = NULL)

                return(NULL)  # Stop further processing if sheets are missing
            }

            samples_df <- as.data.frame(
                suppressWarnings(readxl::read_excel(final_peptide_table_path_1(), sheet = "Samples", col_types = "text"))
            )


            # Select only unique origins to make it selectable
            origins <- unique(samples_df$origin)
            
            if (length(origins) <2) {
                shinyalert(
                        title = "Not enough origins!", 
                        text = paste("Error missing distinct groups in 'origin' only one available: ", paste(origins, collapse = ", ")),
                        type = "error"
                )
                runjs("document.getElementById('status_4').innerText = 'Need at least two distinct groups in origin or short to compare.';")
            
                return(NULL)
            
            } else {
                # Use 'origin' for grouping
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

                    <h4>🗂️ Excel File Format (2-all-metrics.xlsx)</h4>
                    <p>Upload an Excel file containing the following three sheets:</p>
                    <ul>
                    <li><strong>Samples</strong></li>
                    <li><strong>Construct_Metadata</strong></li>
                    <li><strong>Construct_Counts</strong></li>
                    </ul>
                    <p>⚠️ All three sheets must contain only data from a single patient.</p>

                    <h4>🧾 1. Samples Sheet</h4>
                    <p>Must include one unique patient only and at least two different origins. Format example:</p>
                    <pre style='background:#f8f9fa; border:1px solid #dee2e6; padding:10px;'>
            sample_id  patient       rep  origin  barcode  total_reads  mapped_reads  smp        short                  label
            mock1      pat1+common   1    Mock    TGAGTCC  224687       224687        Mock-1     pat1+common Mock-1     pat1+common Mock-1 (mock1)
            screen1    pat1+common   1    Sample  AACCGAC  454355       454355        Sample-1   pat1+common Sample-1   pat1+common Sample-1 (screen1)
                    </pre>

                    <h4>🧬 2. Construct_Metadata Sheet</h4>
                    <p>This sheet must also contain only entries for one patient.</p>
                    <pre style='background:#f8f9fa; border:1px solid #dee2e6; padding:10px;'>
            barcode        bc_type  var_id             gene_name  mut_id         pep_id        pep_type  gene_id        tx_id           n_tiles  BbsI_replaced  tiled  nt  peptide
            AACAACCATCCA   pat1     chr1:46458643_T/C  BKLA1      KLN2A1_EdP     KLN2A1_Ed     ref       KNLG...        MDST...         1        0              ...    93  LEDDAA...
            AACAACCGCATT   pat1     chr1:56458644_C/T  BKLA1      KLN2A1_EdL     KLN2A1_EdL    alt       KNLG...        MDST...         1        0              ...    93  LEDDAA...
                    </pre>

                    <h4>📊 3. Construct_Counts Sheet</h4>
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
           req(final_peptide_table_path_1(), input$ref_group , input$comp_group )

            
            runjs("document.getElementById('status_4').innerText = 'Step 1/8 - Starting analysis...';")

            if (input$ref_group == input$comp_group) {
                runjs("document.getElementById('status_4').innerText = 'Reference group and comparison group must be different.';")
                shinyalert(
                        title = "Can´t compare group with itself.", 
                        text = paste("Error Reference group and comparison group must be different!"),
                        type = "error"
                )
                return(NULL)
            }

            dset <- self$dataHandling$transform_xlsx(final_peptide_table_path_1())
            
            # Run the differential expression analysis
            res_list <- self$run_differential_expression(dset, input$ref_group , input$comp_group )
        
            runjs("document.getElementById('status_4').innerText = 'Step 7/9 - Analysis completed successfully';")
            

            # Extract the comparison name from the configuration
            comparison_name <- paste(input$comp_group, "vs", input$ref_group)

            # Create the 'text' column in res_list$res
            res_list$res[[1]]$text <- with(res_list$res[[1]], sprintf(
                "%s %s (%s)\nFC %.1fx p=%.2g",
                pep_id, pep_type, barcode,
                sign(stat) * 2^abs(log2FoldChange), pvalue
            ))
            # Remove " NA" from the text
            res_list$res[[1]]$text <- gsub(" NA", "", res_list$res[[1]]$text)

            plt <- ggplot(res_list$res[[1]], aes(x = baseMean, y = log2FoldChange, text = text)) +
                geom_point(aes(color = bc_type, shape = pep_type, size = padj < 0.1, alpha = padj < 0.1)) +
                scale_shape_manual(values = c(ref = 1, alt = 19), na.value = 19) +
                scale_size_manual(values = c("TRUE" = 2, "FALSE" = 0.8), na.value = 0.8) +
                scale_alpha_manual(values = c("TRUE" = 0.7, "FALSE" = 0.3), na.value = 0.3) +
                scale_color_brewer(palette = "Set1") +
                scale_x_log10(limits = c(10, NA)) +
                ggtitle("Differential Expression Plot") +
                theme_minimal()

            # Render the interactive Plotly plot
            output$de_plot <- renderPlotly({
                if (exists("plt") && !is.null(plt)) {
                    ggplotly(plt, tooltip = "text")  # Display the interactive plot with tooltips
                } else {
                    warning("WARNING: No valid plots available to render")
                    NULL
                }
            })
    
            # Provide the PDF file for download the plot
            output$download_pdf <- downloadHandler(
                filename = function() {
                    paste0("Differential_Expression_Result_", Sys.Date(), ".pdf")
                },
                content = function(file) {
                    # Create a temporary file for the plot
                    temp_plot <- tempfile(fileext = ".pdf")
                    
                    # Create a static version of the plot using ggplot2
                    pdf(temp_plot, width = 8, height = 6)  # Open PDF device
                    print(plt)  # Print the ggplot object directly
                    dev.off()  # Close PDF device

                    # Copy the plot to the desired file location
                    file.copy(temp_plot, file, overwrite = TRUE)
                }
            )


            # Flatten each comparison and tag with sample name to be able to export proper as table
            #export_result <- self$dataHandling$summarize_full_diff_expr(res_list$res)
            # Drop columns that are completely NA
            cleaned_df <- res_list$res[[comparison_name]] %>%
                select_if(~ !all(is.na(.)))
            
            # Print the cleaned data frame
            print("Cleaned_df for export:")
            print(head(cleaned_df))

            shinyjs::show("export_plot_data")

            # Export table
            output$download_plot_data <- downloadHandler(
            filename = function() {
                paste0("diff_expression_summary_", Sys.Date(), "_", comparison_name, ".xlsx")
            },
            content = function(file) {
                # Create a named list where the name is set using the variable comparison_name
                sheet_list <- list()
                sheet_list[[comparison_name]] <- cleaned_df
                
                # Write the data to an Excel file with the specified sheet name
                writexl::write_xlsx(sheet_list, path = file)
            }
            )
            runjs("document.getElementById('status_4').innerText = 'Step 9/9 - Plotting successful';")
        })
    }
  )
)
