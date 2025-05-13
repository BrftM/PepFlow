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

        res = pepitope::screen_calc(dset, cfg)

        runjs("document.getElementById('status_4').innerText = 'Step 4/9 - Export 3-screen.pdf...';")
                
        return(list(res = res, cfg=cfg))
      },

      ui = function() {
          tabPanel("Co-culture screen",
              sidebarLayout(
                  sidebarPanel(
                      width = 3,
                      fileInput("final_peptide_table_1", "Please select 2-all-metrics.xlsx file"),
                      selectInput("ref_group", "Reference Group", choices = NULL),
                      selectInput("comp_group", "Comparison Group", choices = NULL),
                      actionButton("differential", "Perform Analysis"),
                      verbatimTextOutput("status_4"),
                      div(id= "export_plot_data", style = "display: none;",
                        downloadButton("download_pdf", "Download graphics PDF"),
                        downloadButton("download_plot_data", "Download data as .xlsx"),
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
            
            sheet_names <- readxl::excel_sheets(final_peptide_table_path_1())

            required_sheets <- c("Samples", "Construct_Metadata", "Construct_Counts")
            missing_sheets <- setdiff(required_sheets, sheet_names)

            if (length(missing_sheets) > 0) {
                # Display an error message indicating missing sheets
                runjs(paste("document.getElementById('status_4').innerText = 'Error missing sheet: ", paste(missing_sheets, collapse = ", "), "';", sep = ""))

                # Optionally, disable or hide further controls to prevent the user from continuing with the analysis
                runjs("document.getElementById('status_4').innerText += ', please upload a valid file sheets (Samples, Construct_Metadata, Construct_Counts)';")
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

            tables_list <- lapply(sheet_names, function(sheet) {
                readxl::read_excel(final_peptide_table_path_1(), sheet = sheet)
            })
            names(tables_list) <- sheet_names
            samples_df <- as.data.frame(tables_list[["Samples"]])

            
            # Extract and set dropdown choices
            samples <- as.data.frame(samples_df)
            #origins <- unique(make.names(samples$origin))
            origins <- unique(samples$origin)
            print(origins)

            if (length(origins) <= 2) {
                runjs("document.getElementById('status_4').innerText = 'Need at least two distinct groups in 'origin' to compare.';")
                return(NULL)
            }

            updateSelectInput(inputId = "ref_group", choices = origins, selected = origins[1])
            updateSelectInput(inputId = "comp_group", choices = origins, selected = origins[2])
        })

        observeEvent(input$differential , {
           req(final_peptide_table_path_1(), input$ref_group , input$comp_group )

            
            runjs("document.getElementById('status_4').innerText = 'Step 1/8 - Starting analysis...';")

            if (input$ref_group == input$comp_group) {
                runjs("document.getElementById('status_4').innerText = 'Reference group and comparison group must be different.';")
                return(NULL)
            }

            dset <- self$dataHandling$transform_xlsx(final_peptide_table_path_1())
            
            print("1. Sample metadata -> COL DATA DSET")
            print(head(colData(dset)))
            print("2. Construct metadata -> ROW DATA DSET")
            print(head(rowData(dset)))
            print("3. Construct metadata -> ASSAY COUNTS DSET")
            print(head(assay(dset)))

            print("Comp vs. Comp")
            print(input$comp_group)
            print(input$ref_group)

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
                print("DEBUG: Rendering interactive plot")
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
