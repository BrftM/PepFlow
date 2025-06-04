VariantProcessor <- R6Class("VariantProcessor",
  public = list(
    ens106 = NULL,
    asm = NULL,
    rv = NULL,

    rv_sheet = reactiveValues(sheet_data = NULL, table_data = NULL),

    dataHandling = NULL,  # Field for the helper instance

    initialize = function(dataHandling) {
      self$ens106 <- AnnotationHub::AnnotationHub()[["AH100643"]]
      seqlevelsStyle(self$ens106) <- "UCSC"
      self$asm <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
      self$asm@seqinfo@genome[] <- "GRCh38"
      self$rv <- reactiveValues(report_data = list())
      
      self$dataHandling <- dataHandling
    },

    process_vcf = function(vcf_file) {
      runjs("document.getElementById('status_1').innerText = 'Step 2/8 - Read VcF as VRanges...';")

      # Reading and filtering mutations
      vr1 <- pepitope::readVcfAsVRanges(vcf_file) |>
        pepitope::filter_variants(min_cov=2, min_af=0.05, pass=TRUE, chrs="default")

      runjs("document.getElementById('status_1').innerText = 'Step 3/8 - Annotate coding and Subset context...';")
      # Annotating and subsetting expressed variants
      ann <- pepitope::annotate_coding(vr1, self$ens106, self$asm)  
      subs <- ann |>
        pepitope::subset_context(15)

      runjs("document.getElementById('status_1').innerText = 'Step 4/8 Tiling cDNAs of interest into smaller peptides ...';")
      # Tiling cDNAs of interest into smaller peptides
      tiled <- pepitope::make_peptides(subs) |>
        pepitope::pep_tile() |>
        pepitope::remove_cutsite(BbsI="GAAGAC")


      runjs("document.getElementById('status_1').innerText = 'Step 5/8 - Make report...';")
      # Put all sheets of Report into reactive value to work with globally
      self$rv_sheet$report <- tryCatch({
        pepitope::make_report(vars=ann, subs=subs, tiled=tiled)
        shinyalert(
          title = "Report completed", 
          text = paste("Barcodes can be added in sheet \"93 nt Peptides.\""),
          type = "success"
        )
      }, error = function(e) {
        warning(paste("Error in make_report:", e$message))
          shinyalert(
                title = "Failed to generate report!", 
                text = paste("Error: Failed to generate report: ", paste(e$message)),
                type = "error"
          )
        return(NULL)
      })
    },

    display_table = function(output, input) {
      runjs("document.getElementById('status_1').innerText = 'Step 7/8 - Add barcodes or download';")
      output$dynamic_table <- renderUI({
        req(self$rv_sheet$report)

        tab_list <- lapply(names(self$rv_sheet$report), function(sheet_name) {
         
          is_editable <- identical(sheet_name, "93 nt Peptides")
          ns <- NS(sheet_name)

          tabPanel(
            title = sheet_name,
            tagList(
              if (is_editable) {
                tagList(
                  textAreaInput(ns("barcode_input"), "Barcodes", "", rows = 4),
                  actionButton(ns("add_barcode"), "Add Barcodes")
                )
              },
              DTOutput(ns("datatable"))
            )
          )
        })
        do.call(tabsetPanel, c(Filter(Negate(is.null), tab_list), id = "sheet_tabs"))
      })

      # Render tables and setup barcode handler
      lapply(names(self$rv_sheet$report), function(sheet_name) {
       
        local({
          sheet <- sheet_name
          ns <- NS(sheet)

          is_editable <- identical(sheet_name, "93 nt Peptides")

          output[[ns("datatable")]] <- renderDT({
            df <- self$rv_sheet$report[[sheet]]
            datatable(df, editable = FALSE,
                      options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE),
                      rownames = FALSE)
          })

          if (is_editable) {
            observeEvent(input[[ns("add_barcode")]], {
              req(input[[ns("barcode_input")]])

              barcode_list <- unlist(strsplit(input[[ns("barcode_input")]], "[,\n]+"))
              barcode_list <- trimws(barcode_list)

              # Filter out empty lines so that in the next step when count happens the error can occure.
              barcode_list <- barcode_list[barcode_list != ""]

              df <- self$rv_sheet$report[[sheet]]

              if (length(barcode_list) != nrow(df)) {
                error_msg <- sprintf(
                  "Error - %d barcodes provided, but %d rows are needed.",
                  length(barcode_list),
                  nrow(df)
                )
                runjs(sprintf("document.getElementById('status_1').innerText = '%s';", error_msg))
                return()
              }

              df$barcode <- barcode_list
              print("Add barcodes to sheet:")
              self$rv_sheet$report[["93 nt Peptides"]] <- df
              print(head(self$rv_sheet$report[["93 nt Peptides"]]))

              runjs("document.getElementById('status_1').innerText = 'Step 8/8 - Barcodes added!';")

              output[[ns("datatable")]] <- renderDT({
                datatable(df, editable = FALSE,
                          options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE),
                          rownames = FALSE)
              })
            })
          }
        })
      })
      shinyjs::show("barcode_export")
    },

    ui = function() {
      tabPanel(" Variant calling",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            fileInput("vcf_file", "Upload VCF File (.vcf.gz)", accept = ".vcf.gz"),
            div(id = "run_annotation", style = "display: none;",
                actionButton("process", "Process VCF")
            ),
            
            div(id = "barcode_export", style = "display: none;",
              downloadButton("download_peptide_table", "Download Peptide Table")
            ),
            verbatimTextOutput("status_1")
          ),
          mainPanel(
            width = 9,
            uiOutput("dynamic_table")
          )
        )
      )
    },

    server = function(input, output) {
      output$status_1 <- renderText({ "Waiting for input..." })

       observeEvent(input$vcf_file, {
        req(input$vcf_file)
        shinyjs::show("run_annotation")
       })

      observeEvent(input$process, {
        req(input$vcf_file)

        runjs("document.getElementById('status_1').innerText = 'Step 1/8 - Processing VCF file...';")

        self$process_vcf(input$vcf_file$datapath)

        runjs("document.getElementById('status_1').innerText = 'Step 6/8 - VCF file processed';")

        self$display_table(output, input)
      })
      
      output$download_peptide_table <- downloadHandler(
        filename = function() {
          paste0("peptide_table_full_", Sys.Date(), ".xlsx")
        },
        content = function(file) {
          req(self$rv_sheet$report)

          export_data <- self$rv_sheet$report
          
          writexl::write_xlsx(export_data, path = file)
        }
      )
    }
  )
)