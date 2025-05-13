library(R6)
library(shiny)
library(pepitope)
library(AnnotationHub)
library(writexl)

library(DT)
library(shinyjs)

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
      seqlevelsStyle(self$asm) <- "UCSC"
      self$asm@seqinfo@genome[] <- "GRCh38"
      self$rv <- reactiveValues(report_data = list(), sheet_names = NULL)
      
      self$dataHandling <- dataHandling
    },

    process_vcf = function(vcf_file) {
      runjs("document.getElementById('status_1').innerText = 'Step 1/11 - Read VcF as VRanges...';")


      vr1 <- pepitope::readVcfAsVRanges(vcf_file) |>
        pepitope::filter_variants(min_cov=2, min_af=0.05, pass=TRUE, chrs="default")

      runjs("document.getElementById('status_1').innerText = 'Step 2/11 - Annotate coding...';")
      ann <- pepitope::annotate_coding(vr1, self$ens106, self$asm)  

      runjs("document.getElementById('status_1').innerText = 'Step 4/11 - Subset context...';")
      subs <- ann |>
        pepitope::subset_context(15)

      runjs("document.getElementById('status_1').innerText = 'Step 5/11 Tiling cDNAs of interest into smaller peptides ...';")
      tiled <- pepitope::make_peptides(subs) |>
        pepitope::pep_tile() |>
        pepitope::remove_cutsite(BbsI="GAAGAC")


      runjs("document.getElementById('status_1').innerText = 'Step 6/11 - Make report...';")
      report <- tryCatch({
        pepitope::make_report(vars=ann, subs=subs, tiled=tiled)
      }, error = function(e) {
        warning(paste("Error in make_report:", e$message))
          shinyalert(
                title = "Failed to generate report!", 
                text = paste("Error: Failed to generate report: ", paste(e$message)),
                type = "error"
          )
        return(NULL)
      })

      if (is.null(report)) return("Report processing failed.")

      runjs("document.getElementById('status_1').innerText = 'Step 7/11 - Convert report to tmp XLSX...';")
      temp_file <- self$dataHandling$writeTMPxlsx(report)
      runjs("document.getElementById('status_1').innerText = 'Step 8/11 - Reading tmp XLSX...';")
      # Now read from that temp file
      self$read_xlsx_sheet(temp_file, "93 nt Peptides")


      self$prepare_sheet_data()
      runjs("document.getElementById('status_1').innerText = 'Step 9/11 - Add barcodes to download';")
    },

    read_xlsx_sheet = function(xlsx_file, sheet_name) {
      req(xlsx_file, sheet_name)

      tryCatch({
        available_sheets <- excel_sheets(xlsx_file)
        self$rv_sheet$barcode_sheet_name <- sheet_name

        all_sheets_data <- lapply(available_sheets, function(sheet) {
          read_xlsx(xlsx_file, sheet = sheet)
        })
        names(all_sheets_data) <- available_sheets
        self$rv_sheet$all_sheets <- all_sheets_data

        if (sheet_name %in% available_sheets) {
          self$rv_sheet$barcode_sheet <- all_sheets_data[[sheet_name]]
        } else {
          self$rv_sheet$barcode_sheet <- NULL
          warning(paste("Sheet", sheet_name, "not found."))
        }
      }, error = function(e) {
        warning(paste("Error reading XLSX sheet:", e$message))
      })
    },

    prepare_sheet_data = function() {
      req(self$rv_sheet$barcode_sheet)
      self$rv_sheet$barcode_table_data <- self$rv_sheet$barcode_sheet
      self$rv_sheet$barcode_table_data$barcode <- ""
    },

    display_table = function(output, input) {
      output$dynamic_table <- renderUI({
        req(self$rv_sheet$all_sheets)

        tab_list <- lapply(names(self$rv_sheet$all_sheets), function(sheet_name) {
          is_editable <- sheet_name == self$rv_sheet$barcode_sheet_name
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

        do.call(tabsetPanel, c(tab_list, id = "sheet_tabs"))
      })

      lapply(names(self$rv_sheet$all_sheets), function(sheet_name) {
        local({
          sheet <- sheet_name
          ns <- NS(sheet)
          is_editable <- sheet == self$rv_sheet$barcode_sheet_name

          output[[ns("datatable")]] <- renderDT({
            df <- if (is_editable) self$rv_sheet$barcode_table_data else self$rv_sheet$all_sheets[[sheet]]
            datatable(df, editable = FALSE, options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE), rownames = FALSE)
          })

          if (is_editable) {
            observeEvent(input[[ns("add_barcode")]], {
              req(input[[ns("barcode_input")]], self$rv_sheet$barcode_table_data)

              barcode_list <- unlist(strsplit(input[[ns("barcode_input")]], "[,\n]+"))
              barcode_list <- trimws(barcode_list)
              barcode_list <- barcode_list[barcode_list != ""]

              if (length(barcode_list) != nrow(self$rv_sheet$barcode_table_data)) {
                error_msg <- sprintf(
                  "Error - %d barcodes provided, but %d rows are needed.",
                  length(barcode_list),
                  nrow(self$rv_sheet$barcode_table_data)
                )
                runjs(sprintf("document.getElementById('status_1').innerText = '%s';", error_msg))
                return()
              }

              self$rv_sheet$barcode_table_data$barcode <- barcode_list
              runjs("document.getElementById('status_1').innerText = 'Step 10/11 - Barcodes added!';")

              output[[ns("datatable")]] <- renderDT({
                datatable(self$rv_sheet$barcode_table_data, editable = FALSE, options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE), rownames = FALSE)
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

        runjs("document.getElementById('status_1').innerText = 'Step 11/11 - Processing VCF file...';")
        self$process_vcf(input$vcf_file$datapath)
        shinyalert(
            title = "Count completed", 
            text = paste("Have fun checking the results! "),
            type = "success"
        )
        self$display_table(output, input)
      })
      
      output$download_peptide_table <- downloadHandler(
        filename = function() {
          paste0("peptide_table_full_", Sys.Date(), ".xlsx")
        },
        content = function(file) {
          req(self$rv_sheet$all_sheets, self$rv_sheet$barcode_table_data, self$rv_sheet$barcode_sheet_name)

          export_data <- self$rv_sheet$all_sheets
          export_data[[self$rv_sheet$barcode_sheet_name]] <- self$rv_sheet$barcode_table_data

          writexl::write_xlsx(export_data, path = file)
        }
      )
    }
  )
)