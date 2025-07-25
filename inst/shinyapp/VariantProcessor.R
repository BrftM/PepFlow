library(R6)
library(shiny)
library(pepitope)
library(AnnotationHub)
library(writexl)
library(VariantAnnotation)

library(DT)
library(shinyjs)

VariantProcessor <- R6Class("VariantProcessor",
  public = list(
    ens106 = NULL,
    asm = NULL,
    rv = NULL,
    test_mode_1 = FALSE,

    rv_sheet = reactiveValues(sheet_data = NULL, table_data = NULL),
    dataHandling = NULL,  

    initialize = function(dataHandling, test_mode_1 = FALSE) {
      self$ens106 <- AnnotationHub::AnnotationHub()[["AH100643"]]
      seqlevelsStyle(self$ens106) <- "UCSC"
      self$asm <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
      self$asm@seqinfo@genome[] <- "GRCh38"
      self$rv <- reactiveValues(report_data = list())
      
      self$dataHandling <- dataHandling
      self$test_mode_1 <- test_mode_1
    },

    process_vcf = function(vcf_file = NULL, sample_selected = NULL) {
      runjs("document.getElementById('status_1').innerText = 'Step 2/8 - Read VCF as VRanges...';")

      # Step 1: Load test VCF if needed
      if (self$test_mode_1 || is.null(vcf_file)) {
        vcf_file <- system.file("my_variants.vcf", package = "pepitope")
        sample_selected <- NULL
      }

      # Step 3: Convert to VRanges and filter variants
      vr1 <- tryCatch({
        vr1 <- pepitope::readVcfAsVRanges(vcf_file) |>
                  pepitope::filter_variants(min_cov = 2, min_af = 0.05, pass = TRUE, sample=sample_selected)
        vr1
      }, error = function(e) {
        shinyalert("Error reading VCF or filtering variants", e$message, type = "error")
        return(NULL)
      })

      if (is.null(vr1) || length(vr1) == 0) {
        shinyalert("No variants read", "No variants found for the selected sample.", type = "error")
        return()
      }

      runjs("document.getElementById('status_1').innerText = 'Step 3/8 - Annotate coding and Subset context...';")

      # Extract the list of sequence levels (chromosomes/contigs) from the BSgenome object (e.g., "chr1", "chr2", etc.)
      # This list is fixed and cannot be changed, but can be used to filter other objects.
      genome_seqlevels <- seqlevels(self$asm)

      # ---- Level-fix-part 1: Filter VRanges to match BSgenome ----
      # Keep only the chromosomes that are present in both the VCF-derived VRanges (vr1) and the BSgenome.
      # This ensures that downstream annotations only happen on contigs with known sequence lengths and data.
      vr1 <- keepSeqlevels(
        vr1,
        intersect(seqlevels(vr1), genome_seqlevels),
        pruning.mode = "coarse"  # drops any ranges on excluded chromosomes
      )

      # ---- Level-fix-part 2: Extract gene annotations as GRanges from the EnsDb ----
      # Since you cannot modify EnsDb directly, you pull out the gene models as a GRanges object.
      genes_gr <- GenomicFeatures::genes(self$ens106)

      # Ensure the gene annotations use UCSC-style naming (e.g., "chr1", "chrM" instead of "1", "MT")
      seqlevelsStyle(genes_gr) <- "UCSC"

      # Explicitly set the genome build label to match other objects (GRCh38 in this case)
      genome(genes_gr) <- "GRCh38"

      # Keep only chromosomes present in both the gene annotations and the genome reference
      genes_gr <- keepSeqlevels(
        genes_gr,
        intersect(seqlevels(genes_gr), genome_seqlevels),
        pruning.mode = "coarse"
      )

      # ---- Level-fix-part 3: Ensure vr1 and gene annotations share the exact same chromosomes ----
      # Now that both have been filtered to what's available in the BSgenome,
      # further restrict to chromosomes common to both VRanges and gene annotation data.
      common <- Reduce(intersect, list(seqlevels(vr1), seqlevels(genes_gr)))

      # Apply this reduced set of seqlevels to both objects again
      vr1 <- keepSeqlevels(vr1, common, pruning.mode = "coarse")
      genes_gr <- keepSeqlevels(genes_gr, common, pruning.mode = "coarse")

      # ---- Level-fix-part 4 (Optional but recommended): Trim any variants that exceed known chromosome bounds ----
      # This step removes or shortens any variant positions that lie outside of the valid sequence lengths
      # defined in the Seqinfo object of the genome.
      vr1 <- GenomicRanges::trim(vr1)

      # Step 4: Annotate coding
      ann <- tryCatch({

        ann <- pepitope::annotate_coding(vr1, self$ens106, self$asm)
        if (is.null(ann) || length(ann) == 0) {
          stop("Annotation result is empty.")
        }
        ann
      }, error = function(e) {
        shinyalert("Annotation error", e$message, type = "error")
        return(NULL)
      })

      if (is.null(ann) || length(ann) == 0) return()

      # Step 5: Subset context
      subs <- tryCatch({
            subs = ann |>
              pepitope::subset_context(15)

      }, error = function(e) {
        shinyalert("Subset context error", e$message, type = "error")
        return(NULL)
      })

      if (is.null(subs)) {
        warning("Subset context failed. Halting downstream processing.")
        return()
      }

      runjs("document.getElementById('status_1').innerText = 'Step 4/8 - Tiling cDNAs...';")

      # Step 6: Tiling
      tiled <- tryCatch({
        pepitope::make_peptides(subs) |>
          pepitope::pep_tile() |>
          pepitope::remove_cutsite(BbsI = "GAAGAC")
      }, error = function(e) {
        shinyalert("Tiling error", e$message, type = "error")
        return(NULL)
      })

      if (is.null(tiled)) return()

      runjs("document.getElementById('status_1').innerText = 'Step 5/8 - Make report...';")

      # Step 7: Generate report
      self$rv_sheet$report <- tryCatch({
        pepitope::make_report(vars = ann, subs = subs, tiled = tiled)
      }, error = function(e) {
        shinyalert("Report generation error", e$message, type = "error")
        return(NULL)
      })

      shinyalert(
        title = "Report completed",
        text = "Barcodes can be added in sheet '93 nt Peptides.'",
        type = "success"
      )
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

              # Filter out empty lines 
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
              self$rv_sheet$report[["93 nt Peptides"]] <- df

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
            actionButton("help_btn_1", "Upload info ℹ️", title = "Need help for what to upload?"),
            checkboxInput("use_test_data_1", "Use test data", value = FALSE),
            fileInput("vcf_file", "Upload VCF File (.vcf.gz)", accept = ".vcf.gz"),
            selectInput("sampleNames", "Select sample", choices = NULL),
            actionButton("process", "Process VCF"),
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

        # If test mode was enabled, reset it and show alert
        if (isTRUE(input$use_test_data_1)) {
          updateCheckboxInput(inputId = "use_test_data_1", value = FALSE)

          self$test_mode_1 <- FALSE

          shinyalert(
            title = "Switched to real data",
            text = "A VCF file was uploaded. Test mode is now disabled.",
            type = "info"
          )
        }

        vr <- pepitope::readVcfAsVRanges(input$vcf_file$datapath)
        samples <- levels(sampleNames(vr))
        
        updateSelectInput(inputId = "sampleNames", choices = samples, selected = samples[1])
      })

      observeEvent(input$process, {
         use_test_1 <- isTRUE(input$use_test_data_1)

          if (!use_test_1 && is.null(input$vcf_file)) {
            shinyalert(
              title = "Missing Input",
              text = "Please upload a VCF file or enable 'Use test data'.",
              type = "warning"
            )
            runjs("document.getElementById('status_1').innerText = 'Waiting for input...';")
            return()
          }

        runjs("document.getElementById('status_1').innerText = 'Step 1/8 - Processing VCF file...';")

        if (use_test_1) {
          self$test_mode_1 <- TRUE
          shinyalert(
            title = "Switched to test data",
            text = "Test mode is activated.",
            type = "info"
          )
          self$process_vcf()
        } else {
          req(input$vcf_file)
          self$process_vcf(input$vcf_file$datapath , input$sampleNames)
        }

        runjs("document.getElementById('status_1').innerText = 'Step 6/8 - VCF file processed';")

        self$display_table(output, input)
      })

      observeEvent(input$help_btn_1, {
        showModal(modalDialog(
          title = "📘 Help: Variant Calling Workflow",
          HTML(
            "<div style='line-height: 1.5;'>
              <h4>🔬 1. Upload a VCF File</h4>
              <p>
                Upload a <strong>VCF file (.vcf.gz)</strong> containing patient-specific tumor mutations. 
                This file will be used to extract, annotate, and tile coding variants into peptide sequences.
              </p>
              <p>The VCF should follow standard format conventions (e.g., GRCh38 reference).</p>

              <h4>🧪 2. Use Built-In Test Data</h4>
              <p>
                Alternatively, you can check the <strong>'Use test data'</strong> box to run this workflow 
                on a predefined test VCF file. This is useful for testing the application without uploading your own file.
              </p>

              <h5>Test Data Preview (First 5 entries):</h5>
              <pre style='background:#f8f9fa; border:1px solid #dee2e6; padding:10px; font-size: 0.85em;'>
                ##fileformat=VCFv4.2
                #CHROM  POS        ID           REF  ALT  FORMAT     SAMPLE
                chr1    114713908  NRAS_Q61L    T    A    AD:DP      90,10:100
                chr1    114713908  NRAS_Q61R    T    C    AD:DP      90,10:100
                chr2    177234082  NFE2L2_E79K  G    A    AD:DP      90,10:100
                chr2    198267522  SF3B1_K700E  A    G    AD:DP      90,10:100
                chr2    208248388  IDH1_R132H   G    A    AD:DP      90,10:100
              </pre>

              <h4>⚙️ 3. Trigger the Workflow</h4>
              <ol>
                <li>Either upload a valid <code>.vcf.gz</code> file, or check the 'Use test data' option.</li>
                <li>Click <strong>'Process VCF'</strong> to start the workflow.</li>
                <li>Once complete, you'll be able to view and edit the resulting peptide tables.</li>
                <li>You can add barcodes manually to the <strong>'93 nt Peptides'</strong> sheet before downloading.</li>
              </ol>

              <p><em>Note: All critical steps include status updates and error alerts in case something goes wrong.</em></p>
            </div>"
          ),
          easyClose = TRUE,
          size = "l",
          footer = NULL
        ))
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