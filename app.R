library(shiny)
library(ggplot2)
library(DT)
library(git2r)
library(readr)
library(dplyr)
library(tidyr)
library(reticulate)

use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
print(py_config())
options(shiny.maxRequestSize = 20 * 1024^2)

# --- Paths ---
risk_model_dir <- "/models/tallsorts"

# --- Load Gene Expression Data ---
gene_expr <- read.table("/srv/shiny-server/KFTALL_1026_samples_gene_tpm.GeneNames.txt", header = TRUE, as.is = TRUE)
umap_data <- read.table("/srv/shiny-server/KFTALL_1026_GeneExpression_batchCorrected_vsd.top5K.UMAP.txt", header = TRUE)

gene_names <- colnames(gene_expr)
umap_data <- umap_data[match(gene_names, rownames(umap_data)), ]
stopifnot(all(gene_names == rownames(umap_data)))

merged_data <- as.data.frame(t(gene_expr))
colnames(merged_data) <- rownames(gene_expr)
merged_data$cancer_subtype <- umap_data[, 3]

subtype_means_kftall <- merged_data %>%
  pivot_longer(-cancer_subtype, names_to = "gene", values_to = "expression") %>%
  group_by(gene, cancer_subtype) %>%
  summarize(mean_expr = mean(expression, na.rm = TRUE), .groups = 'drop')

# CCLE
gene_ccle_raw <- read.table("/srv/shiny-server/CCLE_CellLines.log2TPM.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rownames(gene_ccle_raw) <- make.unique(rownames(gene_ccle_raw))
ccle_expr <- gene_ccle_raw

# TCellALL
tcell_expr_raw <- read.table("/srv/shiny-server/TcellALLData_filtered_TPM.txt", header = TRUE, stringsAsFactors = FALSE)
tcell_expr_raw$ensembl_gene_id <- gsub("\\..*", "", rownames(tcell_expr_raw))

tcell_annotations <- readRDS("/srv/shiny-server/tcell_annotations.rds")
tcell_annotations <- subset(tcell_annotations, ensembl_gene_id %in% tcell_expr_raw$ensembl_gene_id & hgnc_symbol != "")
tcell_merged <- merge(tcell_annotations, tcell_expr_raw, by = "ensembl_gene_id")
rownames(tcell_merged) <- make.unique(as.character(tcell_merged$hgnc_symbol))
tcell_expr <- tcell_merged[, !(names(tcell_merged) %in% c("ensembl_gene_id", "hgnc_symbol"))]
tcell_data <- as.data.frame(t(tcell_expr))
rownames(tcell_data) <- colnames(tcell_expr)

sample_info <- read.table("/srv/shiny-server/sample_selected.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
subtype_vector <- setNames(sample_info$Subtype, gsub("\\.txt$", "", sample_info$Sample))
rownames(tcell_data) <- gsub("\\.txt$", "", rownames(tcell_data))
tcell_data$cancer_subtype <- subtype_vector[rownames(tcell_data)]
tcell_data <- subset(tcell_data, !is.na(cancer_subtype))

# --- UI ---
ui <- fluidPage(
  titlePanel("Gene Expression & Model Results"),
  tabsetPanel(
    tabPanel("KFTALL Data", sidebarLayout(
      sidebarPanel(
        selectInput("gene", "Select Gene:", choices = rownames(gene_expr)),
        selectInput("selected_subtype", "Select Cancer Subtype for Mean Line:", choices = unique(merged_data$cancer_subtype)),
        actionButton("plot_button", "Plot Gene Expression")
      ),
      mainPanel(plotOutput("boxPlot"), dataTableOutput("table"))
    )),
    tabPanel("CCLE Data", sidebarLayout(
      sidebarPanel(
        selectizeInput("gene_ccle", "Select Gene:", choices = NULL),
        actionButton("plot_button_ccle", "Plot Gene Expression")
      ),
      mainPanel(plotOutput("boxPlot_ccle"), dataTableOutput("table_ccle"))
    )),
    tabPanel("TCellALL Data", sidebarLayout(
      sidebarPanel(
        selectInput("gene_tcell", "Select Gene:", choices = rownames(tcell_expr)),
        selectInput("selected_subtype_tcell", "Select Cancer Subtype for Mean Line:", choices = unique(na.omit(tcell_data$cancer_subtype))),
        actionButton("plot_button_tcell", "Plot Gene Expression")
      ),
      mainPanel(plotOutput("boxPlot_tcell"), dataTableOutput("table_tcell"))
    )),
    tabPanel("TALLsorts Local Analysis", sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Upload CSV File", accept = ".csv"),
        actionButton("run_analysis", "Run TALLsorts")
      ),
      mainPanel(
        textOutput("status"),
        h4("Predictions"), tableOutput("contents"),
        plotOutput("tallsorts_plot"),
        h4("Probabilities"), dataTableOutput("probabilities_table"),
        h4("Report"), uiOutput("tallsorts_html")
      )
    )),
    tabPanel("Risk Model Analysis", sidebarLayout(
      sidebarPanel(
        fileInput("risk_file", "Upload CSV for Risk Model", accept = ".csv"),
        actionButton("run_risk", "Run Risk Model")
      ),
      mainPanel(
        textOutput("risk_status"),
        h4("Predictions"), tableOutput("risk_contents"),
        plotOutput("risk_plot"),
        h4("Probabilities"), dataTableOutput("risk_probabilities_table"),
        h4("Report"), uiOutput("risk_html")
      )
    ))
  )
)

# --- Server ---
server <- function(input, output, session) {
  updateSelectizeInput(session, "gene_ccle", choices = rownames(ccle_expr), server = TRUE)
  
  observeEvent(input$plot_button, {
    req(input$gene, input$selected_subtype)
    gd <- data.frame(gene_expr = merged_data[[input$gene]], cancer_subtype = merged_data$cancer_subtype)
    mv <- mean(gd$gene_expr[gd$cancer_subtype == input$selected_subtype], na.rm = TRUE)
    
    output$boxPlot <- renderPlot({
      ggplot(gd, aes(cancer_subtype, gene_expr, fill = cancer_subtype)) +
        geom_boxplot() +
        geom_hline(yintercept = mv, color = "red", linetype = "dashed") +
        labs(title = paste("Gene Expression for", input$gene), x = "Cancer Subtype", y = "Expression") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    output$table <- renderDataTable(gd)
  })
  
  observeEvent(input$plot_button_ccle, {
    req(input$gene_ccle)
    ced <- data.frame(cell_line = colnames(ccle_expr), expression = as.numeric(ccle_expr[input$gene_ccle, ]))
    
    output$boxPlot_ccle <- renderPlot({
      ggplot(ced, aes(reorder(cell_line, -expression), expression)) +
        geom_bar(stat = "identity") +
        labs(title = paste("Expression of", input$gene_ccle), x = "Cell Line", y = "Log2 TPM") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
    })
    
    output$table_ccle <- renderDataTable(ced)
  })
  
  observeEvent(input$plot_button_tcell, {
    req(input$gene_tcell, input$selected_subtype_tcell)
    gtd <- data.frame(gene_expr = tcell_data[[input$gene_tcell]], cancer_subtype = tcell_data$cancer_subtype)
    mv2 <- mean(gtd$gene_expr[gtd$cancer_subtype == input$selected_subtype_tcell], na.rm = TRUE)
    
    output$boxPlot_tcell <- renderPlot({
      ggplot(gtd, aes(cancer_subtype, gene_expr, fill = cancer_subtype)) +
        geom_boxplot() +
        geom_hline(yintercept = mv2, color = "red", linetype = "dashed") +
        labs(title = paste("Gene Expression for", input$gene_tcell), x = "Cancer Subtype", y = "log2TPM + 1") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    })
    
    output$table_tcell <- renderDataTable(gtd)
  })
  
  observeEvent(input$run_analysis, {
    req(input$file1)
    output$status <- renderText("Running analysis…")
    
    tmp <- file.path(tempdir(), "tallsorts_tmp")
    dir.create(tmp, showWarnings = FALSE)
    
    inFile <- file.path(tmp, basename(input$file1$name))
    file.copy(input$file1$datapath, inFile, overwrite = TRUE)
    
    outDir <- file.path(tmp, paste0("output_", as.integer(Sys.time())))
    dir.create(outDir)
    
    py_code <- sprintf(
      "import subprocess\nsubprocess.run(['TALLSorts', '--samples', %s, '--destination', %s], check=True)",
      shQuote(inFile), shQuote(outDir)
    )
    
    result <- tryCatch({
      py_run_string(py_code)
      output$status <- renderText("Complete")
      showNotification("TALLSorts results loaded", type = "message")
    }, error = function(e) {
      output$status <- renderText("Error")
      showNotification(paste("TALLSorts error:", e$message), type = "error")
      return(NULL)
    })
    
    preds <- list.files(outDir, "^predictions\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(preds)) {
      pd <- read_csv(preds[1], show_col_types = FALSE)
      output$contents <- renderTable(head(pd, 20))
      output$tallsorts_plot <- renderPlot({
        ggplot(pd, aes(x = Predictions)) + geom_bar(fill = "#67a9cf") +
          labs(title = "Predicted Subtypes", x = "Subtype", y = "Count")
      })
    }
    
    pr <- list.files(outDir, "^probabilities\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(pr)) {
      pb <- read_csv(pr[1], show_col_types = FALSE)
      output$probabilities_table <- renderDataTable(pb)
    }
    
    output$tallsorts_html <- renderUI({
      report_files <- list.files(outDir, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
      if (!length(report_files)) return(NULL)
      lapply(report_files, includeHTML)
    })
  })
  
  observeEvent(input$run_risk, {
    req(input$risk_file)
    output$risk_status <- renderText("Running risk model…")
    
    tmpR <- file.path(tempdir(), "riskmodel_tmp")
    dir.create(tmpR, showWarnings = FALSE)
    
    inR <- file.path(tmpR, basename(input$risk_file$name))
    file.copy(input$risk_file$datapath, inR, overwrite = TRUE)
    
    outR <- file.path(tmpR, paste0("risk_output_", as.integer(Sys.time())))
    dir.create(outR)
    
    py_code <- sprintf(
      "import subprocess\nsubprocess.run(['TALLSorts', '--samples', %s, '--destination', %s, '--model', %s], check=True)",
      shQuote(inR), shQuote(outR), shQuote(file.path(risk_model_dir, "risk_trained_model.pkl.gz"))
    )
    
    result_risk <- tryCatch({
      py_run_string(py_code)
      output$risk_status <- renderText("Complete")
      showNotification("Risk model results loaded", type = "message")
    }, error = function(e) {
      output$risk_status <- renderText("Error")
      showNotification(paste("Risk model error:", e$message), type = "error")
      return(NULL)
    })
    
    rpreds <- list.files(outR, "^predictions\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(rpreds)) {
      rp <- read_csv(rpreds[1], show_col_types = FALSE)
      output$risk_contents <- renderTable(head(rp, 20))
      output$risk_plot <- renderPlot({
        ggplot(rp, aes(x = Predictions)) + geom_bar(fill = "#67a9cf") +
          labs(title = "Risk Model Predictions", x = "Label", y = "Count")
      })
    }
    
    rpr <- list.files(outR, "^probabilities\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(rpr)) {
      rpb <- read_csv(rpr[1], show_col_types = FALSE)
      output$risk_probabilities_table <- renderDataTable(rpb)
    }
    
    output$risk_html <- renderUI({
      risk_files <- list.files(outR, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
      if (!length(risk_files)) return(NULL)
      lapply(risk_files, includeHTML)
    })
  })
}

options(shiny.host = "0.0.0.0", shiny.port = 3838)
shinyApp(ui = ui, server = server)




