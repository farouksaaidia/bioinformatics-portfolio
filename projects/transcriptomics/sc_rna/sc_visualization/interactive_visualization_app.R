#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(shiny)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input Seurat .rds")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
if (is.null(opt$input)) stop("❌ Provide --input Seurat file")

so <- readRDS(opt$input)
meta_cols <- colnames(so@meta.data)

ui <- fluidPage(
  titlePanel("🧬 Interactive Single-Cell Visualization"),
  sidebarLayout(
    sidebarPanel(
      selectInput("meta", "Color by metadata:", choices=meta_cols, selected="cell_type"),
      textInput("gene", "Gene for expression plot:", ""),
      actionButton("refresh", "Refresh")
    ),
    mainPanel(
      plotOutput("umapPlot", height="600px"),
      plotOutput("featurePlot", height="600px")
    )
  )
)

server <- function(input, output) {
  dataReactive <- reactiveVal(so)
  
  output$umapPlot <- renderPlot({
    DimPlot(dataReactive(), group.by=input$meta, label=TRUE) + ggtitle(paste("UMAP -", input$meta))
  })
  
  output$featurePlot <- renderPlot({
    req(input$gene)
    if (input$gene %in% rownames(dataReactive())) {
      FeaturePlot(dataReactive(), features=input$gene)
    } else {
      ggplot() + ggtitle("Gene not found")
    }
  })
}

shinyApp(ui=ui, server=server)
