#!/usr/bin/env Rscript
# Lightweight Shiny app to explore a Seurat spatial object interactively.
# Run: Rscript interactive_spatial_viewer.R --args path/to/object.rds
suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(plotly)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript interactive_spatial_viewer.R path/to/spatial_object.rds")
}
so_path <- args[1]
so <- readRDS(so_path)

ui <- fluidPage(
  titlePanel("Interactive Spatial Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("feature", "Feature (gene or meta):", choices = c(rownames(so), colnames(so@meta.data)), selected = colnames(so@meta.data)[1]),
      numericInput("point_size", "Point size:", value = 1, min = 0.2, step = 0.2),
      checkboxInput("use_image", "Show histology image background (if available)", value = TRUE)
    ),
    mainPanel(
      plotlyOutput("spatialPlot", height = "700px")
    )
  )
)

server <- function(input, output, session) {
  output$spatialPlot <- renderPlotly({
    feat <- input$feature
    p <- if (feat %in% rownames(so)) {
      SpatialFeaturePlot(so, features = feat, image.alpha = ifelse(input$use_image, 1, 0)) + theme(legend.position = "right")
    } else {
      SpatialDimPlot(so, group.by = feat, label = TRUE, image.alpha = ifelse(input$use_image, 1, 0))
    }
    ggplotly(p, tooltip = "text") %>% layout(dragmode = "zoom")
  })
}

# Run the app
shinyApp(ui, server)
