#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(readr)
  library(pheatmap)
})

ui <- fluidPage(
  titlePanel("Long RNA-seq Interactive Visualization Dashboard"),
  sidebarLayout(
    sidebarPanel(
      fileInput("de_file", "Upload DE Results (CSV)"),
      numericInput("padj", "Adjusted p-value cutoff:", 0.05, step=0.01),
      selectInput("plot_type", "Plot Type:", choices=c("Volcano", "MA Plot")),
      actionButton("update", "Generate Plot")
    ),
    mainPanel(plotOutput("mainPlot", height="600px"))
  )
)

server <- function(input, output) {
  observeEvent(input$update, {
    req(input$de_file)
    de <- read.csv(input$de_file$datapath)
    output$mainPlot <- renderPlot({
      if (input$plot_type == "Volcano") {
        ggplot(de, aes(x=log2FoldChange, y=-log10(padj), color=padj < input$padj)) +
          geom_point(alpha=0.5) + theme_minimal() +
          labs(title="Volcano Plot", x="log2 Fold Change", y="-log10(padj)")
      } else {
        ggplot(de, aes(x=baseMean, y=log2FoldChange, color=padj < input$padj)) +
          geom_point(alpha=0.5) + theme_minimal() +
          labs(title="MA Plot", x="Mean Expression", y="log2 Fold Change")
      }
    })
  })
}
shinyApp(ui, server)
