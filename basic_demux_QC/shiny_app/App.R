# based on script: basic_demux_QC.Rmd
# SP@NC; 2023-07-14 v1.0

library("shiny")
library("shinyBS")
library("readr")
library("ggplot2")
library("gridExtra")
library("grid")
library("gridGraphics")
library("DT")

ui <- fluidPage(
  titlePanel("Read Count variability of demultiplexed Samples"),
  
  sidebarLayout(
    sidebarPanel(
      tipify(fileInput("file", "Choose CSV File", accept = ".csv"), "select at Nucleomics Core poolingQC-FinalTable.csv file"),
      br(),
      tipify(radioButtons("center", "Choose Centering Method:",
                          choices = c("Mean", "Median"),
                          selected = "Mean"), "choose the centering method to compute the count variability"),
      br(),
      tipify(textInput("expnum", "Enter Expnum:"), "if more than one project is reported in the CSV file, type the correct one"),
      actionButton("startAnalysis", "Start Analysis")
    ),
    
    mainPanel(
      plotOutput("plot"),
      dataTableOutput("table"),
      downloadButton("downloadPlot", "Download Plot as PNG")
    )
  )
)

calculate_diff <- function(subset_data, center) {
  count_column <- subset_data$Total
  mean_value <- round(mean(count_column), 1)
  median_value <- round(median(count_column), 1)
  center_value <- ifelse(center == "Mean", mean_value, median_value)
  count_diff <- count_column - center_value
  percentage_diff <- (count_diff / center_value) * 100
  
  mean_diff <- round(mean(abs(percentage_diff)), 2)
  
  return(list(count_diff = count_diff, percentage_diff = percentage_diff, mean_diff = mean_diff, mean_value = mean_value, median_value = median_value, center_value = center_value))
}

server <- function(input, output, session) {
  options(device = "quartz")
  
  data <- reactive({
    req(input$file)
    read_csv(input$file$datapath, show_col_types = FALSE)
  })
  
  observeEvent(input$file, {
    filename <- basename(input$file$name)
    expnum <- strsplit(filename, "_")[[1]][1]
    updateTextInput(session, "expnum", value = expnum)
  })
  
  analysis <- eventReactive(input$startAnalysis, {
    req(input$center, input$expnum)
    
    center <- input$center
    expnum <- as.numeric(input$expnum)
    
    subset_data <- subset(data(), Project == expnum)
    
    diff_data <- calculate_diff(subset_data, center)
    count_diff <- diff_data$count_diff
    percentage_diff <- diff_data$percentage_diff
    mean_diff <- diff_data$mean_diff
    center_value <- diff_data$center_value
    mean_value <- diff_data$mean_value
    median_value <- diff_data$median_value
    
    # Create bar plot
    bar_plot <- ggplot(subset_data, aes(x = Sample, y = Total)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(x = NULL, y = "Total read counts") +
      ggtitle("A") +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_flip(ylim = c(0, max(subset_data$Total))) +
      geom_hline(yintercept = mean_value, linetype = "dashed", color = "red") +
      geom_hline(yintercept = median_value, color = "blue")
    
    # Create box plot
    box_plot <- ggplot(subset_data, aes(x = "Sample", y = Total)) +
      geom_boxplot(fill = "lightgray", color = "black") +
      labs(y = NULL, x = "Total read counts") +
      ggtitle("B") +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            panel.background = element_blank(),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
    
    # Create diff box plot
    box_plot_diff <- ggplot(subset_data, aes(x = "Sample", y = percentage_diff)) +
      geom_boxplot(fill = "lightgray", color = "black") +
      labs(y = NULL, x = paste("% difference from the", center)) +
      ggtitle("C") +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            panel.background = element_blank(),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
    
    # Create text grob
    text_grob <- textGrob(
      label = paste(
        "Legend:\n",
        "* A,B Total read counts\n",
        "* C: Boxplot of count difference from the", center, "\n\n",
        "Sample count: ", nrow(subset_data), "\n",
        "Total read count: ", sum(subset_data$Total), "\n",
        center, "read count: ", center_value, "\n",
        "mean %Difference from", center, " (Absolute Value):", mean_diff, "%", "\n",
        "mean: red dashed line", "\n",
        "median: blue line"
      ),
      x = 0,  # Adjust the x position of the text
      y = 1,  # Adjust the y position of the text
      just = c("left", "top"),  # Set alignment to "left" and "top"
      gp = gpar(fontsize = 9)
    )
    
    # Create an empty data frame and blank plot
    blank <- grid.rect(gp = gpar(col = "white"))
    
    # Combine plots using gridExtra
    combined_plot <- grid.arrange(bar_plot, box_plot, box_plot_diff, text_grob, blank, blank, nrow = 2, ncol = 3, widths = c(8, 2, 2), heights = c(6, 3.5))
    
    return(combined_plot)
  })
  
  output$plot <- renderPlot({
    req(analysis())
    analysis()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      plot <- analysis()
      ggsave(file, plot, device = "png", dpi = 300)
    }
  )
  
  output$table <- renderDataTable({
    req(data())
    datatable(data(), options = list(pageLength = 10, columnDefs = list(list(targets = 0, visible = FALSE)))) %>%
      formatStyle(names(data()), fontSize = "10px")
  })
}

shinyApp(ui, server)