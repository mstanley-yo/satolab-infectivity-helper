# Load R packages
library(shiny)
library(shinythemes)
library(tidyverse)
library(ggprism)
library(flextable)

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
  navbarPage(
    "Infectivity calculator",
    tabPanel("Standard curve",
             sidebarPanel(
               tags$h4("RLU of serially diluted standard:"),
               numericInput("std1", "1141.1 ng/mL (stock):", value = NA, min = 0),
               numericInput("std2", "570.6 ng/mL:", value = NA, min = 0),
               numericInput("std3", "285.3 ng/mL:", value = NA, min = 0),
               numericInput("std4", "142.6 ng/mL:", value = NA, min = 0),
               numericInput("std5", "71.3 ng/mL:", value = NA, min = 0),
               numericInput("std6", "0 ng/mL:", value = NA, min = 0)
               
             ), # sidebarPanel
             mainPanel(
               h3("Standard curve"),
               h5("Typical slope is 0.0002 ~ 0.0003"),
               h5("Typical intercept is -100 ~ 0"),
               plotOutput("rlu_plot"),

             ) # mainPanel
    ), # Navbar 1, tabPanel
    tabPanel("HiBiT normalization",
             sidebarPanel(
             textInput("sample_id_input", "Sample ID", value = ""),
             numericInput("rlu_input", "RLU", value = NA, min = 0),
             numericInput("pre_hibit_input", "Dilution before HiBiT assay (times)", value = 4, min = 0),
             actionButton("add_sample", "Add Sample"),
             br(), br(), br(),
             
             h4("Dilution settings"),
             numericInput("target_volume_uL", "Target volume (uL)", value = 300),
             numericInput("target_p24_ng_mL", "Target p24 (ng/mL)", value = 400),
             br(), br(), br(),
             
             actionButton("remove_all_samples", "Remove All Samples")
             
             ), # sidebarPanel
             mainPanel(
             tableOutput("sample_table"),
             
             ) # mainPanel
             
    ), # Navbar, tabPanel2
    tabPanel("Pretty output",
             mainPanel(
             uiOutput("pretty_output")
             )
    ), # Navbar, tabPanel3
    
    tabPanel("About",
             mainPanel(
               h3("About this website"),
               p("Written in R Shiny by Maximilian Stanley Yo."),
               p(
                 "Follow development here: ",
                 tags$a("GitHub Repository", href = "https://github.com/mstanley-yo/satolab-infectivity-helper", target = "_blank")
               )
             )
    ) # Navbar, tabPanel4
    
  ) # navbarPage
) # fluidPage

  
# Define server function  
server <- function(input, output) {
  
  # Create a reactiveValues list to hold slope and intercept and R-squared
    coeffs <- reactiveValues(intercept = NULL, slope = NULL, r2 = NULL)
    
  # --- Standard curve plot ---
  output$rlu_plot <- renderPlot({
    
    # validate that there is at least one RLU value entered
    validate(
      need(
        any(c(input$std1, input$std2, input$std3, 
              input$std4, input$std5, input$std6) != ""),
        'Please input luminescence (RLU) values!'
      )
    )

    # Extract numeric RLU values from inputs
    rlus <- as.numeric(c(input$std1, input$std2, input$std3, 
                         input$std4, input$std5, input$std6))
    
    # Known concentrations
    concentrations <- c(1141.1/2^(0:4), 0)
    
    # Build data frame
    df <- data.frame(
      RLU = rlus,
      Concentration = concentrations) %>%
      drop_na()
    
    # Fit linear model: Concentration as a function of RLU
    model <- lm(Concentration ~ RLU, data = df)
    
    # Extract coefficients and compute R-squared
    coeffs$intercept <- coef(model)[1]
    coeffs$slope <- coef(model)[2]
    coeffs$r2 <- summary(model)$r.squared
    
    # Format equation string and R-squared
    eq <- paste0("y = ", round(coeffs$slope, 6), "x + ", round(coeffs$intercept, 3))
    r2_label <- paste0("R² = ", round(coeffs$r2, 3))
    
    # ggplot
    ggplot(df, aes(x = RLU, y = Concentration)) +
      geom_point(size = 3, color = "blue") +
      geom_smooth(method = "lm", se = FALSE, color = "grey", linetype = "dashed") +
      labs(
        x = "RLU (Relative Light Units)",
        y = "p24 concentration (ng/mL)",
        title = "Standard Curve: RLU vs. Concentration",
        subtitle = paste(eq, "   ", r2_label)
      ) +
      theme_prism()
    
  }) #renderPlot
  
  # --- Dynamic sample table ---
  sample_data <- reactiveVal(data.frame(sample_id = character(),
                                        rlu = numeric(),
                                        p24_concentration_ng_mL = numeric(),
                                        volume_to_dilute = numeric(),
                                        stringsAsFactors = FALSE))
  
  observeEvent(input$add_sample, {
    
    # validate that there is a standard curve. 
    validate(need(coeffs$slope != "" | coeffs$intercept != "", "Please set up the standard curve!"))
    
    current <- sample_data()
    
    p24_conc <- (coeffs$slope * input$rlu_input + coeffs$intercept) * input$pre_hibit_input
    vol_to_dilute = input$target_volume_uL * input$target_p24_ng_mL / p24_conc
    
    new_row <- data.frame(
      sample_id = input$sample_id_input,
      rlu = input$rlu_input,
      p24_concentration_ng_mL = p24_conc,
      volume_to_dilute = vol_to_dilute,
      stringsAsFactors = FALSE
    )
    
    sample_data(rbind(current, new_row))
  }) # observeEvent - add_sample
  
  observeEvent(input$remove_all_samples, {
    sample_data(data.frame(
      sample_id = character(),
      rlu = numeric(),
      p24_concentration_ng_mL = numeric(),
      volume_to_dilute = numeric(),
      stringsAsFactors = FALSE
    ))
  }) # observeEvent - remove_all_samples
  
  output$sample_table <- renderTable({
    
    # validate that there is a standard curve. 
    validate(need(coeffs$slope != "" | coeffs$intercept != "", "Please set up the standard curve!"))
    
    # rename column headers and render table
    sample_data() %>%
      rename(
        ID = sample_id,
        `Measured RLU` = rlu,
        `Estimated p24 (ng/mL)` = p24_concentration_ng_mL,
        `Volume of stock to dilute (uL)` = volume_to_dilute
      )
  }) # renderTable
  
  output$pretty_output <- renderUI({
    
    # validate that there is a standard curve. 
    validate(need(coeffs$slope != "" | coeffs$intercept != "", "Please set up the standard curve!"))
    
    #req(nrow(sample_data()) > 0)  # only show when table is not empty
    
    # Extract input values
    total_volume <- input$target_volume_uL
    target_p24 <- input$target_p24_ng_mL
    
    # Get current data
    df <- sample_data()
    
    # Optional: round numeric columns for display
    df_display <- df %>%
      mutate(
        p24_concentration_ng_mL = round(p24_concentration_ng_mL, 2),
        volume_to_dilute = round(volume_to_dilute, 1)
      ) %>%
      rename(
        `Sample ID` = sample_id,
        `RLU` = rlu,
        `p24 (ng/mL)` = p24_concentration_ng_mL,
        `volume (uL)` = volume_to_dilute
      )
    
    # Build flextable
    ft <- flextable(df_display) %>%
      autofit(add_w = 100) %>%
      add_footer_lines(
        values = as_paragraph(
          paste0("Dilute each sample up to ", total_volume, " uL to reach ",
                 target_p24, " ng/mL p24 concentration.")
        )
      ) %>%
      align(j = c("p24 (ng/mL)", "volume (uL)"), align = "left", part = "all") %>%
      bold(bold = TRUE, part = "header") %>%
      bold(j = "volume (uL)", bold = TRUE, part = "body") %>%
      color(j = "volume (uL)", color = "red", part = "body")
    
    # Render as HTML widget
    tagList(
      h4("Dilution Table"),
      flextable::htmltools_value(ft)
    )
  })
  
} #server

# Create Shiny object
shinyApp(ui = ui, server = server)
