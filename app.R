# Load R packages
library(shiny)
library(bslib)
library(dplyr) # for mutate(), rename()
library(tidyr) # for drop_na()
library(ggplot2) # for ggplot()
library(ggprism)
library(flextable)

github_link <- tags$a(
    shiny::icon("github"), "GitHub",
    href = "https://github.com/mstanley-yo/satolab-infectivity-helper",
    target = "_blank"
)

rm_samples_button <- actionButton(
    "remove_all_samples",
    "Remove All Samples",
    icon = icon("trash")
)

# Define UI #####
ui <- page_navbar(
    title = "Infectivity calculator",
    theme = bs_theme(bootswatch = "flatly"),
    
    # Standard curve tab ----
    nav_panel(
        "Standard curve",
        layout_columns(
            col_widths = c(4, 8),
            card(
                card_header("RLU of serially diluted standard:"),
                numericInput("std1", "1141.1 ng/mL (stock):", NA, min = 0),
                numericInput("std2", "570.6 ng/mL:", NA, min = 0),
                numericInput("std3", "285.3 ng/mL:", NA, min = 0),
                numericInput("std4", "142.6 ng/mL:", NA, min = 0),
                numericInput("std5", "71.3 ng/mL:", NA, min = 0),
                numericInput("std6", "0 ng/mL:", NA, min = 0)
            ),
            card(
                card_header("Standard curve"),
                plotOutput("rlu_plot")
            )
        )
    ),
    
    # HiBiT normalization tab ----
    nav_panel(
        "HiBiT normalization",
        layout_columns(
            col_widths = c(4, 8),
            layout_columns(
                col_widths = c(12, 12),
                navset_card_pill(
                    nav_panel(
                        "Manual input",
                        textInput("sample_id_input", "Sample ID", ""),
                        numericInput("rlu_input", "RLU", NA, min = 0),
                        actionButton(
                            "add_sample", 
                            "Add Sample", 
                            icon = icon("plus")
                        ),
                        rm_samples_button
                    ),
                    nav_panel(
                        "Input from Excel",
                        textInput("excel_input", "Copy from Excel", ""),
                        p("Input as 'sample_id RLU' separated by spaces/tabs."),
                        p(
                            paste(
                                "For example:",
                                "'virus1\t2.57E+04 virus2\t2.92E+04'..."
                            )
                        ),
                        actionButton(
                            "add_sample_excel",
                            "Add Sample from Excel",
                            icon = icon("file-excel")
                        ),
                        rm_samples_button
                    )
                ),
                card(
                    card_header("Dilution settings"),
                    numericInput(
                        "pre_hibit_input",
                        "Dilution before HiBiT assay (times)",
                        value = 4,
                        min = 0
                    ),
                    layout_columns(
                        numericInput(
                            "target_volume_uL", 
                            "Target volume (uL)", 
                            value = 300
                        ),
                        numericInput(
                            "target_p24_ng_mL", 
                            "Target p24 (ng/mL)", 
                            value = 400)
                    )
                )
            ),
            card(
                card_header("Dilution Table"),
                uiOutput("table_output"),
                downloadButton("download_docx", "Download table as .docx"),
                p("Written in R Shiny by Maximilian Stanley Yo."),
                github_link
            )
        )
    )
)

# Define server #####
server <- function(input, output) {
    # Define validate standard curve function
    validate_stdcurve <- function() {
        validate(need(
            coeffs$slope != "" | coeffs$intercept != "",
            "Please set up the standard curve!"
        ))
    }
    
    # Create a reactiveValues list to hold slope and intercept and R-squared
    coeffs <- reactiveValues(intercept = NULL, slope = NULL, r2 = NULL)

    # ggplot - draw to reactive object so you can reuse in .docx output
    plot_obj <- reactive({
        # Extract numeric RLU values from inputs
        rlus <- as.numeric(c(
            input$std1, input$std2, input$std3,
            input$std4, input$std5, input$std6
        ))

        # Known concentrations
        concentrations <- c(1141.1 / 2^(0:4), 0)

        # Build data frame
        df <- data.frame(
            RLU = rlus,
            Concentration = concentrations
        ) %>%
            drop_na()

        # Fit linear model: Concentration as a function of RLU
        model <- lm(Concentration ~ RLU, data = df)

        # Extract coefficients and compute R-squared
        coeffs$intercept <- coef(model)[1]
        coeffs$slope <- coef(model)[2]
        coeffs$r2 <- summary(model)$r.squared

        # Format equation string and R-squared
        eq <- paste0(
            "y = ", 
            round(coeffs$slope, 6), 
            "x + ", 
            round(coeffs$intercept, 3)
        )
        r2_label <- paste0("R² = ", round(coeffs$r2, 3))

        # draw ggplot
        ggplot(df, aes(x = RLU, y = Concentration)) +
            geom_point(size = 3, color = "blue") +
            geom_smooth(
                method = "lm", 
                se = FALSE, 
                color = "grey", 
                linetype = "dashed"
            ) +
            labs(
                x = "RLU (Relative Light Units)",
                y = "p24 concentration (ng/mL)",
                title = "Standard Curve: RLU vs. Concentration",
                subtitle = paste(eq, "     ", r2_label)
            ) +
            theme_prism()
    })

    # render standard curve plot
    output$rlu_plot <- renderPlot({
        # validate that there is at least one RLU value entered
        validate(
            need(
                any(c(
                    input$std1, input$std2, input$std3,
                    input$std4, input$std5, input$std6
                ) != ""),
                "Please input luminescence (RLU) values!"
            )
        )

        # render the ggplot object defined earlier.
        plot_obj()
    }) # renderPlot

    # Create dynamic sample table
    sample_data <- reactiveVal(data.frame(
        sample_id = character(),
        rlu = numeric(),
        p24_concentration_ng_mL = numeric(),
        volume_to_dilute = numeric(),
        stringsAsFactors = FALSE
    ))

    # add sample
    observeEvent(input$add_sample, {
        # validate that there is a standard curve.
        validate_stdcurve()

        # Calculate p24 concentration & volumes
        p24_conc <- (coeffs$slope * input$rlu_input + coeffs$intercept) * 
                    input$pre_hibit_input
        vol_to_dilute <- input$target_volume_uL * input$target_p24_ng_mL / 
                         p24_conc
        
        # Create new row
        new_row <- data.frame(
            sample_id = input$sample_id_input,
            rlu = input$rlu_input,
            p24_concentration_ng_mL = p24_conc,
            volume_to_dilute = vol_to_dilute,
            stringsAsFactors = FALSE
        )

        # Bind to current sample_data
        current <- sample_data()
        sample_data(rbind(current, new_row))
    })

    # add sample from excel
    observeEvent(input$add_sample_excel, {
        # validate that there is a standard curve
        validate_stdcurve()

        # Try to parse samples from excel.
        # Example: "SARS-CoV-2 2.31E+04\tHIV-1,3.89E+03\nMERS-CoV 4.12E+04"
        tryCatch(
            expr = {
                # 1. Normalize all delims 
                # (tabs, commas, multiple spaces → single space)
                excel_input <- input$excel_input
                excel_input <- gsub("[,\\t]+", " ", excel_input)
                excel_input <- gsub("\\s+", " ", excel_input)
                
                # 2. Insert newline after each numeric value
                excel_input <- gsub(
                    "([0-9.Ee+-]+)\\s+(?=[A-Za-z0-9._-])", 
                    "\\1\n", 
                    excel_input, 
                    perl = TRUE
                )
                
                # 3. Replace single space between name and number with a tab
                excel_input <- gsub(
                    "([A-Za-z0-9._-]+)\\s+([0-9.Ee+-]+)", 
                    "\\1\t\\2", 
                    excel_input
                )
                
                # 4. Trim whitespace
                excel_input <- trimws(excel_input)
                
                # 5. Read into a data frame.
                df_in <- read.table(
                    text = excel_input,
                    header = FALSE, sep = "\t", stringsAsFactors = FALSE
                )
                names(df_in) <- c("sample_id", "rlu")
                
                # Calculate p24 concentration & volumes vectorised
                p24_conc <- (coeffs$slope * df_in$rlu + coeffs$intercept) * 
                    input$pre_hibit_input
                vol_to_dilute <- input$target_volume_uL * input$target_p24_ng_mL / 
                    p24_conc
                
                # Create new rows
                new_rows <- data.frame(
                    sample_id = df_in$sample_id,
                    rlu = df_in$rlu,
                    p24_concentration_ng_mL = p24_conc,
                    volume_to_dilute = vol_to_dilute,
                    stringsAsFactors = FALSE
                )
                
                # Bind to current sample_data
                current <- sample_data()
                sample_data(rbind(current, new_rows))
            },
            error = function(e) {
                showNotification(
                    paste(
                        "Error parsing Excel input.",
                        "Please ensure the format is 'sample_id RLU'",
                        "separated by spaces, tabs, or commas."
                    ),
                    type = "error",
                    duration = 8
                )
                
                # return empty
                return(NULL)
            }
        )
    })

    # remove all samples
    observeEvent(input$remove_all_samples, {
        sample_data(data.frame(
            sample_id = character(),
            rlu = numeric(),
            p24_concentration_ng_mL = numeric(),
            volume_to_dilute = numeric(),
            stringsAsFactors = FALSE
        ))
    })
    
    # define round volumes function
    round_volumes <- function(vol) {
        case_when(
            vol < 20 ~ round(vol, 2), # p20
            vol > 200 ~ round(vol, 0), # p1000
            .default = round(vol, 1) # p200
        )
    }

    # Create flextable for output in both pretty output and download
    table_obj <- reactive({
        # Extract input values
        total_volume <- input$target_volume_uL
        target_p24 <- input$target_p24_ng_mL

        # Get current data
        df <- sample_data()

        # Round numeric columns for display and rename
        df_display <- df %>%
            mutate(
                p24_concentration_ng_mL = round(p24_concentration_ng_mL, 2),
                volume_to_dilute = round_volumes(volume_to_dilute)
            ) %>%
            rename(
                `ID` = sample_id,
                `RLU` = rlu,
                `p24\n(ng/mL)` = p24_concentration_ng_mL,
                `volume\n(uL)` = volume_to_dilute
            )

        # create flextable
        flextable(df_display) %>%
            autofit(add_w = 100) %>%
            add_header_lines(
                values = paste0(
                    "Pseudovirus infectivity assay (", 
                    Sys.Date(), 
                    ") dilution table"
                )
            ) %>%
            add_footer_lines(
                values = paste0(
                    "Dilute each sample up to ", total_volume, " uL to reach ",
                    target_p24, " ng/mL p24 concentration."
                )
            ) %>%
            align(
                j = c("p24\n(ng/mL)", "volume\n(uL)"), 
                align = "left", 
                part = "all"
            ) %>%
            bold(bold = TRUE, part = "header") %>%
            bold(j = "volume\n(uL)", bold = TRUE, part = "body") %>%
            color(j = "volume\n(uL)", color = "red", part = "body") %>%
            fontsize(size = 14, part = "all")
    })

    # Render table_obj() as HTML object for pretty output
    output$table_output <- renderUI({
        # validate that there is a standard curve.
        validate_stdcurve()

        # Render as HTML widget
        tagList(
            flextable::htmltools_value(table_obj())
        )
    })

    # Add ggplot in the footer and render as .docx for download.
    output$download_docx <- downloadHandler(
        filename = function() {
            paste0("dilution_table_", Sys.Date(), ".docx")
        },
        content = function(file) {
            # Save ggplot to temporary file to insert later
            tmpfile <- tempfile(fileext = ".png")
            ggsave(
                tmpfile, 
                plot = plot_obj(), 
                width = 7, 
                height = 4, 
                dpi = 150
            )
            
            # Add ggplot as image
            ft <- table_obj() %>%
                add_footer_lines(
                    values = as_paragraph(
                        as_image(src = tmpfile, width = 7, height = 4) 
                    )
                )
            
            # Format for docx.
            ft <- ft %>%
                autofit() %>%
                width(width = dim(.)$widths * 8 / (flextable_dim(.)$widths))

            # Save as docx
            flextable::save_as_docx(ft, path = file)
        }
    )
} # server


# Create Shiny object #####
shinyApp(ui = ui, server = server)
