# Load R packages
library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggprism)
library(flextable)
library(officer)

##### UI Function #####
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

ui <- page_navbar(
    title = "Infectivity calculator",
    theme = bs_theme(bootswatch = "flatly"),
    
    # Standard curve tab
    nav_panel(
        "Standard curve",
        layout_columns(
            col_widths = c(4, 8),
            navset_card_pill(
                nav_panel(
                    "Manual Input",
                    numericInput("std1", "1141.1 ng/mL (stock):", NA, min = 0),
                    numericInput("std2", "570.6 ng/mL:", NA, min = 0),
                    numericInput("std3", "285.3 ng/mL:", NA, min = 0),
                    numericInput("std4", "142.6 ng/mL:", NA, min = 0),
                    numericInput("std5", "71.3 ng/mL:", NA, min = 0),
                    numericInput("std6", "0 ng/mL:", NA, min = 0)
                ),
                nav_panel(
                    "Input From Excel",
                    textInput("std_excel", "Copy from Excel", ""),
                    p("Input as 6 values separated by spaces/tabs."),
                    p("e.g.: '4.10E+06 2.57E+06 ...'")
                )
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
                        "Manual Input",
                        textInput("sample_id_input", "Sample ID", ""),
                        numericInput("rlu_input", "RLU", NA, min = 0),
                        actionButton("add", "Add Sample", icon = icon("plus")),
                        rm_samples_button
                    ),
                    nav_panel(
                        "Input From Excel",
                        textInput("excel_input", "Copy from Excel", ""),
                        p("Input as 'sample RLU' separated by spaces/tabs."),
                        p("e.g.: 'virus1\t2.57E+04 virus2\t2.92E+04'..."),
                        actionButton(
                            "add_excel",
                            "Add Sample from Excel",
                            icon = icon("file-excel")
                        ),
                        rm_samples_button
                    )
                ),
                card(
                    card_header("Dilution settings"),
                    numericInput(
                        "pre_dilution",
                        "Dilution before HiBiT assay (times)",
                        value = 4,
                        min = 0
                    ),
                    layout_columns(
                        numericInput(
                            "target_vol", 
                            "Target volume (uL)", 
                            value = 300
                        ),
                        numericInput(
                            "target_p24", 
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
            cfs$slope != "" | cfs$intercept != "",
            "Please set up the standard curve!"
        ))
    }
    
    # Create a reactiveValues list to hold slope and intercept and R-squared
    cfs <- reactiveValues(intercept = NULL, slope = NULL, r2 = NULL)

    # ggplot - draw to reactive object so you can reuse in .docx output
    plot_obj <- reactive({
        # Extract numeric RLU values from inputs
        rlus <- str_split_1(input$std_excel, "\\s+")
        if (length(rlus) != 6) {
            rlus <- c(
                input$std1, input$std2, input$std3,
                input$std4, input$std5, input$std6
            )
        }
        
        # Build data frame
        known_concs <- c(1141.1 / 2^(0:4), 0)
        df <- data.frame(RLU = as.numeric(rlus), conc = known_concs) %>%
            drop_na()

        # validate that there is at least one RLU value entered
        validate(
            need(nrow(df) > 0, "Please input luminescence (RLU) values!")
        )
        
        # Fit linear model: Concentration as a function of RLU
        model <- lm(conc ~ RLU, data = df)

        # Extract coefficients and compute R-squared
        cfs$intercept <- coef(model)[1]
        cfs$slope <- coef(model)[2]
        cfs$r2 <- summary(model)$r.squared

        # Format equation string and R-squared
        eq <- paste0(
            "y = ", 
            round(cfs$slope, 6), 
            "x + ", 
            round(cfs$intercept, 3)
        )
        r2_label <- paste0("R² = ", round(cfs$r2, 3))

        # draw ggplot
        ggplot(df, aes(x = RLU, y = conc)) +
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
        plot_obj()
    })

    # Create dynamic sample table sample_data()
    empty_df <- data.frame(
        sample_id = character(),
        rlu = numeric(),
        p24_concentration_ng_mL = numeric(),
        volume_to_dilute = numeric()
    )
    
    sample_data <- reactiveVal(empty_df)
    
    # Function to add to sample_data()
    append_sample_data <- function(data) {
        # validate that there is a standard curve.
        validate_stdcurve()
        
        # Calculate p24 concentration & volumes vectorised
        p24_conc <- (cfs$slope * data$rlu + cfs$intercept) * input$pre_dilution
        vol_to_dilute <- input$target_vol * input$target_p24 / p24_conc
        
        # Create new rows
        new_rows <- data.frame(
            sample_id = data$sample_id,
            rlu = data$rlu,
            p24_concentration_ng_mL = p24_conc,
            volume_to_dilute = vol_to_dilute,
            stringsAsFactors = FALSE
        )
        
        # Bind to current sample_data
        current <- sample_data()
        sample_data(rbind(current, new_rows))
    }
    
    # Parse manual inputs for append_sample_data
    observeEvent(input$add, {
        samples_rlu <- data.frame(
            sample_id = input$sample_id_input,
            rlu = input$rlu_input
        )
        
        # append to sample_data() 
        append_sample_data(samples_rlu)
    })

    # Parse inputs from excel for append_sample_data
    observeEvent(input$add_excel, {
        # Try to parse samples from excel.
        tryCatch(
            expr = {
                # 1. Normalize all delims → single space
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
                samples_rlu <- read.table(
                    text = excel_input,
                    header = FALSE, 
                    sep = "\t", 
                    stringsAsFactors = FALSE
                )
                names(samples_rlu) <- c("sample_id", "rlu")
                
                # append to sample_data()
                append_sample_data(samples_rlu)
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
            }
        )
    })

    # Remove all samples
    observeEvent(input$remove_all_samples, {
        sample_data(empty_df)
    })
    
    # Observe changes in target_vol, target_p24, and pre_dilution and recalc
    observe({
        # Require sample_data()
        req(sample_data())
        req(cfs$slope, cfs$intercept)
        
        # Trigger when any of these inputs change
        input$target_vol
        input$target_p24
        input$pre_dilution
        
        # Get current data
        data <- sample_data()
        
        # Skip if no rows
        validate(need(nrow(data) > 0, "No samples to recalculate."))
        
        # Recalculate p24 concentration and volume
        p24_conc <- (cfs$slope * data$rlu + cfs$intercept) * input$pre_dilution
        vol_to_dilute <- input$target_vol * input$target_p24 / p24_conc
        
        # Update dataset
        data$p24_concentration_ng_mL <- p24_conc
        data$volume_to_dilute <- vol_to_dilute
        
        # Save back
        sample_data(data)
    })
    
    # Define round volumes function
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
        total_volume <- input$target_vol
        target_p24 <- input$target_p24

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

        # Create flextable
        set_flextable_defaults(
            font.family = "Helvetica",
            font.size = 14,
            word_wrap = FALSE
        )
        
        df_display %>%
            flextable() %>%
            autofit(add_w = 0.5) %>%
            add_header_lines(
                values = paste0(
                    "Pseudovirus infectivity assay (", 
                    Sys.Date(), 
                    ") dilution table"
                )
            ) %>%
            add_footer_lines(
                values = as_paragraph(
                    paste(
                        "Dilute each sample up to", 
                        total_volume, 
                        "uL to reach",
                        target_p24, 
                        "ng/mL p24 concentration.",
                        "\n\n"
                    ),
                    paste0(
                        "Samples were diluted ", 
                        input$pre_dilution, 
                        "× prior to HiBiT measurement."
                    )
                )
            ) %>%
            align(align = "left", part = "all") %>%
            bold(bold = TRUE, part = "header") %>%
            bold(j = "volume\n(uL)", bold = TRUE, part = "body") %>%
            color(j = "volume\n(uL)", color = "red", part = "body")
    })

    # Render table_obj() as HTML object for pretty output
    output$table_output <- renderUI({
        # validate that there is a standard curve.
        validate_stdcurve()

        # Render as HTML widget
        tags$div(
            style = "margin-left: 0; margin-right: auto; width: fit-content;",
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
                        "\n",
                        as_image(src = tmpfile, width = 7, height = 4) 
                    )
                )
            
            # Format for docx.
            ft <- ft %>%
                autofit() %>%
                width(width = dim(.)$widths * 7.5 / (flextable_dim(.)$widths))

            # Save as docx
            sect_prop <- prop_section(
                page_size = page_size(
                    orient = "portrait",
                    width = 8.27, 
                    height = 11.69
                ),
                type = "continuous",
                page_margins = page_mar(
                    top = 0.5, 
                    right = 0.5, 
                    bottom = 0.5, 
                    left = 0.5, 
                    header = 0.3, 
                    footer = 0.3, 
                    gutter = 0
                )
            )
            
            flextable::save_as_docx(ft, path = file, pr_section = sect_prop)
        }
    )
} # server


# Create Shiny object #####
shinyApp(ui = ui, server = server)
