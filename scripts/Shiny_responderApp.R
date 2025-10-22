library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(openxlsx)

# Load data
influenza_folds_fullid  <- read.xlsx(here::here("processed", "influenza_folds_fullid.xlsx")) %>% 
  mutate(fold_change = case_when(fold_change >= 40 ~ 45,
                                 TRUE ~ as.double(fold_change)))

# Define UI
ui <- fluidPage(
  titlePanel("Interactive Patient Data Filtering"),
  sidebarLayout(
    sidebarPanel(
      # Dropdown to select the classification category type (d30 or d20)
      selectInput("classification_type", "Select Classification Type:",
                  choices = c("classification_d30", "classification_d20"),
                  selected = "classification_d30"),
      
      # Dropdown to select the classification category (based on selected type)
      uiOutput("classification"),
      
      # Buttons to cycle through the patients
      actionButton("previous_patient", "Previous Patient"),
      actionButton("next_patient", "Next Patient"),
      
      # Button to update the plot
      actionButton("update", "Update Plot")
    ),
    mainPanel(
      # First plot: Fold Change Plot
      plotlyOutput("interactivePlot"),
      
      # Second plot: Baseline Titer Plot
      plotlyOutput("baselineTiterPlot"),
      
      textOutput("patientID"),  # Display the current patient ID
      textOutput("totalPatients")  # Display the total number of patients in selected category
    )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  
  # Reactive value to keep track of the current patient index
  current_patient_index <- reactiveVal(1)
  
  # Dynamically generate classification category options based on selected classification type
  output$classification <- renderUI({
    classification_column <- input$classification_type  # Select column based on the input
    
    # Get the list of unique patient IDs from the data, filtered by the selected classification category
    selectInput("classification", "Select Classification Category:",
                choices = unique(influenza_folds_fullid[[classification_column]]),
                selected = unique(influenza_folds_fullid[[classification_column]])[1])
  })
  
  # Get the list of unique patient IDs from the data, filtered by the selected classification category
  filtered_patients <- reactive({
    classification_column <- input$classification_type  # Select column based on input
    patient_ids <- influenza_folds_fullid %>%
      filter(!!sym(classification_column) == input$classification) %>%
      select(Pat_ID) %>%
      distinct() %>%
      pull(Pat_ID)  # Get unique patient IDs for the selected classification
    patient_ids
  })
  
  # Display the total number of patients in the selected classification category
  output$totalPatients <- renderText({
    total_patients <- length(filtered_patients())  # Total number of patients in the selected category
    paste("Total Patients in category '", input$classification, "': ", total_patients)
  })
  
  # Update the patient index when the "Next Patient" button is clicked
  observeEvent(input$next_patient, {
    current_index <- current_patient_index()
    next_index <- ifelse(current_index < length(filtered_patients()), current_index + 1, 1)  # Loop back to first patient after last one
    current_patient_index(next_index)  # Update the reactive value
  })
  
  # Update the patient index when the "Previous Patient" button is clicked
  observeEvent(input$previous_patient, {
    current_index <- current_patient_index()
    prev_index <- ifelse(current_index > 1, current_index - 1, length(filtered_patients()))  # Loop back to last patient if at first
    current_patient_index(prev_index)  # Update the reactive value
  })
  
  # Reactive expression to filter the data based on the selected patient and classification category
  filtered_data <- reactive({
    patient <- filtered_patients()[current_patient_index()]  # Get the current patient based on the index
    classification_column <- input$classification_type  # Select column based on input
    data <- influenza_folds_fullid %>%
      filter(Pat_ID == patient, !!sym(classification_column) == input$classification) %>%
      select(strain, fold_change, !!sym(classification_column), Pat_ID) %>%
      mutate(strain = as.factor(strain))
    
    # Create a new color column based on fold_change (green if fold_change > 3.9)
    data <- data %>%
      mutate(color = ifelse(fold_change > 3.9, "green", "blue"))  # Green for fold_change > 3.9, blue otherwise
    
    data
  })
  
  # Reactive expression to get baseline titer data for the selected patient
  baseline_data <- reactive({
    patient <- filtered_patients()[current_patient_index()]  # Get the current patient based on the index
    classification_column <- input$classification_type  # Select column based on input
    data <- influenza_folds_fullid %>%
      filter(Pat_ID == patient, !!sym(classification_column) == input$classification) %>%
      select(strain, baseline_titer, !!sym(classification_column), Pat_ID) %>%
      mutate(strain = as.factor(strain))
    
    data
  })
  
  # Create the first plot (fold change)
  output$interactivePlot <- renderPlotly({
    req(input$update)  # Wait for the "Update Plot" button to be clicked
    
    # Plotting the filtered data
    p <- filtered_data() %>%
      ggplot(aes(x = strain, y = fold_change, color = color, text = paste("Fold Change:", fold_change, "<br>Strain:", strain))) +
      geom_point(size = 3) +  # Adjust dot size for better visibility
      geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
      geom_hline(yintercept = 2, color = "orange", linetype = "dashed") +
      geom_hline(yintercept = 4, color = "green", linetype = "dashed") +
      coord_cartesian(ylim = c(0, 45)) +  # Fixed y-axis range from 0 to 45
      theme_minimal() +
      scale_color_manual(values = c("green" = "green", "blue" = "blue"))+ # Color mapping for the dots
      theme(legend.position = "none") 
    
    # Convert ggplot to plotly for interactivity
    ggplotly(p, tooltip = "text")
  })
  
  # Create the second plot (baseline titer)
  output$baselineTiterPlot <- renderPlotly({
    req(input$update)  # Wait for the "Update Plot" button to be clicked
    
    # Plotting the baseline titer data
    p2 <- baseline_data() %>%
      ggplot(aes(x = strain, y = baseline_titer, text = paste("Baseline Titer:", baseline_titer, "<br>Strain:", strain))) +
      geom_point(size = 3, color = "blue") +  # Adjust dot size and color for baseline
      coord_cartesian(ylim = c(0, 450)) +  # Fixed y-axis range from 0 to 45
      theme_minimal()
    
    # Convert ggplot to plotly for interactivity
    ggplotly(p2, tooltip = "text")
  })
  
  # Display the current patient's id
  output$patientID <- renderText({
    patient_id <- filtered_patients()[current_patient_index()]
    paste("Current Patient ID:", patient_id)
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
