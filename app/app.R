library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(bslib)
library(DT)
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(shinyjqui)
library(ProtPQN)
library(shinythemes)
source("functions.R")
library(tidyr)
library(plotly)
library(heatmaply)
library(tibble)
#library(utils)
library(stringr)
#load("shiny_mouse_loaded_data.RData")
#---------------------------                        Loading and prepping data                                           ----------------
# Load the data
load("shiny_mouse_loaded_data.RData")
#source("shiny_mouse_loaded_data.RData")

# Prepare long-file formats:
b1_abspqn_long <- pivot_longer(batch1_abspqn_file, cols = -c(sampleID), values_to = "npx", names_to = "assay")
b2_abspqn_long <- pivot_longer(batch2_abspqn_file, cols = -c(sampleID), values_to = "npx", names_to = "assay")

b1_abspqn_long <- b1_abspqn_long |> mutate(assay = toupper(assay))
b2_abspqn_long <- b2_abspqn_long |> mutate(assay = toupper(assay))
#npx_data_table <- npx_data_table |> mutate(assay = toupper(assay))

colnames(olink_t96_p1)[2:93] <- toupper(colnames(olink_t96_p1)[2:93])
colnames(olink_t96_p2)[2:93] <- toupper(colnames(olink_t96_p2)[2:93])
colnames(batch2_abspqn_file)[2:93] <- toupper(colnames(batch2_abspqn_file)[2:93])
colnames(batch1_abspqn_file)[2:93] <- toupper(colnames(batch1_abspqn_file)[2:93])

# Prepare z-score files: (Z-Score tranformation per plate and per protein):
b1_zscore_long <- b1_abspqn_long |>
  filter(!grepl("Day10|Day21|_A|_B", sampleID)) |>
  group_by(assay) |>
  mutate(z.score = (npx - mean(npx)) / sd(npx)) |>
  ungroup() |>
  mutate(assay = toupper(assay))

b2_zscore_long <- b2_abspqn_long |>
  filter(!grepl("23MFT8:3_D14|_A|_B", sampleID)) |>
  group_by(assay) |>
  mutate(z.score = (npx - mean(npx)) / sd(npx)) |>
  ungroup() |>
  mutate(assay = toupper(assay))

# And wide files:
b1_zscore_wide <- b1_zscore_long |>
  select(sampleID, assay, z.score) |>
  pivot_wider( values_from= z.score, names_from = assay) |>
  mutate( batch = "Study 1")
#colnames(b1_zscore_wide)[2:length(b1_zscore_wide)] <- toupper(colnames(b1_zscore_wide)[2:length(b1_zscore_wide)])

b2_zscore_wide <- b2_zscore_long |>
  select(sampleID, assay, z.score) |>
  pivot_wider( values_from= z.score, names_from = assay) |>
  mutate( batch = "Study 2")
#colnames(b2_zscore_wide)[2:length(b2_zscore_wide)] <- toupper(colnames(b2_zscore_wide)[2:length(b2_zscore_wide)])

## Combine sample information file:
batch1_sample_info <- batch1_sample_info |> select( -sample_number_original) |> mutate( study = 1)
batch2_sample_info <- batch2_sample_info |> select( -sample_name_original) |> mutate( study = 2)
sample_info_both <- rbind(batch1_sample_info, batch2_sample_info)

# Get the package and session information
info <- get_package_info()
package_info <- info$package_info
system_info <- info$system_info

# Convert system info to data frame
system_info_df <- data.frame(
  Item = c("R version", "Platform", "Running under"),
  Value = c(system_info$R.version$version.string, system_info$platform, system_info$running),
  stringsAsFactors = FALSE
)


#---------------------------                        UI                                            ----------------
# Define UI for application
ui <- #navbarPage(
  navbarPage(theme = shinytheme("cosmo"),
  #page_navbar(   theme = bs_theme(version = 5, bootswatch = "minty"), 
  title = "Mouse DBS profiling",
  
  ###### Here : insert shinydashboard dependencies ######
  header = tagList(
    #useShinydashboard()
  ),

  #######################################################
  #### About the app tab:  
  tabPanel( title = "Home", icon = icon("home"),
            value = "about",
           verbatimTextOutput("about"),
           h4("Companion app to"), 
           h1("Longitudinal blood microsampling and proteome monitoring facilitate timely intervention in experimental type 1 diabetes"),
           h5("Anirudra Parajuli, Annika Bendes, Fabian Byvald, Virginia M. Stone, Emma E. Ringqvist, Marta Butrym, Emmanouil Angelis, Sophie Kipper, Stefan Bauer, Niclas Roxhed, Jochen M. Schwenk, Malin Flodström-Tullberg"),
           h2("About"),
           HTML("<p>Frequent self-sampling of dried blood spots (DBS) on volumetric devices could aid the detection of early, disease-predictive protein biomarkers. To test this, a mouse model was infected with a type 1 diabetes-associated virus, coxackievirus B3 and 5 µl blood samples were collected every 1-3 days for 14 days post-infection. DBS samples were analyzed for 92 proteins, and these revealed transient proteome changes in virus-infected animals. This application has been designed to provide an interactive platform for browsing and exploring the proteomic data sets generated in this study.
           <p> This app provides complementary visualizations and aims to help with the exploration of Olink data generated from murine DBS samples. The app includes visualizations of the longitudinal protein profiles for each protein, the protein-protein correlations, and the data quality controls.
              <br>
              The <b> data </b> used in this study will be made available upon publication, and <b> code </b> used is found at the Schwenk Lab Github account.
              <br><br>
DBS samples were collected from CVB3- and mock-infected NOD mice onto Capitainer B (Ref#18-01-001) sampling cards.
<br>
Data was generated using the Olink Target 96 Mouse Exploratory panel (Art# 95380)

<br><br>
<hr>
App developed by Annika Bendes
<br>
App version: 0.0.1
<br>
")
  ),
  
  ## Load input data:
  #plate1_abspqn_file <- read.csv("./data/Mouse DBS_Batch 1_abspqn2 data_230127_wo plasma control.csv")
  
#### Protein tab:  
  tabPanel(title = "Protein profile", value = "protein_profile",
          # uiOutput("dynamic_title"),
           fluidRow(    
             column(width = 3,
                    # Select protein
                    fluidRow(
                      
                      sidebarPanel(title = "Select protein",
                          width = 8,
                          status = "primary",
                          solidHeader = TRUE,
                          selectInput("prot_view_select", "Select protein",
                                      choices = sort(unique(b2_abspqn_long$assay)) ,  selected = "CXCL9" ),
                          
                          # Select batch:
                        #  selectInput("batch_view_select", "Select Study",
                        #              choices = c("Batch 1", "Batch 2") ),
                          
                          # Select data type (AbsPQN or Z-scores):
                          selectInput("datatype_view_select", "Select data",
                                      choices = c("ProtPQN", "Z-score") )
                          
                          
                      ) #end sidebar panel
                    ), # end fluidrow
                  
                    
                    ## Show tab information:
                    # Button to toggle additional information
                    actionButton("toggle_button", HTML("Show more information &#9660;"), 
                                 style="background:none; border:none; color:#474747; font-size:14px; padding:0;"),
                    
                    # Additional information section (initially hidden)
                    div(id = "more_info", style = "display:none;",  # Initially hide the info section
                        p("This page shows protein levels as combined (Z-score) or uncombined (ProtPQN) data over sampling time. 
                        The groups (infected/control) are shown in yellow and green. Loess regression was used to illustrate the trend per group with a 95% confidence interval
                        around the smoothed lines.", tags$br(),tags$br(),
                          "Percentage of samples above and below the LOD for the protein are shown in the barplots for each cohort (ProtPQN) or combined (Z-score).",
                          tags$br(), tags$br(),   
                          "Different proteins and data set can be selected in the selection windows above.", tags$br(),
                        "Uniprot ID for the proteins are shown in green the box.", tags$br(),tags$br()

                        )
                    )
                    
             ), # end column
             
             tags$script(HTML("
    $(document).on('click', '#toggle_button', function() {
      $('#more_info').toggle(); // Toggle the visibility of the information section
      var text = $('#toggle_button').html(); // Get the current text of the button
      if (text.includes('Show')) {
        $('#toggle_button').html('Hide information &#9650;'); // Change to 'Hide' with upward triangle
      } else {
        $('#toggle_button').html('Show more information &#9660;'); // Change back to 'Show' with downward triangle
      }
    });
  ")) ,
             
             
             mainPanel(
               
               h1(uiOutput("name_selected_protein")),
               
               fluidRow( #tags$head(tags$style(HTML('.box{-webkit-box-shadow: none; border-top: none; -moz-box-shadow: none;box-shadow: none;}'))),
                 # box(width = 4,background = "black",
                 div(valueBoxOutput("valueBox_uniprot"), style = "color:white;"),style = "background-color:#639c93;"
                 #    valueBoxOutput("valueBox_uniprot", width = 3)
                 #   valueBoxOutput("valueBox2"),
                 #   valueBoxOutput("valueBox3")
                 #)
               ), # End fluid row with value box
               
               
               conditionalPanel(condition = "input.datatype_view_select == 'ProtPQN'",
             #column(width=9,
               
               # 1st fluid row for value boxes
                    # fluidRow( #tags$head(tags$style(HTML('.box{-webkit-box-shadow: none; border-top: none; -moz-box-shadow: none;box-shadow: none;}'))),
                    #         # box(width = 4,background = "black",
                    #   div(valueBoxOutput("valueBox_uniprot"), style = "color:white;"),style = "background-color:#8e9c8a;"
                    #          #    valueBoxOutput("valueBox_uniprot", width = 3)
                    #           #   valueBoxOutput("valueBox2"),
                    #           #   valueBoxOutput("valueBox3")
                    #          #)
                    # ), # End fluid row with value box
                    #br(),
                    hr(),
                    #fluidRow( box(div(style = "height:20px; font-size:30px;",
                    #                  'Study1')) ),
                    br(),
                    #br(),
                    fluidRow(#tags$head(tags$style(HTML('.box{-webkit-box-shadow: none; border-top: none;-moz-box-shadow: none;box-shadow: none;}'))),
                         #    box(   width = 3, height = 2, plotOutput("plot_missingness")
                               #   ),
                              # box(width = 4, DTOutput("cv_table")),
                        #     box(width = 6, plotOutput("trend_plot")),
                             box(title = h2("Study 1", align = "center") , solidHeader = TRUE, width = 6, plotOutput("trend_lod_plot_b1") ),
                             box(title = h2("Study 2", align = "center"), solidHeader = TRUE, width = 6, plotOutput("trend_lod_plot_b2"))
                     ),
                    #fluidRow( box(div(style = "height:20px; font-size:30px;",
                    #                  'Study2')) ),
                    br(),
                    br(),
                    fluidRow(#tags$head(tags$style(HTML('.box{-webkit-box-shadow: none; border-top: none;-moz-box-shadow: none;box-shadow: none;}'))),
                      box( width = 6, plotOutput("lod_plot_b1")),
                      box( widht = 6, plotOutput("lod_plot_b2") )
                    )
                    
             ), # end conditional panel
             
             conditionalPanel(condition = "input.datatype_view_select == 'Z-score'",
                              h2("Combined studies - Z-score transformed data"),
                              box( width = 8, plotOutput("trendline_zscore")), 
                              box( width = 8, plotOutput("lod_both_studies_plot")),
                              
                    
             ) # end conditional panel
             
           ) # end mainpanel
           ) #end fluid row
  ), # end protein view tab panel

# Protien-protein correlation tab:
tabPanel("Protein-protein correlation",
        # h2("Protien-Protien correlation"),
         
         fluidRow(    
           column(width= 3,
         
         fluidRow(
                    
                    sidebarPanel(title = "Select features",
                                 width = 8,
                                 status = "primary",
                                 solidHeader = TRUE,
                                 # selectInput("prot_view_select", "Select protein",
                                 #             choices = sort(unique(b2_abspqn_long$assay)) ),
                                 # Select batch:
                                   selectInput("batch_view_select_hm", "Select Study",
                                               choices = c("Study 1", "Study 2") ),
                                 
                                 # Select data type (AbsPQN or Z-scores):
                                 selectInput("datatype_view_select_hm", "Select data",
                                             choices = c("ProtPQN", "Z-score"),
                                             selected = "Z-score"),# Default value to Z-score
                                 
                                 conditionalPanel( 
                                   condition = "input.datatype_view_select_hm == 'Z-score'",
                                            awesomeCheckbox(
                                              inputId = "combine_data_hm",
                                              label = "Combine data", 
                                              value = TRUE,
                                              status = "info"
                                          )
                                      )
                                 )
                    ), # End fluidRow
         
         ## Show tab information:
         # Button to toggle additional information
         actionButton("toggle_button_corr", HTML("Show more information &#9660;"), 
                      style="background:none; border:none; color:#474747; font-size:14px; padding:0;"),
         
         # Additional information section (initially hidden)
         div(id = "more_info_corr", style = "display:none;",  # Initially hide the info section
             p("Different data sets and studies can be selected in the selection windows above.", tags$br(),tags$br(),
               "This page shows the protein-protein correlations calculated using Spearman's correlation in a heatmap. Red indicates a high correlation and blue indicates an inverse correlation."
             )
         )
         

           ), # end column
         
         
         tags$script(HTML("
    $(document).on('click', '#toggle_button_corr', function() {
      $('#more_info_corr').toggle(); // Toggle the visibility of the information section
      var text = $('#toggle_button_corr').html(); // Get the current text of the button
      if (text.includes('Show')) {
        $('#toggle_button_corr').html('Hide information &#9650;'); // Change to 'Hide' with upward triangle
      } else {
        $('#toggle_button_corr').html('Show more information &#9660;'); // Change back to 'Show' with downward triangle
      }
    });
  ")) ,
         
         mainPanel(
           h2("Protein-protein correlation"),
           box(width = 9, height = 15,
               shinyjqui::jqui_resizable(plotlyOutput("heatmap", height = "800px", width = "1000px" ) ))
           
         )# end mainPanel

        ) # end bigger fluid row

         
), # end protien-protien correlation tab

# QC tab:
tabPanel("Protein QC",
        # h2("QC")
        fluidRow(    
          column(width= 3,
                 
                 fluidRow(
                   
                   sidebarPanel(title="Select protein",
                                width = 8,
                                status = "primary",
                                solidHeader = TRUE,
                                 selectInput("prot_view_select_qc", "Select protein",
                                             choices = sort(unique(b2_abspqn_long$assay)), selected = "CXCL9" ),
                                # Select batch:
                                selectInput("batch_view_select_qc", "Select Study",
                                            choices = c("Study 1", "Study 2") ),
                                
                                # Select data type (AbsPQN or Z-scores):
                                selectInput("datatype_view_select_qc", "Select data",
                                            choices = c("ProtPQN", "Z-score")) ,
                                
                                conditionalPanel( 
                                  condition = "input.datatype_view_select_qc == 'Z-score'",
                                  awesomeCheckbox(
                                    inputId = "combine_data_qc",
                                    label = "Combine data", 
                                    value = FALSE,
                                    status = "info"
                                  )
                                )
                   )
                 ), # End fluidRow
                 
                 ## Show tab information:
                 # Button to toggle additional information
                 actionButton("toggle_button_qc", HTML("Show more information &#9660;"), 
                              style="background:none; border:none; color:#474747; font-size:14px; padding:0;"),
                 
                 # Additional information section (initially hidden)
                 div(id = "more_info_qc", style = "display:none;",  # Initially hide the info section
                     p("Different proteins, study sets, and data sets (ProtPQN or Z-score) can be selected in the selection windows above.", tags$br(),tags$br(),
                       "This page shows precision and variance for each protein target and for the different data sets.", tags$br(),tags$br(),
                       "Precision is determined by CV calculations for replicate samples using the ProtPQN data. The average CV for the replicate samples can be found in the table overview", 
                       tags$br(), tags$br(),
                       "The interquartile range (IQR) is shown in the overview table, and the IQR for the selected target in relation to the other proteins is shown in green in the scatter plot.",
                       tags$br(), tags$br()
                     )
                 )
                 
          ), # end column
          
          tags$script(HTML("
    $(document).on('click', '#toggle_button_qc', function() {
      $('#more_info_qc').toggle(); // Toggle the visibility of the information section
      var text = $('#toggle_button_qc').html(); // Get the current text of the button
      if (text.includes('Show')) {
        $('#toggle_button_qc').html('Hide information &#9650;'); // Change to 'Hide' with upward triangle
      } else {
        $('#toggle_button_qc').html('Show more information &#9660;'); // Change back to 'Show' with downward triangle
      }
    });
  ")) ,
          
          mainPanel(
            h2("Precision and variance"),
            br(),
            box( title = "Precision", width = 6, plotOutput("cv_density_plot")),
            box( title = "Overview", width = 6, DTOutput("cv_iqr_table")),
            br(),
            #br(),
            box(title = "Variance", width = 12, shinyjqui::jqui_resizable(plotOutput("iqr_plot" ) ))
            
          )# end mainPanel
          
        ) # end bigger fluid row
        
), # end QC tab

# 
#   tabPanel("Summary",
#            verbatimTextOutput("summary"),
#            DTOutput("table")
#   ),

tabPanel("Session information",
         #h2("Session inforamtion"),
         h4("System Information:"),
         tableOutput("system_info_table"),
         br(),
         h4("Package Information:"),
         tableOutput("package_table")
         
         )




)

#---------------------------                        SERVER                                            ----------------

# Define server logic
server <- function(input, output) {
  
  #--- Get dynamic title for protein selected:
  output$dynamic_title <- renderUI({
    req(input$prot_view_select, input$batch_view_select)
    title <- paste("Protein:", input$prot_view_select, "Batch:", input$batch_view_select)
    h2(title)
  })
  
  output$name_selected_protein <- renderUI({
    paste(toupper(input$prot_view_select))
  })
  
  
  #-------- Get direction to uniprot:
  observeEvent(input$prot_view_select, {
    selected_protein <- input$prot_view_select
    
    # Get the UniProt ID for the selected protein
    uniprot_id <- olink_t96_p1[olink_t96_p1$Assay == "Uniprot ID", selected_protein, drop = TRUE]
    
    # Construct the UniProt URL
    uniprot_url <- paste0("https://www.uniprot.org/uniprot/", uniprot_id)
    
    # Update the valueBox with the clickable link
    output$valueBox_uniprot <- shinydashboard::renderValueBox({
      valueBox(
        value = tags$a(href = uniprot_url, target = "_blank", uniprot_id),
        subtitle = paste("UniProt ID for", toupper(selected_protein))
        #icon = icon("link"),
      )
    })
    
    output$valueBox2 <- renderValueBox({
      valueBox(10 * 2, "Value Box 2", icon = icon("download"))
    })
    
    output$valueBox3 <- renderValueBox({
      valueBox(10 * 2, "Value Box 3", icon = icon("fa-grill-hot"))
    })
  })
  
  
  #--------- Number of samples below LOD plot:
  output$plot_missingness <- renderPlot({
    # Select the appropriate data frame based on batch_view_select
    df <- if (input$batch_view_select == "Study 1") {
      olink_t96_p1
    } else {
      olink_t96_p2
    }
    
    # Select the specified column
    col_data <- df[[input$prot_view_select]]
    
    # Take out the LOD value:
    reference_value <- col_data[93]
    
    # Compare the values in rows 5:92 (samples) to the value in row 93
    comparison <- ifelse(col_data[5:92] > reference_value, "Above", "Below")
    
    # Calculate percentages
    percentages <- data.frame(
      Category = c("Above", "Below"),
      Percentage = c(sum(comparison == "Above") / length(comparison) * 100,
                     sum(comparison == "Below") / length(comparison) * 100)
    )
    
    # Generate bar plot
    ggplot(percentages, aes(x = Category, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = paste("Samples above and below LOD", input$prot_view_select),
           x = "Category",
           y = "Percentage") +
      ylim(0, 100) + 
      scale_fill_manual(values = c("#639c93", "#FFBF69"), labels = c("> LOD", "< LOD")) +
      theme_bw(14)
  })
  
  #---------------------------------         QC tab                             ---------------
  ### CV plot for replicated sample:
  output$cv_density_plot <- renderPlot({
    
    column_name <- "npx" # if (input$datatype_view_select_qc == "AbsPQN") "npx" else "z.score"
    
   #if (input$datatype_view_select_qc == "ProtPQN") {
    # Select the appropriate data frame based on batch_view_select
    if (input$batch_view_select_qc == "Study 1") {
      abspqn_long <- b1_abspqn_long
      sample_info <- batch1_sample_info
    } else {
      abspqn_long <- b2_abspqn_long
      sample_info <- batch2_sample_info
     } 

    
    selected_protein <- input$prot_view_select_qc
    
    cv_table <- cv_table_fun( npx_data_table = abspqn_long, sample_info_table = sample_info, protein = selected_protein, column_name = column_name)
    
    ggplot(cv_table, aes( x = CV)) +
      geom_histogram( aes(fill = study)) +
      labs(title = paste0("CV for replicated samples - ", toupper(selected_protein)), subtitle = "Calculated on ProtPQN NPX data", x = "CV [%]") +
      xlim(0, 100) +
      theme_classic(14) +
      scale_fill_brewer( palette = "Set2")
    
  })

  #---------------------------------------------------------------------------------------------------
  
  ###               IQR plot                                  #### 
  
  output$iqr_plot <- renderPlot({
    
    column_name <- if (input$datatype_view_select_qc == "ProtPQN") "npx" else "z.score"
    
    if (input$datatype_view_select_qc == "ProtPQN") {
      # Select the appropriate data frame based on batch_view_select
      if (input$batch_view_select_qc == "Study 1") {
        abspqn_long <- b1_abspqn_long
        sample_info <- batch1_sample_info
      } else {
        abspqn_long <- b2_abspqn_long
        sample_info <- batch2_sample_info
      } 
    }  else if ( input$datatype_view_select_qc == "Z-score" ){
      if (input$combine_data_qc) {
        abspqn_long <- bind_rows(b1_zscore_long, b2_zscore_long)
        sample_info <- sample_info_both
      } else { # end if
        if (input$batch_view_select_hm == "Study 1") {
          abspqn_long <- b1_zscore_long
          sample_info <- batch1_sample_info
        } else{
          abspqn_long <- b2_zscore_long
          sample_info <- batch2_sample_info
        }
      }
    } # end else if
    
    selected_protein <- input$prot_view_select_qc
    
    # Remove replicates:
    abspqn_long <- abspqn_long |> filter(!grepl("_A|_B", sampleID))
    # Get the plot:
    iqr_pot_fun( npx_data_table = abspqn_long, sample_info_table = sample_info, protein = selected_protein, column_name = column_name)
    
    # datatable(cv_table, options = list(pageLength = 10))
    
  })
  
  # -------------------------------------------------------------------------------------------------
  ###                           CV IQR TABLE                                              #####
   output$cv_iqr_table <- renderDT({
     
     column_name <- if (input$datatype_view_select_qc == "ProtPQN") "npx" else "z.score"
     
     if (input$datatype_view_select_qc == "ProtPQN") {
       # Select the appropriate data frame based on batch_view_select
       if (input$batch_view_select_qc == "Study 1") {
         abspqn_long <- b1_abspqn_long
         sample_info <- batch1_sample_info
       } else {
         abspqn_long <- b2_abspqn_long
         sample_info <- batch2_sample_info
       } 
     }  else if ( input$datatype_view_select_qc == "Z-score" ){
       if (input$combine_data_qc) {
         abspqn_long <- bind_rows(b1_zscore_long, b2_zscore_long)
         sample_info <- sample_info_both
       } else { # end if
         if (input$batch_view_select_hm == "Study 1") {
           abspqn_long <- b1_zscore_long
           sample_info <- batch1_sample_info
         } else{
           abspqn_long <- b2_zscore_long
           sample_info <- batch2_sample_info
         }
       }
     } # end else if
     
     selected_protein <- input$prot_view_select_qc
     
     ## Obtain CV data frame:
     cv_table <- cv_table_fun( npx_data_table = abspqn_long, sample_info_table = sample_info, protein = selected_protein, column_name = column_name)
     
     # Remove replicates:
     abspqn_long <- abspqn_long |> filter(!grepl("_A|_B", sampleID))
     
     ## Obtain IQR data frame:
    iqr_i <- abspqn_long |> filter(assay == selected_protein) |> select( any_of(column_name)) |> apply( 2, as.numeric) |> IQR( na.rm = TRUE)
    
    # Prepare a summary table:
    df <- data.frame( rown.names = c("Average CV", "Median CV", "IQR"), round(c( mean(cv_table$CV), median(cv_table$CV), iqr_i),  digits = 1))
    colnames(df) <- c("Variable", "Value")

    datatable(df)
    
   }) # end render table
  
  #--------------------------------------------------------------------------------------------
  
  
  
  #####------------                  Protein tab and trendline plots              --------------------
  
  ## Trend line plots grouped by infection status and data points below LOD for Study Batch 2:

  output$trend_lod_plot_b2 <- renderPlot({
      selected_protein <- input$prot_view_select
      sample_info <- batch2_sample_info
 
      
      b2_abspqn_long_temp <-  b2_abspqn_long |>
        filter(!grepl("23MFT8:3_D14|_A|_B", sampleID))
      
    plot2 <- long_plot_fun(b2_abspqn_long_temp, batch2_sample_info, selected_protein)
  
    patchwork::wrap_plots(plot2)
  })
  
  ##              Samples above LOD plot:                       ####
  output$lod_plot_b2 <- renderPlot({
    selected_protein <- input$prot_view_select
    sample_info <- batch2_sample_info
    df <- olink_t96_p2
    
    # Select the specified column
    col_data <- df[[input$prot_view_select]]
    
    # Take out the LOD value:
    reference_value <- col_data[93]
    
    # Compare the values in rows 5:92 (samples) to the value in row 93
    comparison <- ifelse(col_data[5:92] > reference_value, "Above", "Below")
    
    # Calculate percentages
    percentages <- data.frame(
      Category = c("Above", "Below"),
      Percentage = c(sum(comparison == "Above") / length(comparison) * 100,
                     sum(comparison == "Below") / length(comparison) * 100)
    )
    
    # Generate bar plot
    plot1 <- ggplot(percentages, aes(x = Category, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = paste("LOD for ", toupper(input$prot_view_select)),
           x = "Category",
           y = "Percentage of samples [%]") +
      ylim(0, 100) +
      scale_fill_manual(values = c("#639c93", "#FFBF69"), labels = c("> LOD", "< LOD")) +
      theme_bw(15) +
      guides(fill = "none")
    
    patchwork::wrap_plots(plot1)
  })
  

  
  
  output$trend_lod_plot_b1 <- renderPlot({
    selected_protein <- input$prot_view_select
    sample_info <- batch1_sample_info
    selected_protein <- input$prot_view_select
    
    
    b1_abspqn_long_temp <-  b1_abspqn_long |>
      filter(!grepl("Day10|Day21|_A|_B", sampleID))
    
    plot2 <- long_plot_fun(b1_abspqn_long_temp, batch1_sample_info, selected_protein)
    
    patchwork::wrap_plots(plot2)
  })
  
  ### LOD plot Study 1            ###
  
  output$lod_plot_b1 <- renderPlot({
    selected_protein <- input$prot_view_select
    sample_info <- batch1_sample_info
    selected_protein <- input$prot_view_select
    df <- olink_t96_p1
    
    # Select the specified column
    col_data <- df[[input$prot_view_select]]
    
    # Take out the LOD value:
    reference_value <- col_data[93]
    
    # Compare the values in rows 5:92 (samples) to the value in row 93
    comparison <- ifelse(col_data[5:92] > reference_value, "Above", "Below")
    
    # Calculate percentages
    percentages <- data.frame(
      Category = c("Above", "Below"),
      Percentage = c(sum(comparison == "Above") / length(comparison) * 100,
                     sum(comparison == "Below") / length(comparison) * 100)
    )
    
    # Generate bar plot
    plot1 <- ggplot(percentages, aes(x = Category, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = paste("LOD for ", toupper(input$prot_view_select)),
           x = "Category",
           y = "Percentage of samples [%]") +
      ylim(0, 100) +
      scale_fill_manual(values = c("#639c93", "#FFBF69"), labels = c("> LOD", "< LOD")) +
      theme_bw(15) +
      guides(fill = "none")
    
    patchwork::wrap_plots(plot1)
  })
  
  
  #-----------------------------------------------------------------------------------------------
  ## LOD plots for both plates:
  
  output$lod_both_studies_plot <- renderPlot({
    selected_protein <- input$prot_view_select
    sample_info <- sample_info_both
    df <- olink_t96_p1
    
    # Select the specified column
    col_data <- df[[input$prot_view_select]]
    
    # Take out the LOD value:
    reference_value <- col_data[93]
    
    # Compare the values in rows 5:92 (samples) to the value in row 93
    comparison <- ifelse(col_data[5:92] > reference_value, "Above", "Below")
    
    # Calculate percentages
    percentages <- data.frame(
      Category = c("Above", "Below"),
      Percentage = c(sum(comparison == "Above") / length(comparison) * 100,
                     sum(comparison == "Below") / length(comparison) * 100)
    )
    
    # Generate bar plot
    plot1 <- ggplot(percentages, aes(x = Category, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = paste("Study 1 - LOD for ", toupper(input$prot_view_select)),
           x = "Category",
           y = "Percentage of samples [%]") +
      ylim(0, 100) +
      scale_fill_manual(values = c("#639c93", "#FFBF69"), labels = c("> LOD", "< LOD")) +
      theme_bw(15) +
      guides(fill = "none")
    
    df <- olink_t96_p2
    
    # Select the specified column
    col_data <- df[[input$prot_view_select]]
    
    # Take out the LOD value:
    reference_value <- col_data[93]
    
    # Compare the values in rows 5:92 (samples) to the value in row 93
    comparison <- ifelse(col_data[5:92] > reference_value, "Above", "Below")
    
    # Calculate percentages
    percentages <- data.frame(
      Category = c("Above", "Below"),
      Percentage = c(sum(comparison == "Above") / length(comparison) * 100,
                     sum(comparison == "Below") / length(comparison) * 100)
    )
    
    # Generate bar plot
    plot2 <- ggplot(percentages, aes(x = Category, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = paste("Study 2 - LOD for ", toupper(input$prot_view_select)),
           x = "Category",
           y = "Percentage of samples [%]") +
      ylim(0, 100) +
      scale_fill_manual(values = c("#639c93", "#FFBF69"), labels = c("> LOD", "< LOD")) +
      theme_bw(15) +
      guides(fill = "none")
    
    patchwork::wrap_plots(plot1 + plot2)
  })
  
  
  #------------------                               Trendline plots for z-score                -------------------
  ## Trend line plots for Z-scores keep 
  output$trendline_zscore <- renderPlot({
    selected_protein <- toupper(input$prot_view_select)
    sample_info <- sample_info_both
  
    rbind( (b1_zscore_long |> mutate( samp_batch = "Study 1")), (b2_zscore_long |> mutate( samp_batch = "Study 2")) )  |>
      filter( assay == selected_protein) |>
      left_join(sample_info, by = c("sampleID" = "Sample_name"), relationship = "many-to-many") |>
      #filter(!study_subject %in% "Sample pool") |>
      mutate(Infection_status = tolower(Infection_status)) |>
      mutate( Infection_status = case_when(Infection_status == "not infected" ~ "control",
                                           T ~ Infection_status)) |>
      #filter(!str_detect(sample_number_original, "^21")) |>
      ggplot( aes( x = as.numeric(sample_day), y = z.score, group = Infection_status, col = Infection_status )) +
      geom_point(aes(fill=Infection_status), alpha = 0.9, size = 3) +
      labs(title = paste0(selected_protein)) + #,
      xlab("Sample day") + 
      ylab("Z-score") +
      theme_classic(15) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      geom_smooth(method = "loess", aes(fill=Infection_status, col = Infection_status), formula = 'y ~ x'  ) +
      scale_fill_manual( values = c( "infected" = "#FFBF69", "control" = "#639c93")) +
      scale_color_manual(values =  c( "infected" = "#FFBF69", "control" = "#639c93")) +
      guides(fill = "none")
    
  })
  

  
  
  ##----------                    Protein-protein correlation heatmap:                            ------
  # The heatmap itself
  cor_plot <- reactive({
    if (input$datatype_view_select_hm == "ProtPQN") {
      hm_df <- if (input$batch_view_select_hm == "Study 1") {
        batch1_abspqn_file %>%
          filter(!grepl("_A|_B", sampleID)) |>
          column_to_rownames(var = "sampleID")
      } else {
        batch2_abspqn_file |>
          filter(!grepl("_A|_B", sampleID)) |>
          column_to_rownames(var = "sampleID")
      }
    } else if (input$datatype_view_select_hm == "Z-score") {
      if (input$combine_data_hm) {
        combined_zscore_df <- bind_rows(
          b1_zscore_wide %>%
            filter(!grepl("_A|_B", sampleID)) |>
            mutate(batch = "Study 1"),
          b2_zscore_wide %>%
            filter(!grepl("_A|_B", sampleID)) |>
            mutate(batch = "Study 2")
        )
        hm_df <- combined_zscore_df |>
          select(-batch) |>
          column_to_rownames(var = "sampleID")
      } else {
        hm_df <- if (input$batch_view_select_hm == "Study 1") {
          b1_zscore_wide |>
            filter(!grepl("_A|_B", sampleID)) |>
            column_to_rownames(var = "sampleID") |>
            select(!batch)
        } else {
          b2_zscore_wide |>
            filter(!grepl("_A|_B", sampleID)) |>
            column_to_rownames(var = "sampleID") |>
            select(!batch)
        }
      }
    }
    
    return(hm_df)
  })
  
  output$heatmap <- renderPlotly({
    heatmaply_cor(cor(cor_plot(), method = "spearman"), colors = c( "skyblue2", "white" , "sienna2" ), fontsize_col = 6, fontsize_row = 6 )
  })

  
  ### Session information ###
  # Render the system info table
  output$system_info_table <- renderTable({
    system_info_df
  }, sanitize.text.function = function(x) x)
  
  # Render the package table
  output$package_table <- renderTable({
    package_info
  }, sanitize.text.function = function(x) x)
  
  

}

# Run the application 
shinyApp(ui = ui, server = server)

