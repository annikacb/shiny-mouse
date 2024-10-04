## Functions for the shiny app:

# ### CV table for replicated sample:
# cv_table_fun <- function( npx_data_table, sample_info_table, protein){
# # npx_data_table table with sample ID and npx values in long format
# # sample_info_table table with sample inforamtion containing info about technical and biological replicates
# # protein that we want to calculate cv for
#   
#   sample_info <- sample_info_table
#   abspqn_long <- npx_data_table
#   
#   # Filter the rows based on the conditions in sample_info
#   tech_replicated_samples <- sample_info %>%
#     filter(technical_replicate == 1) |>
#     select(Sample_name) |>
#     unlist()
#   
#   bio_replicated_samples <- sample_info %>%
#     filter(biological_replicate_disc == 1) |>
#     select(Sample_name) |>
#     unlist()
# 
# # Select rows from abspqn_long that have assay == input$prot_view_select and match filtered_sample_info
# technical_replicates <- abspqn_long %>%
#     filter(assay == protein & sampleID %in% tech_replicated_samples)
# 
# # Technical replicates:
# tech_replicate_groups <- split(technical_replicates$npx, gsub("_B$|_A$", "", technical_replicates$sampleID))
# 
# # Biological replicates:
# biological_replicates <- abspqn_long %>%
#   filter(assay == protein & sampleID %in% bio_replicated_samples) |>
#   left_join(sample_info |> select(Sample_name, study_subject), by = c("sampleID" = "Sample_name") )
# 
# tech_replicate_groups_names <- names(tech_replicate_groups)
# 
# cv_data <- lapply(seq_along(tech_replicate_groups), function(i) {
#   group <- tech_replicate_groups[[i]]
#   group_name <- tech_replicate_groups_names[i]
#   if (length(group) > 1) {
#     cv <- sd(group) / mean(group) * 100
#     data.frame(Sample_Name = group_name, CV = cv, Num_Samples = length(group))
#   } else {
#     data.frame(Sample_Name = group_name, CV = NA, Num_Samples = length(group))
#   }
# })
# 
# # Combine all CV data into a single data frame
# cv_df <- do.call(rbind, cv_data)
# 
# # Add biological replicate:
# cv_bio <- biological_replicates |>
#   group_by(study_subject) |>
#   summarise( CV = sd(npx)/mean(npx) * 100,
#              Num_Samples = n())
#   #sd(biological_replicates$npx) / mean(biological_replicates$npx) * 100
# 
# bio_sampleID <- biological_replicates |>
#   mutate(length = nchar(sampleID)) |>
#   filter(length == min(length)) |>
#   select(sampleID)|>
#   unlist()
# 
# if (length(bio_sampleID) > 1){
#   bio_sampleID <- bio_sampleID[1]
# }
# 
# cv_bio_df <- data.frame(Sample_Name = bio_sampleID, CV = cv_bio, Num_Samples = length(biological_replicates$sampleID) )
# 
# cv_df <- rbind(cv_df, cv_bio_df)
# rownames(cv_df) <- c(1:nrow(cv_df))
# 
# return(cv_df)
# 
# } 

cv_table_fun <- function(npx_data_table, sample_info_table, protein, column_name) {
  # npx_data_table table with sample ID and npx values in long format
  # sample_info_table table with sample information containing info about technical and biological replicates
  # protein that we want to calculate cv for
  # column_name specifies whether to use "npx" or "z.score"
  
  sample_info <- sample_info_table
  abspqn_long <- npx_data_table
  
  # Filter the rows based on the conditions in sample_info
  tech_replicated_samples <- sample_info %>%
    filter(technical_replicate == 1) %>%
    select(Sample_name) %>%
    unlist()
  
  bio_replicated_samples <- sample_info %>%
    filter(biological_replicate_disc == 1) %>%
    select(Sample_name) %>%
    unlist()
  
  # Select rows from abspqn_long that have assay == protein and match filtered_sample_info
  technical_replicates <- abspqn_long %>%
    filter(assay == protein & sampleID %in% tech_replicated_samples)
  
  # Technical replicates:
  tech_replicate_groups <- split(technical_replicates[[column_name]], gsub("_B$|_A$", "", technical_replicates$sampleID))
  
  # Biological replicates:
  biological_replicates <- abspqn_long %>%
    filter(assay == protein & sampleID %in% bio_replicated_samples) %>%
    left_join(sample_info %>% select(Sample_name, study_subject), by = c("sampleID" = "Sample_name"))
  
  tech_replicate_groups_names <- names(tech_replicate_groups)
  
  cv_data <- lapply(seq_along(tech_replicate_groups), function(i) {
    group <- tech_replicate_groups[[i]]
    group_name <- tech_replicate_groups_names[i]
    if (length(group) > 1) {
      cv <- sd(2^group) / mean(2^group) * 100
      data.frame(Sample_Name = group_name, CV = cv, Num_Samples = length(group))
    } else {
      data.frame(Sample_Name = group_name, CV = NA, Num_Samples = length(group))
    }
  })
  
  # Combine all CV data into a single data frame
  cv_df <- do.call(rbind, cv_data)
  
  # Add biological replicate:
  cv_bio <- biological_replicates %>%
    group_by(study_subject) %>%
   # summarise(CV = sd(.data[[column_name]]) / mean(.data[[column_name]]) * 100,
  #            Num_Samples = n())
  summarise(CV = sd(2^(.data[[column_name]])) / mean(2^(.data[[column_name]])) * 100,
            Num_Samples = n())
  
  # bio_sampleID <- biological_replicates %>%
  #   mutate(length = nchar(sampleID)) %>%
  #   filter(length == min(length)) %>%
  #   select(sampleID) %>%
  #   unlist()
  
  #if (length(bio_sampleID) > 1) {
  #  bio_sampleID <- bio_sampleID[1]
  #}
  
  cv_bio_df <- data.frame(Sample_Name = cv_bio$study_subject, CV = cv_bio$CV, Num_Samples = cv_bio$Num_Samples)
  
  cv_df <- rbind(cv_df, cv_bio_df)
  rownames(cv_df) <- 1:nrow(cv_df)
  
  cv_df <- cv_df |>
    mutate( study = case_when( grepl("F06|F05", Sample_Name ) ~ "Study 1",
                               T ~ "Study 2"))
  
  return(cv_df)
}

#---------------------------------                        Biological replicated CV All proteins         ---------------------------------


cv_table_biological_fun <- function(npx_data_table, sample_info_table) {
  # npx_data_table table with sample ID and npx values in long format
  # sample_info_table table with sample information containing info about technical and biological replicates
  
  sample_info <- sample_info_table
  abspqn_long <- npx_data_table
  
  column_name = "npx"
  
  bio_replicated_samples <- sample_info %>%
    filter(biological_replicate_disc == 1) %>%
    select(Sample_name) %>%
    unlist()

  # Biological replicates:
  biological_replicates <- abspqn_long %>%
    filter( sampleID %in% bio_replicated_samples) %>%
   # group_by( assay) %>%
    left_join(sample_info %>% select(Sample_name, study_subject), by = c("sampleID" = "Sample_name"))
  
  # Add biological replicate:
  cv_bio <- biological_replicates %>%
    group_by(assay, study_subject) %>%
    # summarise(CV = sd(.data[[column_name]]) / mean(.data[[column_name]]) * 100,
    #            Num_Samples = n())
    summarise(CV = sd(2^(.data[[column_name]])) / mean(2^(.data[[column_name]])) * 100,
              Num_Samples = n())
  
  cv_bio_df <- data.frame(cv_bio)
  
  return(cv_df)
}


#---------------------------------                        LOD and longitudinal plot                       --------------------------------

# Plot function for longitudinal and bar plots next to each other:
long_plot_fun <- function(npx_data_table, sample_info_table, protein){
  # npx_data_table table with sample ID and npx values in long format
  # data_wide : Olink output file in wide format
  # sample_info_table table with sample inforamtion containing info about Infection status
  # Specify protein targets
  
  sample_info <- sample_info_table
  selected_protein <- protein # input$prot_view_select
  
  # # Select the specified column
  # col_data <- data_wide[c("Assay", selected_protein)]  #df[[input$prot_view_select]]
  # 
  # # Take out the LOD value:
  # reference_value <- data_wide |>
  #   filter( Assay == c("LOD")) |>
  #   select(any_of(selected_protein))
  # 
  # # Create the template with both categories
  # template <- data.frame(comparison = c("Above", "Below"))
  # 
  # # Compare the values in rows 5:92 (samples) to the value in row 93
  # percentages <- col_data |>
  #   filter(!str_detect(Assay, "Uniprot|OlinkID|Facility|LOD|Missing Data freq.")) |>
  #   mutate( comparison = ifelse( !!sym(selected_protein) > reference_value, "Above", "Below")) |>
  #   group_by(comparison) |>
  #   summarize(count = n(), .groups = 'drop') %>%
  #   right_join(template, by = "comparison") %>%
  #   replace_na(list(count = 0)) %>%
  #   mutate(percentage = count / sum(count) * 100) |>
  #   mutate(Category = comparison)
  # 
  # # Calculate percentages
  # # percentages <- data.frame(
  # #   Category = c("Above", "Below"),
  # #   Percentage = c(sum(comparison == "Above") / length(comparison) * 100,
  # #                  sum(comparison == "Below") / length(comparison) * 100)
  # # )
  # 
  # # Generate bar plot
  # plot1 <- ggplot(percentages, aes(x = Category, y = percentage, fill = Category)) +
  #   geom_bar(stat = "identity") +
  #   theme_minimal() +
  #   labs(title = paste("Samples above and below LOD - ", selected_protein),
  #        x = "Category",
  #        y = "Percentage [%]") +
  #   ylim(0, 100) + 
  #   scale_fill_manual(values = c("#639c93", "#FFBF69"), labels = c("> LOD", "< LOD")) +
  #   theme_bw(14)
  
  plot <- npx_data_table |>
    filter( assay == selected_protein) |>
    left_join(sample_info, by = c("sampleID" = "Sample_name"), relationship = "many-to-many") |>
    #filter(!study_subject %in% "Sample pool") |>
    mutate(Infection_status = tolower(Infection_status),
           Infection_status = case_when( Infection_status == "not infected" ~ "control", T ~ Infection_status )) |>
    #filter(!str_detect(sample_number_original, "^21")) |>
    mutate( category = Infection_status) |>
    ggplot( aes( x = as.numeric(sample_day), y = npx, group = category, col = category )) +
    geom_point(aes(fill=category), alpha = 0.9, size = 3) +
    labs(title = paste0(toupper(selected_protein))) + #,
    xlab("Sampling day") + 
    ylab("ProtPQN NPX") +
    theme_classic(17) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    geom_smooth(method = "loess", aes(fill=category, col = category) ) +
    scale_fill_manual( values = c( "infected" = "#FFBF69", "control" = "#639c93")) +
    scale_color_manual(values =  c( "infected" = "#FFBF69", "control" = "#639c93"))# +
    #guides(fill = "none", color = "none")
  
  return(plot)
} 

##--------------------                        IQR plot                   ----------------------------
# IQR  sorted from lowst to highest for each sample type
#IQR = Q3 – Q1

iqr_pot_fun <- function(npx_data_table, sample_info_table, protein, column_name) {
  # npx_data_table table with sample ID and npx values in long format
  # sample_info_table table with sample information containing info about technical and biological replicates
  # protein that we want to label in the IQ plot
  # column_name specifies whether to use "npx" or "z.score"
  
  sample_info <- sample_info_table
  abspqn_long <- npx_data_table
  
  all_assays <- unique(abspqn_long$assay)
  
  IQR_matrix <- NULL
  for (protein_i in all_assays){
    iqr_i <- abspqn_long |> filter(assay == protein_i) |> select( any_of(column_name)) |> apply( 2, as.numeric) |> IQR( na.rm = TRUE)
    IQR_matrix <- rbind(IQR_matrix,  cbind(protein_i,iqr_i )  )
  }
  
  IQR_matrix <- as.data.frame(IQR_matrix)
  IQR_matrix$iqr_i <- as.numeric(IQR_matrix$iqr_i)
  
  # Add a column to indicate if the protein matches the selected protein
  IQR_matrix$highlight <- ifelse(IQR_matrix$protein_i == protein, "highlight", "normal")
  # Add x position for label
  IQR_matrix$x_pos <- seq_along(IQR_matrix$protein_i)
  
  # Create a plot and highlight the selected protein:
  plot <- ggplot(IQR_matrix, aes(x=reorder(protein_i, iqr_i) , y=iqr_i)) + 
    geom_point(aes(color = highlight), size = 3) +
    scale_color_manual(values = c("highlight" = "darkolivegreen", "normal" = "gray80")) +
    geom_text(aes(label = ifelse(protein_i == protein, toupper(protein), "")), 
              hjust = 0.5, vjust = -1.5, color = "black") +
    labs( title = "IQR plot", x= "Protein", y = "IQR") +
    theme_classic(14) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 9), axis.line = element_line(linewidth = 1)) +
    guides( col = "none")
  
  return(plot)
}


# Function to get package and session information
get_package_info <- function() {
  si <- sessionInfo()
  packages <- c(si$otherPkgs, si$loadedOnly)
  
  format_citation <- function(pkg_name) {
    tryCatch({
      citation_obj <- citation(pkg_name)
      # citation_str <- sapply(citation_obj, function(x) {
      #   author <- if (!is.null(x$author)) paste(unlist(x$author), collapse = ", ") else "No author"
      #   citation_text <- x$textVersion
      #   url <- if (!is.null(x$url)) x$url else ""
      #   citation_formatted <- paste0(author, " - ", citation_text, if (url != "") paste0(" <a href='", url, "' target='_blank'>[link]</a>") else "")
      #   citation_formatted
      # })
      citation_str <-  sapply(pkg_name, function(x) {
        # Get citation in HTML format, edit links to open in new tabs
        citation(x) %>% format(style = "html") %>% str_replace_all("<a href", "<a target='_blank' rel='noopener noreferrer' class='sinfo-link' href") %>%
          unlist() %>% paste(collapse = "")
      })
      paste(citation_str, collapse = "<br>")
    }, error = function(e) {
      return("No citation information available")
    })
  }
  
  packages <- packages[order(names(setNames(packages, packages)))]
  
  package_info <- data.frame(
    Package = names(packages),
    Version = sapply(packages, function(pkg) pkg$Version),
    Citation = sapply(packages, function(pkg) format_citation(pkg$Package)),
    stringsAsFactors = FALSE
  )
  
  return(list(package_info = package_info, system_info = si))
}