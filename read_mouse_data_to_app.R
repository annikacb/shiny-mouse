# Load the different input files:
library(readxl)
library(tidyverse)


## Batch 1 data:
batch1_abspqn_file <- read.csv("./data/Mouse DBS_Batch 1_abspqn2 data_230127_wo plasma control.csv")
batch1_sample_info <-  read_excel("./sample_information/dbs_mouse_sample_information.xlsx")

olink_t96_p1 <- read_csv("./data/Mouse DBS Plate 1_NPX.csv", skip = 3)
#olink_t96_p1 <- olink_t96_p1 |>
#  select( -c(`QC Warning`, ...96 )) |>
#  filter( !Assay %in% c("Uniprot ID", "OlinkID", "LOD", "Missing Data freq.", "Facility Ctrl"))


## Batch 2 data:
batch2_abspqn_file <- read.csv("./data/Mouse DBS_Batch 2_abspqn2 data_230127_wo plasma control.csv")
batch2_sample_info <-  read_excel("./sample_information/dbs_mouse_plate2_sample_information.xlsx")

olink_t96_p2 <- read_csv("./data/Mouse DBS Plate 2_NPX_belowLOD.csv", skip = 3)
#olink_t96_p2 <- olink_t96_p2 |>
#  select( -c(`QC Warning`, `QC Deviation from median...96`, `QC Deviation from median...97`,...98  )) |>
#  filter( !Assay %in% c("Uniprot ID", "MaxLOD", "OlinkID", "LOD", "Missing Data freq.", "Facility Ctrl"))

save.image("shiny_mouse_loaded_data.RData")


#ProtPQN::apply_protpqn()

#olink_t96_p1 %>%
#  filter(!grepl("Uniprot|OlinkID|LOD|Missing|Facility", Assay)) %>%
#  mutate(sample_id = Assay) %>%
#  select(-contains("Assay"), -contains("QC"), -contains("Plate"), -contains("..")) -> test
# test[,1:92] <- apply(test[,1:92], 2, as.numeric)


# ProtPQN::apply_protpqn(test, long_format = FALSE, kitwise = FALSE)