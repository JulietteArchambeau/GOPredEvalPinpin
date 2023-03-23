# Function to extract information from the metadata of ClimateDT

extract_climatedt_metadata <- function(var_clim){
  
  read_excel(here::here("data/ClimaticData/ClimateDT_IndicesMetadata.xlsx"))%>% 
    set_colnames(str_to_lower(colnames(.))) %>% 
    mutate(description = str_to_sentence(description)) %>% 
    mutate(description = case_when(label=="bio12" ~ "Annual precipitation",
                                   label=="MSP" ~ "Summer precipitation",
                                   TRUE ~ description),
           label = if_else(label=="MSP", "SP", label)) %>%  
    mutate(unit_symbol = case_when(grepl("\\(", unit) ~ str_sub(unit, -3,-2),
                                   grepl(" ", unit) ~ str_replace_all(unit, " ", ""),
                                   TRUE ~ unit)) %>% 
    dplyr::filter(label %in% var_clim)
  
}




