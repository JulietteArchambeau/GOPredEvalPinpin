# My function to build kables in reports
# --------------------------------------
kable_mydf <- function(x, 
                       boldfirstcolumn = F, 
                       font_size = 10, 
                       round_number = 2,
                       bootstrap_options = c("hover", "condensed", "responsive")){
  x %>% 
    mutate(across(where(is.numeric), round, round_number)) %>%
    kable() %>%  
    kable_styling(bootstrap_options = bootstrap_options, 
                  full_width = F, 
                  font_size = font_size) %>% 
    {if(boldfirstcolumn == TRUE) column_spec(., 1, bold = T) else .}
}
