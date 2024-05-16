# My function to build kables in reports
# --------------------------------------
kable_mydf <- function(x, 
                       boldfirstcolumn = F, 
                       font_size = 10, 
                       round_number = 2,
                       format_output = "html",
                       bootstrap_options = c("hover", "condensed", "responsive")){
  x %>% 
    mutate(across(where(is.numeric), \(x) round (x, round_number))) %>%
    kable(format = format_output) %>%  
    kable_styling(bootstrap_options = bootstrap_options, 
                  full_width = F, 
                  font_size = font_size) %>% 
    {if(boldfirstcolumn == TRUE) column_spec(., 1, bold = T) else .}
}
