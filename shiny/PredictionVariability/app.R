###################
# Shiny app script
###################

# Visualizing climate differences between current and future climates in Scots pine populations

# Required libraries
library(shiny)
library(reshape2)

# Load genomic offset predictions
df <- readRDS("go_df.rds")

# the colors I will use to color the populations with high GO
my_colors <- c(brewer.pal(n=12, "Paired"),"#FF40EE")

#########
#  UI
#########
ui <- fluidPage(
  titlePanel("Variability in genomic offset predictions"),
  sidebarLayout(
    sidebarPanel(
      width = 2, # Adjust the width here
      selectInput("method", "Method:", choices = unique(df$method), multiple = TRUE, selected = c("GF","GDM","LFMM")),
      selectInput("snpset_names", "SNP set:", choices = unique(df$snpset_names), multiple = TRUE, selected = c("All candidate SNPs", "All SNPs")),
      selectInput("gcm", "GCM:", choices = unique(df$gcm))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("About this app",
                 htmlOutput("about_app_text")),
        tabPanel("Bumpchart of population ranks",
                 plotOutput("bumpchart", height = "900px")),
        tabPanel("Correlations",
                 plotOutput("correlation_plot", height = "800px",width = "1350px")),
      )
    )
  )
)


########
# Server
########
server <- function(input, output) {
  
  
  # About this app
  ################
  output$about_app_text <- renderUI({
    HTML(
      "</br>
      <p>This application is designed to visualize the results from Archambeau et al (2024). 
      Evaluating genomic offset predictions in a forest tree with high population genetic structure. 
      <i>bioRxiv</i> (<a href='https://www.biorxiv.org/content/10.1101/2024.05.17.594631v1.full.pdf'>https://www.biorxiv.org/content/10.1101/2024.05.17.594631v1.full.pdf</a>).</p>
      
      <p><b>A bumpchart</b> shows the <b>population ranks based on their genomic offset predictions</b> for each combination of method and SNP set. 
      Populations with higher genomic offset have a lower value rank. 
      Colored populations are those having genomic offset values in the top three highest values in at least one row 
      (i.e., combination of SNP set and method).</p>
      
      <p><b>Correlation matrices</b> show the correlation among genomic offset predictions from the different methods and SNP sets.</p>
      
      The bumpchart and the correlation matrices can be generated for each Global Climate Model (GCM) or 
      for the mean genomic offset predictions across the five GCMs used in the study (i.e. 'GCMs_average').</p>"
    )
  })
  
  
  # CHESS-SCAPE climatic variables
  ################################
  
  
output$bumpchart <- renderPlot({
  
  plot_df <- df %>% filter(method %in% input$method & snpset_names %in% input$snpset_names & gcm == input$gcm)
  
  high_go_pops <- plot_df %>% filter(rank < 4) %>% pull(pop) %>% unique()    
  
  my_palette <- c(my_colors[1:length(high_go_pops)],"#E8E8E8")
  
  sub <- plot_df %>% 
    mutate(flag = ifelse(pop %in% high_go_pops, TRUE, FALSE),
           pop_col = if_else(flag == TRUE, pop, "Others")) %>% 
    mutate(pop = factor(pop, levels=c(setdiff(unique(plot_df$pop), high_go_pops), high_go_pops)),
           pop_col = factor(pop_col, levels = c(high_go_pops,"Others")))
  
  sub %>% ggplot(aes(x = method_snpset_names, y = rank, group = pop)) +
    coord_flip() +
    scale_y_continuous(breaks = 1:34, minor_breaks = 1:34) +
    geom_point(aes(color = pop_col), size = 2, alpha = 0.9) +
    geom_line(aes(color = pop_col), linewidth = 2, alpha = 0.8) +
    scale_color_manual(values = my_palette) +
    geom_text(data = sub %>% filter(method_snpset_names == last(levels(factor(sub$method_snpset_names)))),
              # to only write the names of the highlighted pops
              #data = sub %>% filter(method_snpset_names == last(levels(factor(sub$method_snpset_names))) & pop %in% high_go_pops), 
              aes(label = pop, x = length(unique(sub$method_snpset_names)) + 0.1), 
              color = "gray20", 
              size = 4, 
              angle = 40) +
    theme_bw() +
    ylab("Population rank") +
    theme(#panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=11),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank()) 
  
    })
  

output$correlation_plot <- renderPlot({

  correlations <- df %>% 
    filter(method %in% input$method & snpset_names %in% input$snpset_names & gcm == input$gcm) %>% 
    pivot_wider(id_cols = pop, values_from = go, names_from = method_snpset_names)  %>% 
    dplyr::select(-pop) %>% 
    cor()
    
  # Order variables alphabetically
  correlations <- correlations[order(rownames(correlations)), order(colnames(correlations))]
    
  # Melt the correlation matrix into a long format
  cor_data <- reshape2::melt(correlations)
  
  # Filter to only show the lower triangle
  cor_data <- cor_data %>% filter(as.numeric(Var1) > as.numeric(Var2))
  
  # Plot
  ggplot(cor_data, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#BB4444", mid = "#FFFFFF", high = "#4477AA", midpoint = 0, limits = c(-1, 1), name="Correlation") +
      geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),
            axis.text.y = element_text(size=15)) +
      labs(x = NULL, y = NULL)

  })

}


######################
# Run the application
######################
shinyApp(ui = ui, server = server)

