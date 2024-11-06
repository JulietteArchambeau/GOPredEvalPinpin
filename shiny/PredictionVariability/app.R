###################
# Shiny app script
###################

# Visualizing climate differences between current and future climates in Scots pine populations

# Required libraries
library(shiny)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)

# Load genomic offset predictions
df <- readRDS("go_df.rds")

# Load the genomic offset predictions and the climatic distances
go_climdist_df <- readRDS("go_climdist_df.rds")

# the colors I will use to color the populations with high GO
my_colors <- c(brewer.pal(n=12, "Paired"),"#FF40EE")

#########
#  UI
#########

ui <- fluidPage(
  titlePanel("Variability in genomic offset predictions"),
  
  tabsetPanel(
    tabPanel("About this app",
             htmlOutput("about_app_text")),
    
    tabPanel("Bumpchart of population ranks",
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 selectInput("method", "Method:", choices = unique(df$method), multiple = TRUE, selected = c("GF","GDM","LFMM")),
                 selectInput("snpset_names", "SNP set:", choices = unique(df$snpset_names), multiple = TRUE, selected = c("All candidate SNPs (380)", "All SNPs (9817)")),
                 selectInput("gcm", "GCM:", choices = unique(df$gcm))
               ),
               mainPanel(
                 plotOutput("bumpchart", height = "900px")
               )
             )
    ),
    
    tabPanel("Correlation matrices",
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 selectInput("method", "Method:", choices = unique(df$method), multiple = TRUE, selected = c("GF","GDM","LFMM")),
                 selectInput("snpset_names", "SNP set:", choices = unique(df$snpset_names), multiple = TRUE, selected = c("All candidate SNPs (380)", "All SNPs (9817)")),
                 selectInput("gcm", "GCM:", choices = unique(df$gcm))
               ),
               mainPanel(
                 plotOutput("correlation_plot", height = "800px", width = "1350px")
               )
             )
    ),
    
    tabPanel("Scatter plots",
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 
                 # Title for the first variable set
                 tags$h4("Variable on the x-axis"),
                 selectInput("method_1", "Method:", choices = unique(go_climdist_df$method_name), selected = "Climatic distance"),
                 selectInput("var_1", "SNP set or climatic variable:", choices = unique(go_climdist_df$var_name), selected = "Euclidean climatic distance"),
                 selectInput("gcm_1", "GCM:", choices = unique(go_climdist_df$gcm), selected = "Average across GCMs"),
                 
                 tags$br(), tags$br(),
                 
                 # Title for the second variable set
                 tags$h4("Variable on the y-axis"),
                 selectInput("method_2", "Method:", choices = unique(go_climdist_df$method_name), selected = "GF"),
                 selectInput("var_2", "SNP set or climatic variable:", choices = unique(go_climdist_df$var_name), selected = "All candidate SNPs (380)"),
                 selectInput("gcm_2", "GCM:", choices = unique(go_climdist_df$gcm), selected = "Average across GCMs"),
               ),
               mainPanel(
                 plotOutput("scatter_plot", height = "800px", width = "1350px")
               )
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
      <p>This application visualizes the variability in genomic offset predictions, based on the findings from Archambeau et al. (2024)
      'Evaluating genomic offset predictions in a forest tree with high population genetic structure.'
      Available on <i>bioRxiv</i> at 
      <a href='https://www.biorxiv.org/content/10.1101/2024.05.17.594631v1.full.pdf'>https://www.biorxiv.org/content/10.1101/2024.05.17.594631v1.full.pdf</a>.</p>
      
      <p><b>Plots available:</b>
      
      <p><b>Bumpcharts of population ranks: </b> these charts display population ranks based on their genomic offset predictions for each combination of method and SNP set. 
      Populations with a lower rank have a higher genomic offset, indicating that they may be more at risk under climate change. 
      Colored populations are those having genomic offset values in the top three highest values in at least one combination of SNP set and method.</p>
      
      <p><b>Correlation matrices: </b> these plots show the correlation among genomic offset predictions from the different methods and SNP sets.</p>
      
      <p><b>Scatter plots: </b> these plots allow users to explore relationships between two selected variables, including genomic offset predictions from 
      different methods and SNP sets, and climatic distances. The climatic distances correspond to the absolute difference between future and reference 
      values for each climatic variable, or the Euclidean climatic distance integrating all the selected climatic variables.</p>
      
      Plots can be generated for each Global Climate Model (GCM) or 
      for the mean genomic offset predictions across the five GCMs used in the study (i.e. 'Average across GCMs').</p>"
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

output$scatter_plot <- renderPlot({
  
  df_1 <- go_climdist_df %>% filter(var_name == input$var_1,
                                    method_name == input$method_1,
                                    gcm == input$gcm_1) %>% 
    dplyr::select(pop, val, main_gp_pop, color_main_gp_pop) %>% 
    dplyr::rename(var_1 = val)
  
  df_2 <-  go_climdist_df %>% filter(var_name == input$var_2,
                                     method_name == input$method_2,
                                     gcm == input$gcm_2) %>% 
    dplyr::select(pop, val, main_gp_pop, color_main_gp_pop) %>% 
    dplyr::rename(var_2 = val) %>% 
    left_join(df_1, by=c("pop", "main_gp_pop", "color_main_gp_pop"))
  
  
  ggplot(df_2, aes(x = var_1, y = var_2)) +
    geom_point(aes(color = main_gp_pop), size = 5) +  # Plot points with color by group
    scale_color_manual(values = setNames(df_2$color_main_gp_pop, df_2$main_gp_pop)) + # Custom colors
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") + # Linear fit line
    stat_poly_eq(
      aes(label = paste(..rr.label.., sep = "~~~")), 
      formula = y ~ x, parse = TRUE, size= 15#label.x.npc = "right", label.y.npc = "top"
    ) +
    labs(x = paste0(input$method_1, " - ", input$var_1, " - ", input$gcm_1), 
         y = paste0(input$method_2, " - ", input$var_2, " - ", input$gcm_2), 
         color = "Main gene pool") +
    geom_text_repel(aes(label = pop), size = 5) +  # Add population names with repel to avoid overlap
    theme_bw() + 
    theme(axis.title = element_text(size=20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20))
  
})



}


######################
# Run the application
######################
shinyApp(ui = ui, server = server)

