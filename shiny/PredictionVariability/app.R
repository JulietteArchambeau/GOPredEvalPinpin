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
    
    tabPanel("Bumpcharts across methods and SNP sets",
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
    
    tabPanel("Bumpcharts across GCMs",
               sidebarLayout(
                 sidebarPanel(
                   width = 2,
                   selectInput("method_bumpchartGCM", "Method:", choices = unique(df$method), selected = "GF"),
                   selectInput("snpset_names_bumpchartGCM", "SNP set:", choices = unique(df$snpset_names), selected = "All candidate SNPs (380)")
                 ),
                 mainPanel(
                   plotOutput("bumpchartGCM", height = "900px")
                 )
               )
    ),
    
    tabPanel("Correlation matrices across methods and SNP sets",
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
    
    
    tabPanel("Correlation matrices across GCMs",
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 selectInput("method_corrGCM", "Method:", choices = unique(df$method), selected = "GF"),
                 selectInput("snpset_names_corrGCM", "SNP set:", choices = unique(df$snpset_names), selected = "All candidate SNPs (380)")
               ),
               mainPanel(
                 plotOutput("correlation_plot_GCM", height = "800px", width = "1350px")
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
      <p>This application visualizes the variability in genomic offset predictions, based on the findings from Archambeau et al. (2025)
      'Evaluating genomic offset predictions in a forest tree with high population genetic structure.'
      Available on <i>bioRxiv</i> at 
      <a href='https://www.biorxiv.org/content/10.1101/2024.05.17.594631v1.full.pdf'>https://www.biorxiv.org/content/10.1101/2024.05.17.594631v1.full.pdf</a>.</p>
      
      
      </br>
      <p> The genomic offset was predicted using five methods, seven sets of genetic markers (single nucleotide polymorphisms, SNPs), and five general circulation models. <br>
      All genomic offset predictions were generated at the locations of the studied populations, using the 1901–1950 climate at these locations as the reference climate. <br>
      For future climates, predictions were based on the 2041–2060 period under the shared socio-economic pathway (SSP) 3.7–0. <br>
      For further details on how the genomic offset predictions were generated, please refer to the manuscript.</p>
      
      <p><b>Genomic offset methods</b> 
      <ul>
      <li>Generalized Dissimilarity Modeling (GDM)
      <li>Gradient Forest (GF)
      <li>Latent Factor Mixed Model (LFMM)
      <li>Redundancy Analysis (RDA)
      <li>Partial Redundancy Analysis (pRDA)
      </ul>
        
      <p><b>SNP sets</b> 
      <ul>
      <li>All candidate SNPs (380 SNPs), i.e., SNPs identified by at least one gene-environment association (GEA) method among RDA, pRDA, LFMM, BayPass and GF.
      <li>Candidate SNPs considering population structure correction (221 SNPs), i.e., SNPs identified by at least one GEA method correcting for population structure (pRDA, LFMM and BayPass).
      <li>Common candidate SNPs (69 SNPs), i.e., SNPs identified by at least two GEA methods.
      <li>Control SNPs unmatching allele frequencies (380 SNPs), i.e., SNPs that were randomly sampled among the non-candidate SNPs and with the same number of SNPs as in the set with all candidate SNPs.
      <li>Control SNPs matching allele frequencies (380 SNPs), i.e., SNPs that were sampled among non-candidate SNPs and that have similar allele frequencies and the same number of SNPs as the set with all candidate SNPs.
      <li>All SNPs, i.e., 9817 SNPs.
      <li>SNPs without any missing data, i.e., 3258 SNPs.
      </ul>
      
      <p><b>General circulation models</b> 
      <ul>
      <li>GFDL-ESM4
      <li>IPSL-CM6A-LR
      <li>MPI-ESM1-2-HR
      <li>MRI-ESM2-0
      <li>UKESM1-0-LL
      </ul>
      
      </br>
      </br>
      <p><b> <u>Plots available in the present app:</u> </b>
      
      <p><b>Bumpcharts: </b> these charts display population ranks based on their genomic offset predictions. <br>
      Populations with a lower rank have a higher genomic offset, indicating that they may be more at risk under climate change. <br>
      Colored populations are those with genomic offset ranks in the top three lowest positions (i.e., highest genomic offset values) in at least one
row, where each row represents either a combination of SNP (single nucleotide polymorphism) set and method (bumpcharts across methods and SNP sets) or a general circulation model (bumpchart across GCMs).

  <ul>
    <li><b>Bumpcharts across methods and SNP sets</b> show the variability across methods and SNP sets for a specific general circulation model (GCM).</li>
    <li><b>Bumpcharts across GCMs</b> show the variability across general circulation models for a particular combination of method and SNP set.</li>
  </ul>
  
      </br>
      <p><b>Correlation matrices: </b> these plots show the Pearson correlation coefficients among genomic offset predictions. 
      
  <ul>
    <li><b>Correlation matrice across methods and SNP sets</b> show the variability across methods and SNP sets for a specific Global Climate Model (GCM).</li>
    <li><b>Correlation matrice across GCMs</b> show the variability across GCMs for a particular combination of method and SNP set.</li>
  </ul>
  
      </br>
      <p><b>Scatter plots: </b> these plots allow users to explore relationships between two selected variables, including genomic offset predictions from 
      different methods and SNP sets, and climatic distances. <br>
      The climatic distances correspond to the absolute difference between future and reference 
      values for each climatic variable, or the Euclidean climatic distance integrating all the selected climatic variables.
      

      
      "
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
    ylab("Population rank (low rank = high genomic offset)") +
    theme(#panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=11),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank()) 
  
    })
  


output$bumpchartGCM <- renderPlot({
  
  plot_df <- df %>% filter(method == input$method_bumpchartGCM & snpset_names == input$snpset_names_bumpchartGCM)
  
  high_go_pops <- plot_df %>% filter(rank < 4) %>% pull(pop) %>% unique()    
  
  my_palette <- c(my_colors[1:length(high_go_pops)],"#E8E8E8")
  
  sub <- plot_df %>% 
    mutate(flag = ifelse(pop %in% high_go_pops, TRUE, FALSE),
           pop_col = if_else(flag == TRUE, pop, "Others")) %>% 
    mutate(pop = factor(pop, levels=c(setdiff(unique(plot_df$pop), high_go_pops), high_go_pops)),
           pop_col = factor(pop_col, levels = c(high_go_pops,"Others")))
  
  sub %>% ggplot(aes(x = gcm, y = rank, group = pop)) +
    coord_flip() +
    scale_y_continuous(breaks = 1:34, minor_breaks = 1:34) +
    geom_point(aes(color = pop_col), size = 2, alpha = 0.9) +
    geom_line(aes(color = pop_col), linewidth = 2, alpha = 0.8) +
    scale_color_manual(values = my_palette) +
    geom_text(data = sub %>% filter(gcm == last(levels(factor(sub$gcm)))),
              # to only write the names of the highlighted pops
              #data = sub %>% filter(method_snpset_names == last(levels(factor(sub$method_snpset_names))) & pop %in% high_go_pops), 
              aes(label = pop, x = length(unique(sub$gcm)) + 0.1), 
              color = "gray20", 
              size = 4, 
              angle = 40) +
    theme_bw() +
    ylab("Population rank (low rank = high genomic offset)") +
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

output$correlation_plot_GCM <- renderPlot({
  
  correlations <- df %>% 
    filter(method == input$method_corrGCM & snpset_names == input$snpset_names_corrGCM) %>% 
    pivot_wider(id_cols = pop, values_from = go, names_from = gcm)  %>% 
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

