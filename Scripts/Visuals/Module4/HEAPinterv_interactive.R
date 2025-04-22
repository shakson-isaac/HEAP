#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
library(plotly)
library(htmlwidgets)
library(qs)

#Load INT structure:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPint <- qread("./Output/HEAPres/HEAPint.qs")

###### ScatterPlot to Compare Effects across UKB and Interventions 

# Function to create plot for each study
create_plot_interact <- function(df, x_col, y_col, eName, effect_name, cor_val, p_val, se_col) {
  
  p <- df %>%
    filter(!is.na(Estimate)) %>%
    mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
    ggplot(aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(text = paste("<br>Prot:", EntrezGeneSymbol,
                                "<br>Exposure:", paste0(!!sym(x_col),"(",signif(!!sym(x_col) - 1.96 * `Std. Error`, digits = 2),
                                                        ",",
                                                        signif(!!sym(x_col) + 1.96 * `Std. Error`,digits = 2),")"),
                                "<br>Intervention:", paste0(!!sym(y_col),"(",signif(!!sym(y_col) - 1.96 * !!sym(se_col), digits = 2),
                                                            ",",
                                                            signif(!!sym(y_col) + 1.96 * !!sym(se_col), digits = 2),")")
    ))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
    # Use dynamic `se_col` for error bars
    geom_linerange(aes(ymin = !!sym(y_col) - 1.96 * !!sym(se_col), ymax = !!sym(y_col) + 1.96 * !!sym(se_col)), 
                   color = "black", size = 0.5, alpha = 0.5) +
    geom_linerange(aes(xmin = Estimate - 1.96 * `Std. Error`, xmax = Estimate + 1.96 * `Std. Error`), 
                   color = "black", size = 0.5, alpha = 0.5) +
    theme_minimal() +
    labs(x = paste0("Beta (", eName, ") - UKB"), y = paste0("Effect Size - ", effect_name)) + 
    geom_text_repel(aes(label = EntrezGeneSymbol), color = "red", size = 3, max.overlaps = 10)
  
  # Get axis limits
  plot_build <- ggplot_build(p)
  x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
  y_range <- plot_build$layout$panel_scales_y[[1]]$range$range
  
  # Dynamically position the correlation label in the top-left corner
  x_pos <- x_range[1] + 0.15 * (x_range[2] - x_range[1])  # 15% from the left
  y_pos <- y_range[2] - 0.05 * (y_range[2] - y_range[1])  # 5% from the top
  
  # Add the annotation
  p1 <- p + annotate("text",
                     x = x_pos, 
                     y = y_pos, 
                     label = paste0("R = ", signif(cor_val, 2), ", p = ", signif(p_val, 2)), 
                     color = "black")
  
  # lapply(c("png", "svg", "pdf"), 
  #        function(fmt) ggsave(paste0("./Figures/HEAP/F5/",
  #                                    make.names(eName),"_", make.names(effect_name), ".", fmt), 
  #                             plot = p1, width = 6, height = 4, dpi = 500))
  
  return(p1)
}

# Save HTML efficiently:
savePlotly <- function(gg, html_file){
  gginteract <- ggplotly(gg, tooltip = "text") 
  #%>% partial_bundle(local = FALSE)
  
  # Save the Plotly widget without embedding resources (Plotly will be fetched from the CDN)
  saveWidget(gginteract, 
             file = html_file,
             selfcontained = FALSE)
  
  # Make sure all the html files refer to static package folder:
  html_content <- readLines(html_file)
  
  fname <- basename(html_file)
  fname <- gsub(".html","_files", fname)
  html_content <- gsub(fname, 'static', html_content)  # Removes fill.css
  
  # Write the modified HTML content back to the file
  writeLines(html_content, html_file)
}

# Function to create and save scatterplots
IntCor_Plotapp <- function(scatList, cList, pvalList, eID, eName, CSpec, sID){
  #Filter the UKBscatter dataframe by eID:
  UKBscatter <- scatList[[CSpec]] %>% filter(ID == eID)
  
  # Prepare the data for plots
  UKBint_cor <- cList[[CSpec]]
  UKBint_p <- pvalList[[CSpec]]
  
  corVal_HERITAGE <- UKBint_cor %>% filter(ID == eID) %>% pull(HERITAGE_effect)
  corVal_GLP1.1 <- UKBint_cor %>% filter(ID == eID) %>% pull(GLP1_effect1)
  corVal_GLP1.2 <- UKBint_cor %>% filter(ID == eID) %>% pull(GLP1_effect2)
  
  pVal_HERITAGE <- UKBint_p %>% filter(ID == eID) %>% pull(HERITAGE_effect)
  pVal_GLP1.1 <- UKBint_p %>% filter(ID == eID) %>% pull(GLP1_effect1)
  pVal_GLP1.2 <- UKBint_p %>% filter(ID == eID) %>% pull(GLP1_effect2)
  
  # Specify standard error column for each study (adapt this based on the dataset)
  se_col_HERITAGE <- "HERITAGE_se"
  se_col_GLP1_STEP1 <- "GLP1_se1"
  se_col_GLP1_STEP2 <- "GLP1_se2"
  
  # Create plots
  gg1 <- create_plot_interact(UKBscatter, "Estimate", "HERITAGE_effect", eName, "HERITAGE", 
                     corVal_HERITAGE, pVal_HERITAGE, se_col_HERITAGE)
  gg2 <- create_plot_interact(UKBscatter, "Estimate", "GLP1_effect1", eName, "GLP1 STEP1", 
                     corVal_GLP1.1, pVal_GLP1.1, se_col_GLP1_STEP1)
  gg3 <- create_plot_interact(UKBscatter, "Estimate", "GLP1_effect2", eName, "GLP1 STEP2", 
                     corVal_GLP1.2, pVal_GLP1.2, se_col_GLP1_STEP2)
  
  savePlotly(gg1, paste0("./Output/App/Interactive/A3/",
                         sID,"_","HERITAGE",".html"))
  savePlotly(gg2, paste0("./Output/App/Interactive/A3/",
                         sID,"_","GLP1_STEP1",".html"))
  savePlotly(gg3, paste0("./Output/App/Interactive/A3/",
                         sID,"_","GLP1_STEP2",".html"))
}


#For loop through these exposures:
ExposureDict <- data.frame(
  origID = rownames(HEAPint@cList$Model6),
  sID =  paste0("E",1:length(rownames(HEAPint@cList$Model6))),
  name = c("White Bread Intake",
           "Exercise: Swimming, Cycling, etc.",
           "Incr. Alcohol Intake (vs. 10 yrs ago)",
           "Cereal Intake",
           "Former Daily Smoking",
           "Fresh Fruit Intake",
           "Smoked in Lifetime (Yes)",
           "Never Smoked",
           "Beef Intake 1x/WK",
           "MET min/WK of Vigorous Activity",
           "# Days/WK of Vigorous Activity",
           "Exercise: Strenuous Sports")
)


fwrite(ExposureDict, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Interactive/A3/exposure_interv.csv")

### Create plots across relevant significant exposure and intervention study pairs:
lapply(1:nrow(ExposureDict), function(x){
  IntCor_Plotapp(HEAPint@sList, HEAPint@cList, HEAPint@pList,
                 eID = ExposureDict[x,"origID"], 
                 eName = ExposureDict[x,"name"], 
                 CSpec = "Model6", 
                 sID = ExposureDict[x,"sID"])
})


