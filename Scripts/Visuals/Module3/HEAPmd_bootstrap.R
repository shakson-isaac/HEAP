#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(pbapply)
library(dplyr)
library(purrr)
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

#Load the results:
# Load Files w/ Progress Bar:
loadMDAssocBoot <- function(covarType,nPairs){
  error_files <- list()  # Track files w/ errors
  stat_all <- list()
  
  # Use pblapply for progress bar and iteration
  for(i in 1:nPairs){
    tryCatch({
      load <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_mediationboot/",
                                covarType,"/MDboot_",i,".txt"))
      
      stat_all[[i]] <- load
      
    }, error = function(e) {
      message(paste("Error reading file:", i))
      error_files <<- append(error_files, i)
    })
  }
  
  
  stat_all <- do.call("rbind", stat_all)
  return(stat_all)
}

# Load results:
Type1 <- loadMDAssocBoot(covarType = "Type1", 48)
Type2 <- loadMDAssocBoot(covarType = "Type2", 48)
Type3 <- loadMDAssocBoot(covarType = "Type3", 48)
Type4 <- loadMDAssocBoot(covarType = "Type4", 48)
Type5 <- loadMDAssocBoot(covarType = "Type5", 48)


# Function to plot mediation results w/ Bootstrap.
mediation_plot <- function(MDboot, dz, protIDs, title){
  bootstrap_res_upd <- MDboot %>% #Specify mediation results: Disease code + proteins
    filter(DZid %in% dz & ID %in% protIDs) %>% 
    select(all_of(c("ID","DZid",
                    "E_HRi_orig","E_HRi_2.5","E_HRi_97.5",
                    "G_HRi_orig","G_HRi_2.5","G_HRi_97.5"))) %>%
    pivot_longer(
      cols = -c(ID,DZid),
      names_to = c("Type",".value"),
      names_pattern = "(E|G)_HRi_(.*)"
    )
  
  gg1 <- bootstrap_res_upd  %>%
    ggplot(aes(x = ID, y = orig, color = Type)) +
    geom_point(position = position_dodge(width = 0.6), size = 1) +  # Point estimate
    geom_errorbar(aes(ymin = `2.5`, ymax = `97.5`), width = 0.2,    # Error bars
                  position = position_dodge(width = 0.6))  +
    geom_hline(yintercept = 1, color = "black", size = 0.8, linetype = "dashed") +  # Reference line at 1
    #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
    scale_color_discrete(labels = c("E", "G")) +
    labs(x = "Proteins", y = "HR: Indirect Effects") +
    coord_flip() +
    theme_minimal() +
    ggtitle(title)
  gg1
  return(gg1)
}
mediation_plot_save <- function(gg1, title){
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/F4/", "Mediation",title,".", fmt), 
                              plot = gg1, width = 5, height = 4, dpi = 500))
  
}


# Type 2 Diabetes:
T2Dplot <- mediation_plot(Type5,
               dz = "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
               protIDs = c("IGFBP1","IGFBP2","CKB","LPL",
                           "FABP4","FGF21","LEP","IL18R1","IGSF9"),
               title = "E11: Type 2 Diabetes")
mediation_plot_save(T2Dplot, title = "E11: Type 2 Diabetes")


# Emphysema:
Eplot <- mediation_plot(Type5,
                          dz = "age_j43_first_reported_emphysema_f131490_0_0",
                          protIDs = c("CXCL17","LAMP3","WFDC2","TNR","GDF15"),
                          title = "J43: Emphysema")
mediation_plot_save(Eplot, title = "J43: Emphysema")


# Mental Disorders - Alcohol Use:
Alcplot <- mediation_plot(Type5,
                          dz = "age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0",
                          protIDs = c("MAMDC4","CEACAM16"),
                          title = "F10: Mental Disorders - Alcohol Use")
mediation_plot_save(Alcplot, title = "F10: Mental Disorders - Alcohol Use")
