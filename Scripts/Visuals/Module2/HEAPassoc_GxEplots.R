#### Get Plots for Significant GxE Associations ####
#Source univarInteract Script:
source("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/runProt_PGS_univarInteract_v2.R")

#'*To Modify Later --> For # of days,etc. that got converted into z-score.*
#'*Discretize the above by days and find the respective quintiles.*

#Run specific associations:
saveGWISplot_app <- function(o,e,c,n,l=c(),p,f){
  gg <- runGWISplot(omic = o,
                    exposure = e,
                    CType = c,
                    ename = n,
                    levels = l,
                    plotname = p)
  ggsave(paste0("./Output/App/Interactive/P1/",paste0(o,"_",f),".png"), 
         plot = gg, width = 5, height = 4, dpi = 1000)
}

#Find top E hits that REPLICATE to show:
GxEappType5 <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec5GxEassocReplicated.txt")
GxEappType6 <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec6GxEassocReplicated.txt")


GxEids <- GxEappType6 %>%
  select(c("Eid_train","omicID"))

IDmatching <- data.frame(
  ID = c("alcohol_intake_frequency_f1558_0_0",
         "plays_computer_games_f2237_0_0",                                                                                     
         "types_of_transport_used_excluding_work_f6162_0_0.multi_Public_transport",                                            
         "alcohol_drinker_status_f20117_0_0_Current",                                                                          
         "alcohol_drinker_status_f20117_0_0_Never",                                                                            
         "bread_intake_f1438_0_0",                                                                                             
         "coffee_intake_f1498_0_0",                                                                                            
         "pork_intake_f1389_0_0",                                                                                              
         "time_spent_driving_f1090_0_0",                                                                                       
         "current_tobacco_smoking_f1239_0_0_No",                                                                               
         "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",                                                         
         "smoking_status_f20116_0_0_Current",                                                                                  
         "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Heavy_DIY_.eg._weeding._lawn_mowing._carpentry._digging.",
         "types_of_transport_used_excluding_work_f6162_0_0.multi_Car.motor_vehicle",                                           
         "water_intake_f1528_0_0",                                                                                             
         "cheese_intake_f1408_0_0"),
  Name = c("Alcohol Intake Freq.",
           "Plays Computer Games",
           "Uses Public Transport",
           "Current Alcohol Drinker",
           "Never Drinks Alcohol",
           "Bread Intake",
           "Coffee Intake",
           "Pork Intake",
           "Time Spent Driving",
           "Non Tobacco Smoker",
           "Current Tobacco Smoker (Most Days)",
           "Current Smoker",
           "Physical Activity: weeding, lawn mower, carpentry",
           "Motor Vehicle Transport",
           "Water Intake",
           "Cheese Intake")
)

GxEmatch <- merge(GxEids, IDmatching, by.x = "Eid_train", by.y = "ID")


### Find the scoring for specific IDs.
IDscores <- data.frame(
  ID = c("alcohol_drinker_status_f20117_0_0_Current",
         "alcohol_drinker_status_f20117_0_0_Never",
         "alcohol_intake_frequency_f1558_0_0",
         "bread_intake_f1438_0_0",
         "cheese_intake_f1408_0_0",
         "coffee_intake_f1498_0_0",
         "current_tobacco_smoking_f1239_0_0_No",
         "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
         "pork_intake_f1389_0_0",
         "plays_computer_games_f2237_0_0",
         "smoking_status_f20116_0_0_Current",
         "time_spent_driving_f1090_0_0",
         "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Heavy_DIY_.eg._weeding._lawn_mowing._carpentry._digging.",
         "types_of_transport_used_excluding_work_f6162_0_0.multi_Car.motor_vehicle",
         "types_of_transport_used_excluding_work_f6162_0_0.multi_Public_transport",
         "water_intake_f1528_0_0"),
  Scores = I(list(c("NO","YES"),
                  c("NO","YES"),
                  c("Never","Special Occasions",
                    "1-3 times monthly","1-2 times weekly",
                    "3-4 times weekly","Daily"),
                  c(),
                  c(),
                  c(),
                  c("NO","YES"),
                  c("NO","YES"),
                  c(),
                  c("Never","Sometimes","Often"),
                  c("NO","YES"),
                  c(),
                  c("NO","YES"),
                  c("NO","YES"),
                  c("NO","YES"),
                  c()))
)

GxEmatchv2 <- merge(GxEmatch, IDscores, by.x = "Eid_train", by.y="ID")

lapply(1:nrow(GxEmatch), function(x){
  oID <- GxEmatchv2[x, ][["omicID"]]
  eID <- GxEmatchv2[x, ][["Eid_train"]]
  nID <- GxEmatchv2[x, ][["Name"]]
  lID <- GxEmatchv2[x,][["Scores"]][[1]]
  
  saveGWISplot_app(o = oID,
                   e = eID,
                   c = "Type6",
                   n = nID,
                   l = lID,
                   p = "Polygenic GxE Interaction",
                   f = paste0(nID,"_GxE"))
}) 



c("Never","Special Occasions",
  "1-3 times monthly","1-2 times weekly",
  "3-4 times weekly","Daily")


saveGWISplot_app(o = "APOF",
                 e = "alcohol_intake_frequency_f1558_0_0",
                 c = "Type6",
                 n = "Alcohol Intake Freq.",
                 l = c("Never","Special Occasions",
                       "1-3 times monthly","1-2 times weekly",
                       "3-4 times weekly","Daily"),
                 p = "Polygenic GxE Interaction",
                 f = "Alcohol_GxE")








