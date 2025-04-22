#'*Option*
#SAVE AS VECTOR GRAPHICS
#Examples:
#library(svglite)
#ggsave(paste0(".svg"),plot = E_ds1, width = 5, height = 3, dpi = 500)
#ggsave(paste0(".pdf"), plot = E_ds1, width = 5, height = 3, dpi = 500)
##Visualization THINGS: TITLES DO I NEED TO SPECIFY THE 'TEST SET' ALL THE TIME

#'*GOALS*
#Things to look at: Distribution of Total R2 from lasso
#Train vs Test R2 across covariate specifications
#Central tendency of R2 for G, E, and GxE terms
#Function types: summarise vs summarize
# purrr::map2 iterate over the list created to generate one dataframe:

#'*Additional Notes*'
#State in paper polygenic GxE is not well powered in this study
#- given analysis of univariate tests.



#'*Runtime Comments on Results from HEAP: using nested CV*
#'General Comments about Runtime of Covariate Specifications:
#'Type 1-3 - works with 1 hour, 10 minutes
#'Type 4-5 - needs more time to run: set to maybe 1 hour, 30 minutes.