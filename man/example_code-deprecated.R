
library(tidyverse)
library(hafez)

features_cellcycle = c('Ki67','pRbS780', 'pH3_s10','CDT1', 'IdU', 'Geminin', 'PLK1', 'DNA',  'CyclinB1', 'PCNA', 'SLBP')%>% make.names()
features_other = c( 'EZH2','HH3','CyclinE','FoxM1','CDC20','RUNX1','CTCF','H3K9ac','H3K27me3','H3K18ac','WGA','MLL1')

# output_files = paste0('/Users/meelad/phd-projects/cell_cycle-trajectory_reconstruction/analysis/scripts/2023-02-01_CellLines/v2_data/analysis-ready/')
# CELL_LINE_TRAJECTORY_INFERENCE_df =  readRDS(paste0(output_files,'CELL_LINE_TRAJECTORY_INFERENCE_df-final_object_v1_2024.rds') )

example_cytof_data= CELL_LINE_TRAJECTORY_INFERENCE_df %>%
     dplyr::filter(!is.na(gate)) %>%
     group_by(cell_line) %>%
     sample_n(5000) %>%
     ungroup() %>%
     dplyr::filter(cell_line == 'JURKAT') %>%
     dplyr::select(all_of(c(features_cellcycle,features_other)), cell.id, cell_line, PC1, PC2, gate ) %>%  ungroup
save(example_cytof_data, file = '/Users/meelad/Rpackages/hafez/data/cytof_example_data.rda')

#Î
## load example data
glm_input=hafez::glm_input

## features of interest
my_features = colnames(glm_input)[-c(1:3)]
covariates_in_model = c('mouse') ## covariates of interest
contrast_variables = c('tissue') ## variable you want to compare

###############################################################################################
## (A) example for running analysis on glm_input features that look like celltype abundance fractions
###############################################################################################
## note: be sure to make covariates as factor variables when it's not a continuous variable. Eg, mouse ids as a covariate for matched analysis
glm_input = glm_input%>%
     mutate(mouse = factor(mouse))
dif_res =hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm',)
head(dif_res)



###############################################################################################
## (B) to run interactions, you would run something like this:
###############################################################################################
contrast_variables = c('tissue','ed')
dif_res =hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm',)

###############################################################################################
## (C) example for running analysis on glm_input features that look like marker expression
###############################################################################################
glm_input =hafez::ct %>% mutate(mouse = factor(mouse)) %>%
     bind_rows(hafez::ct %>% mutate(mouse = factor(mouse)))
contrast_variables = c('tissue')
my_features = c('CD4','CD40','FOXP3') ## make sure these are in the object (no spelling errors)
# my_features %in% colnames(glm_input)
dif_res = hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm',)

###############################################################################################
# (D) example code to split across multiple subsets of interest
###############################################################################################
dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm',
                                                      SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))


## you can also use the tukey method instead of emm
dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'tukey',
                                                      SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))


###############################################################################################
# (E) If your single-cell dataset is very large, we can speed up the process by downsampling.
###############################################################################################
glm_input_sampled = glm_input %>%
     group_by(patient_id, tissue, annotation_lin_sub) %>%
     group_map(~downsampleWith_group_by(., 50),.keep = T) %>% bind_rows()
dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input_sampled ,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = 'patient_id',intercept = TRUE,contrast_method = 'emm',
                                                      SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))











