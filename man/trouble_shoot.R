
library(trajiq)
library(tidyverse)
contrast_variables = c('tissue')
my_features = c('PDL1','CTLA4','PD1')

# scaled_cell_table =readRDS('/Users/meelad/Downloads/scaled_celltable.RDS')
# glm_test = scaled_cell_table %>%
#      # dplyr::filter(annotation_subset  ==  'Naive' &  annotation_lin_sub %in%  c('B Cells_Naive') | annotation_subset == 'DN T Cells' & annotation_lin_sub ==  'DN T Cells_DN T Cells') %>%
#      group_by(patient_id, annotation_lin_sub, annotation_subset,tissue , mouse,ed)  %>%
#      group_map(~downsampleWith_group_by(., 50),.keep = T) %>% bind_rows()
# saveRDS(glm_test,'~/Downloads/sampled_inputkb.rds')
glm_test =readRDS('~/Downloads/sampled_inputkb.rds')


glm_test=glm_test %>%
     mutate(patient_id = factor(patient_id))
dif_res_by_ct = trajiq::differential_analysis_program(glm_input = glm_test ,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = 'patient_id',intercept = TRUE,contrast_method = 'emm',
                                                      SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))
dif_res_by_ct

dif_res_by_ct = trajiq::differential_analysis_program(glm_input = glm_test2 ,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = 'patient_id',intercept = TRUE,contrast_method = 'emm',
                                                      SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))
dif_res_by_ct




dif_res_by_ct = trajiq::differential_analysis_program(glm_input = glm_test ,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = 'patient_id',intercept = TRUE,contrast_method = 'tukey',
                                                      SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))
dif_res_by_ct


###
######
#########



glm_test2 = glm_test  %>%
     dplyr::filter(annotation_subset == 'DN T Cells' & annotation_lin_sub ==  'DN T Cells_DN T Cells')#  %>%
     # bind_rows(.,.)  %>%
     # bind_rows(.,.)

glm_resas = trajiq:::construct_model(glm_input = glm_test2,outcome_features = c('CD4'), contrast_variables = c('tissue'), covariates_in_model = 'patient_id',intercept = TRUE)

glm_test2$group = glm_resas$group
my_group = unique(glm_test2$group)

## filter  to only include contrasts of interest (comparisons)
my_tukeywukey_combos = trajiq:::create_TukeyWukey_combos(levels(my_group))

## filter  to only include contrasts of interest (comparisons)
tukey_res_test = trajiq:::TukeyWukey(myGLM = glm_res$myGLM,combos_df =  my_tukeywukey_combos,contrast_variables = contrast_variables)


# multitest_res <- multcomp::glht(glm_resas$myGLM$CD4, mcp = "Tukey")
# multitest_res <- multcomp::glht(glm_res$myGLM$CD4, linfct = c("`groupLN` - `groupPL`== 0"))



emm_df = pairs(emmeans::emmeans(glm_resas$myGLM[[1]], "group")) %>%
     data.frame() %>%
     dplyr::rename(comparison = contrast)


multcomp::glht(lmod, linfct = c("Agriculture = 0",
                                "Examination = 0",
                                "Education = 0",
                                "Catholic = 0",
                                "Infant.Mortality = 0"))
