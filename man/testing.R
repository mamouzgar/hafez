
library(trajiq)
library(tidyverse)
# contrast_variables = c('tissue')
# my_features = c('PDL1','CTLA4','PD1')

scaled_cell_table =readRDS('/Users/meelad/Downloads/scaled_celltable.RDS')
# glm_test = scaled_cell_table %>%
#      # dplyr::filter(annotation_subset  ==  'Naive' &  annotation_lin_sub %in%  c('B Cells_Naive') | annotation_subset == 'DN T Cells' & annotation_lin_sub ==  'DN T Cells_DN T Cells') %>%
#      group_by(patient_id, annotation_lin_sub, annotation_subset,tissue , mouse,ed)  %>%
#      group_map(~downsampleWith_group_by(., 50),.keep = T) %>% bind_rows()
# saveRDS(glm_test,'~/Downloads/sampled_inputkb.rds')
# glm_test =readRDS('~/Downloads/sampled_inputkb.rds')
#
#
# glm_test=glm_test %>%
#      mutate(patient_id = factor(patient_id))
# dif_res_by_ct = trajiq::differential_analysis_prog

contrast_variables <- c("ed_bin", "tissue")
my_features <- c("CD115", "CD69", "CD25", "CD64", "CD80", "CD11b", "CD40", "CTLA4", "Ly6C", "CD194", "CD62L", "PDL1", "CD44", "PD1", "MHCII", "CD86")
covariates_in_model=c('patient_id')
scaled_cell_table =readRDS('/Users/meelad/Downloads/scaled_celltable.RDS')

# Downsamples to a maximum of 1000 events per features specified in group_by
glm_input_sampled <- scaled_cell_table %>%
     group_by(patient_id, tissue, annotation_lin_sub, ed_bin) %>%
     group_map(~downsampleWith_group_by(., 100),.keep = T) %>% bind_rows()
#
# glm_test <- glm_input_sampled %>%
#      mutate(patient_id = factor(patient_id))

glm_output <- differential_analysis_program(
     glm_input = glm_input_sampled,
     outcome_features = my_features,
     contrast_variables = contrast_variables,
     covariates_in_model = NULL,
     intercept = TRUE,
     contrast_method = 'emm',
     SPLIT_BY_NAMES = c("annotation_lin_sub", "merge1"))
saveRDS(glm_output,'~/Downloads/test_filt.rds')




# #
# filt_output = readRDS('~/Downloads/filteredoutput.RDS')
# until_output = readRDS('~/Downloads/unfilteredoutput_downsampled.RDS')
# glm_tissue_edbin <- filt_output %>%
#      mutate(tissue1 =gsub('\\(|\\)','',tissue1),
#             tissue2 =gsub('\\(|\\)','',tissue2),) %>%
#      filter(tissue1 == tissue2)
#
# glm_tissue_edbin_unfilt <- until_output %>%
#      mutate(tissue1 =gsub('\\(|\\)','',tissue1),
#             tissue2 =gsub('\\(|\\)','',tissue2),) %>%
#      filter(tissue1 == tissue2)
# #
# #
# ggplot(glm_tissue_edbin_unfilt, aes(x= estimate, y =  -1*log10(pval))) +
#      geom_point()
#
#
#
glm_output = readRDS('~/Downloads/test_filt.rds')
glm_output_filt = glm_output %>%
     mutate(tissue1 =gsub('\\(|\\)','',tissue1),
            tissue2 =gsub('\\(|\\)','',tissue2)) %>%
     filter(tissue1 == tissue2)
ggplot(glm_output_filt, aes(x= estimate, y =  -1*log10(pval))) +
     geom_point()
# #
# # distinct(glm_input_sampled, ed_bin, tissue)
# # dim(distinct(glm_input_sampled, ed_bin, tissue))
# # dim(distinct(scaled_cell_table,ed_bin, tissue))
#
#
#
md = distinct(glm_input_sampled, mouse, patient_id, ed_bin, tissue)
# ggplot(md, aes(x = patient_id, y = ed_bin)) + geom_point() +
#      facet_wrap(~tissue)
#
ggplot(md %>% group_by(ed_bin, tissue) %>% summarize(count = n()), aes(x = tissue, y = count, fill = ed_bin)) +
     geom_col()



omg_summary = glm_input_sampled %>%
     group_by(patient_id, mouse , ed_bin, tissue ) %>%
     summarize_at(my_features, median)

my_features

ggplot(omg_summary%>% dplyr::filter(tissue =='PB') , aes(x = ed_bin,  y= CD40)) +
     geom_point()+
     geom_line(aes(color = mouse, group = patient_id))
# pheatmap::



ggplot(md, aes(x=  patient_id, y = tissue ))+
     geom_tile(fill = 'blue') +
     ggpubr::rotate_x_text(90)+
     facet_wrap(~ed_bin, ncol = 1)

