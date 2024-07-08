
library(tidyverse)
library(hafez)
#
# features_cellcycle = c('Ki67','pRbS780', 'pH3_s10','CDT1', 'IdU', 'Geminin', 'PLK1', 'DNA',  'CyclinB1', 'PCNA', 'SLBP')%>% make.names()
# features_other = c( 'EZH2','HH3','CyclinE','FoxM1','CDC20','RUNX1','CTCF','H3K9ac','H3K27me3','H3K18ac','WGA','MLL1')
#
# output_files = paste0('/Users/meelad/phd-projects/cell_cycle-trajectory_reconstruction/analysis/scripts/2023-02-01_CellLines/v2_data/analysis-ready/')
# CELL_LINE_TRAJECTORY_INFERENCE_df =  readRDS(paste0(output_files,'CELL_LINE_TRAJECTORY_INFERENCE_df-final_object_v1_2024.rds') )
# phate_res = readRDS(paste0( '/Users/meelad/phd-projects/T_cell-CellCycle-trajectory_analysis/2023-05-01_Jurkat_drug_sync_expr_NolanCyTOF/data/analysis-ready//2024-4-18_drug_tx_phate_object.rds'))
# #
# features_cellcycle  = c( 'dsDNA', 'pH3_s10','CDT1', 'Geminin', 'PLK1', 'DNA2', 'CyclinB1', 'PCNA', 'Ki67','SLBP','IdU','pRb_S780')%>% unique()
# example_cytof_data= phate_res %>%
#      dplyr::filter(!is.na(gate)) %>%
#      dplyr::filter(timepoint %in% c('0') |( condition == 'WT'  )) %>%
#      mutate(timepoint = 0) %>%
#      group_by(timepoint, condition) %>%
#      sample_n(5000) %>%
#      ungroup() %>%
#      # dplyr::filter(cell_line == 'JURKAT') %>%
#      dplyr::select(all_of(c(features_cellcycle)), cell.id, gate=traj_phase, condition, timepoint, pseudotime =pseudotime_adjusted ) %>%
#      ungroup() %>%
#      mutate(gate = factor(gate, levels = c('G1','S','G2','M')),
#             condition = factor(condition, levels =c('WT','PALBO','HU','NOC')))
# save(example_cytof_data, file = '/Users/meelad/Rpackages/hafez/data/cytof_example_data.rda')
features_cellcycle  = c( 'dsDNA', 'pH3_s10','CDT1', 'Geminin', 'PLK1', 'DNA2', 'CyclinB1', 'PCNA', 'Ki67','SLBP','IdU','pRb_S780')%>% unique()

Phases_gate =  c('G1','S','G2','M')
Phases_gate_Colors = viridis::magma(n=4)
# allPhasesColors[1] = "#000000"
names(Phases_gate_Colors) = Phases_gate

colScale_phases_manualAll = scale_color_manual(values = Phases_gate_Colors,name = 'Cell cycle\nphase')
fillScale_phases_manualAll = scale_fill_manual(values = Phases_gate_Colors, name = 'Cell cycle\nphase')


example_cytof_data = hafez::example_cytof_data %>%
     mutate_at(features_cellcycle, scale)


############################################################
## DBPN: density based pseudotime normalization
############################################################
example_cytof_data = DBPN(dataset = example_cytof_data, dataset.subset_to_use = example_cytof_data %>%dplyr::filter(condition == 'WT'),column_to_normalize = 'pseudotime',adjust.value = 0.5)

ggplot(example_cytof_data, aes(x = PSEUDOTIME_NORMALIZED, y = pseudotime))+
     theme_bw() +
     theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
     geom_point(aes(color = gate)) +
     colScale_phases_manualAll


############################################################
## calculate cell density by group
############################################################
density_by_group_results = calculate_cell_density_by_group(dataset = example_cytof_data, group_splits = c('condition','timepoint'),pseudotime_column_name = 'PSEUDOTIME_NORMALIZED')
myRug = compute_rug(dataset = example_cytof_data,pseudotime_column = 'PSEUDOTIME_NORMALIZED',group_column = 'gate',num_bins = 100,recompute_pseudotime_bins = T) %>%
     mutate(group = factor(group, levels =c('G1','S','G2','M')))

ggplot(density_by_group_results, aes(x = x, y=y, color = condition))+
     theme_bw() +
     theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
     geom_line(size = 1)+
     ggnewscale::new_scale_colour()+
     viridis::scale_color_viridis(option = 'magma', discrete = T) +
     geom_rug(data = myRug, aes(x = pseudotime, color  = group), size =1,inherit.aes = F)


# example_cytof_data = hafez::example_cytof_data
## requires a cell.id column

######################################################################################################
## perform pseudotime mapped mahalanobis distance  without and with landmark cells (WT condition)
######################################################################################################
## perform desired discetization of pseudotime bins
example_cytof_data = computePseudotimeBins(dataset = example_cytof_data,pseudotime_column = 'PSEUDOTIME_NORMALIZED',interval_sequence = seq(0,1,0.05))

pm_mahalanobis_results_all_cells = pseudotime_mapped_mahalanobis_analysis(dataset = example_cytof_data,method = 'all',features = features_cellcycle, reference_group_column = NULL,reference_group_name = NULL,CELL_COUNT_THRESHOLD = 0)
pm_mahalanobis_results_wt_landmarks = pseudotime_mapped_mahalanobis_analysis(dataset = example_cytof_data,method = 'landmark', features = features_cellcycle, reference_group_column = 'condition',reference_group_name = 'WT', CELL_COUNT_THRESHOLD = 0)


example_cytof_data_mahalanobis_results = example_cytof_data %>%
     dplyr::left_join(pm_mahalanobis_results_all_cells %>% dplyr::rename(mahalanobis_distance_all = mahalanobis_distance)) %>%
     dplyr::left_join(pm_mahalanobis_results_wt_landmarks %>% dplyr::rename(mahalanobis_distance_lm = mahalanobis_distance))


##################################################################################
## identify differential features of noncaonical cells ##
##################################################################################
## identify differential features between noncanonical and canonical cells in drug systems.
glm_input = example_cytof_data_mahalanobis_results %>%
     computePseudotimeBins(dataset = .,pseudotime_column = 'PSEUDOTIME_NORMALIZED',interval_sequence = seq(0,1,0.5)) %>%
     mutate(noncanonical = ifelse(mahalanobis_distance_lm > quantile(example_cytof_data_mahalanobis_results$mahalanobis_distance_lm, 0.75,na.rm =T), 'noncanonical', 'canonical'))
dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input ,outcome_features = features_cellcycle, contrast_variables = 'noncanonical',intercept = TRUE,contrast_method = 'emm',
                                                     SPLIT_BY_NAMES = c('condition', 'pseudotime_bins'))
dif_res_by_ct = dif_res_by_ct #%>%
     # dplyr::filter(!is.na(comparison))

library(patchwork)
ggplot(dif_res_by_ct , aes(x = estimate, y = minus_log10padj, color = condition))+
     geom_point()+
     geom_vline(xintercept= c(0.5,-0.5), linetype = 'dashed')+
     geom_vline(xintercept= c(0.5,-0.5), linetype = 'dashed')+
     ggrepel::geom_text_repel(data =dif_res_by_ct %>% dplyr::filter(abs(estimate)>=0.5, padj<=0.1), aes(label = feature), size = 2.5) +
     facet_wrap(~pseudotime_bins)

## differential analysis can also  be used to compare conditions directly
dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input ,outcome_features = features_cellcycle, contrast_variables = 'condition',intercept = TRUE,contrast_method = 'emm',
                                                     SPLIT_BY_NAMES = c( 'pseudotime_bins'))
dif_res_by_ct_filt = dif_res_by_ct %>%
     dplyr::filter(condition1 =='WT'|condition2 == 'WT')
ggplot(dif_res_by_ct_filt , aes(x = estimate, y = minus_log10padj, color = condition_comparison))+
     geom_point()+
     geom_vline(xintercept= c(0.5,-0.5), linetype = 'dashed')+
     geom_vline(xintercept= c(0.5,-0.5), linetype = 'dashed')+
     ggrepel::geom_text_repel(data =dif_res_by_ct_filt %>% dplyr::filter(abs(estimate)>=0.5, padj<=0.1), aes(label = feature), size = 2.5) +
     facet_wrap(~pseudotime_bins)






#
# # ggplot(example_cytof_data_mahalanobis_results, aes(x = PSEUDOTIME_NORMALIZED, y=mahalanobis_distance_all, color = condition))+
# #      theme_bw() +
# #      theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
# #      # geom_point() +
# #      # geom_line(size = 0.5, linetype = 'dashed',aes(y = mahalanobis_distance_lm + mahalanobis_distance_lm_sd)) +
# #      # geom_line(size = 1) +
# #      geom_smooth()+
# #      # geom_line(size = 1)+
# #      ggnewscale::new_scale_colour()+
# #      viridis::scale_color_viridis(option = 'magma', discrete = T) +
# #      geom_rug(data = myRug, aes(x = pseudotime, color  = group), inherit.aes = F)
#
# ggplot(example_cytof_data_mahalanobis_results, aes(x = PSEUDOTIME_NORMALIZED, y=mahalanobis_distance_lm, color = condition))+
#      theme_bw() +
#      theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
#      # geom_point() +
#      # geom_line(size = 0.5, linetype = 'dashed',aes(y = mahalanobis_distance_lm + mahalanobis_distance_lm_sd)) +
#      # geom_line(size = 1) +
#      geom_smooth()+
#      # geom_line(size = 1)+
#      ggnewscale::new_scale_colour()+
#      viridis::scale_color_viridis(option = 'magma', discrete = T) +
#      geom_rug(data = myRug, aes(x = pseudotime, color  = group), inherit.aes = F)


# example_cytof_data_mahalanobis_results_summary = example_cytof_data_mahalanobis_results %>%
#      group_by(pseudotime_bins) %>%
#      mutate(PSEUDOTIME_NORMALIZED  = min(PSEUDOTIME_NORMALIZED)) %>%
#      group_by(PSEUDOTIME_NORMALIZED, condition) %>%
#      summarize(mahalanobis_distance_lm_sd = sd(mahalanobis_distance_lm),
#                mahalanobis_distance_all_sd = sd(mahalanobis_distance_all),
#                mahalanobis_distance_lm=median(mahalanobis_distance_lm),
#                mahalanobis_distance_all = median(mahalanobis_distance_all)) %>%
#      ungroup()
#
# ## all cells
# ggplot(example_cytof_data_mahalanobis_results_summary, aes(x = PSEUDOTIME_NORMALIZED, y=mahalanobis_distance_all, color = condition))+
#      theme_bw() +
#      theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
#      # geom_point() +
#      geom_line(size = 0.5, linetype = 'dashed',aes(y =mahalanobis_distance_all+ mahalanobis_distance_all_sd)) +
#      geom_line(size = 1) +
#      # geom_line(size = 1)+
#      ggnewscale::new_scale_colour()+
#      viridis::scale_color_viridis(option = 'magma', discrete = T) +
#      geom_rug(data = myRug, aes(x = pseudotime, color  = group), inherit.aes = F)
#
# ## landmark approach
# ggplot(example_cytof_data_mahalanobis_results_summary, aes(x = PSEUDOTIME_NORMALIZED, y=mahalanobis_distance_lm, color = condition))+
#      theme_bw() +
#      theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
#      # geom_point() +
#      # geom_line(size = 0.5, linetype = 'dashed',aes(y = mahalanobis_distance_lm + mahalanobis_distance_lm_sd)) +
#      # geom_line(size = 1) +
#      geom_smooths()+
#      # geom_line(size = 1)+
#      ggnewscale::new_scale_colour()+
#      viridis::scale_color_viridis(option = 'magma', discrete = T) +
#      geom_rug(data = myRug, aes(x = pseudotime, color  = group), inherit.aes = F)


################################################
## identify differential features of treatment cells
################################################

glm_input = example_cytof_data_mahalanobis_results %>%
     mutate(noncanonical = ifelse(mahalanobis_distance_lm > quantile(example_cytof_data_mahalanobis_results$mahalanobis_distance_lm, 0.75), 'noncanonical', 'canonical'))

dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input ,outcome_features = features_cellcycle, contrast_variables = 'noncanonical',intercept = TRUE,contrast_method = 'emm',
                                                     SPLIT_BY_NAMES = c('condition'))


library(patchwork)
ggplot(dif_res_by_ct , aes(x = estimate, y = minus_log10padj, color = condition))+
     geom_point()+
     geom_vline(xintercept= c(0.2,-0.2), linetype = 'dashed')+
     geom_vline(xintercept= c(0.2,-0.2), linetype = 'dashed')+
     ggrepel::geom_text_repel(data =dif_res_by_ct %>% dplyr::filter(abs(estimate)>=0.2, padj<=0.1), aes(label = feature), size = 2.5)



##############################
## COSINE ANALYSIS
#############################

####################################################################################################
## split CD4 and CD8 and perform statistical analysis  ####################################################################################################
####################################################################################################
# data.output_filt_cycling_G0G1
data.output_filt_cycling_G0G1_CD4CD8 = data.output_filt_cycling_G0G1 %>%
     dplyr::filter(CD4_CD8 %in% c('CD4+','CD8+'))


cosine_input = data.output_filt_cycling_G0G1_CD4CD8 %>%
     dplyr::filter(timepoint !='0hrs', !is.na(DivisionID), DivisionID != '>8')

# cosine_input = glm_input_test


cosine_input_sampled = example_cytof_data_mahalanobis_results %>%
     ungroup() %>%
     group_by(pseudotime_bins) %>%
     mutate(PSEUDOTIME_NORMALIZED  = min(PSEUDOTIME_NORMALIZED)) %>%
     # computePseudotimeBins(.,pseudotime_column = 'G0G1_cycling_trajectory',interval_sequence = seq(-0.5,1,0.025)) %>%
     # group_by(pseudotime_bins) %>%
     mutate(group =condition)

cosine_input_sampled= cosine_input_sampled%>%
     group_by(PSEUDOTIME_NORMALIZED) %>%
     group_map(~downsampleWith_group_by(., 500),.keep = T) %>%
     bind_rows()


pseudotime_match_df = cosine_input_sampled %>%
     distinct(pseudotime_bins, PSEUDOTIME_NORMALIZED) %>%
     mutate(pseudotime_bins = as.character(pseudotime_bins))

calc_cosine_return_unmelted_matrix = function(ddd, features, group_col){
     if( round(runif(n = 1,0,100))==1){
          cat(1)
     }

     if (nrow(ddd) < 1){
          return(NULL)
     }

     print(unique(ddd[[group_col]]))
     cr = coop::cosine(ddd %>% dplyr::select(all_of(features)) %>% t())
     cr[lower.tri(cr,diag = T)] <- NA
     cr = 1-cr
     cr_df = data.frame(cr, check.names = F, check.rows = F)
     colnames(cr_df) = ddd$cell.id
     rownames(cr_df) = ddd$cell.id
     return(cr_df)
}
calc_cosine_return_melted_matrix = function(ddd, features, group_col){
     if( round(runif(n = 1,0,100))==1){
          cat(1)
     }

     if (nrow(ddd) < 1){
          return(NULL)
     }

     print(unique(ddd[[group_col]]))
     cr = coop::cosine(ddd %>% dplyr::select(all_of(features)) %>% t())
     cr[lower.tri(cr,diag = T)] <- NA
     cr = 1-cr
     cr_df = data.frame(cr, check.names = F, check.rows = F)
     colnames(cr_df) = ddd$cell.id
     rownames(cr_df) = ddd$cell.id

     cr_df = cr_df %>%
          mutate(cell.id = rownames(.)) %>%
          gather(key = cell.id2, value = 'cosine', -cell.id) %>%
          dplyr::filter(cell.id != cell.id2) %>%
          na.omit()
     # if (!is.null(self_bin_analysis)){
     #      cr_df = cr_df %>%
     #           left_join(ddd %>% dplyr::select(all_of(self_bin_analysis), cell.id)) %>%
     #           left_join(ddd %>% dplyr::select(all_of(self_bin_analysis), cell.id2=cell.id) %>% dplyr::rename_at(vars(self_bin_analysis), toupper)) %>%
     #           dplyr::filter((!!sym(self_bin_analysis)) == (!!sym(toupper(self_bin_analysis))) ) %>%
     #           group_by(cell.id, !!sym(self_bin_analysis)) %>%
     #           summarize(cosine = mean(cosine, na.rm =T))
     # }

     return(cr_df)
}

## calculate median cosine value
calc_cosine_dist_mean = function(ddd, features, group_col){
     if( round(runif(n = 1,0,100))==1){
          cat(1)
     }

     if (nrow(ddd) < 100) {
          return(NULL)
     }
     cr = coop::cosine(ddd %>% dplyr::select(all_of(features)) %>% t())
     cr[lower.tri(cr,diag = T)] <- NA
     cr = 1-cr
     cr_mean = cr %>% median(.,na.rm = T)
     names(cr_mean) = unique(ddd[[group_col]])
     return(cr_mean)
}

calc_cosine_dist_sample = function(ddd, features, group_col, self_bin_analysis = NULL){
     if( round(runif(n = 1,0,100))==1){
          cat(1)
     }
     cr = coop::cosine(ddd %>% dplyr::select(all_of(features)) %>% t())
     cr[lower.tri(cr,diag = T)] <- NA
     cr = 1-cr
     cr_df = data.frame(cr, check.names = F, check.rows = F)
     colnames(cr_df) = ddd$cell.id
     rownames(cr_df) = ddd$cell.id

     cr_df = cr_df %>%
          mutate(cell.id = rownames(.)) %>%
          gather(key = cell.id2, value = 'cosine', -cell.id) %>%
          na.omit()
     if (!is.null(self_bin_analysis)){
          cr_df = cr_df %>%
               left_join(ddd %>% dplyr::select(all_of(self_bin_analysis), cell.id)) %>%
               left_join(ddd %>% dplyr::select(all_of(self_bin_analysis), cell.id2=cell.id) %>% dplyr::rename_at(vars(self_bin_analysis), toupper)) %>%
               dplyr::filter((!!sym(self_bin_analysis)) == (!!sym(toupper(self_bin_analysis))) ) %>%
               group_by(cell.id, !!sym(self_bin_analysis)) %>%
               summarize(cosine = mean(cosine, na.rm =T))
     }

     return(cr_df)
}







#######################################################################################
# cosine analysis by mean ##########################################################################################
#######################################################################################
# the appropriate analaysis stratetgy
cosine_input_sampled = example_cytof_data_mahalanobis_results %>%
     ungroup() %>%
     group_by(pseudotime_bins) %>%
     mutate(PSEUDOTIME_NORMALIZED  = min(PSEUDOTIME_NORMALIZED)) %>%
     # computePseudotimeBins(.,pseudotime_column = 'G0G1_cycling_trajectory',interval_sequence = seq(-0.5,1,0.025)) %>%
     # group_by(pseudotime_bins) %>%
     mutate(group = paste(pseudotime_bins, sep = '_-_')) %>%
     ungroup()

# cosine_input_sampled= cosine_input_sampled%>%
#      group_by(PSEUDOTIME_NORMALIZED) %>%
#      # group_map(~downsampleWith_group_by(., 500),.keep = T) %>%
#      bind_rows()

pseudotime_match_df = cosine_input_sampled %>%
     distinct(pseudotime_bins, PSEUDOTIME_NORMALIZED) %>%
     mutate(pseudotime_bins = as.character(pseudotime_bins))

cosine_results_pseudotime = cosine_input_sampled %>%
     ungroup() %>%
     # mutate_at(features_cellcycle, scale) %>%
     group_by(group) %>%
     group_map(~calc_cosine_dist_mean(.,features_cellcycle,group_col = 'group'),.keep = T)

cosine_df_pseudotime = cosine_results_pseudotime %>% unlist() %>%
     data.frame(cosine = . , group = names(.), stringsAsFactors = F) %>%
     separate(col = 'group',sep = '_-_',into = c('pseudotime_bins'), remove = F) %>%
     left_join(pseudotime_match_df,by = 'pseudotime_bins') #%>%
     # mutate(DivisionID = factor(DivisionID, levels = DIVISION_LEVELS),
     #        timepoint = factor(timepoint, levels =TIMEPOINT_LEVELS),
     #        # G0G1_cycling_trajectory = as.numeric(G0G1_cycling_trajectory))
     # )
ggplot(cosine_df_pseudotime , aes(x = PSEUDOTIME_NORMALIZED,  y = cosine, color = condition))+
     geom_line()
     # colScaleDivisionGreens +
     # facet_wrap(~timepoint, scales = 'free_y')



#
#
#
# ## load example data
# glm_input=hafez::glm_input
#
# ## features of interest
# my_features = colnames(glm_input)[-c(1:3)]
# covariates_in_model = c('mouse') ## covariates of interest
# contrast_variables = c('tissue') ## variable you want to compare
#
# ###############################################################################################
# ## (A) example for running analysis on glm_input features that look like celltype abundance fractions
# ###############################################################################################
# ## note: be sure to make covariates as factor variables when it's not a continuous variable. Eg, mouse ids as a covariate for matched analysis
# glm_input = glm_input%>%
#      mutate(mouse = factor(mouse))
# dif_res =hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm',)
# head(dif_res)
#
#
#
# ###############################################################################################
# ## (B) to run interactions, you would run something like this:
# ###############################################################################################
# contrast_variables = c('tissue','ed')
# dif_res =hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm',)
#
# ###############################################################################################
# ## (C) example for running analysis on glm_input features that look like marker expression
# ###############################################################################################
# glm_input =hafez::ct %>% mutate(mouse = factor(mouse)) %>%
#      bind_rows(hafez::ct %>% mutate(mouse = factor(mouse)))
# contrast_variables = c('tissue')
# my_features = c('CD4','CD40','FOXP3') ## make sure these are in the object (no spelling errors)
# # my_features %in% colnames(glm_input)
# dif_res = hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm',)
#
# ###############################################################################################
# # (D) example code to split across multiple subsets of interest
# ###############################################################################################
# dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm',
#                                                       SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))
#
#
# ## you can also use the tukey method instead of emm
# dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'tukey',
#                                                       SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))
#
#
# ###############################################################################################
# # (E) If your single-cell dataset is very large, we can speed up the process by downsampling.
# ###############################################################################################
# glm_input_sampled = glm_input %>%
#      group_by(patient_id, tissue, annotation_lin_sub) %>%
#      group_map(~downsampleWith_group_by(., 50),.keep = T) %>% bind_rows()
# dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input_sampled ,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = 'patient_id',intercept = TRUE,contrast_method = 'emm',
#                                                       SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))
#










