
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
#      dplyr::select(all_of(c(features_cellcycle)), cell.id, gate=traj_phase, condition, timepoint, pseudotime_orig =pseudotime_adjusted ) %>%
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


# example_cytof_data = hafez::example_cytof_data %>%
#      mutate_at(features_cellcycle, scale)

##################################################################################
## perform hafez landmark processing and trajectory inference
##################################################################################
# start=Sys.time()
hafez_input =  hafez::example_cytof_data
hafez_input = hafez_landmark_processing(train = hafez_input %>% dplyr::filter(condition == 'WT'),full_data = hafez_input,features = features_cellcycle,method = 'PCA',return_object = F)
features_pc=paste0('PC',1:6)

## without landmarks
hafez_output = hafez_TI(FULL_DATA = hafez_input,LM_DATA = hafez_input,
                        features = features_pc,NumNodes =5,lambda = 0.01,mu = 0.01,return_pseudotime_only = F,
                        features_for_start_cell_id =  features_cellcycle,branch_type  = 'circle',verbose = F, ## verbose does not currenly work
                        return_node_pos = F )

## with landmarks
hafez_output = hafez_TI(FULL_DATA = hafez_input, LM_DATA = hafez_output %>% sample_n(500) ,
                        features = features_pc,NumNodes =5,lambda = 0.01,mu = 0.01,return_pseudotime_only = F,
                        features_for_start_cell_id =  features_cellcycle,branch_type  = 'curve',verbose = F, ## verbose does not currenly work
                        return_node_pos = F )
hafez_output


############################################################
## DBPN: density based pseudotime normalization
############################################################
## perform pseudotime normalization on the newly computed pseudotime with the toy dataset
hafez_output = hafez_DBPN(dataset = hafez_output, dataset.subset_to_use = hafez_output %>%dplyr::filter(condition == 'WT'),column_to_normalize = 'LM_TI_path1',adjust.value = 0.5)

ggplot(hafez_output, aes(x = PSEUDOTIME_NORMALIZED, y = LM_TI_path1))+
     theme_bw() +
     theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
     geom_point(aes(color = gate)) +
     colScale_phases_manualAll

## The toy dataset contructs a poor pseudotime so we provide the original pseudotime as an example
hafez_output = hafez_DBPN(dataset = hafez_output, dataset.subset_to_use = hafez_output %>%dplyr::filter(condition == 'WT'),column_to_normalize = 'pseudotime_orig',adjust.value = 0.5)

ggplot(hafez_output, aes(x = PSEUDOTIME_NORMALIZED, y = pseudotime_orig))+
     theme_bw() +
     theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
     geom_point(aes(color = gate)) +
     colScale_phases_manualAll

## you can also density normalize by different groups, like different donors.
## it will compute a density curve for each donor, then average them before performing pseudotime normalization.
## this is useful if you want to stabilize donor heterogeneity, and controls for variation in cell numbers so donors with more cells do not bias the normalization
## here we just normalize by the condition column as the group.
hafez_output = hafez_DBPN(dataset = hafez_output, dataset.subset_to_use = hafez_output ,
                          column_to_normalize = 'pseudotime_orig', adjust.value = 0.5,
                          normalize_by_sample_column = 'condition')

## We demonstrate that pseudotime normalization quantifies inhibitor action as a biological concept, but we recommend looking at cell density distributions of different inhibitors/groups along a normal trajectory for interepretation purposes.
## That said, you can reproduce this inhibitor action demonstrated in the manuscript using then following code.
## herre is an example code to efficiently do this using tidyverse
hafez_output_pst_action = hafez_output %>%
     group_by(condition2=condition) %>%
     group_map(~{
          norm_pst_col_name = unique(.$condition)
          # print(norm_pst_col_name)
          hafez_DBPN(dataset=hafez_output, dataset.subset_to_use = .,column_to_normalize = 'pseudotime_orig',
                     new_dbp_name = 'pst_inhib_action'  ) %>%
               mutate(condition_used = norm_pst_col_name)
          }) %>%
     bind_rows()

## then we compute the average phenotype in each bin
hafez_output_pst_action_res = hafez_output_pst_action %>%
     group_by(condition2=condition,condition_used2=condition_used)%>%
     group_map(~{
          condition = unique(.$condition)
          condition_used = unique(.$condition_used)
          compute_rug(.,pseudotime_column = 'pst_inhib_action',group_column = 'gate',mySeq = seq(0,1,0.01),groups_to_keep = c('G1','S','G2','M' )) %>%
                           mutate(condition=condition,
                                  condition_used=condition_used)
     }) %>% bind_rows()
ggplot(hafez_output_pst_action_res %>% dplyr::filter(condition==condition_used), aes(x = pseudotime_bins, y =condition  , fill =group ))+
     geom_tile(color = 'black')+
     fillScale_phases_manualAll


############################################################
## calculate cell density by group
############################################################
density_by_group_results = calculate_cell_density_by_group(dataset = hafez_output, group_splits = c('condition','timepoint'),pseudotime_column_name = 'PSEUDOTIME_NORMALIZED')
myRug = compute_rug(dataset = hafez_output,pseudotime_column = 'PSEUDOTIME_NORMALIZED',group_column = 'gate',num_bins = 100,recompute_pseudotime_bins = T) %>%
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
## The previous landmark preprocessing used a PCA approach for LM normalization. The outputted preprocessed data only returns the landmark constructs PCs even though the features were landmark processed under the hood.
## To perform further landmark analysis directly on the features, we reprocess the features directly using then method arguement.
hafez_output_2 = hafez_landmark_processing(train = hafez_output %>% dplyr::filter(condition == 'WT'),
                                           full_data = hafez_output,features = features_cellcycle,
                                         method = 'raw',return_object = F)

## perform desired discetization of pseudotime bins for post-traj landmark analysis.
hafez_output_2 = computePseudotimeBins(dataset = hafez_output_2,pseudotime_column = 'PSEUDOTIME_NORMALIZED',interval_sequence = seq(0,1,0.05))

pm_mahalanobis_results_all_cells = pseudotime_mapped_mahalanobis_analysis(dataset = hafez_output_2,method = 'all',features = features_cellcycle, reference_group_column = NULL,reference_group_name = NULL,CELL_COUNT_THRESHOLD = 0)
pm_mahalanobis_results_wt_landmarks = pseudotime_mapped_mahalanobis_analysis(dataset = hafez_output_2,method = 'landmark', features = features_cellcycle, reference_group_column = 'condition',reference_group_name = 'WT', CELL_COUNT_THRESHOLD = 0)


hafez_output_2_mahalanobis_results = hafez_output_2 %>%
     dplyr::left_join(pm_mahalanobis_results_all_cells %>% dplyr::rename(mahalanobis_distance_all = mahalanobis_distance)) %>%
     dplyr::left_join(pm_mahalanobis_results_wt_landmarks %>% dplyr::rename(mahalanobis_distance_lm = mahalanobis_distance))

ggplot(hafez_output_2_mahalanobis_results, aes(x = PSEUDOTIME_NORMALIZED, y=mahalanobis_distance_all, color = condition))+
     theme_bw() +
     theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
     geom_smooth()+
     ggnewscale::new_scale_colour()+
     viridis::scale_color_viridis(option = 'magma', discrete = T) +
     geom_rug(data = myRug, aes(x = pseudotime, color  = group), inherit.aes = F)

ggplot(hafez_output_2_mahalanobis_results, aes(x = PSEUDOTIME_NORMALIZED, y=mahalanobis_distance_lm, color = condition))+
     theme_bw() +
     theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
     geom_smooth()+
     ggnewscale::new_scale_colour()+
     viridis::scale_color_viridis(option = 'magma', discrete = T) +
     geom_rug(data = myRug, aes(x = pseudotime, color  = group), inherit.aes = F)


######################################################################################################
## performing landmark, or landmark-pseudotime knn  distance analysis
######################################################################################################

## map to landmark cells
score_res_distance = find_closest_neighbor_distance(full_data = hafez_output %>% ungroup()  %>% mutate_at(features_cellcycle, scale),
                                                    LM_col = 'condition',LM_group = 'WT',features = features_cellcycle,
                                                    k_ave = 10,
                                                    k_solo = 1, return_knn_graph = F)

hafez_output$nn_dist = score_res_distance

ggplot(hafez_output, aes(x = PSEUDOTIME_NORMALIZED, y=nn_dist, color = condition))+
     theme_bw() +
     theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
     geom_smooth()+
     ggnewscale::new_scale_colour()+
     viridis::scale_color_viridis(option = 'magma', discrete = T) +
     geom_rug(data = myRug, aes(x = pseudotime, color  = group), inherit.aes = F)

## map to landmark cells along pseudotime - this utility may be important depending on how your cell ordering was defined.
score_res_distance_pst = hafez_output %>%
     ungroup()  %>%
     mutate_at(features_cellcycle, scale)%>%
     hafez::computePseudotimeBins(pseudotime_column = 'PSEUDOTIME_NORMALIZED',interval_sequence = seq(0,1,0.05))%>%
     group_by(pseudotime_bins) %>%
     group_map(~find_closest_neighbor_distance(full_data = .,
                                               LM_col = 'condition',LM_group = 'WT',features = features_cellcycle,
                                               k_ave = 5,
                                               k_solo = 1, return_knn_graph = F) ) %>%
     unlist()

## make sure to reorder according to cell id
hafez_output$nn_dist_pst = score_res_distance_pst[hafez_output$cell.id]

ggplot(hafez_output, aes(x = PSEUDOTIME_NORMALIZED, y=nn_dist_pst, color = condition))+
     theme_bw() +
     theme(panel.border = element_rect(fill = 'transparent', size = 0.5)) +
     geom_smooth()+
     ggnewscale::new_scale_colour()+
     viridis::scale_color_viridis(option = 'magma', discrete = T) +
     geom_rug(data = myRug, aes(x = pseudotime, color  = group), inherit.aes = F)

################################################################################
## time series visualization and clustering
################################################################################

## using density results
density_by_group_results2=density_by_group_results %>%
     bind_rows(density_by_group_results %>% mutate(timepoint = '1')) ## example pseudo data
ts_res = hafez_tsvz(density_by_group_results2,groups = c('condition','timepoint'),approach = 'density',k = c(2,3),return_dDRclust = F)

ts_res$p.clusterpatterns
ggplot(ts_res$dDR, aes(x= dDR1, y =dDR2, color = dDR_cluster))+
     theme_bw() +
     geom_point()



## compute smoothed expression patterns and perform TSD clustering
exp_patterns = hafez_smoothing(hafez_output,groups = c('condition','timepoint'),features = features_cellcycle,pseudotime_column = 'PSEUDOTIME_NORMALIZED')
ggplot(exp_patterns, aes(x=pst ,y = y, color =feature))+
     geom_line()+
     facet_wrap(~condition)
ts_res_exp = hafez_tsvz(exp_patterns,groups = c('condition','timepoint','feature'), approach = 'density',k = c(2:12),return_dDRclust = F)

ggplot(ts_res_exp$dDR, aes(x=dDR1, y=dDR2, color = dDR_cluster))+
     theme_bw() +
     geom_point()


# ts_res$dDR_clust@cluster
# ts_res$dDR_clust@cldist
# rownames(ts_res@cldist) = ts_res@cluster

# match(rownames(ts_res$dDR_clust@cldist),ts_res$dDR_clust@cluster)

# table(tx2gene_filt$gene_id == rownames(so_d5)) ##  all  should true

# rownames(so) = tx2gene_filt$gene_name ## assign GENE NAME to the transcript name. NOTE: only works for matrices, not dataframes
hafez_tsvz_clust()

##################################################################################
## identify differential features of noncaonical cells ##
##################################################################################
## identify differential features between noncanonical and canonical cells in drug systems.
glm_input = hafez_output_2_mahalanobis_results %>%
     computePseudotimeBins(dataset = .,pseudotime_column = 'PSEUDOTIME_NORMALIZED',interval_sequence = seq(0,1,0.5)) %>%
     mutate(noncanonical = ifelse(mahalanobis_distance_lm > quantile(example_cytof_data_mahalanobis_results$mahalanobis_distance_lm, 0.75,na.rm =T), 'noncanonical', 'canonical'))
dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input ,outcome_features = features_cellcycle, contrast_variables = 'noncanonical',intercept = TRUE,contrast_method = 'emm',
                                                     SPLIT_BY_NAMES = c('condition', 'pseudotime_bins'))
dif_res_by_ct = dif_res_by_ct #%>%
     # dplyr::filter(!is.na(comparison))

library(patchwork)
ggplot(dif_res_by_ct , aes(x = estimate, y = minus_log10padj, color = condition))+
     theme_bw() +
     geom_point()+
     geom_vline(xintercept= c(0.5,-0.5), linetype = 'dashed')+
     geom_vline(xintercept= c(0.5,-0.5), linetype = 'dashed')+
     ggrepel::geom_text_repel(data =dif_res_by_ct %>% dplyr::filter(abs(estimate)>=0.5, padj<=0.1), aes(label = feature), size = 2.5) +
     facet_wrap(~pseudotime_bins)

## omit this in future iteration
# ## differential analysis can also  be used to compare conditions directly
# dif_res_by_ct = hafez::differential_analysis_program(glm_input = glm_input ,outcome_features = features_cellcycle, contrast_variables = 'condition',intercept = TRUE,contrast_method = 'emm',
#                                                      SPLIT_BY_NAMES = c( 'pseudotime_bins'))
# dif_res_by_ct_filt = dif_res_by_ct %>%
#      dplyr::filter(condition1 =='WT'|condition2 == 'WT')
# ggplot(dif_res_by_ct_filt , aes(x = estimate, y = minus_log10padj, color = condition_comparison))+
#      geom_point()+
#      geom_vline(xintercept= c(0.5,-0.5), linetype = 'dashed')+
#      geom_vline(xintercept= c(0.5,-0.5), linetype = 'dashed')+
#      ggrepel::geom_text_repel(data =dif_res_by_ct_filt %>% dplyr::filter(abs(estimate)>=0.5, padj<=0.1), aes(label = feature), size = 2.5) +
#      facet_wrap(~pseudotime_bins)


# hafez_output_ti = hafez_TI(FULL_DATA = hafez_output,
#                               features = features_pc[1:6],NumNodes =5,lambda = 0.01,mu = 0.01,
#                               features_for_start_cell_id =  features_cellcycle,branch_type  = 'curve',
#                               PERFORM_OOS = FALSE, OOS_idx = landmark_dr2 %>% ungroup(), return_node_pos = F )
# end=Sys.time()

# TEST = ElPiGraph.R::computeElasticPrincipalCurve(X = as.matrix(hafez_output %>% dplyr::select(all_of(features_pc))),
#                                                     NumNodes = NumNodes,ProbPoint=ProbPoint, nReps =nReps,
#                                                     MaxNumberOfIterations = MaxNumberOfIterations,
#                                                     verbose=verbose,
#                                                     Lambda = Lambda, Mu = Mu,Do_PCA = F,  drawPCAView = drawPCAView, drawAccuracyComplexity=drawAccuracyComplexity,
#                                                     # verbose = verbose,
#                                                     drawEnergy = drawEnergy)
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










