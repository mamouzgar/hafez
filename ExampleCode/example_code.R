
library(tidyverse)
library(hafez)
library(patchwork)

features_cellcycle  = c( 'dsDNA', 'pH3_s10','CDT1', 'Geminin', 'PLK1', 'DNA2', 'CyclinB1', 'PCNA', 'Ki67','SLBP','IdU','pRb_S780')%>% unique()

Phases_gate =  c('G1','S','G2','M')
Phases_gate_Colors =c("#000004FF", "#721F81FF", "#F1605DFF","#FCFDBFFF")
# allPhasesColors[1] = "#000000"
names(Phases_gate_Colors) = Phases_gate

colScale_phases_manualAll = scale_color_manual(values = Phases_gate_Colors,name = 'Cell cycle\nphase')
fillScale_phases_manualAll = scale_fill_manual(values = Phases_gate_Colors, name = 'Cell cycle\nphase')


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
hafez_output = hafez_DBPN(dataset  = hafez_output,  dataset.subset_to_use = hafez_output %>%dplyr::filter(condition == 'WT'),column_to_normalize = 'LM_TI_path1',adjust.value = 0.5)

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
exp_patterns = hafez_smoothing(hafez_output,groups = c('condition'),features = features_cellcycle,pseudotime_column = 'PSEUDOTIME_NORMALIZED')
ggplot(exp_patterns, aes(x=pst ,y = y, color =feature))+
     geom_line()+
     facet_wrap(~condition)


## detect inflection points from smoothed distributions (requires long-format data)
infl_points=hafez_detect_inflections(dataset = exp_patterns,groups = 'condition')

## errors intentionally to show that it won't work if there are multiple x-y pairs if all the necessary groups to stratify by are not specified.
hafez_detect_inflections(dataset = exp_patterns) ## ERRORS

## plot inflections
infl_points_wt = infl_points %>% dplyr::filter(condition=='WT')
ggplot(infl_points_wt %>% dplyr::filter(feature =='CyclinB1'), aes(x=x, y=y))+
     theme_bw()+
     # geom_point() +
     geom_point(aes(color= factor(slope_change_directionality)))+
     geom_vline(data=infl_points_wt %>% dplyr::filter(feature =='CyclinB1',slope_change==T), aes(xintercept =x))


## TSD embedding and clustering
ts_res_exp = hafez_tsvz(exp_patterns,groups = c('condition','timepoint','feature'), approach = 'density',k = c(2:12),return_dDRclust = F)

ggplot(ts_res_exp$dDR, aes(x=dDR1, y=dDR2, color = dDR_cluster))+
     theme_bw() +
     geom_point()


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









