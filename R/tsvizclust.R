

#' @importFrom dtwclust tsclust
#' @importFrom stats cmdscale
#' @importFrom tidyr gather
#' @importFrom dynutils scale_minmax



# #################################################################################
# ## visualziation and clustering  ##
# #################################################################################

#' @title hafez_tsvz_clust
#' @noRd
#'
hafez_tsvz_clust = function(ts_input,k=5, type='partitional', distance = 'dtw_basic', normalize=TRUE, seed = 0,...){
     dtw_clusters = dtwclust::tsclust(series = ts_input,k = k,type = type,distance = distance,seed = seed,normalize=normalize,...)

     ## extract Cluster sizes with average intra-cluster distance. Select k clusters that minimizes the intra-cluster distance
     if (length(k)> 1){
          cl_dists= lapply(dtw_clusters, function(dtwc){
               # print(dtwc@cluster)
               # print(dtwc@cldist)
               cl_dists = rowsum(dtwc@cldist, (dtwc@cluster))  / as.vector(table(dtwc@cluster))
               cl_dists = data.frame(dist = cl_dists, clust= rownames(cl_dists))  %>%
                    mutate(k_group= as.character(dtwc@k))
          }) %>% bind_rows()
          optim_cluster_k = cl_dists %>% group_by(k_group) %>% summarize(median_icd=median(dist)) %>%
               ungroup() %>%
               dplyr::filter(median_icd == min(median_icd)) %>% .$k_group
          message(paste0('optimal cluster: ' , optim_cluster_k) )
          names(dtw_clusters) = as.character(k)
     } else{
          optim_cluster_k=as.character(k)
          cl_dists = rowsum(dtw_clusters@cldist, (dtw_clusters@cluster))/ as.vector(table(dtw_clusters@cluster))
          cl_dists = data.frame(dist = cl_dists, clust= rownames(cl_dists))  %>%
               mutate(k_group= as.character(dtw_clusters@k))
          dtw_clusters=list(optim_cluster_k=dtw_clusters)
     }
     return(list(cl_dists=cl_dists, dtw_clusters=dtw_clusters, optim_cluster_k=optim_cluster_k))
}

#' @title hafez_tsvz
#' @description Performs time-series distance computation, dDR visualization, and dDR clustering. If multiple clusters_k are provided, it calculates the optimal cluster k based on minimizing intra-cluster distance. Takes as input your processed data (eg, density outputs, smoothed expression, etc). If approach is density or expr, it expects daat structure that is outputed from hafez functions to generate the density or expr values.
#' Otherwise, you have to generate the input data yourself, organized with the same structure. For example,column x:pseudotime that is interpolated to matcha ll timepoints, column y: the corresponding value, and columns separating your groups.
#' To test different clusters, you can input a vector of clusters_k. eg, c(2,3).
#' Cluster patterns are also outputted as a separate object
#' @export
hafez_tsvz=function(ts_input, groups=NULL, approach = c('density','expr','other'), k=5, dims=2, type='partitional', distance = 'dtw_basic', normalize=TRUE, seed = 0, return_dDRclust=FALSE, ...){
     options(dplyr.summarise.inform = FALSE)
     pst_coords=unique(ts_input$x)
     if (approach== 'density' | approach == 'expr'){
          message('using density or expr approach for ts_input constructed using hafez functions...')

          ts_input=ts_input%>%
               group_by(across(all_of(c('x',groups)))) %>%
               summarize(y = mean(y)) %>%
               tidyr::unite("group", groups, sep = "_-_", remove = FALSE) %>%
               ungroup()%>%
               dplyr::select(x,y,group)%>%
               spread(key = x, value = y ) %>%
               column_to_rownames('group')
     }
     message('performing time-series analysis...')
     dtw_clusters_dists = hafez_tsvz_clust(ts_input = ts_input,k = k,type = type,distance = distance,seed = seed,normalize=normalize,...)
     dtw_clusters = dtw_clusters_dists$dtw_clusters[[dtw_clusters_dists$optim_cluster_k]]
     dDR = stats::cmdscale(dtw_clusters@distmat, add = T, k = dims, eig = F ) ## horizontal values should be 0
     colnames(dDR$points) = paste0('dDR',1:dims)
     dDR  = dDR$points %>%
          data.frame(check.rows = F, check.names = F)%>%
          mutate(group = rownames(.)) %>%
          separate(col = 'group',sep = '_-_',into = groups) %>%
          mutate(dDR_cluster = as.character(dtw_clusters@cluster))

     names(dtw_clusters@centroids) = as.character(1:dtw_clusters@k)
     centroid_df = bind_cols(dtw_clusters@centroids)
     centroid_df=centroid_df %>%
          mutate(x= pst_coords) %>%
          gather(key = cluster, value = y ,-x) %>%
          mutate(cluster = as.character(cluster)) %>%
          group_by(cluster)%>%
          mutate(y_01 = dynutils::scale_minmax(y))

     p.clusterpatterns = ggplot(centroid_df, aes(x = x, y = y_01, color =cluster))+
          theme_bw() +
          geom_line()

     if (return_dDRclust==TRUE){
          return(list(dDR=dDR, dDR_centroids = centroid_df, p.clusterpatterns=p.clusterpatterns, dDR_clust=dtw_clusters_dists))
     } else {
         return( list(dDR=dDR, dDR_centroids = centroid_df, p.clusterpatterns=p.clusterpatterns))
     }
}




#' @title hafez_smoothing
#' @description Computes GAM smoothing using a spline.
#' @export
## using smoothed expression patterns
hafez_smoothing = function(data,groups = c(NULL), features, pseudotime_column = 'pst', pst_interval=NULL, bs='cs',...){
     if (is.null(pst_interval)){
          pst_min = min(data[[pseudotime_column]],na.rm=T)
          pst_max = max(data[[pseudotime_column]], na.rm = T)
          pst_interval = seq( pst_min, pst_max, (pst_max-pst_min)/100)
     }

     smoothed_patterns = data %>%
          ungroup()%>%
          mutate(pst = .[[pseudotime_column]]) %>%
          pivot_longer(cols = all_of(features), names_to = 'feature',values_to = 'value') %>%
          group_by(across(all_of(c(groups, 'feature')))) %>%
          group_modify( .f= function(dd,...){
               tryCatch({
                    gr = mgcv::gam(formula = value ~ s(pst, bs = bs), data = dd)
                    gr_pred = predict(gr, data.frame(pst=pst_interval), type = 'response', se.fit = TRUE)
                    result = data.frame(pst=pst_interval) %>%
                         mutate(x=pst,
                                y = as.vector( gr_pred$fit),
                                dev = as.vector(gr_pred$se.fit))
                    return(result)
               }, error = function(e) {
                    # message(unique(dd$feature))
                    # You can log the error message if needed
                    message("An error occurred - unable to fit smoothed curve for a specific group: ", conditionMessage(e))
                    return(data.frame(pst = NA))
               })
          })
     message('complete!')
     return(smoothed_patterns)
}









