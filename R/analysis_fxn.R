
#' @importFrom stats density
#' @importFrom sqldf sqldf



#' @title computePseudotimeBins
#' @description
#' @export
computePseudotimeBins = function(dataset,pseudotime_column, bin_width = 0.01, bin_num = 100, bin_approach = 'find_interval', interval_sequence = seq(0,1,0.01)){
     if (bin_approach == 'find_interval'){
          if (is.null(interval_sequence) == TRUE){
               stop('you must input an arguement for interval_sequence,such as seq(0,1,0.05). The order is seq(from = start pst, to=  end pst, by = increment of sequence ')
          }
          # message('using interval_sequence...' )
          # print(interval_sequence)
          dataset = dataset %>%
               mutate(pseudotime_bins =findInterval(!!sym(pseudotime_column), interval_sequence), ## ~ 180 bins
                      pseudotime_bins = as.numeric(pseudotime_bins))
          # } else if(bin_approach == 'equal_width'){
          #      dataset = dataset %>%
          #           mutate(pseudotime_bins = cut_width(!!sym(pseudotime_column),bin_width,closed = 'left'), ## ~ 180 bins
          #                  pseudotime_bins = as.numeric(pseudotime_bins))
     } else{
          dataset = dataset %>%
               mutate(pseudotime_bins = ggplot2::cut_number(!!sym(pseudotime_column), bin_num), ## ~ 180 bins
                      pseudotime_bins = as.numeric(pseudotime_bins))
     }

     return(dataset)
}


#' @title compute_rug
#' @description outputs a dataset to ggplot a discretized rug along the pseudotime based on categorical labels like cell cycle gates. Each bin is assigned a label based on a label's normalized frequency.
#' @export
compute_rug = function(dataset, pseudotime_column, group_column = 'gate', mySeq=NULL, num_bins=100,groups_to_keep = NULL,recompute_pseudotime_bins=TRUE){


     if (is.null(mySeq)){
          mySeq = seq(0,max(dataset[[pseudotime_column]]),1/num_bins)
     }
     dataset[['group']] = dataset[[group_column]]

     if(is.null(groups_to_keep)){
          groups_to_keep= unique(dataset[['group']])
     }

     RUG_DF = dataset %>%
          dplyr::filter(group %in% groups_to_keep)

     if (recompute_pseudotime_bins == TRUE){
          RUG_DF =  RUG_DF %>%
               computePseudotimeBins(.,pseudotime_column = pseudotime_column,interval_sequence = mySeq)
     }
     RUG_DF=RUG_DF%>%
          group_by(pseudotime_bins, group) %>%
          summarize(count = n(),
                    pseudotime = min(!!sym(pseudotime_column),na.rm = T)) %>%
          group_by(group) %>%
          mutate(percentage = count/sum(count)) %>%
          group_by(pseudotime_bins) %>%
          dplyr::filter(percentage == max(percentage)) %>%
          mutate(group = factor(group, levels = groups_to_keep)) %>%
          ungroup()
     return(RUG_DF)
}




#' @title calculate_cell_density_by_group
#' @description
#' @export
calculate_cell_density_by_group = function(dataset, group_splits = NULL, pseudotime_column_name,adjust.value = 1, n = 512,from=NULL, to =NULL, COMPUTE_CIRCULAR = FALSE ,count_threshold = 100, sep_operator = '_-_', verbose = FALSE, ...){

     if (is.null(from)){
          from = min(dataset[[pseudotime_column_name]])
     }
     if (is.null(to)){
          to = max(dataset[[pseudotime_column_name]])
     }

      if (is.null(group_splits)){
          message('no grouping being performed...')
          myGroups = rep("dataset",nrow(dataset))
          group_splits="dataset"
          dataset[[group_splits]] = "dataset"
     } else {
          myGroups = apply(dataset[group_splits], MARGIN = 1,paste,collapse = sep_operator)

     }
     # if (verbose == TRUE){
     #      message(unique(myGroups))
     # }

     density.results = lapply(split(dataset, myGroups) , function(d.input){
          myGroup = apply(d.input[group_splits], MARGIN = 1,paste,collapse = sep_operator) %>% unique()

          if (nrow(d.input)<=count_threshold){
               return(NULL)
          }
          if (COMPUTE_CIRCULAR==TRUE){
               # message('computing circular...')
               d.input1 = d.input
               d.input2 = d.input
               d.input1[[pseudotime_column_name]] =  1+d.input1[[pseudotime_column_name]]
               d.input2[[pseudotime_column_name]]=  d.input1[[pseudotime_column_name]]-1
               d.input = d.input %>% bind_rows(d.input1) %>% bind_rows(d.input2)
               from = from-to
               to =to+to
          }
          dnsty.res = tryCatch(
               {
                    dnsty.res = density(d.input[[pseudotime_column_name]] %>% na.omit(), adjust = adjust.value, n = n,..., from = from, to =to)
                    # return(dnsty.res)

               },
               error = function(cond) {
                    # warning(cond)
                    message('failed to compute...skipping')
                    return(NULL)
               }
          )
          if(is.null(dnsty.res)){
               myDensityDf=NULL
          } else {
               myDensityDf = data.frame(x=dnsty.res$x, y = dnsty.res$y, group =myGroup, group2 =myGroup) %>%
                    group_by(group2) %>%
                    mutate(y_normalized = dynutils::scale_minmax(y)) %>%
                    ungroup() %>%
                    separate(group2, sep = sep_operator, into = group_splits ) %>%
                    mutate(y_log10 = log10(y+1))
               if (COMPUTE_CIRCULAR==TRUE){
                    myDensityDf = myDensityDf %>% dplyr::filter(x>=from & x <=to) %>%
                         mutate(circularized=TRUE)
               }
          }
          return(myDensityDf)
     }) %>% bind_rows()
     message('calculated cell density by group')
     return(density.results)
}
## example code
## density.res = calculate_cell_density_by_group(data = dataset.dbpn%>% dplyr::filter(auto_gate %in% c('G1','S','G2','M','G1/G2','M',NA)), group_splits = 'condition', pseudotime_column_name='PSEUDOTIME_NORMALIZED', adjust.value = 1, n=1000)

############################################################################################################
## mahalanobis distance  by group ##################################################################################################################################################################
##################################################################################################################################################################
# mahalanobis_by_bins = function(data, features, tol=1e-20 ){
#      if (nrow(data)<20){
#           return(NULL)
#      }
#
#      mahalanobis.input = data %>% ungroup() %>% dplyr::select(all_of(features))# %>% as.matrix()
#      center = colMeans(mahalanobis.input)
#      cov_matrix <- cov(mahalanobis.input)
#      mahalanobis_distance <-  mahalanobis(mahalanobis.input, center, cov_matrix)
#      return(data.frame(cell.id = data$cell.id, mahalanobis_distance=mahalanobis_distance))
#
#      mahalanobis_distance =  try(  mahalanobis_distance <-  mahalanobis(mahalanobis.input, center, cov_matrix), silent = TRUE)
#      if (class(mahalanobis_distance) == 'try-error'){
#           # message('fail')
#           return(data.frame(cell.id = data$cell.id, mahalanobis_distance=NA))
#      } else {
#           return(data.frame(cell.id = data$cell.id, mahalanobis_distance=mahalanobis_distance))
#      }
# }


#' @title mahalanobis_distance_to_reference_group
#' @description
#' @export
mahalanobis_distance_to_reference_group = function(dataset, features, reference_group_column, reference_group_name, referenceSelf = FALSE, sampleInternally = NULL, tol=1e-20,precalculated_mean = NULL,precalculated_cov = NULL,CELL_COUNT_THRESHOLD=NULL){ ## tol=1e-20

     if (referenceSelf == TRUE){
          REFERENCE_mahalanobis.input = dataset %>%
               ungroup() %>%
               dplyr::select(all_of(features))
     } else {
          REFERENCE_mahalanobis.input = dataset %>%
               dplyr::filter(!!sym(reference_group_column) %in% reference_group_name) %>%
               ungroup() %>%
               dplyr::select(all_of(features))



     }
     if (is.numeric(CELL_COUNT_THRESHOLD) ) {
          if (nrow(REFERENCE_mahalanobis.input)<CELL_COUNT_THRESHOLD){
               return(data.frame(cell.id = dataset$cell.id, mahalanobis_distance=NA))
          }
     }
     cov_matrix <- cov(REFERENCE_mahalanobis.input)
     center <- colMeans(REFERENCE_mahalanobis.input)

     mahalanobis.input = dataset %>%
          ungroup() %>%
          dplyr::select(all_of(features))
     if (is.numeric(CELL_COUNT_THRESHOLD) ) {
          if (nrow(mahalanobis.input)<CELL_COUNT_THRESHOLD){
               return(data.frame(cell.id = dataset$cell.id, mahalanobis_distance=NA))
          }
     }
     mahalanobis_distance =  try(mahalanobis_distance <-  mahalanobis(mahalanobis.input, center, cov_matrix), silent = TRUE)
     if (class(mahalanobis_distance) == 'try-error'){
          return(data.frame(cell.id = dataset$cell.id, mahalanobis_distance=NA))
     }
     # mahalanobis_distance <- mahalanobis(mahalanobis.input, center, cov_matrix, tol = tol)
     return(data.frame(cell.id = dataset$cell.id, mahalanobis_distance=mahalanobis_distance))
}

##takes in "landmark" for landmark mapped analysis,  or "all"

#' @title pseudotime_mapped_mahalanobis_analysis
#' @description
#' @export
pseudotime_mapped_mahalanobis_analysis = function(dataset, method = 'landmark_mapped', pseudotime_bin_column='pseudotime_bins', features, reference_group_column=NULL, reference_group_name=NULL,CELL_COUNT_THRESHOLD=NULL){
     dataset$'pseudotime_bins' = dataset[[pseudotime_bin_column]]
     if (method == 'landmark'){
          PM_mahaalanobis_results = dataset %>%
               group_by(pseudotime_bins) %>%
               group_map(~mahalanobis_distance_to_reference_group(., features, reference_group_column, reference_group_name, referenceSelf = FALSE,  tol=tol, CELL_COUNT_THRESHOLD=CELL_COUNT_THRESHOLD)) %>%
               bind_rows()
     } else if (method == 'all'){
          PM_mahaalanobis_results = dataset %>%
               group_by(pseudotime_bins) %>%
               group_map(~mahalanobis_distance_to_reference_group(., features, reference_group_column, reference_group_name, referenceSelf = TRUE,  tol=tol, CELL_COUNT_THRESHOLD=CELL_COUNT_THRESHOLD)) %>%
               bind_rows()
     }
     return(PM_mahaalanobis_results)
}



##' to add:
##' cosine analysis
##'  geometric sketching
##'  density based normalization
# ############################################################
# ## PSEUDOTIME DENSITY NORMALIZATION ##
# ############################################################
## column_to_normalize: name of column to normalize in dataset
#' @title DBPN
#' @description
#' @export
DBPN = function(dataset, density_bins = 1024,  dataset.subset_to_use = NULL ,column_to_normalize=NULL, bandwidth = 'nrd0', adjust.value = 1, CIRCULARIZE = FALSE, RETURN_DENSITY_COORDINATES = FALSE){
     dataset$x = NULL
     dataset$y = NULL
     dataset$y_normalized = NULL
     dataset$pseudotime_bin_01_min = NULL
     dataset$pseudotime_bin_01_max = NULL
     dataset$cumulative.sum = NULL
     dataset$pseudotime_density_min = NULL
     dataset$pseudotime_density_max = NULL
     dataset$pseudotime_bins_density_based = NULL
     dataset$count = NULL
     dataset$percentage_sum = NULL

     if (!is.null(dataset.subset_to_use)){
          use_subset_of_cells = TRUE
          dataset.subset_to_use$PSEUDOTIME_TO_DBPN = dataset.subset_to_use[[column_to_normalize]]
     } else {
          use_subset_of_cells = FALSE
     }
     dataset$PSEUDOTIME_TO_DBPN = dataset[[column_to_normalize]]
     myRange = c(floor(min(dataset$PSEUDOTIME_TO_DBPN , na.rm =T)), ceiling(max(dataset$PSEUDOTIME_TO_DBPN ,na.rm=T)) ) %>% as.numeric()## from and to values for density calculation
     if(is.null(column_to_normalize)){
          stop('must specify the pseudotime column to normalize')
     }

     ## sql hates the a column labeled time so we remove it if it's present
     time_col_present = 'time' %in% colnames(dataset)
     if (time_col_present == TRUE){
          time_df = dataset %>% ungroup() %>% dplyr::select(time, cell.id)
          dataset = dataset %>% dplyr::select(-time)

          if (use_subset_of_cells == TRUE){
               dataset.subset_to_use = dataset.subset_to_use %>% dplyr::select(-time)

          }
     }
     if (use_subset_of_cells == FALSE){
          message('normalizing on all cells')

          cyclingCells = dataset
          density.res = density(cyclingCells$PSEUDOTIME_TO_DBPN, from = myRange[1], to = myRange[2], n = density_bins,adjust = adjust.value)
          dataset.frame(x = density.res$x, y  =density.res$y)
          ## repeat y values on either side to make it connected/circular
          if (CIRCULARIZE == T){
               density.res = density(c(cyclingCells$PSEUDOTIME_TO_DBPN-1, cyclingCells$PSEUDOTIME_TO_DBPN, cyclingCells$PSEUDOTIME_TO_DBPN+1), from = (myRange[1]-myRange[2]), to = ( myRange[2]+ myRange[2]) ,n=density_bins,bw =bandwidth,adjust = adjust.value)
          }
          density.df = data.frame(x = density.res$x, y  =density.res$y)
          # plot(density.df$x, density.df$y)
          density_df_01 = density.df %>% dplyr::filter(x>=0 & x<=1) %>%
               mutate(y_normalized = y/sum(y),
                      pseudotime_bin_01_min =  lag(x) %>% ifelse(is.na(.), -0.01, .),
                      pseudotime_bin_01_max =  c(x[-length(x)], NA)%>% ifelse(is.na(.), 1.01, .),

                      cumulative.sum = cumsum(y_normalized),
                      pseudotime_density_min = lag(cumulative.sum) %>% ifelse(is.na(.), 0, .),
                      pseudotime_density_max = cumulative.sum,
                      pseudotime_bins_density_based = 1:nrow(.))

          ## SQL based join on a range from S/O
          ## https://stackoverflow.com/questions/46795636/r-dplyr-join-by-range-or-virtual-column
          ## 'from' and 'to' are in quotes because it's reserved language for SQL. Here we use pseudotime_bin_01.min and pseudotime_bin_01.max
          # library('sqldf')
          message('combining...')
          cyclingCells_densityBins = cyclingCells # %>% dplyr::select(-time)
          # SQLbins = as_tibble(sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.pseudotime_bins_density_based from cyclingCells_densityBins
          SQLbins = as_tibble(sqldf::sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.* from cyclingCells_densityBins
                join density_df_01 on cyclingCells_densityBins.PSEUDOTIME_TO_DBPN > density_df_01.pseudotime_bin_01_min and
                          cyclingCells_densityBins.PSEUDOTIME_TO_DBPN <= density_df_01.pseudotime_bin_01_max"))
          cyclingCells = cyclingCells %>%
               dplyr::bind_cols(SQLbins %>% dplyr::select(-PSEUDOTIME_TO_DBPN))

          cyclingCells.PST_ADJUSTMENT = cyclingCells %>%
               group_by(pseudotime_bins_density_based) %>%
               summarize(count = n(),
                         min = min(PSEUDOTIME_TO_DBPN),
                         max = max(PSEUDOTIME_TO_DBPN)) %>%
               mutate(percentage = count/sum(count),
                      percentage_sum  = cumsum(percentage))
          cyclingCells= cyclingCells %>%
               left_join(cyclingCells.PST_ADJUSTMENT) %>%
               group_by(pseudotime_bins_density_based) %>%
               mutate(PSEUDOTIME_NORMALIZED = scales::rescale(PSEUDOTIME_TO_DBPN, to = c(unique(pseudotime_density_min), unique(pseudotime_density_max) )))# %>%

          if (time_col_present == TRUE){
               cyclingCells = cyclingCells %>% left_join(time_df, by = 'cell.id')
          }
          if (RETURN_DENSITY_COORDINATES ==TRUE){
               return(list(df.dbpn = cyclingCells %>% ungroup(), density_df_01 = density_df_01))
          }
          return(cyclingCells %>% ungroup())
     } else if (use_subset_of_cells == TRUE){
          message('computing pseudotime normalization on the provided subset of cells...')
          cyclingCells_FULL = dataset
          cyclingCells = dataset.subset_to_use
          density.res = density(cyclingCells$PSEUDOTIME_TO_DBPN, from =  myRange[1], to = myRange[2], n = density_bins,adjust = adjust.value)
          # data.frame(x = density.res$x, y  =density.res$y)
          ## repeat y values on either side to make it connected/circular
          if (CIRCULARIZE == T){
               density.res = density(c(cyclingCells$PSEUDOTIME_TO_DBPN-1, cyclingCells$PSEUDOTIME_TO_DBPN, cyclingCells$PSEUDOTIME_TO_DBPN+1),  from = (myRange[1]-myRange[2]), to = (myRange[2]+myRange[2]) ,n=density_bins,bw =bandwidth,adjust = adjust.value)
          }

          density.df = data.frame(x = density.res$x, y  =density.res$y)
          # plot(density.df$x, density.df$y)
          density_df_01 = density.df %>% dplyr::filter(x>=0 & x<=1) %>%
               mutate(y_normalized = y/sum(y),
                      pseudotime_bin_01_min =  lag(x) %>% ifelse(is.na(.), -0.01, .),
                      pseudotime_bin_01_max =  c(x[-length(x)], NA)%>% ifelse(is.na(.), 1.01, .),

                      cumulative.sum = cumsum(y_normalized),
                      pseudotime_density_min = lag(cumulative.sum) %>% ifelse(is.na(.), 0, .),
                      pseudotime_density_max = cumulative.sum,
                      pseudotime_bins_density_based = 1:nrow(.))

          ## SQL based join on a range from S/O
          ## https://stackoverflow.com/questions/46795636/r-dplyr-join-by-range-or-virtual-column
          ## 'from' and 'to' are in quotes because it's reserved language for SQL. Here we use pseudotime_bin_01.min and pseudotime_bin_01.max
          # library('sqldf')
          message('integrating new pseudotime with full dataset...')
          cyclingCells_densityBins = cyclingCells_FULL # %>% dplyr::select(-time)
          # print(dim(cyclingCells_densityBins))
          # SQLbins = as_tibble(sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.pseudotime_bins_density_based from cyclingCells_densityBins
          SQLbins = as_tibble(sqldf::sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.* from cyclingCells_densityBins
                join density_df_01 on cyclingCells_densityBins.PSEUDOTIME_TO_DBPN > density_df_01.pseudotime_bin_01_min and
                          cyclingCells_densityBins.PSEUDOTIME_TO_DBPN <= density_df_01.pseudotime_bin_01_max"))
          # print(head(SQLbins))
          cyclingCells_FULL = cyclingCells_FULL %>%
               bind_cols(SQLbins %>% dplyr::select(-PSEUDOTIME_TO_DBPN))

          # cyclingCells.PST_ADJUSTMENT = cyclingCells_FULL %>%
          #      group_by(pseudotime_bins_density_based) %>%
          #      summarize(count = n(),
          #                min = min(PSEUDOTIME_TO_DBPN),
          #                max = max(PSEUDOTIME_TO_DBPN)) %>%
          #      mutate(percentage = count/sum(count),
          #             percentage_sum  = cumsum(percentage))
          cyclingCells_FULL= cyclingCells_FULL %>%
               # left_join(cyclingCells.PST_ADJUSTMENT) %>%
               group_by(pseudotime_bins_density_based) %>%
               mutate(PSEUDOTIME_NORMALIZED = scales::rescale(PSEUDOTIME_TO_DBPN, to = c(unique(pseudotime_density_min), unique(pseudotime_density_max) )))# %>%

          cyclingCells_FULL = cyclingCells_FULL[, !colnames(cyclingCells_FULL) %in% c('x','y','y_normalized','pseudotime_bin_01_min','pseudotime_bin_01_max','cumulative.sum','pseudotime_density_min','pseudotime_density_max','pseudotime_bins_density_based') ] ## remove excess columns that users don't need to see
          if (time_col_present == TRUE){
               cyclingCells_FULL = cyclingCells_FULL %>% left_join(time_df, by = 'cell.id')
          }

          if (RETURN_DENSITY_COORDINATES ==TRUE){
               return(list(df.dbpn = cyclingCells_FULL %>% ungroup(), density_df_01 = density_df_01))
          }
          message('density-based pseudotime normalization complete! Column name is now PSEUDOTIME_NORMALIZED')
          return(cyclingCells_FULL %>%ungroup())
     }
}



# # ## performs single value decomposition and geometric sketching for downsampling
# # library(reticulate)
# # library(rsvd)
# # geosketch <- import('geosketch')
# GEOMETRIC_SKETCHING = function(data, features,sketch.size=5000, k_components=5, RSVD=TRUE, scale_center = TRUE){
#      geosketch <- reticulate::import('geosketch')
#      data_input = data %>%  ungroup() %>% dplyr::select(all_of(features))
#      if( RSVD == TRUE){
#           s <- rsvd::rsvd(data_input %>% scale() %>% as.matrix(), k=k_components)
#           X.pcs <- s$u %*% diag(s$d)
#      } else {
#           if(scale_center== TRUE){
#                X.pcs =data_input %>% scale() %>% as.matrix()
#           } else {
#                X.pcs =data_input
#           }
#      }
#      sketch.indices <- geosketch$gs(X.pcs, as.integer(sketch.size), one_indexed = TRUE) %>% unlist() ## For R version: must have one_indexed = TRUE
#      return(data[sketch.indices, ])
# }
