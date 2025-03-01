
#' @importFrom stats density
#' @importFrom sqldf sqldf
#' @importFrom RANN nn2
#' @importFrom scales rescale
#' @importFrom graphics pairs
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats cov
#' @importFrom stats dist
#' @importFrom stats gaussian
#' @importFrom stats lag
#' @importFrom stats mahalanobis
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom stats prcomp
#' @importFrom stats predict
#' @importFrom stats time
#' @importFrom utils capture.output
#' @importFrom utils head
#'

#' @title hafezify
#' @description Takes as input your dataframe/tibble and outputs the same object with necessary columns to run other hafez package functions.
#' @export
hafezify = function(dataset){

     if ('cell.id' %in% colnames(dataset)){
          if (any(duplicated(dataset$cell.id))){
               message('duplicated cell.id names detected, modifying cell.id column to have unique identifiers.')
               dataset$cell.id = make.unique(dataset$cell.id)
               return(dataset)
          } else {
               message('cell.id column meets criteria')
               return(dataset)
          }
     } else {
          message('cell.id column not detected. Creating column with unique identifiers')
          dataset$cell.id = paste0('c',1:nrow(dataset))
          return(dataset)
     }
}



#' @title computePseudotimeBins
#' @description Discretizes pseudotime bins
#' @export
computePseudotimeBins = function(dataset,pseudotime_column, interval_sequence = seq(0,1,0.01)){
     dataset = dataset %>%
          mutate(pseudotime_bins =findInterval(!!sym(pseudotime_column), interval_sequence), ## ~ 180 bins
                 pseudotime_bins = as.numeric(pseudotime_bins))
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
#' @description Calculates kernel density estimate for groups specified in data.
#' @export
calculate_cell_density_by_group = function(dataset, group_splits = NULL, pseudotime_column_name="PSEUDOTIME_NORMALIZED",adjust.value = 1, n = 512,from=NULL, to =NULL, sep_operator = '_-_',  COMPUTE_CIRCULAR = FALSE ,count_threshold = 100,verbose = FALSE, ...){

     if (is.null(from)){
          from = min(dataset[[pseudotime_column_name]],na.rm = T)
     }
     if (is.null(to)){
          to = max(dataset[[pseudotime_column_name]],na.rm = T)
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


#' @title mahalanobis_distance_to_landmark_group
#' @description Calculates mahalanobis distance to a set of landmark cells
#' @noRd
mahalanobis_distance_to_landmark_group = function(dataset, features, reference_group_column, reference_group_name, referenceSelf = FALSE, sampleInternally = NULL, tol=1e-20,precalculated_mean = NULL,precalculated_cov = NULL,CELL_COUNT_THRESHOLD=50){ ## tol=1e-20

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

     mahalanobis_distance = tryCatch({
          mahalanobis_distance <-  mahalanobis(mahalanobis.input, center, cov_matrix)
          return(mahalanobis_distance)
     }, error = function(e) {
          # message(unique(dd$feature))
          # You can log the error message if needed
          message(conditionMessage(e))
          return(data.frame(cell.id = dataset$cell.id, mahalanobis_distance=NA))
     })

     # mahalanobis_distance <- mahalanobis(mahalanobis.input, center, cov_matrix, tol = tol)
     return(data.frame(cell.id = dataset$cell.id, mahalanobis_distance=mahalanobis_distance))
}

#' @title pseudotime_mapped_mahalanobis_analysis
#' @description takes in "landmark" for landmark mapped analysis,  or "all"
#' @export
pseudotime_mapped_mahalanobis_analysis = function(dataset, method = c('landmark_mapped','all'), pseudotime_bin_column='pseudotime_bins', features, reference_group_column=NULL, reference_group_name=NULL,CELL_COUNT_THRESHOLD=50){
     dataset$'pseudotime_bins' = dataset[[pseudotime_bin_column]]
     if (method == 'landmark'){
          PM_mahaalanobis_results = dataset %>%
               group_by(pseudotime_bins) %>%
               group_map(~mahalanobis_distance_to_landmark_group(., features, reference_group_column, reference_group_name, referenceSelf = FALSE,  tol=tol, CELL_COUNT_THRESHOLD=CELL_COUNT_THRESHOLD)) %>%
               bind_rows()
     } else if (method == 'all'){
          PM_mahaalanobis_results = dataset %>%
               group_by(pseudotime_bins) %>%
               group_map(~mahalanobis_distance_to_landmark_group(., features, reference_group_column, reference_group_name, referenceSelf = TRUE,  tol=tol, CELL_COUNT_THRESHOLD=CELL_COUNT_THRESHOLD)) %>%
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
# ############################################################
# OLD VERSION, WORKS
# ############################################################
# ############################################################
#' ## column_to_normalize: name of column to normalize in dataset
#' #' title hafez_DBPN
#' #' description
#' #' use the normalize_by_sample_column to specify the column name that contains the column identifier to calculate density distributions for each subject/sample of choice, then normalize the pseudotime by the density estimate
#' #' export
#' hafez_DBPN = function(dataset, density_bins = 1024, normalize_by_sample_column = NULL, dataset.subset_to_use = NULL ,column_to_normalize=NULL, bandwidth = 'nrd0', adjust.value = 1, CIRCULARIZE = FALSE, RETURN_DENSITY_COORDINATES = FALSE){
#'      dataset$x = NULL
#'      dataset$y = NULL
#'      dataset$y_normalized = NULL
#'      dataset$pseudotime_bin_01_min = NULL
#'      dataset$pseudotime_bin_01_max = NULL
#'      dataset$cumulative.sum = NULL
#'      dataset$pseudotime_density_min = NULL
#'      dataset$pseudotime_density_max = NULL
#'      dataset$pseudotime_bins_density_based = NULL
#'      dataset$count = NULL
#'      dataset$percentage_sum = NULL
#'
#'      if (!is.null(dataset.subset_to_use)){
#'           use_subset_of_cells = TRUE
#'           dataset.subset_to_use$PSEUDOTIME_TO_DBPN = dataset.subset_to_use[[column_to_normalize]]
#'      } else {
#'           use_subset_of_cells = FALSE
#'      }
#'      dataset$PSEUDOTIME_TO_DBPN = dataset[[column_to_normalize]]
#'      myRange = c(floor(min(dataset$PSEUDOTIME_TO_DBPN , na.rm =T)), ceiling(max(dataset$PSEUDOTIME_TO_DBPN ,na.rm=T)) ) %>% as.numeric()## from and to values for density calculation
#'      if(is.null(column_to_normalize)){
#'           stop('must specify the pseudotime column to normalize')
#'      }
#'
#'      ## sql hates the a column labeled time so we remove it if it's present
#'      time_col_present = 'time' %in% colnames(dataset)
#'      if (time_col_present == TRUE){
#'           time_df = dataset %>% ungroup() %>% dplyr::select(time, cell.id)
#'           dataset = dataset %>% dplyr::select(-time)
#'
#'           if (use_subset_of_cells == TRUE){
#'                dataset.subset_to_use = dataset.subset_to_use %>% dplyr::select(-time)
#'
#'           }
#'      }
#'      if (use_subset_of_cells == FALSE){
#'           message('normalizing on all cells')
#'
#'           cyclingCells = dataset
#'           density.res = density(cyclingCells$PSEUDOTIME_TO_DBPN, from = myRange[1], to = myRange[2], n = density_bins,adjust = adjust.value)
#'           dataset.frame(x = density.res$x, y  =density.res$y)
#'           ## repeat y values on either side to make it connected/circular
#'           if (CIRCULARIZE == T){
#'                density.res = density(c(cyclingCells$PSEUDOTIME_TO_DBPN-1, cyclingCells$PSEUDOTIME_TO_DBPN, cyclingCells$PSEUDOTIME_TO_DBPN+1), from = (myRange[1]-myRange[2]), to = ( myRange[2]+ myRange[2]) ,n=density_bins,bw =bandwidth,adjust = adjust.value)
#'           }
#'           density.df = data.frame(x = density.res$x, y  =density.res$y)
#'           # plot(density.df$x, density.df$y)
#'           density_df_01 = density.df %>% dplyr::filter(x>=0 & x<=1) %>%
#'                mutate(y_normalized = y/sum(y),
#'                       pseudotime_bin_01_min =  lag(x) %>% ifelse(is.na(.), -0.01, .),
#'                       pseudotime_bin_01_max =  c(x[-length(x)], NA)%>% ifelse(is.na(.), 1.01, .),
#'
#'                       cumulative.sum = cumsum(y_normalized),
#'                       pseudotime_density_min = lag(cumulative.sum) %>% ifelse(is.na(.), 0, .),
#'                       pseudotime_density_max = cumulative.sum,
#'                       pseudotime_bins_density_based = 1:nrow(.))
#'
#'           ## SQL based join on a range from S/O
#'           ## https://stackoverflow.com/questions/46795636/r-dplyr-join-by-range-or-virtual-column
#'           ## 'from' and 'to' are in quotes because it's reserved language for SQL. Here we use pseudotime_bin_01.min and pseudotime_bin_01.max
#'           # library('sqldf')
#'           message('combining...')
#'           cyclingCells_densityBins = cyclingCells # %>% dplyr::select(-time)
#'           # SQLbins = as_tibble(sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.pseudotime_bins_density_based from cyclingCells_densityBins
#'           SQLbins = as_tibble(sqldf::sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.* from cyclingCells_densityBins
#'                 join density_df_01 on cyclingCells_densityBins.PSEUDOTIME_TO_DBPN > density_df_01.pseudotime_bin_01_min and
#'                           cyclingCells_densityBins.PSEUDOTIME_TO_DBPN <= density_df_01.pseudotime_bin_01_max"))
#'           cyclingCells = cyclingCells %>%
#'                dplyr::bind_cols(SQLbins %>% dplyr::select(-PSEUDOTIME_TO_DBPN))
#'
#'           cyclingCells.PST_ADJUSTMENT = cyclingCells %>%
#'                group_by(pseudotime_bins_density_based) %>%
#'                summarize(count = n(),
#'                          min = min(PSEUDOTIME_TO_DBPN),
#'                          max = max(PSEUDOTIME_TO_DBPN)) %>%
#'                mutate(percentage = count/sum(count),
#'                       percentage_sum  = cumsum(percentage))
#'           cyclingCells= cyclingCells %>%
#'                left_join(cyclingCells.PST_ADJUSTMENT) %>%
#'                group_by(pseudotime_bins_density_based) %>%
#'                mutate(PSEUDOTIME_NORMALIZED = scales::rescale(PSEUDOTIME_TO_DBPN, to = c(unique(pseudotime_density_min), unique(pseudotime_density_max) )))# %>%
#'
#'           if (time_col_present == TRUE){
#'                cyclingCells = cyclingCells %>% left_join(time_df, by = 'cell.id')
#'           }
#'           if (RETURN_DENSITY_COORDINATES ==TRUE){
#'                return(list(df.dbpn = cyclingCells %>% ungroup(), density_df_01 = density_df_01))
#'           }
#'           return(cyclingCells %>% ungroup())
#'
#'      } else if (use_subset_of_cells == TRUE){
#'           message('computing pseudotime normalization on the provided subset of cells...')
#'           cyclingCells_FULL = dataset
#'           cyclingCells = dataset.subset_to_use
#'
#'           if (!is.null(normalize_by_sample_column)){
#'                if (CIRCULARIZE == T){
#'                     density.res = density(c(cyclingCells$PSEUDOTIME_TO_DBPN-1, cyclingCells$PSEUDOTIME_TO_DBPN, cyclingCells$PSEUDOTIME_TO_DBPN+1),  from = (myRange[1]-myRange[2]), to = (myRange[2]+myRange[2]) ,n=density_bins,bw =bandwidth,adjust = adjust.value)
#'                } else {
#'
#'                     density.df = lapply(split(cyclingCells,cyclingCells[[normalize_by_sample_column]]), function(dd){
#'                          density.res = density(cyclingCells$PSEUDOTIME_TO_DBPN, from =  myRange[1], to = myRange[2], n = density_bins,adjust = adjust.value)
#'                          density.df = data.frame(x = density.res$x, y  =density.res$y) %>%
#'                               mutate(normalize_by_sample_column=unique(dd[[normalize_by_sample_column]]))
#'
#'                     }) %>% bind_rows() %>%
#'                          group_by(x)%>%
#'                          summarize(y=mean(y))
#'                }
#'
#'           } else{
#'                if (CIRCULARIZE == T){
#'                     density.res = density(c(cyclingCells$PSEUDOTIME_TO_DBPN-1, cyclingCells$PSEUDOTIME_TO_DBPN, cyclingCells$PSEUDOTIME_TO_DBPN+1),  from = (myRange[1]-myRange[2]), to = (myRange[2]+myRange[2]) ,n=density_bins,bw =bandwidth,adjust = adjust.value)
#'                } else {
#'                     density.res = density(cyclingCells$PSEUDOTIME_TO_DBPN, from =  myRange[1], to = myRange[2], n = density_bins,adjust = adjust.value)
#'                     density.df = data.frame(x = density.res$x, y  =density.res$y)
#'                }
#'           }
#'
#'           density_df_01 = density.df %>%
#'                # dplyr::filter(x>=0 & x<=1) %>%
#'                mutate(y_normalized = y/sum(y),
#'                       pseudotime_bin_01_min =  lag(x) %>% ifelse(is.na(.), -0.01, .),
#'                       pseudotime_bin_01_max =  c(x[-length(x)], NA)%>% ifelse(is.na(.), 1.01, .),
#'
#'                       cumulative.sum = cumsum(y_normalized),
#'                       pseudotime_density_min = lag(cumulative.sum) %>% ifelse(is.na(.), 0, .),
#'                       pseudotime_density_max = cumulative.sum,
#'                       pseudotime_bins_density_based = 1:nrow(.))
#'
#'            # data.frame(x = density.res$x, y  =density.res$y)
#'           ## repeat y values on either side to make it connected/circular
#'
#'
#'           # plot(density.df$x, density.df$y)
#'           # density_df_01 = density.df %>%
#'           #      dplyr::filter(x>=0 & x<=1) %>%
#'           #      mutate(y_normalized = y/sum(y),
#'           #             pseudotime_bin_01_min =  lag(x) %>% ifelse(is.na(.), -0.01, .),
#'           #             pseudotime_bin_01_max =  c(x[-length(x)], NA)%>% ifelse(is.na(.), 1.01, .),
#'           #
#'           #             cumulative.sum = cumsum(y_normalized),
#'           #             pseudotime_density_min = lag(cumulative.sum) %>% ifelse(is.na(.), 0, .),
#'           #             pseudotime_density_max = cumulative.sum,
#'           #             pseudotime_bins_density_based = 1:nrow(.))
#'
#'           ## SQL based join on a range from S/O
#'           ## https://stackoverflow.com/questions/46795636/r-dplyr-join-by-range-or-virtual-column
#'           ## 'from' and 'to' are in quotes because it's reserved language for SQL. Here we use pseudotime_bin_01.min and pseudotime_bin_01.max
#'           # library('sqldf')
#'           message('integrating new pseudotime with full dataset...')
#'           cyclingCells_densityBins = cyclingCells_FULL # %>% dplyr::select(-time)
#'           # print(dim(cyclingCells_densityBins))
#'           # SQLbins = as_tibble(sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.pseudotime_bins_density_based from cyclingCells_densityBins
#'           SQLbins = as_tibble(sqldf::sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.* from cyclingCells_densityBins
#'                 join density_df_01 on cyclingCells_densityBins.PSEUDOTIME_TO_DBPN > density_df_01.pseudotime_bin_01_min and
#'                           cyclingCells_densityBins.PSEUDOTIME_TO_DBPN <= density_df_01.pseudotime_bin_01_max"))
#'           # print(head(SQLbins))
#'           cyclingCells_FULL = cyclingCells_FULL %>%
#'                bind_cols(SQLbins %>% dplyr::select(-PSEUDOTIME_TO_DBPN))
#'
#'           # cyclingCells.PST_ADJUSTMENT = cyclingCells_FULL %>%
#'           #      group_by(pseudotime_bins_density_based) %>%
#'           #      summarize(count = n(),
#'           #                min = min(PSEUDOTIME_TO_DBPN),
#'           #                max = max(PSEUDOTIME_TO_DBPN)) %>%
#'           #      mutate(percentage = count/sum(count),
#'           #             percentage_sum  = cumsum(percentage))
#'           cyclingCells_FULL= cyclingCells_FULL %>%
#'                # left_join(cyclingCells.PST_ADJUSTMENT) %>%
#'                group_by(pseudotime_bins_density_based) %>%
#'                mutate(PSEUDOTIME_NORMALIZED = scales::rescale(PSEUDOTIME_TO_DBPN, to = c(unique(pseudotime_density_min), unique(pseudotime_density_max) )))# %>%
#'
#'           cyclingCells_FULL = cyclingCells_FULL[, !colnames(cyclingCells_FULL) %in% c('x','y','y_normalized','pseudotime_bin_01_min','pseudotime_bin_01_max','cumulative.sum','pseudotime_density_min','pseudotime_density_max','pseudotime_bins_density_based') ] ## remove excess columns that users don't need to see
#'           if (time_col_present == TRUE){
#'                cyclingCells_FULL = cyclingCells_FULL %>% left_join(time_df, by = 'cell.id')
#'           }
#'
#'           if (RETURN_DENSITY_COORDINATES ==TRUE){
#'                return(list(df.dbpn = cyclingCells_FULL %>% ungroup(), density_df_01 = density_df_01))
#'           }
#'           message('density-based pseudotime normalization complete! Column name is now PSEUDOTIME_NORMALIZED')
#'           return(cyclingCells_FULL %>%ungroup())
#'      }
#' }

## column_to_normalize: name of column to normalize in dataset
#' @title hafez_DBPN
#' @description
#' use the normalize_by_sample_column to specify the column name that contains the column identifier to calculate density distributions for each subject/sample of choice, then normalize the pseudotime by the density estimate
#' @export
hafez_DBPN = function(dataset, density_bins = 1024, normalize_by_sample_column = NULL, dataset.subset_to_use = NULL ,column_to_normalize=NULL, bandwidth = 'nrd0', adjust.value = 1, RETURN_DENSITY_COORDINATES = FALSE, new_dbp_name='PSEUDOTIME_NORMALIZED'){
     dataset$x = NULL
     dataset$y = NULL
     ˆdataset$y_normalized = NULL
     dataset$pseudotime_bin_01_min = NULL
     dataset$pseudotime_bin_01_max = NULL
     dataset$cumulative.sum = NULL
     dataset$pseudotime_density_min = NULL
     dataset$pseudotime_density_max = NULL
     dataset$pseudotime_bins_density_based = NULL
     dataset$count = NULL
     dataset$percentage_sum = NULL

     ## must be character
     new_dbp_name=as.character(new_dbp_name)

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

          density.res = density(dataset$PSEUDOTIME_TO_DBPN, from = myRange[1], to = myRange[2], n = density_bins,adjust = adjust.value)
          dataset.frame(x = density.res$x, y  =density.res$y)
          density.df = data.frame(x = density.res$x, y  =density.res$y)
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
          dataset_densityBins = dataset # %>% dplyr::select(-time)
          # SQLbins = as_tibble(sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.pseudotime_bins_density_based from cyclingCells_densityBins
          SQLbins = as_tibble(sqldf::sqldf("select dataset_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.* from dataset_densityBins
                join density_df_01 on dataset_densityBins.PSEUDOTIME_TO_DBPN > density_df_01.pseudotime_bin_01_min and
                          dataset_densityBins.PSEUDOTIME_TO_DBPN <= density_df_01.pseudotime_bin_01_max"))
          dataset = dataset %>%
               dplyr::bind_cols(SQLbins %>% dplyr::select(-PSEUDOTIME_TO_DBPN))

          dataset.PST_ADJUSTMENT = dataset %>%
               group_by(pseudotime_bins_density_based) %>%
               summarize(count = n(),
                         min = min(PSEUDOTIME_TO_DBPN),
                         max = max(PSEUDOTIME_TO_DBPN)) %>%
               mutate(percentage = count/sum(count),
                      percentage_sum  = cumsum(percentage))
          dataset= dataset %>%
               left_join(cyclingCells.PST_ADJUSTMENT) %>%
               group_by(pseudotime_bins_density_based) %>%
               mutate(!!sym(new_dbp_name) := scales::rescale(PSEUDOTIME_TO_DBPN, to = c(unique(pseudotime_density_min), unique(pseudotime_density_max) )))# %>%

          if (time_col_present == TRUE){
               dataset = dataset %>% left_join(time_df, by = 'cell.id')
          }
          if (RETURN_DENSITY_COORDINATES ==TRUE){
               return(list(df.dbpn = dataset %>% ungroup(), density_df_01 = density_df_01))
          }
          return(dataset %>% ungroup())

     } else if (use_subset_of_cells == TRUE){
          message('computing pseudotime normalization on the provided subset of cells...')
          dataset_FULL = dataset
          # dataset_train = dataset.subset_to_use

          if (!is.null(normalize_by_sample_column)){

               density.df = lapply(split(dataset.subset_to_use,dataset.subset_to_use[[normalize_by_sample_column]]), function(dd){
                    density.res = density(dd$PSEUDOTIME_TO_DBPN, from =  myRange[1], to = myRange[2], n = density_bins,adjust = adjust.value)
                    density.df = data.frame(x = density.res$x, y  =density.res$y) %>%
                         mutate(normalize_by_sample_column=unique(dd[[normalize_by_sample_column]]))

               }) %>% bind_rows() %>%
                    group_by(x)%>%
                    summarize(y=mean(y))


          } else{
               density.res = density(dataset.subset_to_use$PSEUDOTIME_TO_DBPN, from =  myRange[1], to = myRange[2], n = density_bins,adjust = adjust.value)
               density.df = data.frame(x = density.res$x, y  =density.res$y)
          }

          density_df_01 = density.df %>%
               # dplyr::filter(x>=0 & x<=1) %>%
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
          # cyclingCells_densityBins = cyclingCells_FULL # %>% dplyr::select(-time)
          # print(dim(cyclingCells_densityBins))
          # SQLbins = as_tibble(sqldf("select cyclingCells_densityBins.PSEUDOTIME_TO_DBPN, density_df_01.pseudotime_bins_density_based from cyclingCells_densityBins
          SQLbins = as_tibble(sqldf::sqldf("select dataset_FULL.PSEUDOTIME_TO_DBPN, density_df_01.* from dataset_FULL
                join density_df_01 on dataset_FULL.PSEUDOTIME_TO_DBPN > density_df_01.pseudotime_bin_01_min and
                          dataset_FULL.PSEUDOTIME_TO_DBPN <= density_df_01.pseudotime_bin_01_max"))
          # print(head(SQLbins))
          dataset_FULL = dataset_FULL %>%
               bind_cols(SQLbins %>% dplyr::select(-PSEUDOTIME_TO_DBPN))

          # cyclingCells.PST_ADJUSTMENT = cyclingCells_FULL %>%
          #      group_by(pseudotime_bins_density_based) %>%
          #      summarize(count = n(),
          #                min = min(PSEUDOTIME_TO_DBPN),
          #                max = max(PSEUDOTIME_TO_DBPN)) %>%
          #      mutate(percentage = count/sum(count),
          #             percentage_sum  = cumsum(percentage))
          dataset_FULL= dataset_FULL %>%
               # left_join(cyclingCells.PST_ADJUSTMENT) %>%
               group_by(pseudotime_bins_density_based) %>%
               mutate(!!sym(new_dbp_name) := scales::rescale(PSEUDOTIME_TO_DBPN, to = c(unique(pseudotime_density_min), unique(pseudotime_density_max) )))# %>%

          dataset_FULL = dataset_FULL[, !colnames(dataset_FULL) %in% c('x','y','y_normalized','pseudotime_bin_01_min','pseudotime_bin_01_max','cumulative.sum','pseudotime_density_min','pseudotime_density_max','pseudotime_bins_density_based') ] ## remove excess columns that users don't need to see
          if (time_col_present == TRUE){
               dataset_FULL = dataset_FULL %>% left_join(time_df, by = 'cell.id')
          }

          if (RETURN_DENSITY_COORDINATES ==TRUE){
               return(list(df.dbpn = dataset_FULL %>% ungroup(), density_df_01 = density_df_01))
          }
          message(paste0('density-based pseudotime normalization complete! Column name is now ',new_dbp_name))
          return(dataset_FULL %>%ungroup())
     }}

##################################################
## landmark  preprocessing and normalization ##
##################################################
#' @title hafez_landmark_processing
#' @description Performs landmark preprocessing
#' @export
hafez_landmark_processing = function(train, full_data, features, method = 'PCA', return_object = FALSE ){

     if (method == 'raw'){
          train.scale = train %>% dplyr::select(all_of(features)) %>% scale()
          scale_attributes = attributes(train.scale)
          data.scaled = full_data %>%
               dplyr::select(all_of(features)) %>%
               scale(., center  = scale_attributes$`scaled:center`, scale = scale_attributes$`scaled:scale`)
          landmark.dr.df = full_data %>%
               dplyr::select(-all_of(features)) %>%
               bind_cols(data.scaled)
          message('Hafez landmark processing complete!')
          return(landmark.dr.df)
     }


     if (method == 'PCA'){
          train.scale = train %>% dplyr::select(all_of(features)) %>% scale()
          pca.res = prcomp(train.scale , scale = F, center = F)
          scale_attributes = attributes(train.scale)
          data.scaled = full_data %>%
               dplyr::select(all_of(features)) %>%
               scale(., center  = scale_attributes$`scaled:center`, scale = scale_attributes$`scaled:scale`)
          my_new_dr = colnames(pca.res$rotation)
          landmark.dr.df = data.scaled %*% pca.res$rotation %>%
               data.frame() %>%
               bind_cols(full_data %>% dplyr::select(-any_of(my_new_dr)),.)
          if (return_object == TRUE){
               message('Hafez landmark processing complete!')
               return(list(data = landmark.dr.df, dr_object = pca.res, scale_attributes=scale_attributes))
          } else{
               message('Hafez landmark processing complete!')
               return(landmark.dr.df)
          }
     }
}



#' @title wrangle_data_for_knn_index_resolution
#' @keywords internal
#' @noRd
## RANN::nn2 has weird behavior with the query search with the indexes of the query dataset (full) and landmark dataset.
## to fix this, the index of the query (fill) and landmark datasets must match the first rows. All remaining indexes won't matter
wrangle_data_for_knn_index_resolution = function(full_data,LM_col=NULL,LM_group=NULL, LM_data=NULL){
     if(!'cell.id' %in% colnames(full_data) ){
          stop('input data must have a cell.id column')
     }
     if( !is.null(LM_data)){
          if(!'cell.id' %in% colnames(LM_data) ){
               stop('landmark data must have a cell.id column')
          }
          QUERY_DF = LM_data %>%
               bind_rows(full_data %>% dplyr::filter(!cell.id %in%  LM_data$cell.id))
     } else {

          if(!LM_col %in% colnames(full_data) ){
               stop(paste0('input data does not contain the provided column for the "LM_col" arguement: ',LM_col))
          }

          if(!LM_group %in% unique(full_data[[LM_col]]) ){
               stop(paste0(LM_col, ' does not contain ',LM_group))
          }
          # LM_col = 'condition'
          # LM_group = 'WT'
          LM_df = full_data %>% .[.[[LM_col]]==LM_group,]
          PROJECT_DF = full_data %>% .[!.[[LM_col]]==LM_group,]
          QUERY_DF = bind_rows(LM_df,PROJECT_DF)
     }

     return(QUERY_DF)
}


#' @title find_closest_neighbor_distance
#' @description Calculates distance to landmark cells.
#' @export
##k_ave is the # of neighbors to use in the distance output. Anything more than 1 will take the average of each index it is proximal to.
## a value of -1 for k_ave finds the closest non-self neighbor, which is equivalently the  2nd closest neighbor.
##k_ave is the # of neighbors to use in the distance output. Anything more than 1 will take the average of each index it is proximal to.
## a value of -1 for k_ave finds the closest non-self neighbor, which is equivalently the  2nd closest neighbor.
find_closest_neighbor_distance = function(full_data, LM_col=NULL, LM_group, features, k_ave =2, k_solo = NULL, return_knn_graph=FALSE, return_nearest_cellid=FALSE){
     if (is.null(LM_col)){
          message('performing search  on full dataset...')
          df_landmark_features_only=full_data%>% dplyr::select(all_of(features))
          full_data_features_only=full_data%>% dplyr::select(all_of(features))
          restructured_data=full_data
     } else {
          message('restructuring data for knn search...')
          restructured_data = wrangle_data_for_knn_index_resolution(full_data, LM_col, LM_group)
          df_landmark_features_only = restructured_data %>% .[.[[LM_col]]==LM_group,] %>% dplyr::select(all_of(features))
          full_data_features_only = restructured_data %>% dplyr::select(all_of(features))
     }


     message('performing knn search...')
     if (is.null(LM_col)){
          knn_res_full = RANN::nn2(data = df_landmark_features_only, searchtype = 'standard',eps=0,treetyp = 'kd',
                                   # radius = 100, ## for exact search
                                   k = k_ave+1)
     } else {
          knn_res_full = RANN::nn2(data = df_landmark_features_only, query = full_data_features_only,searchtype = 'standard',eps=0,treetyp = 'kd',
                                   # radius = 100, ## for exact search
                                   k = k_ave+1)
     }

     if(return_knn_graph ==T){
          return(knn_res_full)
     }

     if (return_nearest_cellid==T){
          score_res_project=knn_res_full$nn.idx[,1:k_ave] %>% data.frame(check.rows = F, check.names = T)
          score_res_project$cell.id = restructured_data$cell.id
          df_landmark_features_only$cell.id = restructured_data %>% .[.[[LM_col]]==LM_group,] %>%.$cell.id
          return( list(df_landmark_features_only=df_landmark_features_only, score_res_project=score_res_project) )
     }

     message('extracting closest cell and returning scores with original indexing...')
     if (!is.null(k_solo)){
          if (k_solo==1  ){
               message('returning distance of nearest landmark neighbor')
               knn_res_full$nn.dists[,1] = ifelse(knn_res_full$nn.idx[,1] == 1:nrow(knn_res_full$nn.idx), knn_res_full$nn.dists[,2], knn_res_full$nn.dists[,1])
               score_res_project= knn_res_full$nn.dists[,1] %>% unlist()
          }
     } else if (k_ave > 1 ){
          knn_res_full$nn.dists[,1] = ifelse(knn_res_full$nn.idx[,1] == 1:nrow(knn_res_full$nn.idx), knn_res_full$nn.dists[,(k_ave+1)], knn_res_full$nn.dists[,1])
          # return(knn_res_full$nn.dists)
          score_res_project= rowMeans(knn_res_full$nn.dists[,1:k_ave])
     } else {
          score_res_project = knn_res_full$nn.dists[,1] %>% unlist()  #%>% data.frame(dist_from_nearest_cell = .)
     }
     names(score_res_project) = restructured_data$cell.id
     # table(score_res_project$cell.id[match(full_data$cell.id, score_res_project$cell.id)] ==full_data$cell.id)
     score_res_project=score_res_project[match(full_data$cell.id, names(score_res_project))]
     return(score_res_project)
}




# MEELAD_WRAPPER = function(FULL_DATA, features, features_for_start_cell_id, PERFORM_OOS=FALSE, OOS_idx=NULL,return_pseudotime_only=TRUE,NumNodes=5,lambda = 0.01, mu = 0.01, branch_type = c('tree','curve','circle'),TRAINING_VARIABLE='WT', return_node_pos = FALSE, use_start_label=NULL, start_label_column_category = c('celltype','HSC')){
#      FULL_DATA=FULL_DATA%>% ungroup()
#      if (PERFORM_OOS==FALSE){
#           FULL_DATA_TRAIN=FULL_DATA
#
#      } else if(PERFORM_OOS==TRUE){
#           if (is.vector(OOS_idx)) {
#                FULL_DATA_TRAIN = FULL_DATA[OOS_idx, ]
#           } else {
#                ## if providing a pre-filtered training data
#                FULL_DATA_TRAIN=OOS_idx %>% ungroup()
#           }
#      }
#
#      if (branch_type == 'circle'){
#
#           ELPPI_OUTPUT = ELASTIC_PRINCIPAL_GRAPH_FUNCTION(CONTROL_DATA = FULL_DATA_TRAIN,
#                                                           FULL_DATA = FULL_DATA,
#                                                           FEATURES = features,Mu = 1.5, Lambda = 0.01, ProbPoint = 0.8,
#                                                           ## changing probPoint
#                                                           LABELS = NULL,  NumNodes = 20, nReps =30, verbose = T )
#
#           proj.data = data.frame(ELPPI_OUTPUT$projection_struture$X_projected)
#           colnames(proj.data) = features_LMDPDR[1:ncol(proj.data)]
#           node.pos = data.frame(ELPPI_OUTPUT$projection_struture$NodePositions)
#           colnames(node.pos) = features_LMDPDR[1:ncol(node.pos)]
#           #
#           # ggplot(ELPPI_OUTPUT$data,
#           #        aes(x = LMDPDR1, y = LMDPDR2, color = gate)) +
#           #      geom_point(size = 0.5) +
#           #      geom_point(data = node.pos,inherit.aes = F, aes(x=LMDPDR1, y=LMDPDR2),size= 4) +
#           #      viridis::scale_color_viridis(option ='magma', discrete = T)
#
#           pst=ELPPI_OUTPUT %>% dplyr::select(cell.id, pseudotime)
#           colnames(pst) = c('cell.id','MA_OOS_path1')
#
#           FULL_DATA = FULL_DATA %>%
#                # dplyr::select(-any_of(   colnames(pst))) %>%
#                left_join(pst)
#      } else {
#
#
#           ELPIGRAPH_RES <- COMPUTE_G0G1M_TRAJECTORY(CONTROL_DATA = FULL_DATA_TRAIN,
#                                                     FULL_DATA = FULL_DATA,
#                                                     FEATURES = features,
#                                                     CCphase_column = 'gate',## used to label the example plot with node locations
#                                                     NumNodes = NumNodes, Lambda = lambda, Mu = mu,Do_PCA = F,nReps = 30,verbose = FALSE,
#                                                     MaxNumberOfIterations=30,ProbPoint=1,branch_type = branch_type,
#                                                     # verbose = T,
#                                                     drawAccuracyComplexity = FALSE, drawEnergy = FALSE)
#
#
#           ## closest node
#           CLOSEST_CELL_IDX = FULL_DATA_TRAIN %>% ungroup() %>%dplyr::select(any_of(features_for_start_cell_id)) %>% apply(.,1, mean) %>% which.min(.)
#           CLOSEST_CELL_ID = FULL_DATA_TRAIN$cell.id[CLOSEST_CELL_IDX]
#
#           colnames(ELPIGRAPH_RES$node.df) = c(features, 'node','path1')
#           node.pos= ELPIGRAPH_RES$node.df
#
#
#           if (!is.null(use_start_label)){
#
#
#                my_start_label = FULL_DATA_TRAIN %>% dplyr::filter(!!sym(start_label_column_category[1]) == start_label_column_category[2]) %>% summarize_at(features, median)
#                ELPIGRAPH_RES$node.df = ELPIGRAPH_RES$node.df %>% bind_rows(my_start_label)
#           } else {
#                ELPIGRAPH_RES$node.df = ELPIGRAPH_RES$node.df %>% bind_rows(FULL_DATA_TRAIN %>% dplyr::filter(cell.id == CLOSEST_CELL_ID) %>% dplyr::select(any_of(features)))
#           }
#
#
#           START_NODE_ID = dist(ELPIGRAPH_RES$node.df %>% dplyr::select(any_of(features)),method = 'euclidean') %>% as.matrix() %>% .[1:nrow(.),ncol(.)] %>% .[.!=0] %>% .[which.min(.)] %>% names() %>% as.numeric()
#           ## select start node at S-phase
#           pst = hafez_lineages_from_root(ELPIGRAPH_RES, OOS_DATA = FULL_DATA, root = START_NODE_ID,features = features)
#           # print(head(pst))
#           if (PERFORM_OOS==TRUE){
#                colnames(pst) = paste('MA_OOS',colnames(pst),sep = '_')
#
#           } else {
#                colnames(pst) = paste('MA',colnames(pst),sep = '_')
#           }
#
#           if (return_pseudotime_only==TRUE){
#                if(return_node_pos==TRUE){
#                     message('returning node positions...')
#                     return(list(pst = data.frame(pst),
#                                 node_df = node.pos))
#                } else {
#                     return(data.frame(pst))
#
#                }
#
#           }
#
#
#           FULL_DATA = FULL_DATA %>%
#                dplyr::select(-any_of(   colnames(pst))) %>%
#                bind_cols(pst)
#      }
#
#      if(return_node_pos==TRUE){
#           return(list(FULL_DATA = FULL_DATA,
#                       node_df = node.pos))
#      }
#
#      return(FULL_DATA)
#
# }


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
