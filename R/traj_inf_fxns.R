
#' @importFrom stats density
#' @importFrom sqldf sqldf
#' @importFrom RANN nn2
#' @importFrom scales rescale
#' @import ElPiGraph.R
#' @import sqldf
#' @import igraph


##################################################
# TRAJECTORY INFERENCE FUNCTIONS ##
##################################################

#' Pick the consensus graph out of an ElPiGraph result list
#'
#' When nReps > 1, ElPiGraph fits nReps bootstrap replicates (each to a random
#' ProbPoint fraction of the data) and then appends ONE extra graph built from
#' all of them -- the consensus. It is the LAST element and is tagged
#' ReplicaID == 0 and ProbPoint == 1. When nReps == 1 no consensus is appended
#' and the single replicate carries ReplicaID == 1.
#'
#' Taking [[1]] therefore returns a subsample fit rather than the consensus
#' whenever nReps > 1, which is the default. Select on ReplicaID == 0 rather
#' than on position, so this stays correct if the return order ever changes,
#' and fall back to the last element if the tag is absent.
#'
#' @keywords internal
#' @noRd
epg_consensus <- function(EPG_LIST) {
     if (is.null(EPG_LIST) || !length(EPG_LIST)) {
          stop('empty ElPiGraph result: no graph to select')
     }
     rid <- vapply(EPG_LIST, function(g) {
          v <- g$ReplicaID
          if (is.null(v) || !length(v)) NA_real_ else as.numeric(v)[1]
     }, numeric(1))
     idx <- if (any(rid == 0, na.rm = TRUE)) which(rid == 0)[1] else length(EPG_LIST)
     EPG_LIST[[idx]]
}

## TI_CIRCLE
#' @description TI_CIRCLE
#' @keywords internal
#' @noRd
hafez_TI_circle = function(LM_DATA, FULL_DATA, FEATURES = paste0('PC',1:10), LABELS = NULL,
                           NumNodes = 25, nReps = 5, ProbPoint = 0.6,
                           Lambda = 0.01, Mu = 0.1, Do_PCA = FALSE,
                           RETURN_PROJECTION = FALSE,
                           MaxNumberOfIterations = 20,
                           DRAW_ACCURACY_COMPLEXITY = FALSE, DRAW_ENERGY = FALSE, DRAW_PCA_VIEW = FALSE,
                           verbose = FALSE){

     TRAINING_DATA = LM_DATA %>% dplyr::select(any_of(FEATURES)) %>% as.matrix()
     OOS_DATA = FULL_DATA %>% dplyr::select(any_of(FEATURES)) %>% as.matrix()

     CircleEPG.Boot <- ElPiGraph.R::computeElasticPrincipalCircle(
          X = TRAINING_DATA,
          NumNodes = NumNodes, nReps = nReps, ProbPoint = ProbPoint,
          Lambda = Lambda, Mu = Mu, Do_PCA = Do_PCA,
          RETURN_PROJECTION = FALSE,
          MaxNumberOfIterations = MaxNumberOfIterations,
          drawAccuracyComplexity = DRAW_ACCURACY_COMPLEXITY, drawEnergy = DRAW_ENERGY, drawPCAView = DRAW_PCA_VIEW,
          verbose = verbose
     )

     ## same selection rule as hafez_TI_LINEAR_BRANCH; previously this function
     ## used [[length()]] while that one used [[1]], so the two disagreed.
     CIRCLE_EPG <- epg_consensus(CircleEPG.Boot)
     PartStruct <- ElPiGraph.R::PartitionData(X = OOS_DATA, NodePositions = CIRCLE_EPG$NodePositions)
     ProjStruct <- ElPiGraph.R::project_point_onto_graph(
          X = OOS_DATA,
          NodePositions = CIRCLE_EPG$NodePositions,
          Edges = CIRCLE_EPG$Edges$Edges,
          Partition = PartStruct$Partition
     )
     Circle_Graph <- ElPiGraph.R::ConstructGraph(CIRCLE_EPG)
     Circle_e2e <- ElPiGraph.R::GetSubGraph(Net = Circle_Graph, Structure = 'circle', Circular = TRUE)
     Root <- 1
     SelPaths <- Circle_e2e[sapply(Circle_e2e, function(x){any(x[c(1, length(x))] == Root)})]

     AllPt <- lapply(SelPaths, function(x){
          ElPiGraph.R::getPseudotime(ProjStruct = ProjStruct, NodeSeq = names(x))
     })

     PointsPT <- apply(sapply(AllPt, "[[", "Pt"), 1, function(x){unique(x[!is.na(x)])})
     myPseudotime = PointsPT %>% dynutils::scale_minmax(.)
     FULL_DATA$pseudotime = myPseudotime

     if (RETURN_PROJECTION == TRUE){
          return(list(projection_struture = ProjStruct, data = FULL_DATA))
     } else {
          return(FULL_DATA)
     }
}


# library(igraph)
# library(ElPiGraph.R)
## TI_LINEAR_BRANCH
#' @description TI_LINEAR_BRANCH
#' @keywords internal
#' @noRd
hafez_TI_LINEAR_BRANCH = function(LM_DATA, FULL_DATA, FEATURES = paste0('PC',1:10), LABELS = NULL,
                                  FEATURES_PLOT = NULL,
                                  NumNodes = 10, nReps = 5, ProbPoint = 0.6,
                                  Lambda = 0.01, Mu = 0.1, Do_PCA = FALSE,
                                  MaxNumberOfIterations = 5, CC_PHASE_COLUMN = NULL,
                                  BRANCH_TYPE = 'tree',
                                  DRAW_ACCURACY_COMPLEXITY = FALSE, DRAW_ENERGY = FALSE, DRAW_PCA_VIEW = FALSE,
                                  verbose = FALSE){

     if (is.null(FEATURES_PLOT)) {
          FEATURES_PLOT = FEATURES
     }

     TRAIN_DATA = LM_DATA %>% ungroup() %>% dplyr::select(any_of(FEATURES))

     if (BRANCH_TYPE == 'tree'){
          TreeEPG <- ElPiGraph.R::computeElasticPrincipalTree(
               X = as.matrix(TRAIN_DATA),
               NumNodes = NumNodes, ProbPoint = ProbPoint, nReps = nReps,
               MaxNumberOfIterations = MaxNumberOfIterations, verbose = verbose,
               Lambda = Lambda, Mu = Mu, Do_PCA = Do_PCA,
               drawPCAView = DRAW_PCA_VIEW, drawAccuracyComplexity = DRAW_ACCURACY_COMPLEXITY,
               drawEnergy = DRAW_ENERGY
          )
     } else if (BRANCH_TYPE == 'curve'){
          TreeEPG <- ElPiGraph.R::computeElasticPrincipalCurve(
               X = as.matrix(TRAIN_DATA),
               NumNodes = NumNodes, ProbPoint = ProbPoint, nReps = nReps,
               MaxNumberOfIterations = MaxNumberOfIterations, verbose = verbose,
               Lambda = Lambda, Mu = Mu, Do_PCA = Do_PCA,
               drawPCAView = DRAW_PCA_VIEW, drawAccuracyComplexity = DRAW_ACCURACY_COMPLEXITY,
               drawEnergy = DRAW_ENERGY
          )
     }

     ## With the default nReps = 5, TreeEPG holds 5 bootstrap replicates plus the
     ## consensus. Indexing [[1]] here returned replicate 1, fitted to a random
     ## ProbPoint (default 0.6) fraction of the landmarks -- a different graph
     ## every call, and not the one the bootstrap was run to produce.
     TREE_EPG <- epg_consensus(TreeEPG)

     node.df = TREE_EPG$NodePositions %>% data.frame() %>%
          mutate(node = paste0(1:nrow(.)), path1 = factor(node))

     if (is.null(CC_PHASE_COLUMN)) {
          p.lineages = ggplot(LM_DATA, aes(x = !!sym(FEATURES_PLOT[1]), y = !!sym(FEATURES_PLOT[2]))) +
               theme_minimal() +
               geom_point(aes(fill = auto_annotation), size = 0.75, shape = 21, stroke = 0.1, color = 'black') +
               geom_point(data = node.df, aes(x = X1, y = X2), color = 'red') +
               geom_text(data = node.df, aes(x = X1, y = X2, label = node), color = 'black') +
               guides(fill = guide_legend(override.aes = list(size = 5)))
     } else {
          p.lineages = ggplot(LM_DATA, aes(x = !!sym(FEATURES_PLOT[1]), y = !!sym(FEATURES_PLOT[2]))) +
               theme_minimal() +
               geom_point(aes(fill = !!sym(CC_PHASE_COLUMN)), size = 0.75, shape = 21, stroke = 0.1, color = 'black') +
               geom_point(data = node.df, aes(x = X1, y = X2), color = 'red') +
               geom_text(data = node.df, aes(x = X1, y = X2, label = node), color = 'black') +
               guides(fill = guide_legend(override.aes = list(size = 5)))
     }

     Tree_Graph <- ElPiGraph.R::ConstructGraph(TREE_EPG)
     Tree_e2e <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = 'end2end')
     NodeLabs <- 1:nrow(TREE_EPG$NodePositions)

     return(list(TreeEPG = TreeEPG, Tree_e2e = Tree_e2e, Tree_Graph = Tree_Graph,
                 NodeLabs = NodeLabs, LM_DATA = LM_DATA, node.df = node.df, plot = p.lineages))
}


#' @description hafez_lineages_from_root
#' @keywords internal
#' @noRd
hafez_lineages_from_root = function(COMPUTE_TI_OUTPUT, OOS_DATA, FEATURES, ROOT, RETURN_OBJECTS = FALSE) {
     Tree_e2e = COMPUTE_TI_OUTPUT[['Tree_e2e']]
     TreeEPG = COMPUTE_TI_OUTPUT[['TreeEPG']]
     Tree_Graph = COMPUTE_TI_OUTPUT[['Tree_Graph']]
     NodeLabs = COMPUTE_TI_OUTPUT[['NodeLabs']]
     LM_DATA = COMPUTE_TI_OUTPUT[['LM_DATA']]

     LM_DATA_features = LM_DATA %>% ungroup() %>% dplyr::select(any_of(FEATURES)) %>% as.matrix()
     OOS_DATA_features = OOS_DATA %>% ungroup() %>% dplyr::select(any_of(FEATURES)) %>% as.matrix()

     node.df = COMPUTE_TI_OUTPUT[['node.df']]
     SelPaths <- Tree_e2e[sapply(Tree_e2e, function(x){any(x[c(1, length(x))] == ROOT)})]
     ## Only paths that START or END at ROOT are kept, so an interior ROOT
     ## leaves SelPaths empty and bind_cols() below silently returns a 0x0
     ## data frame. Fail loudly instead, and say what the valid roots are.
     if (length(SelPaths) == 0) {
          .eps <- unique(unlist(lapply(Tree_e2e, function(p) as.numeric(names(p))[c(1, length(p))])))
          stop("ROOT (", ROOT, ") is not an endpoint of any end2end path. ",
               "Valid endpoints are: ", paste(sort(.eps), collapse = ", "))
     }
     SelPaths <- lapply(SelPaths, function(x){
          if (x[1] == ROOT) return(x) else return(rev(x))
     })

     ## Must be the SAME graph that Tree_Graph / Tree_e2e / node.df were built
     ## from in hafez_TI_LINEAR_BRANCH, otherwise the projection is done against
     ## different node positions than the path structure it is indexed by.
     TREE_EPG <- epg_consensus(TreeEPG)
     PartStruct <- ElPiGraph.R::PartitionData(X = OOS_DATA_features, NodePositions = TREE_EPG$NodePositions)
     ProjStruct <- ElPiGraph.R::project_point_onto_graph(
          X = OOS_DATA_features,
          NodePositions = TREE_EPG$NodePositions,
          Edges = TREE_EPG$Edges$Edges,
          Partition = PartStruct$Partition
     )

     AllPt <- lapply(SelPaths, function(x){
          ElPiGraph.R::getPseudotime(ProjStruct = ProjStruct, NodeSeq = names(x))
     })

     myLineageColumns = lapply(AllPt, function(PT){
          PT = PT$Pt %>% dynutils::scale_minmax()
          return(PT)
     }) %>% bind_cols()

     colnames(myLineageColumns) = paste0('path', 1:ncol(myLineageColumns))

     if (RETURN_OBJECTS == TRUE){
          return(list(SelPaths = SelPaths, myLineageColumns = myLineageColumns,
                      PartStruct = PartStruct, ProjStruct = ProjStruct))
     }

     return(myLineageColumns)
}



#' @title hafez_TI
#' @description Perform linear, branching,or cyclical landmark trajectory inference and projects out-of-sample data.
#' @param FULL_DATA Full dataset (dataframe)
#' @param LM_DATA  Landmark dataset to use(dataframe)
#' @param FEATURES vector of features to train landmark trajectory (genes, PCs, etc)
#' @param features_for_start_cell_id Features to use to find start node. If Null, will return object and plot for manual node initialization.
#' @param return_pseudotime_only Will only return pseudotime estimates
#' @param NumNodes Number of nodes to compute. More nodes typically allows for more flexibility
#' @param Lambda  Tuning parameter, typically good to be about 5x to 10x larger than Mu.
#' @param Mu  Tuning parameter.
#' @param nReps Number of times to repeat construction
#' @param ProbPoint Probability of including a point for computation. Value between 0 and 1.
#' @param MaxNumberOfIterations Number of times to include node.
#' @param branch_type Trajectory topology. Options are either 'curve','tree','circle'.
#' @export

hafez_TI = function(FULL_DATA, LM_DATA=NULL, FEATURES, features_for_start_cell_id=NULL, return_pseudotime_only=TRUE,NumNodes=5,Lambda = 0.01, Mu = 0.01, nReps=30, ProbPoint = 1,MaxNumberOfIterations =30, branch_type = c('curve','tree','circle'),return_node_pos = FALSE, use_start_label=NULL, start_label_column_category = c(NULL), verbose =F){
     ## branch_type defaults to the full choice vector; without match.arg() the
     ## `if (branch_type == 'circle')` below sees a length-3 condition, which is
     ## a hard error in R >= 4.2. match.arg() also gives partial matching and a
     ## clear message for an invalid value.
     branch_type <- match.arg(branch_type)
     START_TIME = Sys.time()
     FULL_DATA=FULL_DATA%>% ungroup()
     if (!is.null(LM_DATA)){
          PERFORM_OOS=TRUE
     } else {
          PERFORM_OOS=FALSE
     }
     if (PERFORM_OOS==FALSE){
          FULL_DATA_TRAIN=FULL_DATA

     } else if(PERFORM_OOS==TRUE){

          # LM data can be either a vector  of indexes or a dataframe
          if (is.vector(LM_DATA)) {
               FULL_DATA_TRAIN = FULL_DATA[LM_DATA, , drop = FALSE]
          } else {
               ## if providing a pre-filtered training data
               FULL_DATA_TRAIN=LM_DATA %>% ungroup()
          }
          ## check the RESOLVED training set, not LM_DATA: nrow() of an index
          ## vector is NULL, so `if (nrow(LM_DATA) == 0)` was `if (logical(0))`
          ## -- "argument is of length zero" -- on the documented vector path.
          if (nrow(FULL_DATA_TRAIN) == 0){
               message('no landmarks found...check inputted index vector or dataframe. Returning NA')
               return(NA)
          }
     }

     if (branch_type == 'circle'){
          message('constructing circular graph... ')

          # invisible(capture.output(ELPPI_OUTPUT <- hafez_TI_circle(LM_DATA = FULL_DATA_TRAIN,
          #                                                          FULL_DATA = FULL_DATA,
          #                                                          FEATURES = features,Mu = Mu, Lambda = Lambda, ProbPoint = ProbPoint,
          #                                                          ## changing probPoint
          #                                                          LABELS = NULL,  NumNodes = NumNodes, nReps =nReps, verbose = verbose )
          # ))
          invisible(capture.output(ELPPI_OUTPUT <- hafez_TI_circle(LM_DATA = FULL_DATA_TRAIN,
                                                                   FULL_DATA = FULL_DATA,
                                                                   FEATURES = FEATURES,Mu = Mu, Lambda = Lambda, ProbPoint = ProbPoint,
                                                                   ## changing probPoint
                                                                   LABELS = NULL,  NumNodes = NumNodes, nReps =nReps, verbose = verbose )
          ))


          message('graph complete... ')
          return(ELPPI_OUTPUT)
          proj.data = data.frame(ELPPI_OUTPUT$projection_struture$X_projected)
          colnames(proj.data) = features_LMDPDR[1:ncol(proj.data)]
          node.pos = data.frame(ELPPI_OUTPUT$projection_struture$NodePositions)
          colnames(node.pos) = features_LMDPDR[1:ncol(node.pos)]
          #
          # ggplot(ELPPI_OUTPUT$data,
          #        aes(x = LMDPDR1, y = LMDPDR2, color = gate)) +
          #      geom_point(size = 0.5) +
          #      geom_point(data = node.pos,inherit.aes = F, aes(x=LMDPDR1, y=LMDPDR2),size= 4) +
          #      viridis::scale_color_viridis(option ='magma', discrete = T)

          pst=ELPPI_OUTPUT %>% dplyr::select(cell.id, pseudotime)
          colnames(pst) = c('cell.id','MA_OOS_path1')

          FULL_DATA = FULL_DATA %>%
               # dplyr::select(-any_of(   colnames(pst))) %>%
               left_join(pst)
     } else if(branch_type %in% c('tree','curve')) {

          message('constructing graph... ')
          # invisible(capture.output(ELPIGRAPH_RES <- hafez_TI_LINEAR_BRANCH(LM_DATA = FULL_DATA_TRAIN,
          #                                                                  FULL_DATA = FULL_DATA,
          #                                                                  FEATURES = features,
          #                                                                  CCphase_column = 'gate',## used to label the example plot with node locations
          #                                                                  NumNodes = NumNodes, Lambda = Lambda, Mu = Mu,Do_PCA = F,nReps = nReps,verbose = verbose,
          #                                                                  MaxNumberOfIterations=MaxNumberOfIterations,ProbPoint=ProbPoint,branch_type = branch_type,
          #                                                                  # verbose = T,
          #                                                                  drawAccuracyComplexity = FALSE, drawEnergy = FALSE)
          # ))

          invisible(capture.output(ELPIGRAPH_RES <- hafez_TI_LINEAR_BRANCH(
               LM_DATA = FULL_DATA_TRAIN,
               FULL_DATA = FULL_DATA,
               FEATURES = FEATURES,
               CC_PHASE_COLUMN = 'gate',  # Used for plot labeling
               NumNodes = NumNodes,
               Lambda = Lambda,
               Mu = Mu,
               Do_PCA = FALSE,
               nReps = nReps,
               verbose = verbose,
               MaxNumberOfIterations = MaxNumberOfIterations,
               ProbPoint = ProbPoint,
               BRANCH_TYPE = branch_type,
               DRAW_ACCURACY_COMPLEXITY = FALSE,
               DRAW_ENERGY = FALSE
          )))
          message('graph complete... ')

          # ELPIGRAPH_RES <- hafez_TI_LINEAR_BRANCH(LM_DATA = FULL_DATA_TRAIN,
          #                                                                  FULL_DATA = FULL_DATA,
          #                                                                  FEATURES = features,
          #                                                                  CCphase_column = 'gate',## used to label the example plot with node locations
          #                                                                  NumNodes = NumNodes, Lambda = Lambda, Mu = Mu,Do_PCA = F,nReps = nReps,verbose = verbose,
          #                                                                  MaxNumberOfIterations=MaxNumberOfIterations,ProbPoint=ProbPoint,branch_type = branch_type,
          #                                                                  # verbose = T,
          #                                                                  drawAccuracyComplexity = FALSE, drawEnergy = FALSE)


          colnames(ELPIGRAPH_RES$node.df) = c(FEATURES, 'node','path1')
          node.pos= ELPIGRAPH_RES$node.df


          if (!is.null(use_start_label)){
               my_start_label = FULL_DATA_TRAIN %>% dplyr::filter(!!sym(start_label_column_category[1]) == start_label_column_category[2]) %>% summarize_at(FEATURES, median)
               ELPIGRAPH_RES$node.df = ELPIGRAPH_RES$node.df %>% bind_rows(my_start_label)
          } else {
               ## Pick the root from the cell with the lowest mean over
               ## features_for_start_cell_id. This block used to run
               ## unconditionally ABOVE the use_start_label branch; with the
               ## default features_for_start_cell_id = NULL, any_of(NULL)
               ## selects zero columns and CLOSEST_CELL_ID became character(0),
               ## which then errored inside the filter() below.
               if (is.null(features_for_start_cell_id)) {
                    stop("no way to choose a root: supply `use_start_label` (with ",
                         "`start_label_column_category`) or `features_for_start_cell_id`.")
               }
               missing_f <- setdiff(features_for_start_cell_id, colnames(FULL_DATA_TRAIN))
               if (length(missing_f)) {
                    stop("features_for_start_cell_id not found in the landmark data: ",
                         paste(missing_f, collapse = ", "))
               }
               CLOSEST_CELL_IDX = FULL_DATA_TRAIN %>% ungroup() %>% dplyr::select(any_of(features_for_start_cell_id)) %>% apply(.,1, mean) %>% which.min(.)
               CLOSEST_CELL_ID = FULL_DATA_TRAIN$cell.id[CLOSEST_CELL_IDX]
               ELPIGRAPH_RES$node.df = ELPIGRAPH_RES$node.df %>% bind_rows(FULL_DATA_TRAIN %>% dplyr::filter(cell.id == CLOSEST_CELL_ID) %>% dplyr::select(any_of(FEATURES)))
          }


          ## The ROOT must be an END-POINT of an end2end path:
          ## hafez_lineages_from_root() keeps only paths whose first or last
          ## node equals ROOT, so an interior node selects no paths at all and
          ## silently yields a 0x0 frame. Snap to the nearest ENDPOINT rather
          ## than the nearest node of any kind. The appended start point is the
          ## last row of node.df, hence the final column of the distance matrix.
          endpoint_ids <- unique(unlist(lapply(ELPIGRAPH_RES$Tree_e2e,
                                               function(p) as.numeric(names(p))[c(1, length(p))])))
          endpoint_ids <- endpoint_ids[!is.na(endpoint_ids)]
          if (length(endpoint_ids) == 0) {
               stop("the fitted graph has no end2end paths, so no root endpoint can be chosen")
          }
          .dmat <- ELPIGRAPH_RES$node.df %>% dplyr::select(any_of(FEATURES)) %>%
               dist(method = 'euclidean') %>% as.matrix()
          START_NODE_ID <- endpoint_ids[which.min(.dmat[endpoint_ids, ncol(.dmat)])]


          ## select start node at S-phase
          message('projecting to landmarks...')
          pst = hafez_lineages_from_root(ELPIGRAPH_RES, OOS_DATA = FULL_DATA, ROOT = START_NODE_ID,FEATURES = FEATURES)
          # print(head(pst))
          if (PERFORM_OOS==TRUE){
               colnames(pst) = paste('LM_TI',colnames(pst),sep = '_')

          } else {
               colnames(pst) = paste('TI',colnames(pst),sep = '_')
          }

          if (return_pseudotime_only==TRUE){
               if(return_node_pos==TRUE){
                    message('returning node positions...')
                    return(list(pst = data.frame(pst),
                                node_df = node.pos))
               } else {
                    return(data.frame(pst))

               }

          }


          if (nrow(FULL_DATA)!= nrow(pst)){
               message('mismatched data and pseudotime dataframes')
               return(NULL)
          }
          FULL_DATA = FULL_DATA %>%
               dplyr::select(-any_of(   colnames(pst))) %>%
               bind_cols(pst)
     } else {
          message('branch_type Must be tree, curve, or circle.')
          stopifnot(branch_type %in% c('tree','curve','circle'))
     }

     if(return_node_pos==TRUE){
          return(list(FULL_DATA = FULL_DATA,
                      node_df = node.pos))
     }
     END_TIME = Sys.time()
     message('Hafez TI complete! Total Time: ', END_TIME-START_TIME)
     return(FULL_DATA)
}
