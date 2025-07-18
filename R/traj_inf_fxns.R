
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

     PartStruct <- ElPiGraph.R::PartitionData(X = OOS_DATA, NodePositions = CircleEPG.Boot[[length(CircleEPG.Boot)]]$NodePositions)
     ProjStruct <- ElPiGraph.R::project_point_onto_graph(
          X = OOS_DATA,
          NodePositions = CircleEPG.Boot[[length(CircleEPG.Boot)]]$NodePositions,
          Edges = CircleEPG.Boot[[length(CircleEPG.Boot)]]$Edges$Edges,
          Partition = PartStruct$Partition
     )
     Circle_Graph <- ElPiGraph.R::ConstructGraph(CircleEPG.Boot[[length(CircleEPG.Boot)]])
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
                                  MaxNumberOfIterations = 20, CC_PHASE_COLUMN = NULL,
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

     node.df = TreeEPG[[1]]$NodePositions %>% data.frame() %>%
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

     Tree_Graph <- ElPiGraph.R::ConstructGraph(TreeEPG[[1]])
     Tree_e2e <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = 'end2end')
     NodeLabs <- 1:nrow(TreeEPG[[1]]$NodePositions)

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
     SelPaths <- lapply(SelPaths, function(x){
          if (x[1] == ROOT) return(x) else return(rev(x))
     })

     PartStruct <- ElPiGraph.R::PartitionData(X = OOS_DATA_features, NodePositions = TreeEPG[[1]]$NodePositions)
     ProjStruct <- ElPiGraph.R::project_point_onto_graph(
          X = OOS_DATA_features,
          NodePositions = TreeEPG[[1]]$NodePositions,
          Edges = TreeEPG[[1]]$Edges$Edges,
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
#' @param features vector of features to train landmark trajectory (genes, PCs, etc)
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

hafez_TI = function(FULL_DATA, LM_DATA=NULL, features, features_for_start_cell_id=NULL, return_pseudotime_only=TRUE,NumNodes=5,Lambda = 0.01, Mu = 0.01, nReps=30, ProbPoint = 1,MaxNumberOfIterations =30, branch_type = c('curve','tree','circle'),return_node_pos = FALSE, use_start_label=NULL, start_label_column_category = c(NULL), verbose =F){
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
               FULL_DATA_TRAIN = FULL_DATA[LM_DATA, ]
          } else {
               ## if providing a pre-filtered training data
               FULL_DATA_TRAIN=LM_DATA %>% ungroup()
          }
          if (nrow(LM_DATA) == 0){
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
                                                                   FEATURES = features,Mu = Mu, Lambda = Lambda, ProbPoint = ProbPoint,
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
               FEATURES = features,
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


          ## closest node
          CLOSEST_CELL_IDX = FULL_DATA_TRAIN %>% ungroup() %>%dplyr::select(any_of(features_for_start_cell_id)) %>% apply(.,1, mean) %>% which.min(.)
          CLOSEST_CELL_ID = FULL_DATA_TRAIN$cell.id[CLOSEST_CELL_IDX]

          colnames(ELPIGRAPH_RES$node.df) = c(features, 'node','path1')
          node.pos= ELPIGRAPH_RES$node.df


          if (!is.null(use_start_label)){
               my_start_label = FULL_DATA_TRAIN %>% dplyr::filter(!!sym(start_label_column_category[1]) == start_label_column_category[2]) %>% summarize_at(features, median)
               ELPIGRAPH_RES$node.df = ELPIGRAPH_RES$node.df %>% bind_rows(my_start_label)
          } else {
               ELPIGRAPH_RES$node.df = ELPIGRAPH_RES$node.df %>% bind_rows(FULL_DATA_TRAIN %>% dplyr::filter(cell.id == CLOSEST_CELL_ID) %>% dplyr::select(any_of(features)))
          }


          START_NODE_ID = dist(ELPIGRAPH_RES$node.df %>% dplyr::select(any_of(features)),method = 'euclidean') %>% as.matrix() %>% .[1:nrow(.),ncol(.)] %>% .[.!=0] %>% .[which.min(.)] %>% names() %>% as.numeric()


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
