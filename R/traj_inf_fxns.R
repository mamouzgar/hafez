
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
hafez_TI_circle = function(CONTROL_DATA, FULL_DATA, FEATURES = paste0('comp_',1:10), LABELS = NULL,
                           NumNodes = 25, nReps = 5, ProbPoint = 0.6 ,
                           Lambda = 0.01, Mu=0.1, Do_PCA = F,
                           RETURN_PROJECTION = FALSE,
                           MaxNumberOfIterations  = 20,
                           drawAccuracyComplexity = FALSE,drawEnergy = FALSE,  drawPCAView = F, verbose = FALSE){
     TRAINING_DATA = CONTROL_DATA %>% dplyr::select(any_of(FEATURES)) %>% as.matrix()
     OOS_DATA = FULL_DATA %>% dplyr::select(any_of(FEATURES)) %>% as.matrix()





     CircleEPG.Boot_controlOnly <- ElPiGraph.R::computeElasticPrincipalCircle(X = TRAINING_DATA,
                                                                              NumNodes = NumNodes, nReps = nReps, ProbPoint = ProbPoint ,
                                                                              Lambda = Lambda, Mu=Mu, Do_PCA = Do_PCA,
                                                                              RETURN_PROJECTION = FALSE,
                                                                              # ReduceDimension = c(1:8),
                                                                              MaxNumberOfIterations  = MaxNumberOfIterations,
                                                                              drawAccuracyComplexity = drawAccuracyComplexity, drawEnergy = drawEnergy, drawPCAView = drawPCAView,
                                                                              verbose = verbose)

     # PlotPG(X = TRAINING_DATA, TargetPG = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]], Main = "A tree", DimToPlot = 1:2,GroupsLab = train.df.filt$flowsom_phase)
     # PlotPG(X = TRAINING_DATA, TargetPG = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]],
     #        BootPG = CircleEPG.Boot_controlOnly[1:(length(CircleEPG.Boot_controlOnly)-1)],
     #        Main = "A bootstrapped circle", DimToPlot = 1:2)

     PartStruct_controlOnly <- ElPiGraph.R::PartitionData(X = OOS_DATA, NodePositions = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]]$NodePositions)
     ProjStruct_controlOnly <- ElPiGraph.R::project_point_onto_graph(X = OOS_DATA,
                                                                     NodePositions = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]]$NodePositions,
                                                                     Edges = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]]$Edges$Edges,
                                                                     Partition = PartStruct_controlOnly$Partition)
     Circle_Graph_controlOnly <- ElPiGraph.R::ConstructGraph(CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]])

     Circle_e2e_controlOnly <- ElPiGraph.R::GetSubGraph(Net = Circle_Graph_controlOnly, Structure = 'circle', Circular = T)
     Root <- 1
     SelPaths_controlOnly <- Circle_e2e_controlOnly[sapply(Circle_e2e_controlOnly, function(x){any(x[c(1, length(x))] == Root)})]

     AllPt_controlOnly <- lapply(SelPaths_controlOnly, function(x){
          ElPiGraph.R::getPseudotime(ProjStruct = ProjStruct_controlOnly, NodeSeq = names(x))
     })

     PointsPT_controlOnly <- apply(sapply(AllPt_controlOnly, "[[", "Pt"), 1, function(x){unique(x[!is.na(x)])})


     myPseudotime = PointsPT_controlOnly %>% dynutils::scale_minmax(.)
     FULL_DATA$pseudotime = myPseudotime
     if(RETURN_PROJECTION== TRUE){
          # colnames(ProjStruct_controlOnly) = FEATURES
          return(list(projection_struture = ProjStruct_controlOnly, data = FULL_DATA))
     } else {
          return(FULL_DATA)
     }

     # } else if (TOPOLOGY == 'linear') {
     #      LinearEPG.Boot_controlOnly <- ElPiGraph.R::computeElasticPrincipalCurve(X = TRAINING_DATA,
     #                                                                              NumNodes = NumNodes, nReps = nReps, ProbPoint = ProbPoint ,
     #                                                                              Lambda = Lambda, Mu=Mu, Do_PCA = Do_PCA,
     #                                                                              # ReduceDimension = c(1:8),
     #                                                                              MaxNumberOfIterations  = MaxNumberOfIterations,
     #                                                                              drawAccuracyComplexity = drawAccuracyComplexity, drawEnergy = drawEnergy, drawPCAView = drawPCAView,
     #                                                                              verbose = T)
     #      return(LinearEPG.Boot_controlOnly)
     #      # PlotPG(X = TRAINING_DATA, TargetPG = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]], Main = "A tree", DimToPlot = 1:2,GroupsLab = train.df.filt$flowsom_phase)
     #      # PlotPG(X = TRAINING_DATA, TargetPG = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]],
     #      #        BootPG = CircleEPG.Boot_controlOnly[1:(length(CircleEPG.Boot_controlOnly)-1)],
     #      #        Main = "A bootstrapped circle", DimToPlot = 1:2)
     #      #
     #      # PartStruct_controlOnly <- PartitionData(X = OOS_DATA, NodePositions = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]]$NodePositions)
     #      # ProjStruct_controlOnly <- project_point_onto_graph(X = OOS_DATA,
     #      #                                                    NodePositions = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]]$NodePositions,
     #      #                                                    Edges = CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]]$Edges$Edges,
     #      #                                                    Partition = PartStruct_controlOnly$Partition)
     #      # Circle_Graph_controlOnly <- ConstructGraph(CircleEPG.Boot_controlOnly[[length(CircleEPG.Boot_controlOnly)]])
     #      #
     #      # Circle_e2e_controlOnly <- GetSubGraph(Net = Circle_Graph_controlOnly, Structure = 'branches', Circular = T)
     #      # Root <- 1
     #      # SelPaths_controlOnly <- Circle_e2e_controlOnly[sapply(Circle_e2e_controlOnly, function(x){any(x[c(1, length(x))] == Root)})]
     #      #
     #      # AllPt_controlOnly <- lapply(SelPaths_controlOnly, function(x){
     #      #      getPseudotime(ProjStruct = ProjStruct_controlOnly, NodeSeq = names(x))
     #      # })
     #      #
     #      # PointsPT_controlOnly <- apply(sapply(AllPt_controlOnly, "[[", "Pt"), 1, function(x){unique(x[!is.na(x)])})
     #      #
     #      #
     #      # myPseudotime = PointsPT_controlOnly %>% dynutils::scale_minmax(.)
     #      # FULL_DATA$pseudotime = myPseudotime
     #      # return(FULL_DATA)
     # }

}


# library(igraph)
# library(ElPiGraph.R)
## TI_LINEAR_BRANCH
#' @description TI_LINEAR_BRANCH
#' @keywords internal
#' @noRd
hafez_TI_LINEAR_BRANCH = function(CONTROL_DATA, FULL_DATA, FEATURES = paste0('PC',1:10), LABELS = NULL,
                                  FEATURES_PLOT = NULL,
                                  NumNodes = 10, nReps = 5, ProbPoint = 0.6 ,
                                  Lambda = 0.01, Mu=0.1, Do_PCA = F,
                                  MaxNumberOfIterations  = 20, CCphase_column = NULL,
                                  branch_type = 'tree',
                                  drawAccuracyComplexity = FALSE,drawEnergy = FALSE,  drawPCAView = F, verbose = FALSE){
     if (is.null(FEATURES_PLOT) ) {
          FEATURES_PLOT = FEATURES
     }

     TRAIN_DATA =  CONTROL_DATA %>%ungroup() %>%dplyr::select(any_of(FEATURES))


     if (branch_type == 'tree'){
          TreeEPG <- ElPiGraph.R::computeElasticPrincipalTree(X = as.matrix(TRAIN_DATA),
                                                              NumNodes = NumNodes,ProbPoint=ProbPoint, nReps =nReps,
                                                              MaxNumberOfIterations = MaxNumberOfIterations,
                                                              verbose = verbose,
                                                              Lambda = Lambda, Mu = Mu,Do_PCA = F,  drawPCAView = drawPCAView, drawAccuracyComplexity=drawAccuracyComplexity,
                                                              # verbose = verbose,
                                                              drawEnergy = drawEnergy)
     } else if (branch_type == 'curve'){
          TreeEPG = ElPiGraph.R::computeElasticPrincipalCurve(X = as.matrix(TRAIN_DATA),
                                                              NumNodes = NumNodes,ProbPoint=ProbPoint, nReps =nReps,
                                                              MaxNumberOfIterations = MaxNumberOfIterations,
                                                              verbose=verbose,
                                                              Lambda = Lambda, Mu = Mu,Do_PCA = F,  drawPCAView = drawPCAView, drawAccuracyComplexity=drawAccuracyComplexity,
                                                              # verbose = verbose,
                                                              drawEnergy = drawEnergy)
     }

     # PlotPG(X = tree_data_features, TargetPG = TreeEPG[[1]],Do_PCA = F,
     #        NodeLabels = 1:nrow(TreeEPG[[1]]$NodePositions),
     #        LabMult = 2.5, PointSize = NA, p.alpha = .1)
     # PlotPG(X = tree_data_features, TargetPG = TreeEPG[[1]],Do_PCA = F, GroupsLab = tree_data$auto_annotation,
     #        NodeLabels = 1:nrow(TreeEPG[[1]]$NodePositions),
     #        # PlotProjections = 'onEdges',
     #        LabMult = 2.5, PointSize = NA, p.alpha = 0.1)
     nodeLabels = TreeEPG[[1]]$NodePositions
     # dist(nodeLabels)
     # pca.res = prcomp(tree_data_features)
     # pca.df = CONTROL_DATA  %>% data.frame() #%>%
     # bind_cols(CONTROL_DATA %>% dplyr::select(cell.id, auto_annotation))
     #
     # FIND_CLOSEST_NODE = function(similarity.df){
     #      ## most similar pairs of people (hacky way to jokingly do this)
     #      # similarity.dist = mds.df %>% dplyr::select(X, MDS2) %>% dist() %>% as.matrix() %>% data.frame()
     #      similarity.df_rows=similarity.df %>%
     #           mutate(node = paste0('X',1:nrow(.)))
     #      colnames(similarity.df) = paste0('X',1:ncol(similarity.df))
     #      rownames(similarity.df) =  paste0('X',1:nrow(similarity.df))
     #      similarity.dist = dist(similarity.df) %>% as.matrix() %>% data.frame()
     #
     #      similar.person = apply(similarity.dist, MARGIN = 1, function(dists) {dists.filt = dists[dists>0] ; min.pos = which.min(dists.filt) ; names(dists.filt[min.pos])})
     #      similar.person_dist = apply(similarity.dist, MARGIN = 1, function(dists) {dists.filt = dists[dists>0] ; min.pos = which.min(dists.filt) ; dists.filt[min.pos]})
     #
     #      similar.person.df = similar.person %>% data.frame(n1 = names(.), n2 = .) %>%
     #           mutate(distance_to_eachother = similar.person_dist) %>%
     #           # left_join(type_primary %>% dplyr::select(person1=name, type_primary1 = type_primary)) %>%
     #           # left_join(type_primary %>% dplyr::select(person2=name, type_primary2 =type_primary)) %>%
     #           left_join(similarity.df_rows   %>% dplyr::select(n1 = node, X1_n1 = X1, X2_n1 = X2)  ) %>%
     #           left_join(similarity.df_rows   %>% dplyr::select(n2 = node, X1_n2 = X1, X2_n2 = X2)  )
     #
     #      node_mapping = data.frame(you = names(similar.person), next_node = similar.person, distance = round(similar.person_dist,3))
     #
     # }
     node.df = nodeLabels %>% data.frame() #%>% # %*% pca.res$rotation
     # rownames(node.df) = 1:nrow(node.df)
     # colnames(node.df) = 1:ncol(node.df)
     # node.dist = dist(node.df)
     # hclust.res = hclust(node.dist)
     # which(as.matrix(node.dist)==min(node.dist),arr.ind=TRUE)

     node.df = node.df %>%
          mutate(node = paste0(1:nrow(.)),
                 path1 = factor(node, ))# %>%
     print(head(CONTROL_DATA))

     if (is.null(CCphase_column)){
          p.lineages = ggplot(CONTROL_DATA , aes(x =!!sym(FEATURES_PLOT[1]), y = !!sym(FEATURES_PLOT[2]))) +
               theme_minimal() +
               geom_point(aes(fill = auto_annotation), size = 0.75, shape = 21, stroke = 0.1, color= 'black')+
               geom_point(data = node.df , aes(x = X1, y = X2), color = 'red') +
               geom_text(data = node.df , aes(x = X1, y = X2, label = node), color = 'black') +
               guides(fill=guide_legend(override.aes=list( size = 5)))
     } else if(!is.null(CCphase_column)){
          p.lineages = ggplot(CONTROL_DATA , aes(x =!!sym(FEATURES_PLOT[1]), y = !!sym(FEATURES_PLOT[2]))) +
               theme_minimal() +
               geom_point(aes(fill = !!sym(CCphase_column)), size = 0.75, shape = 21, stroke = 0.1, color= 'black')+
               geom_point(data = node.df , aes(x = X1, y = X2), color = 'red') +
               geom_text(data = node.df , aes(x = X1, y = X2, label = node), color = 'black') +
               guides(fill=guide_legend(override.aes=list( size = 5)))
     }

     # library(igraph)
     NodeLabs <- 1:nrow(TreeEPG[[1]]$NodePositions)
     # NodeLabs[degree(ConstructGraph(TreeEPG[[1]])) != 1] <- NA

     # PlotPG(X = tree_data_features, TargetPG = TreeEPG[[1]],Do_PCA = F,
     #        NodeLabels = NodeLabs,GroupsLab = tree_data$auto_annotation,
     #        LabMult = 5, PointSize = NA, p.alpha = .1)

     ### get substructure of interest
     Tree_Graph <- ElPiGraph.R::ConstructGraph(TreeEPG[[1]])
     Tree_e2e <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = 'end2end')

     return(list(TreeEPG=TreeEPG, Tree_e2e=Tree_e2e,Tree_Graph=Tree_Graph,NodeLabs=NodeLabs,CONTROL_DATA=CONTROL_DATA,node.df=node.df,plot = p.lineages ))
     # Root <- 4
     # SelPaths <- Tree_e2e[sapply(Tree_e2e, function(x){any(x[c(1, length(x))] == Root)})]
     #
     # SelPaths <- lapply(SelPaths, function(x){
     #      if(x[1] == Root){
     #           return(x)
     #      } else {
     #           return(rev(x))
     #      }
     # })
     # SelPaths

}



#' @description hafez_lineages_from_root
#' @keywords internal
#' @noRd
hafez_lineages_from_root = function(COMPUTE_G0G1M_TRAJECTORY_output, OOS_DATA, features, root, return_objects =FALSE){
     Tree_e2e =  COMPUTE_G0G1M_TRAJECTORY_output[['Tree_e2e']]
     TreeEPG = COMPUTE_G0G1M_TRAJECTORY_output[['TreeEPG']]
     Tree_Graph = COMPUTE_G0G1M_TRAJECTORY_output[['Tree_Graph']]
     NodeLabs = COMPUTE_G0G1M_TRAJECTORY_output[['NodeLabs']]
     CONTROL_DATA = COMPUTE_G0G1M_TRAJECTORY_output[['CONTROL_DATA']]

     CONTROL_DATA_features = CONTROL_DATA %>% ungroup() %>% dplyr::select(any_of(features)) %>% as.matrix()
     OOS_DATA_features =  OOS_DATA %>% ungroup() %>% dplyr::select(any_of(features)) %>% as.matrix()

     node.df = COMPUTE_G0G1M_TRAJECTORY_output[['node.df']]
     Root <- root
     SelPaths <- Tree_e2e[sapply(Tree_e2e, function(x){any(x[c(1, length(x))] == Root)})]

     SelPaths <- lapply(SelPaths, function(x){
          if(x[1] == Root){
               return(x)
          } else {
               return(rev(x))
          }
     })
     PartStruct <- ElPiGraph.R::PartitionData(X = OOS_DATA_features, NodePositions = TreeEPG[[1]]$NodePositions)

     ProjStruct <- ElPiGraph.R::project_point_onto_graph(X = OOS_DATA_features,
                                                         NodePositions = TreeEPG[[1]]$NodePositions,
                                                         Edges = TreeEPG[[1]]$Edges$Edges,
                                                         Partition = PartStruct$Partition)

     AllPt <- lapply(SelPaths, function(x){
          ElPiGraph.R::getPseudotime(ProjStruct = ProjStruct, NodeSeq = names(x))
     })

     myLineageColumns = lapply(AllPt, function(PT){
          PT = PT$Pt %>% dynutils::scale_minmax()
          return(PT)
     }) %>% bind_cols()

     colnames(myLineageColumns) = paste0('path',1:ncol(myLineageColumns))

     if( return_objects == TRUE){
          return(list(SelPaths=SelPaths,myLineageColumns=myLineageColumns,PartStruct=PartStruct,ProjStruct=ProjStruct))
     }

     return(myLineageColumns)
}



#' @title hafez_TI
#' @description hafez_TI
#' @export
hafez_TI = function(FULL_DATA, features, features_for_start_cell_id, LM_DATA=NULL,return_pseudotime_only=TRUE,NumNodes=5,lambda = 0.01, mu = 0.01, nReps=30, ProbPoint = 0.8,MaxNumberOfIterations =30, branch_type = c('tree','curve','circle'),TRAINING_VARIABLE='WT', return_node_pos = FALSE, use_start_label=NULL, start_label_column_category = c('celltype','HSC'), verbose =F){
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

          invisible(capture.output(ELPPI_OUTPUT <- hafez_TI_circle(CONTROL_DATA = FULL_DATA_TRAIN,
                                                                   FULL_DATA = FULL_DATA,
                                                                   FEATURES = features,Mu = mu, Lambda = lambda, ProbPoint = ProbPoint,
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
          invisible(capture.output(ELPIGRAPH_RES <- hafez_TI_LINEAR_BRANCH(CONTROL_DATA = FULL_DATA_TRAIN,
                                                                           FULL_DATA = FULL_DATA,
                                                                           FEATURES = features,
                                                                           CCphase_column = 'gate',## used to label the example plot with node locations
                                                                           NumNodes = NumNodes, Lambda = lambda, Mu = mu,Do_PCA = F,nReps = nReps,verbose = verbose,
                                                                           MaxNumberOfIterations=MaxNumberOfIterations,ProbPoint=ProbPoint,branch_type = branch_type,
                                                                           # verbose = T,
                                                                           drawAccuracyComplexity = FALSE, drawEnergy = FALSE)
          ))
          message('graph complete... ')

          # ELPIGRAPH_RES <- hafez_TI_LINEAR_BRANCH(CONTROL_DATA = FULL_DATA_TRAIN,
          #                                                                  FULL_DATA = FULL_DATA,
          #                                                                  FEATURES = features,
          #                                                                  CCphase_column = 'gate',## used to label the example plot with node locations
          #                                                                  NumNodes = NumNodes, Lambda = lambda, Mu = mu,Do_PCA = F,nReps = nReps,verbose = verbose,
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
          pst = hafez_lineages_from_root(ELPIGRAPH_RES, OOS_DATA = FULL_DATA, root = START_NODE_ID,features = features)
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
          message('branch_type must be tree, curve, or circle.')
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
