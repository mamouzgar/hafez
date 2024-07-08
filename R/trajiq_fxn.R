##' Authors: Meelad Amouzgar
##' @author Meelad Amouzgar
##' date: June 13th, 2024
##' Ctrl+Shift+D namespace
##' usethis::use_pipe()
##' devtools::document()
##'
##'
##'
##'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr rename
#' @importFrom dplyr summarize
#' @importFrom dplyr n
#' @importFrom dplyr %>%
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr rowwise
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr left_join
#' @importFrom dplyr right_join
#' @importFrom dplyr full_join
#' @importFrom dplyr case_when
#' @importFrom dplyr sample_n
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom tidyr separate
#' @importFrom stats setNames
#' @importFrom stats glm
#' @importFrom utils combn
#' @importFrom stats p.adjust
#' @importFrom multcomp glht
#' @importFrom emmeans emmeans
#' @import ElPiGraph.R
#' @import sqldf
#' @import igraph

#' @title downsampleWith_group_by
#' @description downsampleWith_group_by will randomly sample from each group in your dataset that you've specified using group_by. If the # of observations is less than the desired # to sample, then it will sample all the cells in that group.
#' @param df dataset
#' @param sampleN the max # of observations (eg,cells) to sample per group. If
#' @export
downsampleWith_group_by = function(df, sampleN){
     if (nrow(df) >= sampleN){
          sample.ind = sample(1:nrow(df), sampleN)
          return(tibble::tibble(df[sample.ind, ]))
     } else {
          return(tibble::tibble(df))
     }
}
# DATA_SET = DATA_SET %>%
#      group_by(patient_id, annotation_lin_sub, annotation_subset,tissue )  %>%
#      group_map(~sample_value_or_minimum(., 200),.keep = T) %>% bind_rows()

## constructs simple glm
#' @description construct_model
#' @keywords internal
#' @noRd
construct_model = function(glm_input, outcome_features, contrast_variables=NULL, covariates_in_model,polynomial_covariate, intercept = FALSE){
     glm_input = as_tibble(glm_input)
     # glm_input$group = paste(glm_input[[]])
     glm_input$group = apply(glm_input[,contrast_variables], 1, paste, collapse="_-_") %>% factor(.)

     message('constructing models...')
     myGLM = lapply(outcome_features, function(ro){
          # message(ro)
          # myFormula = paste0(ro , '~0 + treatment*DivisionID', '+poly(PSEUDOTIME_NORMALIZED, degree =5 ) + timepoint +donor ') %>% as.formula()
          # if(is.null(group_name)){
          #      group_name='group'
          # }
          if(intercept == TRUE){
               intercept_choice = paste0('~0+', 'group')
          } else {
               intercept_choice = paste0('~', 'group')
          }

          if (!is.null(covariates_in_model)){

               myFormula = paste0(ro , intercept_choice ,'+', paste0( covariates_in_model,collapse = '+')) %>% as.formula()

               if (!is.null(polynomial_covariate)){
                    myFormula = paste0(ro , intercept_choice ,'+', paste0( covariates_in_model,collapse = '+'), '+', paste0('poly(',polynomial_covariate,', degree = 6)')) %>% as.formula()

               }

          } else {
               if (!is.null(polynomial_covariate)){
                    myFormula = paste0(ro , intercept_choice,'+', paste0('poly(',polynomial_covariate,', degree = 6)')) %>% as.formula()

               } else {
                    myFormula = paste0(ro , intercept_choice) %>% as.formula()
               }
          }

          glm_res =  tryCatch(
               {
                    glm_res =glm(formula = myFormula, data = glm_input,family = gaussian)
               },
               error = function(cond) {
                    # warning(cond)
                    message('failed to construct model...', myFormula)
                    glm_res=NA
                    return(glm_res)
               }
          )
          # print(glm_res)
          return(glm_res)
     })
     names(myGLM) = outcome_features
     # TEST <<-myGLM

     return(list(myGLM=myGLM, group =  glm_input$group))
}


#' @description compute_simple_coef
#' @keywords internal
#' @noRd
compute_simple_coef = function(myGLM){
     myCoef =  lapply(names(myGLM), function(glm_res_name){
          glm_res = myGLM[[glm_res_name]]
          glm_df = data.frame(summary(glm_res)$coef, check.rows = F, check.names = F)%>%
               # bind_cols(data.frame(anova(glm_res), check.rows = F, check.names = F)) %>%
               mutate(coef = rownames(.))%>%
               mutate(feature = glm_res_name)
          return(glm_df)
     }) %>%
          bind_rows() %>%
          dplyr::filter(!grepl('Intercept', coef))
     return(myCoef)
}


#' @description example_combos
#' @keywords internal
#' @noRd
example_combos = function(){

     ## createa all compbinations of interest
     myCombos = unique(c('A','B','C'))  #%>% .[!grepl('donor|PSEUDOTIME',.)]
     myCombos = setNames(as.data.frame(t(combn(myCombos, 2))), c("group1", "group2"))
     myCombos_df =  myCombos %>%
          # unique(names(myGLM$CD45RO$coefficients)) %>% .[!grepl('donor|PSEUDOTIME',.)] %>%
          # expand.grid(group1=., group2=.) %>%
          # separate(col = 'group1',sep = '_-_',into = c('treatment1','division1','timepoint1'), remove = F)  %>%
          # separate(col = 'group2',sep = '_-_',into = c('treatment2','division2', 'timepoint2'), remove = F)  %>%
          # mutate(treatment1 = gsub('treatment_DivisionID_timepoint','',treatment1),
          #        treatment2 = gsub('treatment_DivisionID_timepoint','',treatment2)) %>%
          # dplyr::filter((treatment1==treatment2 & division1==division2 & timepoint1 != timepoint2 ) | ## keeps timepoint comparisons
          #                    (treatment1 != treatment2 & division1==division2 & timepoint1 == timepoint2) | ## keeps  treatment comparisons
          #                    (treatment1 ==treatment2 & division1!=division2 & timepoint1 == timepoint2)) %>%  ## keeps division comparisons
          mutate(COMPARISON_DO_NOT_GENERATE_YOURSELF = paste0(paste(paste0("`",group1,"`"), paste0("`",group2,"`"), sep = ' - '), '== 0'))
}
# omg=example_combos()
# omg

## constructs tukey wukey combos.
#' @description create_TukeyWukey_combos
#' @keywords internal
#' @noRd
create_TukeyWukey_combos = function(unique_group_values){
     message('constructing TukeyWukey combos...')

     myCombos = unique(unique_group_values) #%>% .[!grepl('donor|PSEUDOTIME',.)]
     myCombos = setNames(as.data.frame(t(combn(myCombos, 2))), c("group1", "group2"))

     group_df=myCombos  %>%
          mutate(COMPARISON = paste0(paste(paste0("`",'group',group1,"`"), paste0("`",'group',group2,"`"), sep = ' - '), '== 0'))
     return(group_df)
}

## runs tukey analysis
#' @description TukeyWukey
#' @keywords internal
#' @noRd
TukeyWukey = function(myGLM, combos_df, contrast_variables){
     message('TukeyWukey time...')
     multitest_pairstest = lapply(names(myGLM), function(glm_res_name){
          # message(glm_res_name)
          glm_res = myGLM[[glm_res_name]]
          if(any(is.na(glm_res))){
               return(data.frame(estimate = NA,
                                 pval =NA,
                                 comparison = NA,
                                 feature =glm_res_name ) )
          }
          # glm_df = data.frame(summary(glm_res)$coef, check.rows = F, check.names = F)%>%
          #      # bind_cols(data.frame(anova(glm_res), check.rows = F, check.names = F)) %>%
          #      mutate(coef = rownames(.))%>%
          #      mutate(feature = glm_res_name)
          # multitest_res <- glht(glm_res, linfct = mcp(treatment_timepoint = "Tukey"))
          # multitest_res <- multcomp::glht(glm_res, linfct = combos_df$COMPARISON)

          multitest_res =  tryCatch(
               {
                    multitest_res <- multcomp::glht(glm_res, linfct = combos_df$COMPARISON)
               },
               error = function(cond) {
                    # warning(cond)
                    message('failed to compute...', 'contrast: ',  contrast_variables,' - feature: ',glm_res_name,'...skipping')
                    multitest_res=NA
               }
          )
          if(any(is.na(multitest_res))) {
               return( data.frame(estimate = NA,
                                  pval =NA,
                                  comparison = NA,
                                  feature =glm_res_name ) )
          } else {

               multitest_res_summary = summary(multitest_res)

               multitest_res_summary_df = data.frame(estimate = multitest_res_summary$test$coefficients,
                                                     pval = multitest_res_summary$test$pvalues,
                                                     comparison = names(multitest_res_summary$test$coefficients)) %>%
                    mutate(feature = glm_res_name)
          }
          return(multitest_res_summary_df)
     }) %>% bind_rows() %>%
          mutate(comparison =gsub('`|group','',comparison)) %>%
          separate(col = 'comparison',sep = ' - ',into = c('group1','group2'), remove = F)  %>%
          separate(col = 'group1',sep = '_-_',into = paste0(contrast_variables,'1'), remove = F)  %>%
          separate(col = 'group2',sep = '_-_',into = paste0(contrast_variables,'2'), remove = F)  #%>%
     comparison_contrasts = sapply(contrast_variables, function(dd){
          multitest_pairstest[[paste0(dd,'_comparison')]] = paste(multitest_pairstest[[paste0(dd,'1')]], multitest_pairstest[[paste0(dd,'2')]], sep = '-vs-')
     })  %>% data.frame()%>%
          dplyr::rename_all(~paste0(.x,'_comparison'))

     multitest_pairstest=multitest_pairstest  %>%
          bind_cols(comparison_contrasts) %>%
          mutate(padj = p.adjust(pval, method  = 'BH'),
                 minus_log10padj  = -1*log10(padj))
     return(multitest_pairstest)
}


## runs estimated marginal means analysis
#' @description TukeyWukey
#' @keywords internal
#' @noRd
emmy = function(myGLM, contrast_variables){
     multitest_pairstest = lapply(names(myGLM), function(glm_res_name){
          glm_res = myGLM[[glm_res_name]]

          if(any(is.na(glm_res))) {
               # print(glm_res)
               message('failed to compute...', 'contrast: ',  contrast_variables,' - feature: ',glm_res_name,'...skipping')
               emm_df = data.frame(estimate = NA,
                                   pval =NA,
                                   comparison = NA,
                                   feature =glm_res_name )
          } else {

               emm_df = pairs(emmeans::emmeans(glm_res, "group")) %>%
                    data.frame() %>%
                    mutate(feature = glm_res_name) %>%
                    dplyr::select(estimate, pval = p.value, comparison = contrast, feature)
               # print('returning emmdf')
          }
          return(emm_df)
     }) %>% bind_rows() %>%
          mutate(comparison =gsub('`|group','',comparison)) %>%
          mutate(comparison =gsub('\\(|\\)','',comparison)) %>%
          separate(col = 'comparison',sep = ' - ',into = c('group1','group2'), remove = F)  %>%
          separate(col = 'group1',sep = '_-_',into = paste0(contrast_variables,'1'), remove = F)  %>%
          separate(col = 'group2',sep = '_-_',into = paste0(contrast_variables,'2'), remove = F)  #%>%
     comparison_contrasts = sapply(contrast_variables, function(dd){
          multitest_pairstest[[paste0(dd,'_comparison')]] = paste(multitest_pairstest[[paste0(dd,'1')]], multitest_pairstest[[paste0(dd,'2')]], sep = '-vs-')
     })  %>% data.frame()%>%
          dplyr::rename_all(~paste0(.x,'_comparison'))

     multitest_pairstest=multitest_pairstest  %>%
          bind_cols(comparison_contrasts) %>%
          mutate(padj = p.adjust(pval, method  = 'BH'),
                 minus_log10padj  = -1*log10(padj))
     return(multitest_pairstest)
}

## wrapper for everything
# ##' glm input: a dataframe of rows = cells or samples, and columns as your metadata and features of interest. These features should be outcome variables. Make sure column names are in R-friendly format.
# ##' FYI: + and - is not handled well i R. use gsub('\\+','pos',column_name) and gsub('\\+','neg',column_name)  to fix names as R-friendly before using make.names().
# ##' outcome_features: a vector of your outcome features - must have R-friendly format (use make.names() function). eg, c('B_cell_activated')
# ##' contrast_variables: a vector of your contrast variables of interest. You likely just have just 1 value here,  unless you want to study interaction effects. eg, c('tissue')
# ##' covariates_in_model: a vector of your covariates that you want to control for. Also used for a matched analysis. eg, c("mouse_id"),
# ##' intercept: whether you want to set the intercept to 0 or not. Default is TRUE. Leave as-is unless you know to remove intercept.
## runs tukey analysis
#' @description differential_analysis_program_simple
#' @keywords internal
#' @noRd
differential_analysis_program_simple = function(glm_input, outcome_features,contrast_variables,covariates_in_model, polynomial_covariate , intercept=TRUE,contrast_method = 'emm'){
     glm_res = construct_model(glm_input = glm_input,outcome_features = outcome_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model, polynomial_covariate=polynomial_covariate, intercept = intercept)
     # myCoef = compute_simple_coef(glm_res$myGLM)
     glm_input$group = glm_res$group
     my_group = unique(glm_input$group)

     if (contrast_method=='emm'){
          DIF_RES = emmy(glm_res$myGLM, contrast_variables)
          message('emmy complete!')
     } else if (contrast_method == 'tukey') {

          if(is.null(levels(my_group))){
               my_tukeywukey_combos = create_TukeyWukey_combos(my_group)

          } else {
               my_tukeywukey_combos = create_TukeyWukey_combos(levels(my_group))
          }
          DIF_RES = TukeyWukey(myGLM = glm_res$myGLM,combos_df =  my_tukeywukey_combos,contrast_variables = contrast_variables)
          message('TukeyWukey complete!')
     } else {
          message('did not provide valid contrast method of: emm or tukey')
     }

     return(DIF_RES)
}

#' @title differential_analysis_program
#' @description wrapper function for differential analysis
#' @param glm_input a dataframe of rows = cells or samples, and columns as your metadata and features of interest. These features should be outcome variables. Make sure column names are in R-friendly format. Be aware that + and - is not handled well in R. Clean your data by using gsub('\\+','pos',column_name) and gsub('\\+','neg',column_name) before using make.names(). Generally best to avoid any empty spaces.
#' @param outcome_features a vector of your outcome features - must have R-friendly format (use make.names() function). eg, c('B_cell_activated')
#' @param contrast_variables a vector of your contrast variables of interest. You likely just have just 1 value here,  unless you want to study interaction effects. eg, c('tissue')
#' @param covariates_in_model a vector of your covariates that you want to control for. Also used for a matched analysis. eg, c("mouse_id"),
#' @param intercept whether you want to set the intercept to 0 or not. Default is TRUE. Leave as-is unless you know to remove intercept.
#' @param SPLIT_BY_NAMES A vector of groups you want to restrict each analysis within. Use this if you want to subset the analyis within distinct groups (eg, within celltype), you can use this arguement to split the analysis. Eg ,c("celltype","tissue").
#' @param contrast_method Either estimated marginal means or tukey ('emm','tukey'). defaulst to emm
#' @export
differential_analysis_program = function(glm_input, outcome_features,contrast_variables,covariates_in_model=NULL,polynomial_covariate=NULL, intercept=TRUE, contrast_method = 'emm', SPLIT_BY_NAMES=NULL){
     if (!is.null(SPLIT_BY_NAMES)){
          DIF_RES =  glm_input %>%
               group_by_at(SPLIT_BY_NAMES)%>%
               group_modify(~differential_analysis_program_simple(glm_input = .,outcome_features = outcome_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model, polynomial_covariate=polynomial_covariate, intercept = intercept,contrast_method = contrast_method) ,
                            .keep = T
               )
     } else {
          DIF_RES = differential_analysis_program_simple(glm_input = glm_input,outcome_features = outcome_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,polynomial_covariate=polynomial_covariate,intercept = intercept,contrast_method = contrast_method)
     }
     return(DIF_RES)
}
# ##' glm input: a dataframe of rows = cells or samples, and columns as your metadata and features of interest. These features should be outcome variables. Make sure column names are in R-friendly format.
# ##' FYI: + and - is not handled well i R. use gsub('\\+','pos',column_name) and gsub('\\+','neg',column_name)  to fix names as R-friendly before using make.names().
# ##' outcome_features: a vector of your outcome features - must have R-friendly format (use make.names() function). eg, c('B_cell_activated')
# ##' contrast_variables: a vector of your contrast variables of interest. You likely just have just 1 value here,  unless you want to study interaction effects. eg, c('tissue')
# ##' covariates_in_model: a vector of your covariates that you want to control for. Also used for a matched analysis. eg, c("mouse_id"),
# ##' intercept: whether you want to set the intercept to 0 or not. Default is TRUE. Leave as-is unless you know to remove intercept.
# dif_res =differential_analysis_program(glm_input = glm_input,outcome_features = outcome_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE)

# differential_analysis_program but broken down with example code
# glm_res = construct_model(glm_input = glm_input,outcome_features = outcome_features, contrast_variables = c('tissue'), covariates_in_model = covariates_in_model,intercept = TRUE)
# myCoef = compute_simple_coef(glm_res$myGLM)
# glm_input$group = glm_res$group
# my_group = unique(glm_input$group)
#
# ## filter  to only include contrasts of interest (comparisons)
# my_tukeywukey_combos = create_TukeyWukey_combos(levels(my_group))
#
# ## filter  to only include contrasts of interest (comparisons)
# tukey_res = TukeyWukey(myGLM = glm_res$myGLM,combos_df =  my_tukeywukey_combos,contrast_variables = contrast_variables)
#
#

