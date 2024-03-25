#' ggplot2 extension for a MF1_single or MF2_multiple object
#'
#' \code{MFggplot} provides graphical BEF relationships based on the output from the function \code{MF1_single} or \code{MF2_multiple}.
#'
#' @param output the output obtained from \code{MF1_single} or \code{MF2_multiple}. \cr
#' For output obtained from \code{MF1_single}, if BEF relationships are desired within each category specified \code{by_group}, the \code{by_group} column must be included in the input.
#' @param model specifying the fitting model, \code{model = "lm"} for linear model; \code{model = "LMM.intercept"},
#' \code{"LMM.slope"} and \code{"LMM.both"} for linear mixed models with random effects for intercepts, slopes, and both, respectively.
#' Default is \code{model = "LMM.both"}.
#' @param by_group the column name of the stratifying variable that is used to group data for model fitting. For example, if \code{by_group = "country"}, then model will be fitted within each country. Default is \code{NULL}. \cr
#' It is required if a linear mixed model is selected in the \code{model}. \cr
#' If \code{output} is obtained from \code{MF2_multiple}, the \code{by_group} setting must be the same as that set in \code{MF2_multiple}. 
#' @param caption caption that will be shown in the BEF plots; \code{caption = "slope"} to show the estimated slopes in each plot,
#' or \code{caption = "R.squared"} to show the ordinary R-squared for linear models or estimated marginal and conditional R-squared for linear mixed models in each plot.
#' Default is \code{caption = "slope"}.
#'
#' @import tidyverse
#' @import dplyr
#' @import broom
#' @import patchwork
#' @import ggpubr
#' @import grid
#' @import purrr
#' @importFrom dplyr %>%
#' @importFrom lme4 fixef VarCorr
#' @importFrom grDevices hcl
#'
#' @return For an \code{MF1_single} object, this function returns a figure that plots the BEF relationship between multifunctionality of
#' order q (= 0, 1 and 2) and species diversity of the same order q for two cases (i) correlations between functions are not corrected for, and (ii) correlations between
#' functions are corrected. The fitted lines for the chosen model are also shown in the figure. \cr
#' 
#' For an \code{MF2_multiple} object of given individual function weights, this function returns a list of two objects (\code{$corr_uncorrected} and \code{$corr_corrected})
#' that respectively for two cases: (i) correlations between functions are not corrected for, and (ii) correlations between functions are corrected for. \cr
#' 
#' Each object consists of four figures: "\code{$ALL}" returns a figure that depicts the BEF relationship between alpha/beta/gamma multifunctionality of
#' order q (= 0, 1 and 2) and the corresponding species diversity of the same order q. The fitted lines for the chosen model are also shown in the figure.
#' "\code{$Gamma}" returns only the gamma part of "\code{$ALL}", "\code{$Alpha}" returns only the alpha part of "\code{$ALL}", and "\code{$Beta}" returns only the beta part of "\code{$ALL}". 
#' 
#' @examples
#' 
#' library(dplyr)
#' 
#' ### Use data from six countries
#' 
#' \donttest{
#' ## single ecosystem
#' data("forest_function_data_normalized")
#' data("forest_biodiversity_data")
#' output1 <- MF1_single(func_data = forest_function_data_normalized[,6:31], weight = 1,
#'                       species_data = forest_biodiversity_data)
#'   
#' ## Display fitted line of linear mixed model with random slopes and random intercepts  
#' output1 <- data.frame(output1, country=rep(forest_function_data_normalized$country, each = 6))
#' MFggplot(output1, model = "LMM.both", by_group="country", caption = "slope")
#' 
#' 
#' ### Use data from five countries (data in Finland are excluded)
#'  
#' ## multiple ecosystems
#' data("forest_function_data_normalized")
#' data("forest_biodiversity_data")
#' forest_function_data_normalized <- filter(forest_function_data_normalized, country != "FIN")
#' forest_biodiversity_data <- forest_biodiversity_data[-(1:48),]
#' output2 <- MF2_multiple(func_data = forest_function_data_normalized[,6:32],
#'                         species_data = forest_biodiversity_data,
#'                         weight = 1,
#'                         by_group = "country")
#' 
#' ## Display fitted line of linear mixed model with random slopes and random intercepts  
#' figure_LMM <- MFggplot(output2, model = "LMM.both", by_group = "country", 
#'                        caption = "slope")
#' figure_LMM$corr_uncorrected$ALL
#' figure_LMM$corr_corrected$ALL
#' }
#' 
#' ### Use partial data to quickly obtain output  
#' ### (Take the first 18 plots in Germany and the last 18 plots in Italy)
#' 
#' ## single ecosystem
#' data("forest_function_data_raw")
#' data("forest_biodiversity_data")
#' GER_ITA_forest_function_raw <- filter(forest_function_data_raw, 
#'                                       country=="GER"|country=="ITA")[c(1:18,57:74),]
#' GER_ITA_forest_function_normalized <- function_normalization(data = GER_ITA_forest_function_raw,
#'                                                              fun_cols = 6:31, 
#'                                                              negative = c("soil_cn_ff_10","wue"),
#'                                                              by_group = "country")
#' GER_ITA_forest_biodiversity <- forest_biodiversity_data[c(49:82,181:229),]
#' output3 <- MF1_single(func_data = GER_ITA_forest_function_normalized[,6:31], weight = 1,
#'                       species_data = GER_ITA_forest_biodiversity)
#' 
#' 
#' ## Display fitted line of linear mixed model with random slopes and random intercepts  
#' output3 <- data.frame(output3, country=rep(GER_ITA_forest_function_normalized$country, each = 6))
#' MFggplot(output3, model = "LMM.both", by_group="country", caption = "slope")
#' 
#' \donttest{
#' ## multiple ecosystems
#' data("forest_function_data_raw")
#' data("forest_biodiversity_data")
#' GER_ITA_forest_function_raw <- filter(forest_function_data_raw, 
#'                                       country=="GER"|country=="ITA")[c(1:18,57:74),]
#' GER_ITA_forest_function_normalized <- function_normalization(data = GER_ITA_forest_function_raw,
#'                                                              fun_cols = 6:31, 
#'                                                              negative = c("soil_cn_ff_10","wue"),
#'                                                              by_group = "country")
#' GER_ITA_forest_biodiversity <- forest_biodiversity_data[c(49:82,181:229),]
#' output4 <- MF2_multiple(func_data = GER_ITA_forest_function_normalized[,6:32],
#'                         species_data = GER_ITA_forest_biodiversity,
#'                         weight = 1,
#'                         by_group = "country")
#' 
#' 
#' ## Display fitted line of linear mixed model with random slopes and random intercepts  
#' figure_LMM_GER_ITA <- MFggplot(output4, model = "LMM.both", by_group = "country", 
#'                                caption = "slope")
#' figure_LMM_GER_ITA$corr_uncorrected$ALL
#' figure_LMM_GER_ITA$corr_corrected$ALL
#' }
#' 
#' @export

MFggplot <- function(output, model = "LMM.both", caption = "slope", by_group = NULL){
  
  facets_scale = 'fixed'
  fit<-model
  text<-caption
  
  if(!(fit %in% c("lm","LMM.intercept","LMM.slope","LMM.both")))
    stop("Error: the argument fit should be `lm`, `LMM.intercept`, `LMM.slope` or `LMM.both`.")
  else if(fit != "lm" & is.null(by_group))
    stop("Error: linear mixed model must contain the argument by_group.")
  else if(fit != "lm")
    rand_eff <- strsplit(fit,split = "[.]") %>% unlist() %>% .[2]
  
  if(!is.null(by_group)){
    if(length(by_group)!=1) stop("Error: The number of the group variable should not be more than 1.")
    else if(!(by_group %in% names(output))) stop("Error: The group variable is not included in the given output data.")
    
  }
  
  if(!(text %in% c("slope","R.squared")))
    stop("Error: the argument text should be `slope` or `R.squared`.")
  
  # stdPalette <- c("blue", "#00AAAA", "darkorange","gray55","#FF88C2", "purple2")
  
  # Check if the number of unique 'Country' is 6 or less
  if (length(unique(output$country)) <= 6){
    stdPalette <- c("blue", "#00AAAA", "darkorange","gray55","#FF88C2", "purple2")
  }else{
    # If there are more than 6 countries, start with the same predefined color palette
    # Then extend the palette by generating additional colors using the 'ggplotColors' function
    stdPalette <- c("blue", "#00AAAA", "darkorange","gray55","#FF88C2", "purple2")
    stdPalette <- c(stdPalette, ggplotColors(length(unique(output$country))-6))
  }
  abc<-scale_linetype_manual(values = c( "Significant slope (P < 0.05)" = "solid","Nonsignificant slope" = "dashed"), name = NULL, drop = FALSE)
  
  
  #"steelblue1"
  #cbPalette <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
  
  if(!all(c("Type","Order.q","qMF","Species.diversity") %in% names(output)))
    stop("Error: For ggMF function, you should offer the MF_single output data including at least the columns `Type`, `Order.q`, `qMF` and `Species.diversity`, \n
            or the MF_multiple output data including at least the columns `Type`, `Scale`, `Order.q`, `qMF` and `Species.diversity`.")
  
  # For MF_single output
  else if(!("Scale" %in% names(output))){
    output$Type <- factor(output$Type, levels = c("corr_uncorrected", "corr_corrected"))
    
    if(is.null(by_group)){
      ##########################
      ##########11/2 revised
      lm_data <- output %>% group_by(Type, Order.q) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .)))
      papa =  output %>% group_by(Type,Order.q) %>%do(broom::glance(lm(qMF ~ Species.diversity, .)))%>%select(r.squared,Order.q,Type)
      lm_data <- suppressMessages(lm_data %>% dplyr::full_join(papa))
      lm_data <- mutate(lm_data,
                        Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Nonsignificant slope"),
                                            levels = c("Significant slope (P < 0.05)", "Nonsignificant slope")),
                        group="Linear model"
      ) %>%
        filter(term=="Species.diversity") %>% dplyr::select(-c(term, std.error, statistic, p.value))
      
      output <- suppressMessages(output %>% dplyr::left_join(lm_data))
      
      ###################################################
      ########################11/2 revised
      
      if(text=="slope"){
        plot_output <- ggplot(data = output, aes(x = Species.diversity, y = qMF))+
          facet_grid(Type ~ Order.q, scales = facets_scale) +
          geom_point(size=0.7)+
          geom_smooth(aes(lty = Significance, col=group), method = "lm", se = F, size=1.9, formula = y ~ x)+
          geom_text(data = lm_data, aes(x = -Inf, y = Inf, label=paste0("Slope = ", round(estimate, 4)), col=group),
                    hjust= -0.1, vjust= 2, size=3, key_glyph = draw_key_path)+
          scale_color_manual(values = "red")+
          theme_bw() +
          theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,-6,0,-6),legend.text = element_text(size=12, margin = margin(r = 1, unit = 'cm')),
                strip.text = element_text(size=12),axis.text = element_text(size=8),axis.title=element_text(size=18),
                legend.title = element_text(size=12))+
          guides(
            linetype = guide_legend(override.aes = list(col = "#000000",size=0.7,linewidth = 0.7))
          )+
          labs(x = "Species diversity", y = "Multifunctionality")+abc
      }
      else{
        
        lm_text <- lm_data %>% select(Type,Order.q,r.squared) %>%
          mutate(Label=paste("R^2==",round(r.squared,3),sep=" "),
                 h=-0.1,
                 v=2)
        plot_output <- ggplot(data = output, aes(x = Species.diversity, y = qMF))+
          facet_grid(Type ~ Order.q, scales = facets_scale) +
          geom_point(size=0.7)+
          geom_smooth(aes(lty = Significance, col=group), method = "lm", se = F, size=1.9, formula = y ~ x)+
          geom_text(data = lm_text, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=3,key_glyph = draw_key_path, parse = T,col="red")+
          scale_color_manual(values = "red")+
          theme_bw() +
          theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,-6,0,-6),legend.text = element_text(size=12, margin = margin(r = 1, unit = 'cm')),
                strip.text = element_text(size=12),axis.text = element_text(size=8),axis.title=element_text(size=18),
                legend.title = element_text(size=12))+
          guides(
            linetype = guide_legend(override.aes = list(col = "#000000",size=0.7,linewidth = 0.7))
          )+
          labs(x = "Species diversity", y = "Multifunctionality")+abc
      }
      
      
      
    }
    else{
      output$group <- output %>% dplyr::select(by_group) %>% unlist()
      
      if(fit == "lm"){
        ###########################
        ############# 11/2 revised
        lm_all <- output %>% group_by(Type,Order.q) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .))) %>% 
          mutate(group="Linear model")
        papa =  output %>% group_by(Type,Order.q) %>%do(broom::glance(lm(qMF ~ Species.diversity, .)))%>%select(r.squared,Order.q,Type)
        lm_all <- suppressMessages(lm_all %>% dplyr::full_join(papa))
      }
      else lm_all <- output %>% group_by(Type, Order.q) %>% do(Lmm_fit(.,r_effect = rand_eff)) %>% suppressMessages %>%
          mutate(group="Linear mixed model")
      
      lm_overall <- mutate(lm_all,
                           Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Nonsignificant slope"),
                                               levels = c("Significant slope (P < 0.05)", "Nonsignificant slope"))
      ) %>% dplyr::select(-c(std.error, statistic, p.value)) %>%
        pivot_wider(names_from = term,values_from = estimate) %>%
        group_by(Type,Order.q,group) %>%
        summarise(Significance=Significance[!is.na(Species.diversity)],
                  Intercept=sum(`(Intercept)`, na.rm = T),
                  estimate=sum(Species.diversity, na.rm=T)) %>% suppressMessages
      
      
      if(fit %in% c("lm","LMM.intercept")){
        lm_data <- output %>% group_by(Type, Order.q, group) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .))) %>%
          mutate(Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Nonsignificant slope"),
                                     levels = c("Significant slope (P < 0.05)", "Nonsignificant slope"))
          ) %>%
          filter(term=="Species.diversity") %>% dplyr::select(-c(term, std.error, statistic, p.value))
        
        output <- suppressMessages(output %>% dplyr::full_join(lm_data))
        lm_data <- lm_data %>% ungroup() %>%
          dplyr::add_row(lm_overall %>% dplyr::select(-Intercept)) %>%
          mutate(group=factor(group,levels = unique(group)))
        
        output$group <- factor(output$group,levels = levels(lm_data$group))
        plot_output <- ggplot(data = output, aes(x = Species.diversity, y = qMF,col = group))+
          facet_grid(Type ~ Order.q, scales = facets_scale) +
          geom_point(size=0.9,alpha=0.5)+
          geom_smooth(aes(lty = Significance), method = "lm", se = F, size=0.5, formula = y ~ x)+
          geom_abline(data = lm_overall, aes(slope=estimate, intercept=Intercept, lty = Significance), size=1.9, col="red",key_glyph = draw_key_path)
      }
      else{
        lm_data <- output %>% group_by(Type, Order.q) %>%
          do(Lmm_fit(.,r_effect = rand_eff,each_group = T)) %>% suppressMessages
        rangeinfo <- output %>% group_by(Type,Order.q,group) %>%
          summarise(x0 = min(Species.diversity),x1=max(Species.diversity))
        
        output <- suppressMessages(output %>% dplyr::full_join(left_join(lm_data,rangeinfo)))
        lm_data <- lm_data %>% ungroup() %>% select(-Intercept) %>%
          rename(estimate = Slope) %>%
          dplyr::add_row(lm_overall %>% dplyr::select(-c(Significance,Intercept))) %>%
          mutate(group=factor(group,levels = unique(group)))
        
        output$group <- factor(output$group,levels = levels(lm_data$group))
        plot_output <- ggplot(data = output, aes(x = Species.diversity, y = qMF,col = group))+
          facet_grid(Type ~ Order.q, scales = facets_scale) +
          geom_point(size=0.9,alpha=0.5)+
          geom_segment(aes(x=x0,xend=x1,y=Intercept+Slope*x0,yend=Intercept+Slope*x1),size=0.5)+
          geom_abline(data = lm_overall, aes(slope=estimate, intercept=Intercept, lty = Significance), size=1.9, col="red",key_glyph = draw_key_path)
      }
      
      
      lm_data <- lm_data %>% group_by(Type,Order.q,group) %>%
        summarise(Label=ifelse(round(estimate,2)==round(estimate,1),
                               paste0(round(estimate,2),"0"),
                               as.character(round(estimate,2))) %>%
                    paste("Slope = ",.,sep = ifelse(round(estimate,2)<0,""," "))) %>%
        mutate(h=rep(c(-0.1,-1.5),(length(unique(group))+1)/2)[1:length(unique(group))],
               v=rep(2+(0:((length(unique(group))+1)/2-1))*1.5,each=2)[1:length(unique(group))]) %>% suppressMessages
      
      col_manual <- c(stdPalette[1:(length(levels(output$group))-1)],"red") %>%
        `names<-`(levels(output$group))
      #####新增#####
      num_group <- length(unique(output$group))
      
      if(text == "slope"){
        
        plot_output <- plot_output +
          geom_text(data = lm_data, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=3, key_glyph = draw_key_path)+
          scale_colour_manual(values = col_manual, drop = F) +
          theme_bw() +
          guides(linetype = guide_legend(title = "",order = 1,override.aes = list(col = "#000000",size=0.6,linewidth = 0.7)),
                 col = guide_legend(title = by_group, override.aes = list(linewidth=c(rep(1,num_group),2))))+ #####新增#####
        theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,-6,0,-6),legend.text = element_text(size=12, margin = margin(r = 1, unit = 'cm')),
              strip.text = element_text(size=12),axis.text = element_text(size=8),axis.title=element_text(size=18),
              legend.title = element_text(size=12))+
          labs(x = "Species diversity", y = "Multifunctionality")+abc
      }
      else{
        if(fit == "lm"){
          #######################
          #### 11/2 revised
          lm_text <- lm_all %>% select(Type,Order.q,r.squared) %>%
            mutate(Label=paste("R^2==",round(r.squared,3),sep=" "),
                   h=-0.3,
                   v=2)
        }
        else{
          lm_text <- lm_all %>% select(Type, Order.q,R2_C,R2_M) %>%
            pivot_longer(cols = c(R2_C,R2_M),names_to = "R2_type",values_to = "Label") %>%
            mutate(Label=paste(ifelse(R2_type=="R2_C","R[C]^2==","R[M]^2=="),
                               round(Label,3),sep=" "),
                   h=ifelse(R2_type=="R2_C",-0.3,-1.8),
                   v=2)
        }
        
        
        plot_output <- plot_output +
          geom_text(data = lm_text, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=3,key_glyph = draw_key_path, parse = T,col="red")+
          scale_colour_manual(values = col_manual, drop = F) +
          theme_bw() +
          guides(linetype = guide_legend(title = "",order = 1,override.aes = list(col = "red",size=0.6,linewidth = 0.7)),
                 col = guide_legend(title = by_group, override.aes = list(linewidth=c(rep(1,num_group),2))))+ #####新增#####
        theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,-6,0,-6),legend.text = element_text(size=12, margin = margin(r = 1, unit = 'cm')),
              strip.text = element_text(size=12),axis.text = element_text(size=8),axis.title=element_text(size=18),
              legend.title = element_text(size=12))+
          labs(x = "Species diversity", y = "Multifunctionality")+abc
      }
    }
  }
  else{
    multiple_plots <- function(data){
      
      scale_plot <- function(out, xlab, ylab, digit=2, h_j=-1.5){
        if(is.null(by_group)){
          #########################################
          ##################11/2revised
          lm_data <- out %>% group_by(Order.q) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .)))
          papa =  out %>% group_by(Order.q) %>% do(broom::glance(lm(qMF ~ Species.diversity, .)))%>%select(r.squared,Order.q)
          lm_data <- suppressMessages(lm_data %>% dplyr::full_join(papa))
          lm_data <- mutate(lm_data,
                            Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Nonsignificant slope"),
                                                levels = c("Significant slope (P < 0.05)", "Nonsignificant slope")),
                            group="Linear model"
          ) %>%
            filter(term=="Species.diversity") %>% dplyr::select(-c(term, std.error, statistic, p.value))
          out <- suppressMessages(out %>% dplyr::left_join(lm_data))
          #########################################
          #########################################
          #########################################
          #########################################
          ##################11/2revised
          if(text=="slope"){
            plot_output <- ggplot(data = out, aes(x = Species.diversity, y = qMF))+
              facet_grid(~Order.q, scales = facets_scale) +
              geom_point(size=0.7,alpha=0.5)+
              geom_smooth(aes(lty = Significance, col=group), method = "lm", se = F, size=0.5, formula = y ~ x)+
              geom_text(data = lm_data, aes(x = -Inf, y = Inf, label=paste0("Slope = ", round(estimate, 4)), col=group),
                        hjust= -0.1, vjust= 2, size=2.7, key_glyph = draw_key_path)+
              scale_color_manual(values = "red")+
              theme_bw() +
              theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,-6,0,-6), legend.title = element_blank(),legend.text = element_text(size=12, margin = margin(r = 1, unit = 'cm')),
                    strip.text = element_text(size=12),axis.text = element_text(size=8),axis.title=element_text(size=18))+
              guides(
                linetype = guide_legend(override.aes = list(col = "#000000",size=0.7,linewidth = 0.7))
              )+
              labs(x = xlab, y = ylab)+abc
          }
          else{
            lm_text <- lm_data %>% select(Order.q,r.squared) %>%
              mutate(Label=paste("R^2==",round(r.squared,3),sep=" "),
                     h=-0.3,
                     v=2)
            plot_output <- ggplot(data = out, aes(x = Species.diversity, y = qMF))+
              facet_grid(~ Order.q, scales = facets_scale) +
              geom_point(size=0.7,alpha=0.5)+
              geom_smooth(aes(lty = Significance, col=group), method = "lm", se = F, size=0.5, formula = y ~ x)+
              geom_text(data = lm_text, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=2.7,key_glyph = draw_key_path, parse = T,col="red")+
              scale_color_manual(values = "red")+
              theme_bw() +
              theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,-6,0,-6), legend.title = element_blank(),legend.text = element_text(size=12, margin = margin(r = 1, unit = 'cm')),
                    strip.text = element_text(size=12),axis.text = element_text(size=8),axis.title=element_text(size=18))+
              guides(
                linetype = guide_legend(override.aes = list(col = "#000000",size=0.7,linewidth = 0.7))
              )+
              labs(x = xlab, y = ylab)+abc
          }
          
          
          
          
        }
        else{
          out$group <- out %>% dplyr::select(by_group) %>% unlist()
          
          if(fit == "lm") {
            #######################
            #### 11/2 revised
            lm_all <- out %>% group_by(Order.q) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .))) %>% 
              mutate(group="Linear model")
            papa =  out %>% group_by(Order.q) %>%do(broom::glance(lm(qMF ~ Species.diversity, .)))%>%select(r.squared,Order.q)
            lm_all <- suppressMessages(lm_all %>% dplyr::full_join(papa))
            
          }
          else lm_all <- out %>% group_by(Order.q) %>% do(Lmm_fit(.,r_effect = rand_eff)) %>% suppressMessages %>%
              mutate(group="Linear mixed model")
          
          lm_overall <- mutate(lm_all,
                               Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Nonsignificant slope"),
                                                   levels = c("Significant slope (P < 0.05)", "Nonsignificant slope"))
          ) %>% dplyr::select(-c(std.error, statistic, p.value)) %>%
            pivot_wider(names_from = term,values_from = estimate) %>%
            group_by(Order.q,group) %>%
            summarise(Significance=Significance[!is.na(Species.diversity)],
                      Intercept=sum(`(Intercept)`, na.rm = T),
                      estimate=sum(Species.diversity, na.rm=T)) %>% suppressMessages
          
          
          if(fit %in% c("lm","LMM.intercept")){
            lm_data <- out %>% group_by(Order.q, group) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .))) %>%
              mutate(Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Nonsignificant slope"),
                                         levels = c("Significant slope (P < 0.05)", "Nonsignificant slope"))
              ) %>%
              filter(term=="Species.diversity") %>% dplyr::select(-c(term, std.error, statistic, p.value))
            
            out <- suppressMessages(out %>% dplyr::full_join(lm_data))
            lm_data <- lm_data %>% ungroup() %>%
              dplyr::add_row(lm_overall %>% dplyr::select(-Intercept)) %>%
              mutate(group=factor(group,levels = unique(group)))
            
            out$group <- factor(out$group,levels = levels(lm_data$group))
            plot_output <- ggplot(data = out, aes(x = Species.diversity, y = qMF,col = group))+
              facet_grid( ~ Order.q, scales = facets_scale) +
              geom_point(size=0.7,alpha=0.1)+
              geom_smooth(aes(lty = Significance), method = "lm", se = F, size=0.5, formula = y ~ x)+
              geom_abline(data = lm_overall, aes(slope=estimate, intercept=Intercept, lty = Significance), size=1.9, col="red",key_glyph = draw_key_path)
          }
          else{
            lm_data <- out %>% group_by(Order.q) %>%
              do(Lmm_fit(.,r_effect = rand_eff,each_group = T)) %>% suppressMessages
            rangeinfo <- out %>% group_by(Order.q,group) %>%
              summarise(x0 = min(Species.diversity),x1=max(Species.diversity))
            
            out <- suppressMessages(out %>% dplyr::full_join(left_join(lm_data,rangeinfo)))
            lm_data <- lm_data %>% ungroup() %>% select(-Intercept) %>%
              rename(estimate = Slope) %>%
              dplyr::add_row(lm_overall %>% dplyr::select(-c(Significance,Intercept))) %>%
              mutate(group=factor(group,levels = unique(group)))
            
            out$group <- factor(out$group,levels = levels(lm_data$group))
            plot_output <- ggplot(data = out, aes(x = Species.diversity, y = qMF,col = group))+
              facet_grid(~Order.q, scales = facets_scale) +
              geom_point(size=0.7,alpha=0.1)+
              geom_segment(aes(x=x0,xend=x1,y=Intercept+Slope*x0,yend=Intercept+Slope*x1),size=0.5)+
              geom_abline(data = lm_overall, aes(slope=estimate, intercept=Intercept, lty = Significance), size=1.9, col="red",key_glyph = draw_key_path)
          }
          
          
          lm_data <- lm_data %>% group_by(Order.q,group) %>%
            summarise(Label=ifelse(round(estimate,digit)==0,paste0("0.",paste0(rep(0,digit),collapse = "")),
                                   ifelse(round(estimate,digit)==round(estimate,digit-1),
                                          paste0(round(estimate,digit),"0"),
                                          as.character(round(estimate,digit)))) %>%
                        paste("Slope = ",.,sep = ifelse(round(estimate,digit)<0,""," "))) %>%
            mutate(h=rep(c(-0.1,h_j),(length(unique(group))+1)/2)[1:length(unique(group))],
                   v=rep(2+(0:((length(unique(group))+1)/2-1))*1.5,each=2)[1:length(unique(group))]) %>% suppressMessages
          
          
          col_manual <- c(stdPalette[1:(length(levels(out$group))-1)],"red") %>%
            `names<-`(levels(out$group))
          #####新增#####
          num_group <- length(unique(out$group))
          
          if(text == "slope"){
            plot_output <- plot_output +
              geom_text(data = lm_data, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=2.7,key_glyph = draw_key_path)+
              scale_colour_manual(values = col_manual, drop = F) +
              theme_bw() +
              guides(linetype = guide_legend(title = "",order = 1,override.aes = list(col = "#000000",size=0.6,linewidth = 0.7)),
                     col = guide_legend(title = by_group, override.aes = list(linewidth=c(rep(1,num_group),2))))+ #####新增#####
            theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,-6,0,-6),legend.text = element_text(size=12, margin = margin(r = 1, unit = 'cm')),
                  strip.text = element_text(size=12),axis.text = element_text(size=8),axis.title=element_text(size=18),
                  legend.title = element_text(size=12))+
              labs(x = xlab, y = ylab)+abc
          }
          else{
            if(fit == "lm"){
              #######################
              #### 11/2 revised
              lm_text <- lm_all %>% select(Order.q,r.squared) %>%
                mutate(Label=paste("R^2==",round(r.squared,3),sep=" "),
                       h=-0.3,
                       v=2)
            }
            else{
              lm_text <- lm_all %>% select(Order.q,R2_C,R2_M) %>%
                pivot_longer(cols = c(R2_C,R2_M),names_to = "R2_type",values_to = "Label") %>%
                mutate(Label=paste(ifelse(R2_type=="R2_C","R[C]^2==","R[M]^2=="),
                                   round(Label,3),sep=" "),
                       h=ifelse(R2_type=="R2_C",-0.3,-1.8),
                       v=2)
            }
            
            
            
            
            
            plot_output <- plot_output +
              geom_text(data = lm_text, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=2.7,key_glyph = draw_key_path, parse = T,col="red")+
              scale_colour_manual(values = col_manual, drop = F) +
              theme_bw() +
              guides(linetype = guide_legend(title = "",order = 1,override.aes = list(col = "red",size=0.6,linewidth = 0.7)),
                     col = guide_legend(title = by_group, override.aes = list(linewidth=c(rep(1,num_group),2))))+ #####新增#####
            theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,-6,0,-6),legend.text = element_text(size=12, margin = margin(r = 1, unit = 'cm')),
                  strip.text = element_text(size=12),axis.text = element_text(size=8),axis.title=element_text(size=18),
                  legend.title = element_text(size=12))+
              labs(x = xlab, y = ylab)+abc
          }
          
          plot_output
        }
      }
      
      
      gamma_data <- data %>% filter(Scale=="Gamma")
      gamma <- scale_plot(gamma_data, "Gamma species diversity", "Gamma Multifunctionality")
      
      alpha_data <- data %>% filter(Scale=="Alpha")
      alpha <- scale_plot(alpha_data, "Alpha species diversity", "Alpha Multifunctionality")
      
      beta_data0 <- data %>% filter(Scale=="Beta", Order.q == "q = 0")
      beta0 <- scale_plot(beta_data0, NULL, "Beta Multifunctionality", 3,h_j=-1.25)
      
      beta_data12 <- data %>% filter(Scale=="Beta", Order.q != "q = 0")
      beta12 <- scale_plot(beta_data12, NULL, NULL,h_j=-1.25)
      
      
      beta <- ggpubr::ggarrange(beta0,beta12,nrow = 1,widths=c(1.5,2.4),legend = "none") %>%
        ggpubr::annotate_figure(bottom = textGrob("Beta species diversity \n", gp = gpar(cex = 1.5),hjust = 0.35))
      
      combine_plots <- ggpubr::ggarrange(gamma,alpha,beta,ncol = 1,heights = c(1,1,1.12),legend = "bottom",common.legend = T)
      
      out <- list("ALL"=combine_plots,
                  "Gamma"=gamma,
                  "Alpha"=alpha,
                  "Beta"=ggpubr::ggarrange(alpha,beta,ncol = 1,heights = c(0,1),legend = "bottom",common.legend = T))
      return(out)
    }
    
    plot_output <- lapply(unique(output$Type),function(x) multiple_plots(output %>% filter(Type==x))) %>%
      `names<-`(unique(output$Type))
    
  }
  return(plot_output)
}

# Get the summary table of the linear mixed model fitted.
#
# \code{LMM_fit} Get the summary table of the linear mixed model fitted.
#
# @param data a data.frame contains 'qMF', 'Species.diversity' and 'group'.
# @param r_effect decide the random effect is "intercept", "slope" or "both".
# @each_group
# @return a summary table for the fitted model.

Lmm_fit <- function(data, r_effect = "intercept", each_group = FALSE){
  if(r_effect == "intercept") x <- lmerTest::lmer(qMF~Species.diversity+(1|group),data) %>% suppressWarnings()
  else if(r_effect == "slope") x <- lmerTest::lmer(qMF~Species.diversity+(0+Species.diversity|group),data) %>% suppressWarnings()
  else if(r_effect == "both") x <- lmerTest::lmer(qMF~Species.diversity+(1+Species.diversity|group),data) %>% suppressWarnings()
  
  if(each_group){
    ret <- as_tibble(data.frame(group=rownames(coef(x)$group),coef(x)$group))
    colnames(ret) <- c("group","Intercept","Slope")
  }
  
  else{
    johnson_r2 <- function(model){
      X <- model.matrix(model)
      n <- nrow(X)
      Beta <- lme4::fixef(model)
      Sf <- var(X %*% Beta)
      Sigma.list <- lme4::VarCorr(model)
      Sl <-
        sum(
          sapply(Sigma.list,
                 function(Sigma)
                 {
                   Z <-X[,rownames(Sigma)]
                   sum(diag(Z %*% Sigma %*% t(Z)))/n
                 }))
      Se <- attr(Sigma.list, "sc")^2
      Sd <- 0
      total.var <- Sf + Sl + Se + Sd
      Rsq.m <- Sf / total.var
      Rsq.c <- (Sf + Sl) / total.var
      
      return(c(Rsq.c,Rsq.m))
    }
    r2_j <- johnson_r2(x)
    ret <- as_tibble(cbind(summary(x)$coefficients[,-3],"R2_C"=r2_j[1],"R2_M"=r2_j[2]), rownames = "term")
    colnames(ret) <- c("term", "estimate", "std.error","statistic","p.value","R2_C","R2_M")
  }
  
  ret
}

# Generate Color Palette for ggplot2
#
# This function creates a color palette suitable for ggplot2 visualizations by evenly spacing colors in the HCL color space. The function ensures that the colors are well-distributed and visually distinct, making it ideal for categorical data where each category needs to be represented by a different color.
#
# @param g An integer indicating the number of distinct colors to generate. This value should be a positive integer, with higher values resulting in a broader range of colors.
# @return A vector of color codes in hexadecimal format, suitable for use in ggplot2 charts and plots. The length of the vector will match the input parameter `g`.
# @examples
# # Generate a palette of 5 distinct colors
# ggplotColors(5)
#
# # Use the generated colors in a ggplot2 chart
# library(ggplot2)
# df <- data.frame(x = 1:5, y = rnorm(5), group = factor(1:5))
# ggplot(df, aes(x, y, color = group)) +
#   geom_point() +
#   scale_color_manual(values = ggplotColors(5))
#
ggplotColors <- function(g){
  d <- 360/g # Calculate the distance between colors in HCL color space
  h <- cumsum(c(15, rep(d,g - 1))) # Create cumulative sums to define hue values
  hcl(h = h, c = 100, l = 65) # Convert HCL values to hexadecimal color codes
}


## ========== no visible global function definition for R CMD check ========== ##
utils::globalVariables(c("coef", ".", "Order.q", "plotID", "species", "abundance", 
                         "MF_g", "Tau","MF_a", "MF_b", 
                         "L_g", "R_g", "L_a","R_a", "L_b", "R_b",
                         "corr_corrected_Gamma", "corr_corrected_Beta", "x1","x2", 
                         "Species_Gamma", "Species_Alpha", "Species_Beta", "r.squared", 
                         "p.value", "term", "std.error", "statistic", "Type", 
                         "Species.diversity", "Significance", "group", "estimate", "Label",
                         "h", "v", "(Intercept)", "Intercept", "Slope", "x0", "R2_C", "R2_M", "R2_type", "."
                         )) 



