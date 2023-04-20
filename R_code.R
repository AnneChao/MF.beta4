library(lme4)
library(dplyr)
library(haven)
library(tidyr)
library(ggpubr)
library(forcats)
library(ggplot2)
library(lmerTest)
library(magrittr)
library(parallel)
library(reshape2)
library(multifunc)
library(tidyverse)
library(missForest)
library(RColorBrewer)



source("Source R code.txt")


## ========================= Load data ==================================== ##
Ratcliffe = read.csv("Ratcliffe_et_al_ELE_12849_ecosystem_function_variables.csv")

species_com = read.csv("FunDivEUROPE_BetaFor_species_composition.csv", sep = '"')[,c('block', 'plot', 'full_species_original', 'X.6')]
colnames(species_com)[4] = "basal_area"
species_com$block = unlist(strsplit(species_com$block, split = " "))
species_com$basal_area = as.numeric(species_com$basal_area)

tmp = species_com[species_com$block == "Germany",]
for (i in 1:3) {
  tmp[,i] = factor(tmp[,i])
}

species_com[species_com$block == "Germany", "basal_area"] = missForest(tmp[,-1])$ximp$basal_area
rm(tmp)



for (x in c("FIN02", "FIN03", "FIN04", "FIN07", "FIN08", "FIN12", 
            "FIN13", "FIN15", "FIN20", "FIN24", "FIN25", "FIN26", "FIN28")) {
  
  species_com[species_com$plot == x & species_com$full_species_original == "Betula.pendula", "basal_area"] =
    sum(species_com[species_com$plot == x & species_com$full_species_original %in% c("Betula.pendula", "Betula.pubescens"), "basal_area"])
  
  species_com = species_com[ -which(species_com$plot == x & species_com$full_species_original == "Betula.pubescens"),]
  
}

species_com[species_com$plot == "GER05" & species_com$full_species_original == "Quercus.petraea", "basal_area"] =
  sum(species_com[species_com$plot == "GER05" & species_com$full_species_original %in% c("Quercus.petraea", "Quercus.robur"), "basal_area"])

species_com = species_com[ -which(species_com$plot == "GER05" & species_com$full_species_original == "Quercus.robur"),]

species_com$full_species_original[species_com$full_species_original == "Betula.pubescens"] = "Betula.pendula"
species_com$full_species_original[species_com$full_species_original == "Quercus.robur"] = "Quercus.petraea"


##
variables = colnames(Ratcliffe)[-(1:5)]

Ratcliffe = standardize.Ratcliffe(Ratcliffe, n_func = 26)     ## select 26 functions 
# Ratcliffe = standardize.Ratcliffe(Ratcliffe, n_func = 24)   ## select 24 functions (remove functions "wue" and "sapling_growth")

variables.std <- paste0(variables, ".std")


correlation = Ratcliffe[,variables.std]
correlation[correlation == -1] = 0
correlation = cor(correlation)

distM = sqrt(1 - abs(correlation))
distM.Byrnes = (1 - correlation) / 2



## ========================= Plot Figure 2 and Appendices S4 ==================================== ##
Ratcliffe <- Ratcliffe %>%
  mutate(mf_Chao_0 = apply(Ratcliffe[,variables.std], 1, function(x) MF.uncor(x, rep(1, length(x)), 0)$qF),
         mf_Chao_1 = apply(Ratcliffe[,variables.std], 1, function(x) MF.uncor(x, rep(1, length(x)), 1)$qF),
         mf_Chao_2 = apply(Ratcliffe[,variables.std], 1, function(x) MF.uncor(x, rep(1, length(x)), 2)$qF),
         
         mf_eff_0 = apply(Ratcliffe[,variables.std], 1, function(x) Byrnes.uncor(x, 0)$qF),
         mf_eff_1 = apply(Ratcliffe[,variables.std], 1, function(x) Byrnes.uncor(x, 1)$qF),
         mf_eff_2 = apply(Ratcliffe[,variables.std], 1, function(x) Byrnes.uncor(x, 2)$qF),
         
         mf_Chao_AUC_0 = apply(Ratcliffe[,variables.std], 1, function(x) MF.cor(x, rep(1, length(x)), distM, q = 0) %>% filter(tau == 'AUC') %>% select(qF) %>% as.numeric),
         mf_Chao_AUC_1 = apply(Ratcliffe[,variables.std], 1, function(x) MF.cor(x, rep(1, length(x)), distM, q = 1) %>% filter(tau == 'AUC') %>% select(qF) %>% as.numeric),
         mf_Chao_AUC_2 = apply(Ratcliffe[,variables.std], 1, function(x) MF.cor(x, rep(1, length(x)), distM, q = 2) %>% filter(tau == 'AUC') %>% select(qF) %>% as.numeric),
         
         mf_eff_AUC_0 = apply(Ratcliffe[,variables.std], 1, function(x) Byrnes.cor(x, distM.Byrnes, q = 0) %>% filter(tau == 'AUC') %>% select(qF) %>% as.numeric),
         mf_eff_AUC_1 = apply(Ratcliffe[,variables.std], 1, function(x) Byrnes.cor(x, distM.Byrnes, q = 1) %>% filter(tau == 'AUC') %>% select(qF) %>% as.numeric),
         mf_eff_AUC_2 = apply(Ratcliffe[,variables.std], 1, function(x) Byrnes.cor(x, distM.Byrnes, q = 2) %>% filter(tau == 'AUC') %>% select(qF) %>% as.numeric)
         )

mf = Ratcliffe %>%
  select(plotid, target_species_richness, mf_Chao_0:mf_eff_AUC_2) %>%
  pivot_longer(cols = c(mf_Chao_0:mf_eff_AUC_2)) %>%
  mutate(method = fct_inorder(name))

## Computer species diversity
mf = mf %>% mutate(target_species_richness = sapply(unique(species_com$plot), function(x) 
  rep( qD( (species_com %>% filter(plot == x))$basal_area, q = c(0,1,2)), 4)) %>% as.vector)

mf = mf %>% mutate(name = paste(plotid, name, sep = '_'),
                   Order.q = rep(c(0, 1, 2), 4*nrow(Ratcliffe)),
                   Type = rep(rep(c('Proposed MF', 'Byrnes et al. (2023)'), each = 3), 2*nrow(Ratcliffe)),
                   Div = rep(rep(c('Uncorrelated', 'Correlated'), each = 6), nrow(Ratcliffe)))

anovas = data.frame('name' = unique(mf$name),
                    'p_value' = sapply(1:length(unique(mf$name)), function(j) {
                      
                      myout_ <- mf %>% filter(name == unique(mf$name)[j])
                      
                      ( lm(formula = value ~ target_species_richness, data = myout_) %>% summary )$coefficients[2, 'Pr(>|t|)']
                      
                    }))

slope = lapply(1:length(unique(mf$name)), function(j) {
  
  myout_ <- mf %>% filter(name == unique(mf$name)[j])
  
  slope = ( lm(formula = value ~ target_species_richness, data = myout_) %>% summary )$coefficients[2, 'Estimate'] 
  slope = ifelse(round(abs(slope), 3) < 0.001, slope %>% round(., 4) %>% paste("Slope = ", ., sep = ""), 
                 ifelse(round(abs(slope), 2) < 0.01, slope %>% round(., 3) %>% paste("Slope = ", ., sep = ""),
                        format(slope %>% round(., 2), nsmall = 2) %>% paste("Slope = ", ., sep = "")))
  
  tmp = myout_[ !duplicated(myout_[, c('name', 'Order.q', 'Type', 'Div')]),] %>% select(c('name', 'Order.q', 'Type', 'Div'))
  
  cbind(tmp, plotid = unique(myout_$plotid), x = max(myout_$target_species_richness), y = max(myout_$value), Slope = slope)
  
}) %>% do.call(rbind,.)

## Fit linear model
mf = mf %>% arrange(name) %>% group_by(name) %>% 
  do(lm(formula = value ~ target_species_richness, data = . ) %>% predict %>% tibble(fit = .)) %>% 
  ungroup %>% select(fit) %>% bind_cols(mf)

mf = mf %>% left_join(., 
                      anovas %>% mutate('Sig' = ifelse(p_value < 0.05, 'Significant slope (P < 0.05)', 'Insignificant slope')),
                      by = 'name') %>% select(-name)

mf$Order.q = as.factor(mf$Order.q)


## Add linear mixed model
slope = rbind(slope,
              
              lapply(unique(mf$method), function(x) {
                
                myout_ <- mf %>% filter(method == x)
                slope = ( lmer(formula = value ~ 1 + target_species_richness + (1 | plotid), 
                               data    = myout_) %>% summary )$coefficients[2, 'Estimate']
                slope = ifelse(round(abs(slope), 3) < 0.001, slope %>% round(., 4) %>% paste("Slope = ", ., sep = ""), 
                               ifelse(round(abs(slope), 2) < 0.01, slope %>% round(., 3) %>% paste("Slope = ", ., sep = ""),
                                      format(slope %>% round(., 2), nsmall = 2) %>% paste("Slope = ", ., sep = "")))
                
                tmp = myout_[ !duplicated(myout_[, c('Order.q', 'Type', 'Div')]),] %>% select(c('Order.q', 'Type', 'Div'))
                
                cbind(name = unique(myout_$method), tmp, plotid = 'Linear mixed', x = max(myout_$target_species_richness), y = max(myout_$value), Slope = slope)
                
              }) %>% do.call(rbind,.)
              )

mf = rbind(mf,
           
           lapply(unique(mf$method), function(x) {
             lmm.data = mf %>% filter(method == x)
             model <- lmer(formula = value ~ 1 + target_species_richness + (1 | plotid), 
                           data    = lmm.data)
             
             data.frame(fit = predict(model, re.form = NA), 
                        plotid = 'Linear mixed', 
                        target_species_richness = lmm.data$target_species_richness, 
                        value = 0, 
                        method = x, 
                        Order.q = unique(lmm.data$Order.q), 
                        Type = unique(lmm.data$Type), 
                        Div = unique(lmm.data$Div),
                        p_value = summary(model)$coefficients[2, 'Pr(>|t|)'],
                        Sig = ifelse( summary(model)$coefficients[2, 'Pr(>|t|)'] < 0.05, 
                                      'Significant slope (P < 0.05)', 'Insignificant slope' ))
           }) %>% do.call(rbind,.)
           )


fig.theme = theme(legend.position = "bottom", 
                  legend.box = "vertical", 
                  legend.key.width = unit(1.2, "cm"), 
                  text = element_text(size = 16), 
                  plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
                  strip.text = element_text(size = 12, face = "bold"))

guide = guides(colour = guide_legend(title = "Plot ID", override.aes = list(linewidth = 1.5)),
               lty = guide_legend(override.aes = list(linewidth = 1)),
               size = "none")

manual = 
  scale_colour_manual(breaks = c("FIN", "GER", "ITA", "POL", 
                                 "ROM", "SPA", "Linear mixed"),
                      label = c("Finland (3)", "Germany (5)", "Italy (5)",  
                                "Poland (5)", "Romania (4)", "Spain (4)", "Linear mixed"),
                      values = c("FIN" = "black", "GER" = "purple2", "ITA" = "darkorange", "POL" = "steelblue1", 
                                 "ROM" = "blue", "SPA" = "gray55", "Linear mixed" = "red")) +
  scale_size_manual(values = c("FIN" = 0.8, "GER" = 0.8, "ITA" = 0.8, "POL" = 0.8, 
                               "ROM" = 0.8, "SPA" = 0.8, "Linear mixed" = 1.9)) +
  scale_linetype_manual(values = c("Insignificant slope" = "dashed", "Significant slope (P < 0.05)" = "solid"), name = NULL)


fig.proposed = ggplot() +
  geom_line(data = mf %>% filter(Div == 'Uncorrelated'),
            aes(x = target_species_richness, y = fit, col = plotid, lty = Sig, size = plotid)) +
  geom_point(data = mf %>% filter(Div == 'Uncorrelated', plotid != 'Linear mixed'),
             aes(x = target_species_richness, y = value, color = plotid), alpha = 0.2) +
  geom_text(data = slope %>% filter(Div == 'Uncorrelated') %>% 
              mutate(x = max(x) - rep(c(3, 1, 3, 1, 3, 1, 3), each = 6), 
                     y = max(y) + c(rep(0.8, 12), rep(0.2, 12), rep(-0.4, 12), rep(-1, 6))),
            aes(x = x, y = y, label = Slope, color = plotid), size = 5, show.legend = FALSE) +
  theme_bw() +
  facet_grid(Type ~ Order.q, scale = "fixed", labeller = labeller(Order.q = c(`0` = "q = 0", `1` = "q = 1", `2` = "q = 2"))) +
  labs(x = "Species diversity", y = "Multifunctionality") +
  manual +
  fig.theme +
  guide +
  coord_cartesian(ylim = c(6, 15.4))


fig.uncor = ggplot() +
  geom_line(data = mf %>% filter(Type == "Proposed MF" & Div != "Functional (dmean)") %>% 
              mutate(Div = fct_inorder(Div)),
            aes(x = target_species_richness, y = fit, col = plotid, lty = Sig, size = plotid)) +
  geom_point(data = mf %>% filter(Type == "Proposed MF" & Div != "Functional (dmean)", plotid != 'Linear mixed') %>% 
               mutate(Div = fct_inorder(Div)),
             aes(x = target_species_richness, y = value, color = plotid), alpha = 0.2) +
  geom_text(data = slope %>% filter(Type == "Proposed MF" & Div != "Functional (dmean)") %>% 
              mutate(Div = fct_inorder(Div),
                     x = max(x) - rep(c(3, 1, 3, 1, 3, 1, 3), each = 6), 
                     y = max(y) + c(rep(0.7, 12), rep(0.2, 12), rep(-0.3, 12), rep(-0.8, 6))),
            aes(x = x, y = y, label = Slope, color = plotid), size = 5, show.legend = FALSE) +
  theme_bw() +
  facet_grid(Div ~ Order.q, scale = "fixed", 
             labeller = labeller(Order.q = c(`0` = "q = 0", `1` = "q = 1", `2` = "q = 2"),
                                 Div = c(`Uncorrelated` = "Uncorrelated functions", 
                                         `Correlated` = "Correlated functions"))) +
  labs(x = "Species diversity", y = "Multifunctionality") +
  manual +
  fig.theme +
  guide +
  coord_cartesian(ylim = c(7.3, 15.2))


fig.cor = ggplot() +
  geom_line(data = mf %>% filter( (Div == "Correlated" & Type == 'Byrnes et al. (2023)' & (value >= 3.3 | value == 0)) | 
                                    (Div == "Correlated" & Type == 'Proposed MF' & (value >= 7.9 | value == 0)) ),
            aes(x = target_species_richness, y = fit, col = plotid, lty = Sig, size = plotid)) +
  geom_point(data = mf %>% filter( plotid != 'Linear mixed' & Div == "Correlated" &
                                     ((Type == 'Byrnes et al. (2023)' & Div == "Correlated" & (value >= 3.3 | value == 0)) | 
                                        (Type == 'Proposed MF' & Div == "Correlated" & (value >= 7.9 | value == 0)))
  ), aes(x = target_species_richness, y = value, color = plotid), alpha = 0.2) +
  geom_text(data = slope %>% filter(Div == "Correlated") %>% 
              mutate(x = max(x) - rep(c(3, 1, 3, 1, 3, 1, 3), each = 6), 
                     y = max(y) + c(rep(rep(c(0, -6), each = 3), 2), rep(rep(c(-0.5, -6.3), each = 3), 2), 
                                    rep(rep(c(-1, -6.6), each = 3), 2), rep(c(-1.5, -6.9), each = 3))),
            aes(x = x, y = y, label = Slope, color = plotid), size = 5, show.legend = FALSE) +
  theme_bw() +
  facet_grid(Type ~ Order.q, scale = "free_y", 
             labeller = labeller(Order.q = c(`0` = "q = 0", `1` = "q = 1", `2` = "q = 2"))) +
  labs(x = "Species diversity", y = "Multifunctionality") +
  manual +
  fig.theme +
  guide 


ggsave("fig.proposed.png", plot = fig.proposed, width = 10, height = 10)
ggsave("fig.uncor.png",    plot = fig.uncor,    width = 10, height = 10)
ggsave("fig.cor.png",      plot = fig.cor,      width = 10, height = 10)



## ====================== Plot Figure 4, 5, and Appendices S3 ==================================== ##
cpu.cores <- detectCores()-1
cl <- makeCluster(cpu.cores)
clusterExport(cl, varlist = c("Ratcliffe", "variables.std", "distM", "qD", "MF.uncor", "MF.cor", "beta.MF.uncor", "beta.MF.cor", "species_com"), 
              envir = environment())
clusterEvalQ(cl, c(library(dplyr), library(tidyr), library(reshape2)))


beta.result.uncor = parLapply(cl, unique(Ratcliffe$plotid), function(x) {
  
  country.data = Ratcliffe %>% filter(plotid == x)
  comb = combn(1:nrow(country.data), 2)
  
  species.data = species_com[substr(species_com$plot, 1, 3) == x,]
  
  lapply(1:ncol(comb), function(i) {
    
    index = comb[,i]
    
    data = country.data[index, ]
    data = data %>% mutate(Species = strsplit(composition, split = '.', fixed = TRUE))
    
    N = length(index)
    
    Species = species.data %>% filter(plot %in% unique(species.data$plot)[index]) %>% acast(., full_species_original ~ plot, value.var = "basal_area")
    Species[is.na(Species)] = 0
    
    Gamma.div = qD(rowSums(Species), q = c(0,1,2))
    Alpha.div = qD(as.vector(Species), q = c(0,1,2)) / N
    Beta.div = Gamma.div / Alpha.div
    
    out = beta.MF.uncor(data[,variables.std] %>% t, q = c(0,1,2))
    
    out %>% filter(Type %in% c('Gamma','Alpha','Beta')) %>% 
      mutate(Gamma = rep(Gamma.div, 3), Alpha = rep(Alpha.div, 3), Beta = rep(Beta.div, 3)) %>%
      pivot_longer(cols = c(Gamma:Beta), names_to = 'Dissimilarity')
    
  }) %>% do.call(rbind,.) %>% mutate(plotid = x)
  
}) %>% do.call(rbind,.)


beta.result.cor = parLapply(cl, unique(Ratcliffe$plotid), function(x) {
  
  country.data = Ratcliffe %>% filter(plotid == x)
  comb = combn(1:nrow(country.data), 2)
  
  species.data = species_com[substr(species_com$plot, 1, 3) == x,]
  
  lapply(1:ncol(comb), function(i) {
    
    index = comb[,i]
    
    data = country.data[index, ]
    data = data %>% mutate(Species = strsplit(composition, split = '.', fixed = TRUE))
    
    N = length(index)
    
    Species = species.data %>% filter(plot %in% unique(species.data$plot)[index]) %>% acast(., full_species_original ~ plot, value.var = "basal_area")
    Species[is.na(Species)] = 0
    
    Gamma.div = qD(rowSums(Species), q = c(0,1,2))
    Alpha.div = qD(as.vector(Species), q = c(0,1,2)) / N
    Beta.div = Gamma.div / Alpha.div
    
    out = beta.MF.cor(data[,variables.std] %>% t, distM = distM, q = c(0,1,2)) %>% filter(tau == 'AUC') %>% select(-tau)
    
    out %>% filter(Type %in% c('Gamma','Alpha','Beta')) %>% 
      mutate(Gamma = rep(Gamma.div, 3), Alpha = rep(Alpha.div, 3), Beta = rep(Beta.div, 3)) %>%
      pivot_longer(cols = c(Gamma:Beta), names_to = 'Dissimilarity')
    
  }) %>% do.call(rbind,.) %>% mutate(plotid = x)
  
}) %>% do.call(rbind,.)


stopCluster(cl)


## Fit linear model
result = beta.result.uncor %>% mutate(name = paste(q, Type, Dissimilarity, plotid, sep = '.'))   ## Select uncorrelated beta result
# result = beta.result.cor %>% mutate(name = paste(q, Type, Dissimilarity, plotid, sep = '.'))   ## Select correlated beta result

anovas = data.frame('name' = unique(result$name),
                    'p_value' = sapply(1:length(unique(result$name)), function(j) {
                      
                      myout_ <- result %>% filter(name == unique(result$name)[j])
                      ( lm(formula = qF ~ value, data = myout_) %>% summary )$coefficients[2, 'Pr(>|t|)']
                      
                    }))

slope = lapply(1:length(unique(result$name)), function(j) {
  myout_ <- result %>% filter(name == unique(result$name)[j])
  slope = ( lm(formula = qF ~ value, data = myout_) %>% summary )$coefficients[2, 'Estimate'] 
  slope = ifelse(round(abs(slope), 3) < 0.001, slope %>% round(., 4) %>% paste("Slope = ", ., sep = ""), 
                 ifelse(round(abs(slope), 2) < 0.01, slope %>% round(., 3) %>% paste("Slope = ", ., sep = ""),
                        format(slope %>% round(., 2), nsmall = 2)  %>% paste("Slope = ", ., sep = "")))
  
  tmp = myout_[ !duplicated(myout_[, c('q', 'Type', 'Dissimilarity', 'plotid')]),] %>% select(c('q', 'Type', 'Dissimilarity', 'plotid'))
  
  cbind(tmp, x = max(myout_$value), y = max(myout_$qF), Slope = slope)
}) %>% do.call(rbind,.)

result = result %>% arrange(name) %>% group_by(name) %>% 
  do(lm(formula = qF ~ value, data = . ) %>% predict %>% tibble(fit = .)) %>% 
  ungroup %>% select(fit) %>% bind_cols(result)

result = result %>% left_join(., 
                              anovas %>% mutate('Sig' = ifelse(p_value < 0.05, 'Significant slope (P < 0.05)', 'Insignificant slope')),
                              by = 'name') %>% select(-name)


## Add Linear mixed model
slope = rbind(slope,
              lapply(0:2, function(i) {
                
                lapply(c('Gamma', 'Alpha', 'Beta'), function(x) {
                  
                  myout_ <- result %>% filter(Type == x, Dissimilarity == x, plotid != 'FIN', q == i)
                  
                  slope = ( lmer(formula = qF ~ 1 + value + (1 | plotid), 
                                 data    = myout_) %>% summary )$coefficients[2, 'Estimate']
                  slope = ifelse(round(abs(slope), 3) < 0.001, slope %>% round(., 4) %>% paste("Slope = ", ., sep = ""), 
                                 ifelse(round(abs(slope), 2) < 0.01, slope %>% round(., 3) %>% paste("Slope = ", ., sep = ""),
                                        format(slope %>% round(., 2), nsmall = 2) %>% paste("Slope = ", ., sep = "")))
                  
                  tmp = myout_[ !duplicated(myout_[, c('q', 'Type', 'Dissimilarity')]),] %>% select(c('q', 'Type', 'Dissimilarity'))
                  
                  cbind(tmp, plotid = 'Linear mixed', x = max(myout_$value), y = max(myout_$qF), Slope = slope)
                  
                }) %>% do.call(rbind,.)
                
              }) %>% do.call(rbind,.)
)

result = rbind(result,
               
               lapply(c(0,1,2), function(i) {
                 
                 lapply(c('Gamma', 'Alpha', 'Beta'), function(x) {
                   lmm.data = result %>% filter(Type == x, Dissimilarity == x, plotid != 'FIN', q == i)
                   
                   model <- lmer(formula = qF ~ 1 + value + (1 | plotid), 
                                 data    = lmm.data)
                   data.frame(fit = predict(model, re.form = NA), 
                              qF = 0, 
                              q = i, 
                              Type = x, 
                              Dissimilarity = x, 
                              value = lmm.data$value, 
                              Species = 'LLM',
                              plotid = 'Linear mixed', 
                              p_value = summary(model)$coefficients[2, 'Pr(>|t|)'],
                              Sig = ifelse( summary(model)$coefficients[2, 'Pr(>|t|)'] < 0.05, 
                                            'Significant slope (P < 0.05)', 'Insignificant slope'))
                 }) %>% do.call(rbind,.)
                 
               }) %>% do.call(rbind,.)
               
)



manual.beta = 
  scale_colour_manual(breaks = c("GER", "ITA", "POL", 
                                 "ROM", "SPA", "Linear mixed"),
                      label = c("Germany (5)", "Italy (5)", "Poland (5)", 
                                "Romania (4)", "Spain (4)", "Linear mixed"),
                      values = c("GER" = "purple2", "ITA" = "darkorange", "POL" = "steelblue1", 
                                 "ROM" = "blue", "SPA" = "gray55", "Linear mixed" = "red")) +
  scale_size_manual(values = c("GER" = 0.8, "ITA" = 0.8, "POL" = 0.8, 
                               "ROM" = 0.8, "SPA" = 0.8, "Linear mixed" = 1.9)) +
  scale_linetype_manual(values = c("Insignificant slope" = "dashed", "Significant slope (P < 0.05)" = "solid"), name = NULL)



## Beta multifunctionality v.s. Beta species diversity
fig.beta.q0 = ggplot() +
  geom_point(data = result %>% filter(Type == 'Beta', Dissimilarity == 'Beta', !plotid %in% c('Linear mixed', 'FIN'), q == 0),
             aes(x = value, y = qF, col = plotid), alpha = 0.05) +
  geom_line(data = result %>% filter(Type == 'Beta', Dissimilarity == 'Beta', plotid != 'FIN', q == 0),
            aes(x = value, y = fit, col = plotid, lty = Sig, size = plotid)) +
  geom_text(data = slope %>% filter(Type == 'Beta', Dissimilarity == 'Beta', plotid != 'FIN', q == 0) %>% 
              mutate(x = max(x) - c(0.8, 0.2, 0.8, 0.2, 0.8, 0.2),
                     y = max(y) + c(rep(-0.037, 2), rep(-0.0375, 2), rep(-0.038, 2))),       ## for uncorrelated (26 functions)
                     # y = max(y) + c(rep(-0.0337, 2), rep(-0.0342, 2), rep(-0.0347, 2))),   ## for correlated (26 functions)
                     # y = max(y) + c(rep(-0.042, 2), rep(-0.0425, 2), rep(-0.043, 2))),     ## for uncorrelated (24 functions)
                     # y = max(y) + c(rep(-0.0388, 2), rep(-0.0393, 2), rep(-0.0398, 2))),   ## for correlated (24 functions)
            aes(x = x, y = y, label = Slope, col = plotid), size = 4, show.legend = FALSE) +
  theme_bw() +
  facet_wrap(q ~ ., labeller = labeller(q = c(`0` = "q = 0"))) +
  labs(x = NULL, y = "Beta multifunctionality") +
  manual.beta +
  fig.theme +
  guide +
  coord_cartesian(ylim = c(1, 1.01))       ## for uncorrelated (26 functions)
  # coord_cartesian(ylim = c(1, 1.009))    ## for correlated (26 functions)
  # coord_cartesian(ylim = c(1, 1.01))     ## for uncorrelated (24 functions)
  # coord_cartesian(ylim = c(1, 1.0085))   ## for correlated (24 functions)


fig.beta.q12 = ggplot() +
  geom_point(data = result %>% filter(Type == 'Beta', Dissimilarity == 'Beta', !plotid %in% c('Linear mixed', 'FIN'), q %in% c(1,2)),
             aes(x = value, y = qF, col = plotid), alpha = 0.05) +
  geom_line(data = result %>% filter(Type == 'Beta', Dissimilarity == 'Beta', plotid != 'FIN', q %in% c(1,2)),
            aes(x = value, y = fit, col = plotid, lty = Sig, size = plotid)) +
  geom_text(data = slope %>% filter(Type == 'Beta', Dissimilarity == 'Beta', plotid != 'FIN', q %in% c(1,2)) %>% 
              mutate(x = max(x) - rep(c(0.8, 0.2, 0.8, 0.2, 0.8, 0.2), each = 2),
                     y = max(y) + c(rep(-0.16, 4), rep(-0.167, 4), rep(-0.174, 4))),     ## (26 functions)
                     # y = max(y) + c(rep(-0.16, 4), rep(-0.167, 4), rep(-0.174, 4))),   ## for uncorrelated (24 functions)
                     # y = max(y) + c(rep(-0.167, 4), rep(-0.174, 4), rep(-0.181, 4))),  ## for correlated (24 functions)
            aes(x = x, y = y, label = Slope, col = plotid), size = 4, show.legend = FALSE) +
  theme_bw() +
  facet_wrap(q ~ ., labeller = labeller(q = c(`1` = "q = 1", `2` = "q = 2"))) +
  labs(x = "Beta species diversity", y = NULL) +
  fig.theme +
  theme(axis.title.x = element_text(hjust = 0)) +
  manual.beta +
  guide +
  coord_cartesian(ylim = c(1.02, 1.15))      ## for uncorrelated (26 functions)
  # coord_cartesian(ylim = c(1.02, 1.145))   ## for correlated (26 functions)
  # coord_cartesian(ylim = c(1.02, 1.15))    ## for uncorrelated (24 functions)
  # coord_cartesian(ylim = c(1.025, 1.14))   ## for correlated (24 functions)

fig.beta = annotate_figure(ggarrange(fig.beta.q0, fig.beta.q12, ncol = 2, nrow = 1, common.legend = TRUE, 
                                 align = 'h', legend = 'bottom', widths = c(1.2, 2))) 



## Gamma multifunctionality v.s. Gamma species diversity
fig.gamma = ggplot() +
  geom_point(data = result %>% filter(Type == 'Gamma', Dissimilarity == 'Gamma', !plotid %in% c('Linear mixed', 'FIN')),
             aes(x = value, y = qF, col = plotid), alpha = 0.05) +
  geom_line(data = result %>% filter(Type == 'Gamma', Dissimilarity == 'Gamma', plotid != 'FIN'),
            aes(x = value, y = fit, col = plotid, lty = Sig, size = plotid)) +
  geom_text(data = slope %>% filter(Type == 'Gamma', Dissimilarity == 'Gamma', plotid != 'FIN') %>% 
              mutate(x = max(x) - rep(c(3, 1, 3, 1, 3, 1), each = 3),
                     y = max(y) + c(rep(1.5, 6), rep(1, 6), rep(0.5, 6))),         ## (26 functions)
                     # y = max(y) + c(rep(0.8, 6), rep(0.3, 6), rep(-0.2, 6))),    ## (24 functions)
            aes(x = x, y = y, label = Slope, col = plotid), size = 5, show.legend = FALSE) +
  theme_bw() +
  facet_wrap(q ~ ., scale = "fixed", labeller = labeller(q = c(`0` = "q = 0", `1` = "q = 1", `2` = "q = 2"))) +
  labs(x = "Gamma species diversity", y = "Gamma multifunctionality") +
  manual.beta +
  fig.theme +
  guide +
  coord_cartesian(ylim = c(7.8, 15.6))    ## for uncorrelated (26 functions)
  # coord_cartesian(ylim = c(7.5, 15.2))  ## for correlated (26 functions)
  # coord_cartesian(ylim = c(7.2, 14.6))  ## for uncorrelated (24 functions)
  # coord_cartesian(ylim = c(6.8, 14))    ## for correlated (24 functions)



## Alpha multifunctionality v.s. Alpha species diversity
fig.alpha = ggplot() +
  geom_point(data = result %>% filter(Type == 'Alpha', Dissimilarity == 'Alpha', !plotid %in% c('Linear mixed', 'FIN')),
             aes(x = value, y = qF, col = plotid), alpha = 0.05) +
  geom_line(data = result %>% filter(Type == 'Alpha', Dissimilarity == 'Alpha', plotid != 'FIN'),
            aes(x = value, y = fit, col = plotid, lty = Sig, size = plotid)) +
  geom_text(data = slope %>% filter(Type == 'Alpha', Dissimilarity == 'Alpha', plotid != 'FIN') %>% 
              mutate(x = max(x) - rep(c(3.3, 1, 3.3, 1, 3.3, 1), each = 3),
                     y = max(y) + c(rep(1.5, 6), rep(1, 6), rep(0.5, 6))),        ## for uncorrelated (26 functions)
                     # y = max(y) + c(rep(1.9, 6), rep(1.2, 6), rep(0.5, 6))),    ## for correlated (26 functions)
                     # y = max(y) + c(rep(1, 6), rep(0.5, 6), rep(0, 6))),        ## (24 functions)
            aes(x = x, y = y, label = Slope, col = plotid), size = 5, show.legend = FALSE) +
  theme_bw() +
  facet_wrap(q ~ ., scale = "fixed", labeller = labeller(q = c(`0` = "q = 0", `1` = "q = 1", `2` = "q = 2"))) +
  labs(x = "Alpha species diversity", y = "Alpha multifunctionality") + 
  manual.beta +
  fig.theme +
  guide +
  coord_cartesian(ylim = c(7, 15.7))      ## for uncorrelated (26 functions)
  # coord_cartesian(ylim = c(6.8, 15.7))  ## for correlated (26 functions)
  # coord_cartesian(ylim = c(6.3, 14.6))  ## for uncorrelated (24 functions)
  # coord_cartesian(ylim = c(6, 14.2))    ## for correlated (24 functions)


ggsave("fig.beta.png",  plot = fig.beta,  width = 10, height = 7)
ggsave("fig.gamma.png", plot = fig.gamma, width = 9,  height = 7)
ggsave("fig.alpha.png", plot = fig.alpha, width = 9,  height = 7)



