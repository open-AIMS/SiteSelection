# Functions for the analyses of the manuscript 
# "Testing a site-prioritisation framework to maximise survival of seeded corals"
# Journal of Applied Ecology
# 2026-05

# Specify the working directory to the folder that contains the folder "Analyses"

# Specify the path to the folder "Analyses"
Path = paste(getwd(),"Analyses", sep = "/")

check_range = function(data, species) {
  
  # Function to check the range of the different environmental variables for
  # each species. It returns the variables for which there is variability.
  
  # data is the data frame that includes environmental variables
  # species is the name of the targeted species ("A. millepora","A. muricata",
  # "M. aequituberculata","A. tersa","A. kenti","A. digitifera","A. spathulata",
  # "A. loripes", "G. fascicularis", or "M. turtlensis")
  
  dat = data[data$Species == species, ]
  
  vars =   c("Water.velocity", 
             "Slope",
             "Rugosity",
             "Sedimentation",
             "Consolidated.substrate",
             "Free.space",
             "Macroalgae", 
             "Hard.coral")
  
  range = c()
  
  for (var in vars) {
    
    Range = unique(dat[, var])
    Range = Range[!is.na(Range)]
    range = c(range, length(Range))
  }
  
  return(vars[which(range > 1)])
}

group_mean_clustering = function(data, site_info) {
  
  # Function to compute mean clustering by reef
  # Returns each variable as a standard deviation of the mean for that reef
  
  reef_env = site_info %>% group_by(Reef) %>%
    summarise("Water.velocity.cm"= mean(Water.velocity, na.rm = TRUE),
              "Slope.cm" = mean(Slope, na.rm = TRUE),
              "Consolidated.substrate.cm" = mean(Consolidated.substrate, na.rm = TRUE),
              "Free.space.cm" = mean(Free.space, na.rm = TRUE),
              "Hard.coral.cm" = mean(Hard.coral, na.rm = TRUE),
              "Macroalgae.cm" = mean(Macroalgae, na.rm = TRUE),
              
              "Water.velocity.sd"= sd(Water.velocity, na.rm = TRUE),
              "Slope.sd" = sd(Slope, na.rm = TRUE),
              "Consolidated.substrate.sd" = sd(Consolidated.substrate, na.rm = TRUE),
              "Free.space.sd" = sd(Free.space, na.rm = TRUE),
              "Hard.coral.sd" = sd(Hard.coral, na.rm = TRUE),
              "Macroalgae.sd" = sd(Macroalgae, na.rm = TRUE),
    )
  
  data = data %>% left_join(reef_env)
  data$Water.velocity = (data$Water.velocity - data$Water.velocity.cm) /data$Water.velocity.sd
  data$Slope = (data$Slope - data$Slope.cm)/ data$Slope.sd
  data$Consolidated.substrate = (data$Consolidated.substrate - data$Consolidated.substrate.cm)/ data$Consolidated.substrate.sd
  data$Free.space = (data$Free.space - data$Free.space.cm)/data$Free.space.sd
  data$Hard.coral = (data$Hard.coral - data$Hard.coral.cm)/data$Hard.coral.sd
  data$Macroalgae = (data$Macroalgae - data$Macroalgae.cm)/data$Macroalgae.sd
  return(data)
}


fit_models <- function(data, species) {
  
  # Function to do model selection and identify key environmental variables
  # for each species.
  # It returns the best-fit model, the AICc comparisons, model variables,
  # the direction of the effect for each variable, p-values,
  # Akaike weights, and the best-fit model's R squared.
  
  
  
  vars = check_range(data, species)
  
  # create a matrix to save the AICc values
  M = matrix(NA, nrow = length(vars) + 1, ncol = 1)
  row.names(M) = c("Previous model", vars)
  
  M2 = M # matrix to save the coefficient direction
  M3 = M # matrix to save the p-values
  M4 = M # matrix to save the Akaike weights
  
  data$n_failures = data$n_trial - data$n_surv
  
  # Specify baseline model depending on the species
  if (species %in% c("A. tersa", "A. millepora")) {
    Form = "cbind(n_surv, n_failures) ~   (1|Site) + Reef "
    model_vars = c("Reef")
  } else {
    Form = "cbind(n_surv, n_failures) ~    (1|Site)  "
    model_vars = c()
  }
  if (species %in% c("A. millepora", "A. muricata", "M. aequituberculata",
                     "A. tersa", "A. digitifera")) {
    Form = paste(Form, "+ offset(log(Deployment.duration)) ", sep = "")
  }
  
  # eliminate rows with missing values in the variables to be selected from
  dat = data[data$Species == species, ]
  dat = dat[complete.cases(dat[ ,which(colnames(dat) %in% vars)]), ]
  dat = dat[complete.cases(dat$Deployment.duration), ]
  
  # specify family distribution 
  if (species %in% c("A. spathulata", "A. digitifera")) {
    # Overdispersion issues in A. spathulata and A. digitifera survival
    Family = "betabinomial"
  } else {
    Family = "binomial"
  }
  
  
  # null model
  m0 = glmmTMB(as.formula(Form), data = dat,
               family = Family)
  
  AIC_m0 = AICc(m0)
  
  
  # Comparison of models with only one environmental variables
  j = 1
  
  imp_vars = c()
  rows = 2:nrow(M)
  M[1,1] = AIC_m0
  
  for (i in rows) {
    Form_i = paste(Form, row.names(M)[i], sep = "+")
    m =  glmmTMB(as.formula(Form_i), data = dat,
                 family = Family)
    summ = summary(m)
    M[i ,j] = AICc(m)
    
    if (species %in% c( "A. tersa", "A. millepora")) {
      start_i = 3
    } else {
      start_i = 2
    }
    M2[i,j] = ifelse(fixef(m)$cond[start_i] > 0, "positive", "negative")
    M3[i,j] = min(summ$coefficients$cond[start_i:nrow(summ$coefficients$cond), 4])
    
    if (row.names(M)[i]  == "Sedimentation") {
      
      CE = as.data.frame(predict_response(m,
                                          terms = list("Sedimentation" =  unique(dat$Sedimentation)),
                                          condition = c(Deployment.duration = 365), margin = "marginalmeans"))
      CE$x = factor(CE$x, levels = c("There are no sediments deposited",
                                     "Thin layer easy to resuspend",
                                     "Moderate layer that offers resistance to be resuspended",
                                     "Thick layer deep layer that is not possible to resuspend" 
      ))
      CE = CE[order(CE$x, decreasing = FALSE), ]
      dir = c()
      
      for (k in 2:nrow(CE)) {
        dir = c(dir, ifelse(CE$predicted[k] > CE$predicted[k-1], "increasing", "decreasing"))
        
      }
      
      if (length(unique(dir)) > 1) {
        M2[i,j] = "mixed"
      } else {
        if (unique(dir) == "decreasing") {
          M2[i,j] = "negative"
        } else {
          M2[i,j] = "positive"
        }
      }
      
    }
    if (row.names(M)[i] == "Rugosity") {
      CE = as.data.frame(predict_response(m,
                                          terms = list("Rugosity" =  unique(dat$Rugosity)),
                                          condition = c(Deployment.duration = 365), margin = "marginalmeans"))
      CE$x = factor(CE$x, levels = c("Low", "Medium", "High"))
      CE = CE[order(CE$x, decreasing = FALSE), ]
      
      dir = c()
      
      for (k in 2:nrow(CE)) {
        dir = c(dir, ifelse(CE$predicted[k] > CE$predicted[k-1], "increasing", "decreasing"))
        
      }
      
      if (length(unique(dir)) > 1) {
        M2[i,j] = "mixed"
      } else {
        if (unique(dir) == "decreasing") {
          M2[i,j] = "negative"
        } else {
          M2[i,j] = "positive"
        }
      }
      
    }
  }
  
  aw = akaike.weights(M[ ,j])
  
  # Identify best-fit model from initial comparison
  min_AIC = min(M[ ,j])
  n = which(aw$weights == max(aw$weights))
  list_n = n
  M4[ ,j] = aw$weights
  
  
  # Add one environmental variable at a time to the best-fit model, and update
  # the best fit model following a model comparison. Repeat this process until
  # adding environmental variables does not improve model fit
  
  if (n != 1) {
    imp_vars = c(imp_vars, names(n))
    Form = paste(Form, names(n), sep = "+")
    M = cbind(M, matrix(NA, nrow = length(vars) +1, ncol = 1))
    M4 = cbind(M4, matrix(NA, nrow = length(vars) +1, ncol = 1))
    rows = rows[rows != which(row.names(M) == names(n ))]
    
    model_vars = c(model_vars, names(n))
    j = j +1
    
    M[1,j] = min_AIC
    
    for (i in rows) {
      Form_i = paste(Form, row.names(M)[i], sep = "+")
      m =  glmmTMB(as.formula(Form_i), data = dat,
                   family = Family)
      M[i,j] = AICc(m)
    }
    
    aw = akaike.weights( M[ ,j])
    M4[-list_n,j] = aw$weights
    
    n = which(M4[,j] == max(M4[,j], na.rm = TRUE))
    list_n = c(list_n, n)
  }
  
  
  
  
  if (species %in% c("A. digitifera", "A. hyacinthus", "A. tersa",
                     "A. kenti", "A. loripes", "A. muricata", "A. spathulata",
                     "M. aequituberculata", "M. turtlensis")) {
    # Data from G. fascicularis and A. millepora spat do not allow for complex
    # models. Additional variables prevent model convergence or result in 
    # model fit issues.
    if( n != 1) {
      while (n != 1 ) {
        
        min_AIC = min( M[ ,j], na.rm = TRUE)
        Form = paste(Form, names(n), sep = "+")
        M = cbind(M, matrix(NA, nrow = length(vars) +1, ncol = 1))
        M4 = cbind(M4, matrix(NA, nrow = length(vars) +1, ncol = 1))
        model_vars = c(model_vars, names(n))
        rows = rows[rows != which(row.names(M) == names(n ))]
        
        j = j +1
        M[1,j] = min_AIC
        for (i in rows) {
          Form_i = paste(Form, row.names(M)[i], sep = "+")
          m =  glmmTMB(as.formula(Form_i), data = dat,
                       family = Family)
          M[i,j] = AICc(m)
          
        }
        aw = akaike.weights( M[ ,j])
        M4[-list_n,j] = aw$weights
        n = which(M4[,j] == max(M4[,j], na.rm = TRUE))
        list_n = c(list_n, n)
        imp_vars = c(imp_vars, names(n))
        
      }
    }
  } 
  
  if (species == "G. fascicularis") {
    
    if( n != 1) {
      while (n != 1 & j < 3) {
        
        min_AIC = min( M[ ,j], na.rm = TRUE)
        Form = paste(Form, names(n), sep = "+")
        M = cbind(M, matrix(NA, nrow = length(vars) +1, ncol = 1))
        M4 = cbind(M4, matrix(NA, nrow = length(vars) +1, ncol = 1))
        model_vars = c(model_vars, names(n))
        rows = rows[rows != which(row.names(M) == names(n ))]
        
        j = j +1
        M[1,j] = min_AIC
        for (i in rows) {
          Form_i = paste(Form, row.names(M)[i], sep = "+")
          m =  glmmTMB(as.formula(Form_i), data = dat,
                       family = Family)
          M[i,j] = AICc(m)
          
        }
        aw = akaike.weights( M[ ,j])
        M4[-list_n,j] = aw$weights
        n = which(M4[,j] == max(M4[,j], na.rm = TRUE))
        
        list_n = c(list_n, n)        
        Form = paste(Form, names(n), sep = "+")
        if (species == "G. fascicularis" & j == 2) { n <- 1
        }
        
        
      }
    }
  }
  
  
  
  
  m =  glmmTMB(as.formula(Form), data = dat,
               family = Family)
  
  
  return(list("model" = m, "M" = M,
              "vars" = model_vars,
              "coef" = M2, "pvals" = M3,
              "aw" = M4,
              "R2" = performance::r2(m)))
  
  
  
}


plot_DHARMa = function(model) {
  # Function to check model fit using the package DHARMa
  require(DHARMa)
  simulationOutput = simulateResiduals(fittedModel = model, quantreg = TRUE)
  print(testOutliers(simulationOutput , type  = "bootstrap"))
  print(testDispersion(simulationOutput, alternative = "two.sided"))
  plot(simulationOutput)
  
}

get_direction = function(data, species) {
  
  # Extract the direction of the relationship between survival and each 
  # environmental variable
  
  vars = check_range(data, species)
  M = matrix(NA, nrow = length(vars), ncol = 1)
  row.names(M) = c( vars)
  
  M2 = M
  
  j = 1
  
  data$n_failures = data$n_trial - data$n_surv
  n_reefs = length(unique(data[data$Species == species, ]$Reef))
  
  if (n_reefs > 1) {
    Form = "cbind(n_surv, n_failures) ~   (1|Site) + Reef "
    model_vars = c("Reef")
    
  } else {
    Form = "cbind(n_surv, n_failures) ~   (1|Site) "
    model_vars = c()
  }
  
  
  if (species %in% c("A. millepora", "A. muricata", "M. aequituberculata",
                     "A. tersa", "A. digitifera")) {
    Form = paste(Form, "+ offset(log(Deployment.duration)) ", sep = "")
  }
  
  dat = data[data$Species == species, ]
  dat = dat[complete.cases(dat[ ,which(colnames(dat) %in% vars)]), ]
  dat = dat[complete.cases(dat$Deployment.duration), ]
  
  
  
  
  if (species %in% c("A. spathulata", "A. digitifera")) {
    # Overdispersion issues in A. spathulata survival
    Family = "betabinomial"
  } else {
    Family = "binomial"
  }
  
  
  for (i in 1:nrow(M)) {
    Form_i = paste(Form, row.names(M)[i], sep = "+")
    m =  glmmTMB(as.formula(Form_i), data = dat,
                 family = Family)
    summ = summary(m)
    
    if (n_reefs > 1) {
      start_i = 3
    } else {
      start_i = 2
    }
    M2[i,j] = ifelse(fixef(m)$cond[start_i] > 0, "positive", "negative")
    
    if (names(M[i,])  == "Sedimentation") {
      
      CE = as.data.frame(predict_response(m,
                                          terms = list("Sedimentation" =  unique(dat$Sedimentation)),
                                          condition = c(Deployment.duration = 365), margin = "marginalmeans"))
      CE$x = factor(CE$x, levels = c("There are no sediments deposited",
                                     "Thin layer easy to resuspend",
                                     "Moderate layer that offers resistance to be resuspended",
                                     "Thick layer deep layer that is not possible to resuspend" 
      ))
      CE = CE[order(CE$x, decreasing = FALSE), ]
      dir = c()
      
      for (k in 2:nrow(CE)) {
        dir = c(dir, ifelse(CE$predicted[k] > CE$predicted[k-1], "increasing", "decreasing"))
        
      }
      
      if (length(unique(dir)) > 1) {
        M2[i,j] = "mixed"
      }
      
    }
    if (names(M[i,])  == "Rugosity") {
      CE = as.data.frame(predict_response(m,
                                          terms = list("Rugosity" =  unique(dat$Rugosity)),
                                          condition = c(Deployment.duration = 365), margin = "marginalmeans"))
      CE$x = factor(CE$x, levels = c("Low", "Medium", "High", "Very high"))
      CE = CE[order(CE$x, decreasing = FALSE), ]
      
      dir = c()
      
      for (k in 2:nrow(CE)) {
        dir = c(dir, ifelse(CE$predicted[k] > CE$predicted[k-1], "increasing", "decreasing"))
        
      }
      
      if (length(unique(dir)) > 1) {
        M2[i,j] = "mixed"
      }
      
    }
    
    
    
  }
  return(M2)
}


plot_models = function(data, species) {
  
  path_name = paste(Path, "outputs", species, sep = "/")
  out = fit_models(data, species)
  
  
  vars = check_range(data, species)
  
  dat = data[data$Species == species, ]
  dat = dat[complete.cases(dat[ ,which(colnames(dat) %in% vars)]), ]
  dat = dat[complete.cases(dat$Deployment.duration), ]
  
  plot_list = list()
  
  
  if (length(out$vars) < 1) {
    require(arm)
    S = summary(out$model)
    
    CE = as.data.frame(predict(out$model, dat[1, ], re.form = NA, type = "response", se.fit = TRUE))
    CE$Reef = unique(dat$Reef)
    
    p1 = ggplot() + 
      geom_jitter(data = dat, aes( x = Reef, y = n_surv/n_trial), width = 0.2, 
                  height = 0.025, pch = 21, col = "black", fill = "grey", alpha = 0.5) +
      geom_errorbar(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), col = "black", width = 0.1, linewidth = 1.1) +
      geom_point(data = CE, aes(x = x, y = predicted), fill = "black", pch = 21, col = "black",
                 size = 3) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      xlab("") + 
      ylab("Proportion of tabs with surviving spat") +
      theme_bw()
    print(p1)
    
    
  }
  
  if ("Reef" %in% out$vars) {
    CE = predict_response(out$model,
                          terms = list("Reef" =  unique(dat$Reef)),
                          condition = c(Deployment.duration = 365), margin = "marginalmeans")
    CE = as.data.frame(CE)
    
    p1 = ggplot() + 
      geom_jitter(data = dat, aes( x = Reef, y = n_surv/n_trial), width = 0.2, 
                  height = 0.025, pch = 21, col = "black", fill = "grey", alpha = 0.5) +
      geom_errorbar(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), col = "black", width = 0.1, linewidth = 1.1) +
      geom_point(data = CE, aes(x = x, y = predicted), fill = "black", pch = 21, col = "black",
                 size = 3) +
      xlab("") + 
      ylab("Proportion of tabs with surviving spat") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      theme_bw()
    print(p1)
    saveRDS(p1, paste(path_name, "reef.rds", sep = "_"))
    
  }
  
  if ("Slope" %in% out$vars) {
    CE = predict_response(out$model,
                          terms = list("Slope" =  seq(min(dat$Slope), max(dat$Slope), l = 100)),
                          condition = c(Deployment.duration = 365), margin = "marginalmeans")
    CE = as.data.frame(CE)
    
    p2 = ggplot() + 
      geom_jitter(data = dat, aes( x = Slope, y = n_surv/n_trial), width = 0.01, 
                  height = 0.025, pch = 21, col = "black", fill = "#CC79A7", alpha = 0.5) +
      geom_ribbon(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), fill = "#CC79A7", alpha = 0.5) +
      geom_line(data = CE, aes(x = x, y = predicted), col = "#CC79A7", linewidth = 1.1) +
      xlab("Slope (S.D.)") + 
      ylab("Proportion of tabs with surviving spat") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      theme_bw()
    print(p2)
    saveRDS(p2, paste(path_name, "slope.rds", sep = "_"))
  }
  
  if ("Rugosity" %in% out$vars) {
    CE = predict_response(out$model,
                          terms = list("Rugosity" =  unique(dat$Rugosity)),
                          condition = c(Deployment.duration = 365), margin = "marginalmeans")
    CE = as.data.frame(CE)
    limits = c("Low", "Medium", "High", "Very high")
    limits = limits[which(limits %in% unique(dat$Rugosity))]
    p3 = ggplot() + 
      geom_jitter(data = dat, aes( x = Rugosity, y = n_surv/n_trial), width = 0.2, 
                  height = 0.025, pch = 21, col = "black", fill = "#009373", alpha = 0.5) +
      geom_errorbar(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), col = "black", width = 0.1, linewidth = 1.1) +
      geom_point(data = CE, aes(x = x, y = predicted), fill = "#009373", pch = 21, col = "black",
                 size = 3) +
      xlab("Rugosity") + 
      ylab("Proportion of tabs with surviving spat") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      scale_x_discrete(limits = limits) +
      theme_bw()
    print(p3)
    saveRDS(p3, paste(path_name, "rugosity.rds", sep = "_"))
  }
  
  
  if ("Sedimentation" %in% out$vars) {
    CE = predict_response(out$model,
                          terms = list("Sedimentation" =  unique(dat$Sedimentation)),
                          condition = c(Deployment.duration = 365), margin = "marginalmeans")
    CE = as.data.frame(CE)
    
    limits = c( "There are no sediments deposited",
                "Thin layer easy to resuspend",
                "Moderate layer that offers resistance to be resuspended",
                "Thick layer deep layer that is not possible to resuspend"
    )
    labels = c("No sediments", "Thin layer", "Moderate layer", "Thick layer")
    labels = labels[which(limits %in% unique(dat$Sedimentation))]
    limits = limits[which(limits %in% unique(dat$Sedimentation))]
    
    p4 = ggplot() + 
      geom_jitter(data = dat, aes( x = Sedimentation, y = n_surv/n_trial), width = 0.2, 
                  height = 0.025, pch = 21, col = "black", fill = "#F0E442", alpha = 0.5) +
      geom_errorbar(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), col = "black", width = 0.1, linewidth = 1.1) +
      geom_point(data = CE, aes(x = x, y = predicted), fill = "#F0E442", pch = 21, col = "black",
                 size = 3) +
      xlab("Sedimentation") + 
      ylab("Proportion of tabs with surviving spat") +
      scale_x_discrete(limits = limits, labels = labels) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      theme_bw()
    print(p4)
    saveRDS(p4, paste(path_name, "sedimentation.rds", sep = "_"))
  }
  
  
  if ("Consolidated.substrate" %in% out$vars) {
    CE = predict_response(out$model,
                          terms = list("Consolidated.substrate" =  seq(min(dat$Consolidated.substrate, na.rm = TRUE),
                                                                       max(dat$Consolidated.substrate, na.rm = TRUE), l = 100)),
                          condition = c(Deployment.duration = 365), margin = "marginalmeans")
    CE = as.data.frame(CE)
    
    
    p5 = ggplot() + 
      geom_jitter(data = dat, aes( x = Consolidated.substrate, y = n_surv/n_trial), width = 0.2, 
                  height = 0.025, pch = 21, col = "black", fill = "#56B4E9", alpha = 0.5) +
      geom_ribbon(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), fill = "#56B4E9", alpha = 0.5) +
      geom_line(data = CE, aes(x = x, y = predicted), col = "#56B4E9", linewidth = 1.1) +
      xlab("Percentage of consolidated substrate substrate (S.D.)") + 
      ylab("Proportion of tabs with surviving spat") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      theme_bw()
    
    print(p5)
    saveRDS(p5, paste(path_name, "consolidated.rds", sep = "_"))
  }
  
  if ("Free.space" %in% out$vars) {
    CE = predict_response(out$model,
                          terms = list("Free.space" =  seq(min(dat$Free.space, na.rm = TRUE), max(dat$Free.space, na.rm = TRUE), l = 100)),
                          condition = c(Deployment.duration = 365), margin = "marginalmeans")
    CE = as.data.frame(CE)
    
    
    
    p6 = ggplot() + 
      geom_jitter(data = dat, aes( x = Free.space, y = n_surv/n_trial), width = 0.2, 
                  height = 0.025, pch = 21, col = "black", fill = "#0072B2", alpha = 0.5) +
      geom_ribbon(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), fill = "#0072B2", alpha = 0.5) +
      geom_line(data = CE, aes(x = x, y = predicted), col = "#0072B2", linewidth = 1.1) +
      xlab("Percentage of free consolidated substrate (S.D.)") + 
      ylab("Proportion of tabs with surviving spat") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      theme_bw()
    print(p6)
    saveRDS(p6, paste(path_name, "free.rds", sep = "_"))
  }
  
  if ("Macroalgae" %in% out$vars) {
    CE = predict_response(out$model,
                          terms = list("Macroalgae" =  seq(min(dat$Macroalgae), max(dat$Macroalgae), l = 100)),
                          condition = c(Deployment.duration = 365), margin = "marginalmeans")
    CE = as.data.frame(CE)
    
    p7 = ggplot() + 
      geom_jitter(data = dat, aes( x = Macroalgae, y = n_surv/n_trial), width = 0.01, 
                  height = 0.025, pch = 21, col = "black", fill = "#D55E00", alpha = 0.5) +
      geom_ribbon(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), fill = "#D55E00", alpha = 0.5) +
      geom_line(data = CE, aes(x = x, y = predicted), col = "#D55E00", linewidth = 1.1) +
      xlab("Percentage of macroalgal cover (S.D.)") + 
      ylab("Proportion of tabs with surviving spat") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      theme_bw()
    print(p7)
    saveRDS(p7, paste(path_name, "macroalgae.rds", sep = "_"))
  }
  
  if ("Hard.coral" %in% out$vars) {
    CE = predict_response(out$model,
                          terms = list("Hard.coral" =  seq(min(dat$Hard.coral), max(dat$Hard.coral), l = 100)),
                          condition = c(Deployment.duration = 365), margin = "marginalmeans")
    CE = as.data.frame(CE)
    
    p8 = ggplot() + 
      geom_jitter(data = dat, aes( x = Hard.coral, y = n_surv/n_trial), width = 0.01, 
                  height = 0.025, pch = 21, col = "black", fill = "#E69F00", alpha = 0.5) +
      geom_ribbon(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), fill = "#E69F00", alpha = 0.5) +
      geom_line(data = CE, aes(x = x, y = predicted), col = "#E69F00", linewidth = 1.1) +
      xlab("Percentage of hard coral cover (S.D.)") + 
      ylab("Proportion of tabs with surviving spat") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      theme_bw()
    print(p8)
    saveRDS(p8, paste(path_name, "hardcoral.rds", sep = "_"))
  }
  
  
  if ("Water.velocity" %in% out$vars) {
    CE = predict_response(out$model,
                          terms = list("Water.velocity" =  seq(min(dat$Water.velocity), max(dat$Water.velocity), l = 100)),
                          condition = c(Deployment.duration = 365), margin = "marginalmeans")
    CE = as.data.frame(CE)
    
    p9 = ggplot() + 
      geom_jitter(data = dat, aes( x = Water.velocity, y = n_surv/n_trial), width = 0.01, 
                  height = 0.025, pch = 21, col = "black", fill = "#999999", alpha = 0.5) +
      geom_ribbon(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high), fill = "#999999", alpha = 0.5) +
      geom_line(data = CE, aes(x = x, y = predicted), col = "#999999", linewidth = 1.1) +
      xlab("Water velocity (S.D.)") + 
      ylab("Proportion of tabs with surviving spat") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      theme_bw()
    print(p9)
    saveRDS(p9, paste(path_name, "watervel.rds", sep = "_"))
  }
  
  
  return(plot_list)
  
}



################################################################################
############################### Score functions ################################
################################################################################

# Score environmental variables according to Humanes et al. 2025. The direction
# component indicates the direction in which the environmental variable is 
# empirically related to survival. If the direction does not match the direction 
# assumed by Humanes et al, the score calculation is modified to match the 
# direction indicated by the empirical data.

# Humanes, A., Fabricius, K. E., Ferrari, R., & Ortiz, J. C. (2025). 
# Ecological quantitative criteria for reef site prioritisation to maximise
# survivorship and growth of outplanted corals. Journal of Environmental 
# Management, 392, 126585.

# Water velocity
water_score = function(water_vel, direction) {
  
  if (direction == "positive") {
    
    score = ifelse(water_vel < 0.1, 0,
                   ifelse(water_vel < 0.2, 1,
                          ifelse(water_vel < 0.3, 2, 
                                 ifelse(water_vel < 0.4, 3, 4))))
  }
  
  if (direction == "negative") {
    
    score = ifelse(water_vel < 0.1, 4,
                   ifelse(water_vel < 0.2, 3,
                          ifelse(water_vel < 0.3, 2, 
                                 ifelse(water_vel < 0.4, 1, 0))))
  }
  
  if (direction == "original") {
    score = ifelse(water_vel < 0.1, 1,
                   ifelse(water_vel < 0.2, 3,
                          ifelse(water_vel < 0.3, 4, 
                                 ifelse(water_vel < 0.4, 2, 0))))
  }
  return(score)
}

# Reef slope
slope_score = function(slope_val, direction) {
  
  if (direction == "positive") {
    
    score = ifelse(slope_val <= 10, 0,
                   ifelse(slope_val <= 20, 1,
                          ifelse(slope_val <= 30, 2, 
                                 ifelse(slope_val <= 40, 3, 4))))
  }
  
  if (direction == "negative") {
    
    score = ifelse(slope_val <= 10, 4,
                   ifelse(slope_val <= 20, 3,
                          ifelse(slope_val <= 30, 2, 
                                 ifelse(slope_val <= 40, 1, 0))))
  }
  
  if (direction == "original") {
    
    score = ifelse(slope_val <= 10, 2,
                   ifelse(slope_val <= 20, 4,
                          ifelse(slope_val <= 30, 3, 
                                 ifelse(slope_val <= 40, 1, 0))))
  }
  return(score)
}

# Rugosity
rugosity_score = function(rugosity_val, direction) {
  
  if (direction == "positive") {
    score = ifelse(rugosity_val == "Low", 2,
                   ifelse(rugosity_val == "Medium", 3, 4))
  }
  
  if (direction == "negative") {
    score = ifelse(rugosity_val == "Low", 4,
                   ifelse(rugosity_val == "Medium", 3, 2))
  }
  
  if (direction == "mixed") {
    score = ifelse(rugosity_val == "Low", 2,
                   ifelse(rugosity_val == "Medium", 4, 3))
  }
  
  if (direction == "original") {
    score = ifelse(rugosity_val == "Low", 2,
                   ifelse(rugosity_val == "Medium", 4, 3))
  }
  return(score)
  
}

# Sedimentation
sedimentation_score = function(sedimentation_val, direction) {
  
  if (direction == "negative") {
    score = ifelse(sedimentation_val == "There are no sediments deposited", 4,
                   ifelse(sedimentation_val == "Thin layer easy to resuspend", 3,
                          ifelse(sedimentation_val == "Moderate layer that offers resistance to be resuspended", 1, 0)))
  }
  
  if (direction == "mixed") {
    score = ifelse(sedimentation_val == "There are no sediments deposited", 3,
                   ifelse(sedimentation_val == "Thin layer easy to resuspend", 3,
                          ifelse(sedimentation_val == "Moderate layer that offers resistance to be resuspended", 3, 0)))
  }
  
  if (direction == "positive") {
    score = ifelse(sedimentation_val == "There are no sediments deposited", 0,
                   ifelse(sedimentation_val == "Thin layer easy to resuspend", 1,
                          ifelse(sedimentation_val == "Moderate layer that offers resistance to be resuspended", 2, 3)))
  }
  
  if (direction == "original") {
    score = ifelse(sedimentation_val == "There are no sediments deposited", 4,
                   ifelse(sedimentation_val == "Thin layer easy to resuspend", 3,
                          ifelse(sedimentation_val == "Moderate layer that offers resistance to be resuspended", 1, 0)))
  }
  return(score)
  
}

consolidated_score = function(consolidated_val, direction) {
  
  
  if (direction == "positive") {
    
    score = ifelse(consolidated_val <= 40, 0,
                   ifelse(consolidated_val <= 55, 1,
                          ifelse(consolidated_val <= 70, 2, 
                                 ifelse(consolidated_val <= 85, 3, 4))))
  }
  
  if (direction == "negative") {
    
    score = ifelse(consolidated_val <= 40, 4,
                   ifelse(consolidated_val <= 55, 3,
                          ifelse(consolidated_val <= 70, 2, 
                                 ifelse(consolidated_val <= 85, 1, 0))))
  }
  
  if (direction == "original") {
    
    score = ifelse(consolidated_val <= 40, 0,
                   ifelse(consolidated_val <= 55, 1,
                          ifelse(consolidated_val <= 70, 2, 
                                 ifelse(consolidated_val <= 85, 3, 4))))
  }
  return(score)
  
}


# Free available space
free_score = function(free_val, direction) {
  
  
  if (direction == "positive") {
    
    score = ifelse(free_val <= 30, 0,
                   ifelse(free_val <= 45, 1,
                          ifelse(free_val <= 60, 2, 
                                 ifelse(free_val <= 75, 3, 4))))
  }
  
  if (direction == "negative") {
    
    score = ifelse(free_val <= 30, 4,
                   ifelse(free_val <= 45, 3,
                          ifelse(free_val <= 60, 2, 
                                 ifelse(free_val <= 75, 1, 0))))
  }
  
  
  if (direction == "original") {
    
    score = ifelse(free_val <= 30, 0,
                   ifelse(free_val <= 45, 1,
                          ifelse(free_val <= 60, 2, 
                                 ifelse(free_val <= 75, 3, 4))))
  }
  return(score)
  
}

# Macroalgal cover
macroalgae_score = function(macroalgae_val, direction) {
  
  
  if (direction == "positive") {
    
    score = ifelse(macroalgae_val == 0, 0,
                   ifelse(macroalgae_val <= 10, 1,
                          ifelse(macroalgae_val <= 20, 2, 
                                 ifelse(macroalgae_val <= 30, 3, 4))))
  }
  
  if (direction == "negative") {
    
    score = ifelse(macroalgae_val == 0, 4,
                   ifelse(macroalgae_val <= 10, 3,
                          ifelse(macroalgae_val <= 20, 2, 
                                 ifelse(macroalgae_val <= 30, 1, 0))))
  }
  
  if (direction == "original") {
    
    score = ifelse(macroalgae_val == 0, 4,
                   ifelse(macroalgae_val <= 10, 3,
                          ifelse(macroalgae_val <= 20, 2, 
                                 ifelse(macroalgae_val <= 30, 1, 0))))
  }
  return(score)
  
}


# Hard coral cover
hard.coral_score = function(hard.coral_val, direction) {
  
  
  if (direction == "positive") {
    
    score = ifelse(hard.coral_val <= 5, 0,
                   ifelse(hard.coral_val <= 15, 1,
                          ifelse(hard.coral_val <= 20, 2, 
                                 ifelse(hard.coral_val <= 30, 3, 4))))
  }
  
  if (direction == "negative") {
    
    score = ifelse(hard.coral_val <= 5, 4,
                   ifelse(hard.coral_val <= 15, 3,
                          ifelse(hard.coral_val <= 20, 2, 
                                 ifelse(hard.coral_val <= 30, 1, 0))))
  }
  
  
  if (direction == "original") {
    
    score = ifelse(hard.coral_val <= 5, 1,
                   ifelse(hard.coral_val <= 15, 4,
                          ifelse(hard.coral_val <= 20, 3, 
                                 ifelse(hard.coral_val <= 30, 2, 0))))
  }
  return(score)
  
}
################################################################################


abiotic.biotic_score = function(data, model, weighting) {
  
  # Function to calculate the overall site score.
  # An "even" weighting gives the same weight to all environmental variables
  # An "uneven" weighting (not used in the analyses) gives a higher weight
  # to environmental variables included in the best-fit models for each species.
  # Returns three vectors: one for the abiotic scores, one for the biotic scores,
  # and an overall site score.
  
  variables = model$vars
  variables = variables[!variables == "Reef"]
  variables_abiotic = variables[which(variables %in% c("Water.velocity",
                                                       "Slope",
                                                       "Rugosity", "Sedimentation",
                                                       "Consolidated.substrate"))]
  variables_biotic = variables[which(variables %in% c("Free.space", "Macroalgae",
                                                      "Hard.coral"))]
  
  data$score.all = NA
  
  
  if (weighting == "uneven") {
    
    if (length(variables_abiotic) > 0) {
      high_weight_abiotic = (1 - 0.10 * (5 - length(variables_abiotic))) / length(variables_abiotic)
      low_weight_abiotic = 0.10
    } else {
      low_weight_abiotic = 1 / 5
      high_weight_abiotic = NA
    }
    
    if (length(variables_biotic) > 0) {
      high_weight_biotic = (1 - 0.165 * (3 - length(variables_biotic))) / length(variables_biotic)
      low_weight_biotic = 0.165
    } else {
      low_weight_biotic = 1 / 3
      high_weight_biotic = NA
    }
    
    if (length(variables) > 0) {
      low_weight = (1/8) * 0.5
      high_weight = (1 -  (1/8) * 0.5 * (8 - length(variables))) / length(variables)
    } else {
      low_weight = 1/8
    }
    weights = c(rep(low_weight_abiotic, 5), rep(low_weight_biotic, 3))
    weights_all = rep(low_weight, 8)
    names(weights) = c("Water.velocity", "Slope", "Rugosity", "Sedimentation",
                       "Consolidated.substrate", "Free.space", "Macroalgae",
                       "Hard.coral")
    names(weights_all) = c("Water.velocity", "Slope", "Rugosity", "Sedimentation",
                           "Consolidated.substrate", "Free.space", "Macroalgae",
                           "Hard.coral")
    weights[which(names(weights) %in% variables_abiotic )] = high_weight_abiotic
    weights[which(names(weights) %in% variables_biotic )] = high_weight_biotic
    weights_all[which(names(weights_all) %in% variables)] = high_weight
  } else {
    weights_all  = rep(1/8, 8)
    weights = rep(1/8, 8)
  }
  
  
  
  
  for (i in 1:nrow(data)) {
    data$abiotic.score[i] = data$score.water.velocity[i] * weights[1] +
      data$score.slope[i] * weights[2] + 
      data$score.rugosity[i] * weights[3] + 
      data$score.sedimentation[i] * weights[4] +
      data$score.consolidated[i] * weights[5]
    
    data$biotic.score[i] = data$score.free.space[i] * weights[6] +
      data$score.macroalgae[i] * weights[7] + 
      data$score.hard.coral[i] * weights[8] 
    data$score.all[i] = data$score.water.velocity[i] * weights_all[1] +
      data$score.slope[i] * weights_all[2] + 
      data$score.rugosity[i] * weights_all[3] + 
      data$score.sedimentation[i] * weights_all[4] +
      data$score.consolidated[i] * weights_all[5] +
      data$score.free.space[i] * weights_all[6] +
      data$score.macroalgae[i] * weights_all[7] + 
      data$score.hard.coral[i] * weights_all[8] 
    
  }
  
  return(list("abiotic.score" = data$abiotic.score,
              "biotic.score" = data$biotic.score,
              "score.all" = data$score.all))
  
}

################################################################################
####### Function to modify scores to match direction of empirical data #########
################################################################################

modify_scores = function(data1, data, species) {
  
  # Function to modify the scores to match the direction of the relationship
  # between spat survival and each environmental variable from the empirical
  # data.
  
  # data1 is in z-scores
  # data is in raw scale
  
  # It returns the input data frame with three additional columns: an abiotic
  # score, a biotic score, and an overall score (assuming even weights).
  
  model = fit_models(data1, species)
  
  data = data[data$Species == species, ]
  data$score.all = NA
  for (i in 1:nrow(data)) {
    data$score.all[i] = mean(c(data$score.water.velocity[i],
                               data$score.slope[i],
                               data$score.rugosity[i],
                               data$score.sedimentation[i],
                               data$score.consolidated[i],
                               data$score.free.space[i],
                               data$score.macroalgae[i],
                               data$score.hard.coral[i]))
  }
  
  if ("Water.velocity" %in% model$vars) {
    data$score.water.velocity = sapply(data$Water.velocity, water_score,
                                       direction = model$coef[which(row.names(model$coef) == "Water.velocity")])
  }
  
  if ("Slope" %in% model$vars) {
    data$score.slope = sapply(data$score.slope, slope_score,
                              direction = model$coef[which(row.names(model$coef) == "Slope")])
  }
  
  if ("Rugosity" %in% model$vars) {
    data$score.rugosity = sapply(data$Rugosity, rugosity_score,
                                 direction = model$coef[which(row.names(model$coef) == "Rugosity")])
  }
  
  if ("Sedimentation" %in% model$vars) {
    data$score.sedimentation = sapply(data$Sedimentation, sedimentation_score,
                                      direction = model$coef[which(row.names(model$coef) == "Sedimentation")])
  }
  if ("Consolidated" %in% model$vars) {
    data$score.consolidated = sapply(data$Consolidated.substrate, consolidated_score,
                                     direction = model$coef[which(row.names(model$coef) == "Consolidated.substrate")])
  }
  if ("Free.space" %in% model$vars) {
    data$score.free.space = sapply(data$Free.space, free_score,
                                   direction = model$coef[which(row.names(model$coef) == "Free.space")])
  }
  if ("Macroalgae" %in% model$vars) {
    data$score.macroalgae = sapply(data$Macroalgae, macroalgae_score,
                                   direction = model$coef[which(row.names(model$coef) == "Macroalgae")])
  }
  
  if ("Hard.coral" %in% model$vars) {
    data$score.hard.coral = sapply(data$Hard.coral, hard.coral_score,
                                   direction = model$coef[which(row.names(model$coef) == "Hard.coral")])
  }
  
  # Get the updated overall site score
  abio.bio = abiotic.biotic_score(data, model, "even")
  data$abiotic.score = abio.bio$abiotic.score
  data$biotic.score = abio.bio$biotic.score
  data$score.all.up = abio.bio$score.all
  
  return(data)
  
}

get_overallscore = function(data, weights_all) {
  
  # Function to get a site score with specific weighting for each variable. 
  # The weights_all vector must sum to 1 and be in the following order:
  # weight for water velocity, slope, rugosity, sedimentation, total consolidated space,
  # free consolidate space, macroalgal cover, and hard coral cover.
  
  data$score.all = NA
  for (i in 1:nrow(data)){
    data$score.all[i] = data$score.water.velocity[i] * weights_all[1] +
      data$score.slope[i] * weights_all[2] + 
      data$score.rugosity[i] * weights_all[3] + 
      data$score.sedimentation[i] * weights_all[4] +
      data$score.consolidated[i] * weights_all[5] +
      data$score.free.space[i] * weights_all[6] +
      data$score.macroalgae[i] * weights_all[7] + 
      data$score.hard.coral[i] * weights_all[8] 
  }
  return(data$score.all)
  
}



modify_scores_sites = function(data, species, Eff, weights) {
  # Function to modify scores according to the direction of the effect in the
  # empirical data. It gets the direction from the function "get_direction"
  
  model = Eff[Eff$Species == species, ]
  colnames(model)[which(colnames(model) == "direction")] = "coef"
  
  
  data$score.water.velocity = sapply(data$Water.velocity, water_score,
                                     direction = model[model$var == "Water.velocity", ]$coef)
  
  
  data$score.slope = sapply(data$score.slope, slope_score,
                            direction = model[model$var == "Slope", ]$coef)
  
  
  data$score.rugosity = sapply(data$Rugosity, rugosity_score,
                               direction = model[model$var == "Rugosity", ]$coef)
  
  
  data$score.sedimentation = sapply(data$Sedimentation, sedimentation_score,
                                    direction = model[model$var == "Sedimentation", ]$coef)
  
  
  data$score.consolidated = sapply(data$Consolidated.substrate, consolidated_score,
                                   direction = model[model$var == "Consolidated.substrate", ]$coef)
  
  data$score.free.space = sapply(data$Free.space, free_score,
                                 direction =model[model$var == "Free.space", ]$coef)
  
  data$score.macroalgae = sapply(data$Macroalgae, macroalgae_score,
                                 direction = model[model$var == "Macroalgae", ]$coef)
  
  data$score.hard.coral = sapply(data$Hard.coral, hard.coral_score,
                                 direction = model[model$var == "Hard.coral", ]$coef)
  
  
  # Get the updated overall site score
  abio.bio = get_overallscore(data, weights[weights$Species == species, ]$weights)
  
  return(abio.bio)
  
}


get_overallscore_original = function(data) {
  
  # Function to get the original site score with even weighting for all variables. 
  data$score.all = NA
  for (i in 1:nrow(data)){
    data$score.all[i] = (data$score.water.velocity[i] +
                           data$score.slope[i]+ 
                           data$score.rugosity[i]  + 
                           data$score.sedimentation[i] +
                           data$score.consolidated[i]  +
                           data$score.free.space[i]  +
                           data$score.macroalgae[i]  + 
                           data$score.hard.coral[i]) /8 
    
    min_score = NA
    min_score = min(c(data$score.water.velocity[i],
                      data$score.slope[i],
                      data$score.rugosity[i], 
                      data$score.sedimentation[i],
                      data$score.consolidated[i],
                      data$score.free.space[i],
                      data$score.macroalgae[i], 
                      data$score.hard.coral[i]))
    if (min_score == 0) {
      data$score.all[i] = 0
    }
  }
  return(data$score.all)
  
}



sensitivity <- function(data1, data, species_list) {
  
  # Function to investigate the sensitivity of the model fit to the weighting
  # of the environmental variables using a 'one factor at a time' approach (OAT).
  # Returns a data frame that contains species, environmental variable, weight assigned,
  # and R squared for a model predicting survival as a function of the resulting
  # site score.
  
  # data1 is in z-scores
  # data is in raw scale
  
  scores = c(  "score.water.velocity",
               "score.slope",
               "score.rugosity",
               "score.sedimentation",
               "score.consolidated",
               "score.free.space",
               "score.macroalgae",
               "score.hard.coral")
  
  # Weights assigned to the target variable.
  weights_v = seq(0, 1, by = 0.1)
  
  df = data.frame("Species" = factor(), "var" = factor(),
                  "weight" = numeric(), "R2" = numeric())
  
  data1$n_failures = data1$n_trial - data1$n_surv
  
  for (species in species_list) {
    abio.bio = modify_scores(data1, data, species)
    data3 = data1[data1$Species == species, ]
    
    # weight order: water.velocity, slope, rugosity, sedimentation, consolidated,
    # free.space, macroalgae, hard.coral
    data3$score.all = get_overallscore(abio.bio, rep(1/8, 8) )
    
    if (species %in% c("A. tersa", "A. millepora")) {
      Form = "cbind(n_surv, n_failures) ~  (1|Site) + Reef + score.all"
    } else {
      Form = "cbind(n_surv, n_failures) ~ (1|Site) + score.all"
    }
    
    if (species %in% c("A. millepora", "A. muricata", "M. aequituberculata",
                       "A. tersa", "A. digitifera")) {
      Form = paste(Form, "+ offset(log(Deployment.duration)) ", sep = "")
    }
    
    Family = ifelse(species %in% c("A. digitifera", "A. spathulata", "A. millepora"),
                    "betabinomial", "binomial")
    
    m =  glmmTMB(as.formula(Form), data = data3,
                 family = Family)
    
    df = rbind(df, data.frame("Species" = species, "var" = "original", 
                              "weight" = 1/8, "R2" = round(as.numeric(performance::r2(m)[2]), 3)))
    
    for(scorei in scores) {
      
      for (weight_n in weights_v) {
        
        # weights assigned to the non-targeted environmental variable
        other_weight = (1 - weight_n) / 7
        weights_list = rep(other_weight, 8)
        weights_list[which(scores == scorei)] = weight_n
        data3$score.all = get_overallscore(abio.bio, weights_list)
        
        m =  glmmTMB(as.formula(Form), data = data3,
                     family = "binomial")
        
        df = rbind(df, data.frame("Species" = species, "var" = unlist(strsplit(scorei, "score."))[2], 
                                  "weight" = weight_n, "R2" = round(as.numeric(performance::r2(m)[2]), 3)))
        
      }
      
      
    }
    
  }
  
  return(df)
  
  
}


optim_scores2 = function(param, data, species) {
  
  # Function to optimise the weighting of the environmental factors when
  # calculating the site score for each species. 
  # It returns 1 - R square of the model fit (model predicting survival as 
  # a function of site score) so that the optimising function can minimise this
  # value.
  
  data = data[data$Species == species, ]
  param2 = param
  
  n_reefs = length(unique(data$Reef))
  data$score.all = get_overallscore(data, param2)
  
  if (n_reefs >1 ){
    Form = "cbind(n_surv, n_failures) ~  (1|Site) + Reef + score.all"
  } else {
    Form = "cbind(n_surv, n_failures) ~  (1|Site) + score.all"
  }
  
  
  if (species %in% c("A. millepora", "A. muricata", "M. aequituberculata",
                     "A. tersa", "A. digitifera")) {
    Form = paste(Form, "+ offset(log(Deployment.duration)) ", sep = "")
  }
  
  m =  glmmTMB(as.formula(Form), data = data,
               family = "binomial")
  R2 = as.numeric(performance::r2(m)[2])
  sm = NULL
  sm = summary(m)$coefficients$cond
  
  if (is.null(sm)) {
    R2 = 0
  } else {
    if (sm[nrow(sm), 1] < 0) {
      R2 = 0
    }
  }
  
  
  
  return(1 - R2)
  
}


fit_score2 = function(data1, data, species, weights) {
  
  
  # Function to fit a regression of survival as a function of site score.
  # It returns the plotted regressions, the R squared, the slope (coefficient) of the regression,
  # and the p-value.
  
  data1 = data1[data1$Species == species, ]
  data = data[data$Species == species, ]
  data1$n_failures = data1$n_trial - data1$n_surv
  data$n_failures = data$n_trial - data$n_surv
  
  if (species %in% c("A. tersa", "A. millepora")) {
    Form3 = "cbind(n_surv, n_failures) ~  (1|Site) + Reef + score.all"
    Form4 = "cbind(n_surv, n_failures) ~  (1|Site) + Reef + score.all.up"
  } else {
    Form3 = "cbind(n_surv, n_failures) ~  (1|Site) + score.all"
    Form4 = "cbind(n_surv, n_failures) ~ (1|Site) + score.all.up"
  }
  
  if (species %in% c("A. millepora", "A. muricata", "M. aequituberculata",
                     "A. tersa", "A. digitifera")) {
    Form3 = paste(Form3, "+ offset(log(Deployment.duration)) ", sep = "")
    Form4 = paste(Form4, "+ offset(log(Deployment.duration)) ", sep = "")
  }
  
  dat = data1
  
  abio.bio = modify_scores(data1, data, species)
  
  dat$score.all.up = get_overallscore(abio.bio[abio.bio$Species == species, ], weights)
  dat$score.all = get_overallscore(abio.bio[abio.bio$Species == species, ], rep(1/8, 8))
  
  m3 =  glmmTMB(as.formula(Form3), data = dat,
                family = "binomial")
  
  m4 =  glmmTMB(as.formula(Form4), data = dat,
                family = "binomial")
  
  
  CE = predict_response(m3,
                        terms = list("score.all" =  seq(min(dat$score.all), max(dat$score.all), l = 100)),
                        condition = c(Deployment.duration = 365), margin = "marginalmeans")
  CE = as.data.frame(CE)
  CE$Species = species
  
  CE2 = predict_response(m4,
                         terms = list("score.all.up" =  seq(min(dat$score.all.up), max(dat$score.all.up), l = 100)),
                         condition = c(Deployment.duration = 365), margin = "marginalmeans")
  CE2 = as.data.frame(CE2)
  CE2$Species = species
  
  
  p3 = ggplot() + 
    geom_ribbon(data = CE, aes(x  = x, ymin = conf.low, ymax = conf.high, col = Species), alpha = 0.5, fill = NA) +
    geom_line(data = CE, aes(x = x, y = predicted, col = Species),linewidth = 1.1) +
    
    geom_ribbon(data = CE2, aes(x  = x, ymin = conf.low, ymax = conf.high, col = Species, fill = Species), alpha = 0.2, linetype = "dashed") +
    geom_line(data = CE2, aes(x = x, y = predicted, col = Species),linewidth = 1.1, linetype = "dashed") +
    xlab("Site score") + 
    ylab("Proportion of tabs with surviving spat") +
    theme_classic() +
    scale_colour_manual(values = colour_p) +
    scale_fill_manual(values = colour_p) +
    scale_x_continuous(limits = c(1, 4.3), breaks = seq(1,4, by = 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = "none", text = element_text(size = 10))
  
  
  sm3 = summary(m3)$coefficients$cond
  sm4 = summary(m4)$coefficients$cond
  
  return(list("score.plot" = p3,
              "R2_all1" = round(as.numeric(performance::r2(m3)[2]), 3), "R2_all2" = round(as.numeric(performance::r2(m4)[2]), 3),
              "R2_allc1" = round(as.numeric(performance::r2(m3)[1]), 3), "R2_allc2" = round(as.numeric(performance::r2(m4)[1]), 3),
              "all1" = sm3[nrow(sm3), 1:2], "all2" = sm4[nrow(sm4), 1:2],
              "p_all1" = sm3[nrow(sm3), 4],
              "p_all2" = sm4[nrow(sm4), 4]
              
  ))
}



fit_score3 = function(data1, data, species, weight_list) {
  
  # Function to fit a regression of survival as a function of site score.
  # It returns the plotted regressions, the R squared, the slope (coefficient) of the regression,
  # and the p-value.
  
  # To be used in the bootstrapping function. It avoids plotting and returns 
  # less outputs so the bootstrapping is faster.
  
  data1 = data1[data1$Species == species, ]
  data = data[data$Species == species, ]
  data1$n_failures = data1$n_trial - data1$n_surv
  data$n_failures = data$n_trial - data$n_surv
  
  if (species %in% c( "A. tersa", "A. millepora")) {
    
    Form4 = "cbind(n_surv, n_failures) ~  (1|Site) + Reef + score.all.up"
  } else {
    
    Form4 = "cbind(n_surv, n_failures) ~ (1|Site) + score.all.up"
  }
  
  if (species %in% c("A. millepora", "A. muricata", "M. aequituberculata",
                     "A. tersa", "A. digitifera")) {
    Form4 = paste(Form4, "+ offset(log(Deployment.duration)) ", sep = "")
  }
  
  
  
  
  abio.bio = modify_scores(data1, data, species)
  
  R2_all2 = c()
  R2_allc2 = c()
  all2 = c()
  for (k in 1:nrow(weight_list)) {
    abio.bio$score.all.up = get_overallscore(abio.bio[abio.bio$Species == species, ], weight_list[k, ])
    
    
    
    m4 =  glmmTMB(as.formula(Form4), data = abio.bio,
                  family = "binomial")
    
    
    
    sm4 = summary(m4)$coefficients$cond
    R2_all2[k] = round(as.numeric(performance::r2(m4)[2]), 3)
    R2_allc2[k] = round(as.numeric(performance::r2(m4)[1]), 3)
    all2[k] = sm4[nrow(sm4), 1]
    
  }
  
  
  return(list("R2_all2" = R2_all2,
              "R2_allc2" = R2_allc2,
              "all2" = all2
              
  ))
}




fit_score4 = function(data1, data, species, weight_list) {
  
  # Function to fit a regression of survival as a function of site score.
  # It returns the plotted regressions, the R squared, the slope (coefficient) of the regression,
  # and the p-value.
  
  # To be used in the score_predict() function. It returns min and max scores.
  
  data1 = data1[data1$Species == species, ]
  data = data[data$Species == species, ]
  data1$n_failures = data1$n_trial - data1$n_surv
  data$n_failures = data$n_trial - data$n_surv
  
  if (species %in% c( "A. tersa", "A. millepora")) {
    
    Form4 = "cbind(n_surv, n_failures) ~  score.all.up"
  } else {
    
    Form4 = "cbind(n_surv, n_failures) ~  score.all.up"
  }
  
  
  
  
  
  
  data1$score.all.up = get_overallscore(data1[data1$Species == species, ], weight_list)
  
  
  
  m4 =  glmmTMB(as.formula(Form4), data = data1,
                family = "binomial")
  
  
  
  sm4 = summary(m4)$coefficients$cond
  R2_all2 = round(as.numeric(performance::r2(m4)[2]), 3)
  R2_allc2 = round(as.numeric(performance::r2(m4)[1]), 3)
  all2 = sm4[nrow(sm4), 1]
  
  
  
  
  return(list("R2_all2" = R2_all2,
              "R2_allc2" = R2_allc2,
              "all2" = all2,
              "intercept" = sm4[1, 1],
              "slope" = sm4[2,1],
              "min_score" = min(data1$score.all.up, na.rm = TRUE),
              "max_score" = max(data1$score.all.up, na.rm = TRUE)
              
  ))
}


# Predictability
score_predict = function(data2, dat, nsim) {
  
  # Function to test how well can site scores predict spat survival in new 
  # data. It splits the data in three, using two thirds to fit the models and
  # one third to check out of sample prediction. 
  
  # nsim is the number of times the script randomly splits the data in three and
  # checks the predictions. The models do not always converge so the for loop
  # requires updates depending on the seed and number of iterations. e.g., run it
  # first with "for (i in 1:nsim) " but may need to be updated to e.g.,
  # "for (i in 70:nsim) {" if the 70th iteration is not successful at first. If this
  # is the case, do not rewrite R2_df, just add rows in each iteration.
  
  R2_df = data.frame("Species" = factor(),
                     "R2" = numeric(),
                     "R2_conditional" = numeric(),
                     "intercept" = numeric(),
                     "slope" = numeric(),
                     "min_score" = numeric(),
                     "max_score" = numeric(),
                     "convergence" = integer()
  )
  
  set.seed(12)
  empty_mat = matrix(NA, nrow = 8, ncol = nsim)
  row.names(empty_mat) = c("Water.velocity", "Slope", "Rugosity",
                           "Sedimentation", "Consolidated", "Free.space", "Macroalgae",
                           "Hard.coral")
  
  
  sp_weights = list(empty_mat, empty_mat, empty_mat, empty_mat, empty_mat,
                    empty_mat, empty_mat, empty_mat, empty_mat, empty_mat)
  
  for (i in 1:nsim) {
    data2$data.use = NA
    dat$data.use = NA
    
    
    for (species in unique(data2$Species)) {
      n0 = which(data2$Species == species)
      n1 = sample(n0, round(length(n0)*(2/3), 0), replace = FALSE )
      n2 = setdiff(n0, n1)
      data2[n1, ]$data.use = "fit"
      data2[n2, ]$data.use = "test"
      dat[n1, ]$data.use = "fit"
      dat[n2, ]$data.use = "test"
      
      MinLS <- nlminb(start = rep(1/8, 8), 
                      objective = optim_scores2, 
                      data = modify_scores2(data2[data2$data.use == "fit" & data2$Species == species, ],
                                            species),
                      species = species,
                      control= list(abs.tol = 1e-1, rel.tol = 1e-1, step.min = 5, step.max = 10),
                      lower = 0)
      params  = MinLS$par
      params = params/sum(params)
      
      sp_weights[[which(unique(data2$Species) == species)]][ ,i] = params
      
      Fit = fit_score4(data2[data2$data.use == "test" & data2$Species == species, ],
                       dat[dat$data.use == "test" & dat$Species == species, ],
                       species, params)
      R2_df = rbind(R2_df,
                    data.frame("Species" = species,
                               "R2" = Fit$R2_all2,
                               "R2_conditional" = Fit$R2_allc2,
                               "intercept" = Fit$intercept,
                               "slope" = Fit$slope,
                               "min_score" = Fit$min_score,
                               "max_score" = Fit$max_score,
                               "convergence" = MinLS$convergence))
      print(species)
      
    }
    print(i)
    
  }
  
  weights = data.frame("Species" = factor(),
                       "var" = factor(),
                       "mean" = numeric(),
                       "se" = numeric())
  
  sp_i = 1
  
  for (species in unique(data2$Species)) {
    weights_mat = sp_weights[[sp_i]]
    
    for (j in 1:nrow(weights_mat)) {
      weights = rbind(weights,
                      data.frame("Species" = species,
                                 "var" = row.names(weights_mat)[j],
                                 "mean" = mean(weights_mat[j, ]),
                                 "se" = sd(weights_mat[j, ]/sqrt(nsim))) )
    }
    sp_i = sp_i + 1
  }
  
  return(R2_df)
}






score_predict_sites = function(data2, nsim, weights, Eff) {
  
  # Function to quantify the improvement in mean survival when using score-based
  # site selection
  
  
  R2_df = data.frame("Species" = factor(),
                     "Reef" = factor(),
                     "Year" = factor(),
                     "survival_random" = numeric(),
                     "survival_original" = numeric(),
                     "survival_random2" = numeric(),
                     "survival_modified" = numeric(),
                     "survival_max" = numeric(),
                     "survival_min" = numeric(),
                     "n_sites_chosen" = numeric(),
                     "n_sites_total" = numeric()
  )
  
  set.seed(12)
  
  
  
  
  
  
  for (i in 1:nsim) {
    
    
    for (species in unique(data2$Species)) {
      
      Reefs = unique(data2[data2$Species == species, ]$Reef) 
      
      for (reef in Reefs) {
        sub_data2 = data2[data2$Species == species & data2$Reef == reef, ]
        Years = unique(sub_data2$Deployment.year)
        
        for (year in Years) {
          
          sub_data3 = sub_data2[sub_data2$Deployment.year == year, ]
          sub_data3 = sub_data3[sample(1:nrow(sub_data3), nrow(sub_data3), replace = TRUE), ]
          sites = unique(sub_data3[sub_data3$Species == species, ]$Site)
          
          sub_data3$score.all = get_overallscore_original(sub_data3)
          
          sub_data3$score.all.up = modify_scores_sites(sub_data3, species, Eff, weights )
          
          
          
          
          site_survival = sub_data3 %>% 
            group_by(Site, score.all, score.all.up) %>%
            summarise("mean_survival" = mean(n_surv / n_trial, na.rm = TRUE),
                      "n" = n())
          site_survival.d.n.d = site_survival[!site_survival$score.all == 0, ]
          
          survival_random = NULL; survival_original = NULL; survival_modified = NULL; survival_random2 = NULL
          survival_random = sample(site_survival$mean_survival, round((1/2)*length(sites), 0), replace = FALSE)
          survival_random2 = sample(site_survival$mean_survival, round((1/2)*length(sites), 0), replace = FALSE)
          
          if (nrow(site_survival.d.n.d) >= round((1/2)*length(sites), 0)) {
            survival_original = mean(site_survival.d.n.d[order(site_survival.d.n.d$score.all, decreasing = TRUE), ][1:round((1/2)*length(sites), 0), ]$mean_survival)
          } else {
            if (nrow(site_survival.d.n.d) > 0) {
              survival_original = mean(site_survival.d.n.d$mean_survival)
            } 
          }
          
          if (nrow(site_survival.d.n.d) == 0) {
            survival_original = NA
          }
          
          survival_modified = site_survival[order(site_survival$score.all.up, decreasing = TRUE), ][1:round((1/2)*length(sites), 0), ]$mean_survival
          survival_max = site_survival[order(site_survival$mean_survival, decreasing = TRUE), ][1:round((1/2)*length(sites), 0), ]$mean_survival
          survival_min = site_survival[order(site_survival$mean_survival, decreasing = FALSE), ][1:round((1/2)*length(sites), 0), ]$mean_survival
          
          
          R2_df = rbind(R2_df,
                        data.frame("Species" = species,
                                   "Reef" = reef,
                                   "Year" = year,
                                   "survival_random" = mean(survival_random),
                                   "survival_original" = survival_original,
                                   "survival_random2" = mean(survival_random2),
                                   "survival_modified" = mean(survival_modified),
                                   "survival_max" = mean(survival_max),
                                   "survival_min" = mean(survival_min),
                                   "n_sites_chosen" = round((1/2)*length(sites), 0),
                                   "n_sites_total" = length(sites)))
          
          
          write.csv(R2_df, paste(Path,"outputs/R2_df_sites_modified.csv", sep = "/" ), row.names = FALSE)
        }
        
        
        
      }
      
    }
    print(i)
  }
  
  
  
  
  
  return(R2_df)
}


