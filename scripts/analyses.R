# Statistical analyses for the manuscript:
# "Testing a site-prioritisation framework to maximise survival of seeded corals"
# Journal of Applied Ecology
# 2026-05


# Load libraries
library(readxl)
library(qpcR)
library(dplyr)
library(glmmTMB)
library(ggplot2)
library(ggeffects)
library(AICcmodavg)
library(DHARMa)
library(viridis)
library(patchwork)
library(ggstance)
library(ggpubr)

set.seed(9)

## Specify the working directory to the folder that contains the folder "Analyses"

# Specify the path to the folder "Analyses"
Path = paste(getwd(),"Analyses", sep = "/")

# Load required functions
source(paste(Path, "scripts", "functions.R", sep = "/"))

# Read survival and environmental data
survival = read_xlsx(paste(Path, "data", "survival.xlsx", sep = "/"))
site_info = read_xlsx(paste(Path, "data", "environment.xlsx", sep = "/"))

# Join datasets
dat = as.data.frame(survival[ , -which(colnames(survival) %in% c("Deployment.duration"))] %>% left_join(site_info[ , -which(colnames(site_info) == "Species")]))

# Group categories when required for analyses
dat[complete.cases(dat$Rugosity) & dat$Rugosity == "Very high", ]$Rugosity = "High"

# Get number of tabs with no living spat
dat$n_failures = dat$n_trial - dat$n_surv

# Group deployment years when they are separated by a few days (e.g., Dec 2022 and Jan 2023)
dat[dat$Species == "A. loripes" & dat$Reef == "Davies Reef", ]$Deployment.year = "2022 & 2023"
dat[dat$Species == "A. tersa" & dat$Reef == "Davies Reef", ]$Deployment.year = "2022 & 2023"
dat[dat$Species == "M. turtlensis" & dat$Reef == "Davies Reef", ]$Deployment.year = "2022 & 2023"


# Set colour palette for plotting
colour_p = viridis(10, option = "H")
names(colour_p) = c("A. millepora", "A. muricata", "M. aequituberculata", "A. loripes", "M. turtlensis",
                    "G. fascicularis", "A. tersa", "A. kenti", "A. digitifera","A. spathulata")   

# Get site score as proposed by Humanes et al. (2025)
dat$overall.score = get_overallscore(dat, rep(1/8, 8))

# Get mean group clustering
data2 = group_mean_clustering(dat, site_info)

# Fit models
Maeq_m = fit_models(data2, "M. aequituberculata") # OKAY
Aspa_m = fit_models(data2, "A. spathulata") # needs betabinomial 
Amil_m = fit_models(data2, "A. millepora") # can only handle one variable
Adig_m = fit_models(data2, "A. digitifera") #  needs betabinomial 
Ater_m = fit_models(data2, "A. tersa") # GOOD
Amur_m = fit_models(data2, "A. muricata") # KS deviation
Alor_m = fit_models(data2, "A. loripes") # GOOD
Aken_m = fit_models(data2, "A. kenti") # OKAY
Gfas_m = fit_models(data2, "G. fascicularis") # GOOD
Mtur_m = fit_models(data2, "M. turtlensis") # GOOD


# Check model fit
plot_DHARMa(Maeq_m$model)
plot_DHARMa(Aspa_m$model)
plot_DHARMa(Amil_m$model)
plot_DHARMa(Adig_m$model)
plot_DHARMa(Ater_m$model)
plot_DHARMa(Amur_m$model)
plot_DHARMa(Alor_m$model)
plot_DHARMa(Aken_m$model)
plot_DHARMa(Gfas_m$model)
plot_DHARMa(Mtur_m$model)



# Plot fitted models 
plot_models(data2, "M. aequituberculata")
plot_models(data2, "A. spathulata")
plot_models(data2, "A. millepora")
plot_models(data2, "A. digitifera")
plot_models(data2, "A. tersa")
plot_models(data2, "A. muricata")
plot_models(data2, "A. loripes")
plot_models(data2, "A. kenti")
plot_models(data2, "G. fascicularis")
plot_models(data2, "M. turtlensis")

A1 = readRDS(paste(Path, "outputs", "A. tersa_slope.rds", sep = "/"))
B1 = readRDS(paste(Path, "outputs", "G. fascicularis_slope.rds", sep = "/")) 
C1 = readRDS(paste(Path, "outputs", "A. loripes_slope.rds", sep = "/"))  
A2 = readRDS(paste(Path, "outputs", "A. spathulata_rugosity.rds", sep = "/")) 
B2 = readRDS(paste(Path, "outputs", "A. muricata_rugosity.rds", sep = "/")) 
C2 = readRDS(paste(Path, "outputs", "A. digitifera_rugosity.rds", sep = "/")) 
A3 = readRDS(paste(Path, "outputs", "A. spathulata_watervel.rds", sep = "/")) 
B3 = readRDS(paste(Path, "outputs", "A. muricata_watervel.rds", sep = "/")) 
C3 = readRDS(paste(Path, "outputs", "A. spathulata_sedimentation.rds", sep = "/"))
A4 = readRDS(paste(Path, "outputs", "A. digitifera_consolidated.rds", sep = "/"))
B4 = readRDS(paste(Path, "outputs", "A. tersa_consolidated.rds", sep = "/"))


A5 = readRDS(paste(Path, "outputs", "A. spathulata_macroalgae.rds", sep = "/")) 
B5 = readRDS(paste(Path, "outputs", "A. muricata_macroalgae.rds", sep = "/")) 
C5 = readRDS(paste(Path, "outputs", "M. aequituberculata_macroalgae.rds", sep = "/")) 
A6 = readRDS(paste(Path, "outputs", "A. kenti_hardcoral.rds", sep = "/")) 
B6 = readRDS(paste(Path, "outputs", "G. fascicularis_hardcoral.rds", sep = "/")) 
C6 = readRDS(paste(Path, "outputs", "M. turtlensis_hardcoral.rds", sep = "/")) 
A7 = readRDS(paste(Path, "outputs", "A. spathulata_free.rds", sep = "/")) 
B7 = readRDS(paste(Path, "outputs", "A. millepora_free.rds", sep = "/")) 


abiotic_plot = A1 + ylab("") + xlab("") + B1 + ylab("") + xlab("")  + C1 + ylab("") + xlab("")  +
  A2 + ylab("") + xlab("") + B2 + ylab("") + xlab("") + C2 + ylab("") + xlab("")  + 
  A3 + ylab("") + xlab("") + B3 + ylab("") + xlab("")  + C3 + ylab("") + xlab("")  +
  A4 + ylab("") + xlab("") + B4 + ylab("") + xlab("") + plot_layout(ncol = 3)

#ggsave(paste(Path, "outputs", "Abiotic_vars.png", sep = "/"), width = 18, 
#       height = 15, units = "cm")

biotic_plot = A5 + ylab("") + xlab("") + B5 + ylab("") + xlab("")  + C5 + ylab("") + xlab("")  +
  A6 + ylab("") + xlab("") + B6 + ylab("") + xlab("") + C6 + ylab("") + xlab("")  + 
  A7 + ylab("") + xlab("") + B7 + ylab("") + xlab("") + plot_layout(ncol = 3)
#ggsave(paste(Path, "outputs", "Biotic_vars.png", sep = "/"), width = 18, 
#       height = 12, units = "cm")



# Extract and plot the direction of the relationship between each environmental
# variable and species' survival
Eff = expand.grid("Species" = unique(data2$Species),
                  "var" = c("Water.velocity", "Slope", "Rugosity", "Sedimentation",
                            "Consolidated.substrate", "Free.space", "Macroalgae", 
                            "Hard.coral"))
Eff$direction = NA
Eff$bf = 0
Eff$p_val = NA

for (sp in unique(data2$Species)) {
  
  model = fit_models(data2, sp)
  if (length(model$vars) > 0 ) {
    Eff[Eff$Species == sp & Eff$var %in% model$vars, ]$bf = 1
  }
  
  
  for (var in unique(Eff$var)) {
    Eff[Eff$Species == sp & Eff$var == var, ]$direction = model$coef[which(row.names(model$coef) == var)]
    Eff[Eff$Species == sp & Eff$var == var, ]$p_val = model$pvals[which(row.names(model$pvals) == var)]
  }
  exclude = setdiff(unique(Eff$var), check_range(data2, sp))
  
}

# Change direction of sedimentation for A. spathulata. When sedimentation is the only
# explanatory variable accounted for, the relationship is positive. However, when 
# accounting for other variables in its best-fit model, the relationship is negative.
Eff[Eff$Species == "A. spathulata" & Eff$var == "Sedimentation", ]$direction <- "negative"


Fig2 = ggplot() +
  geom_tile(data = Eff, aes(y = Species, x = var, fill = direction), col = "black", alpha = 0.4) +
  geom_tile(data = Eff[Eff$bf ==1, ], aes(y = Species, x = var, fill = direction), col = "black", alpha = 1) +
  geom_point(data = Eff[Eff$bf == 1, ], aes(y = Species, x = var), pch = 23,
             size = 6) +
  scale_fill_manual(values = c("#0072B2", "#666666", "#D55E00"), name = "",
                    limits =c("positive", "mixed", "negative"),
                    labels = c("positive", "nonlinear", "negative")) +
  scale_x_discrete(limits = c("Water.velocity", "Slope", "Rugosity", "Sedimentation",
                              "Consolidated.substrate", "Free.space", "Macroalgae", "Hard.coral"),
                   labels = c("water velocity", "slope", "rugosity", "sedimentation", 
                              "% consolidated", "% free space", "% macroalgal cover",
                              "% hard coral cover"))+
  scale_y_discrete(limits = c("M. turtlensis", "M. aequituberculata", "G. fascicularis",
                              "A. spathulata", "A. muricata", "A. millepora",
                              "A. loripes", "A. kenti", "A. tersa", 
                              "A. digitifera")) +
  
  xlab("") + 
  ylab("") + xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text( face = "italic"),
        text = element_text(size = 12)) 
Fig2

#ggsave(paste(Path, "outputs", "Figure2.png", sep = "/"), Fig1, dpi = 300,
# width = 13, height = 10, units = "cm")



# Get sensitivity of R2 to how weight is assigned across environmental variables
# using a 'one factor at a time' approach (OAT).
sens_df = sensitivity(data2, dat, unique(dat$Species))

# Extract the R2 for an even (1/8) weighting for all variables and add it
# as a column
or_df = sens_df[sens_df$var == "original", ]
or_df$original_r2 = or_df$R2
sens_df = sens_df[!sens_df$var == "original", ]
sens_df = sens_df %>% left_join(or_df[ , c("Species", "original_r2")])

# Get maximimum R2
max_df = sens_df %>% group_by(Species) %>%
  summarise("max_R2" = max(R2))
sens_df = sens_df %>% left_join(max_df)


# Optimise the weighting for each species such that the R2 is maximised
weights_df = data.frame("Species" = factor(),
                        "var" = factor(),
                        "weights" = numeric())

for (species in unique(sens_df$Species)) {
  
  MinLS <- nlminb(start = rep(1/8, 8), 
                  objective = optim_scores2, 
                  data = modify_scores(data2, dat, species),
                  species = species,
                  control= list(abs.tol = 1e-1, rel.tol = 1e-1, step.min = 5, step.max = 10),
                  lower = 0)
  print(MinLS)
  params  = MinLS$par
  params = params/sum(params)
  
  weights_df = rbind(weights_df,
                     data.frame("Species" = species,
                                "var" =  c("water.velocity",
                                           "slope",
                                           "rugosity",
                                           "sedimentation",
                                           "consolidated",
                                           "free.space",
                                           "macroalgae",
                                           "hard.coral"),
                                "weights" = params))
  print(species)
}


#write.table(weights_df, paste(Path, "outputs", "weights_m.csv", sep = "/"),
#            sep = ",", row.names = FALSE)


# Plot the sensitivity analysis
FigS4 = ggplot(sens_df, aes(x = weight, y = var, fill = R2/max_R2))+
  geom_tile(col = "white") +
  facet_wrap(~ Species, ncol = 2) + 
  scale_fill_viridis(option = "A", name = expression(paste("Relative ", R^2, sep = " "))) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_discrete(limits = c("water.velocity", "slope", "rugosity", "sedimentation",
                              "consolidated", "free.space", "macroalgae", "hard.coral"),
                   labels = c("water velocity", "slope", "rugosity", "sedimentation", 
                              "% consolidated", "% free space", "% macroalgal cover",
                              "% hard coral cover")) +
  xlab("Score weight") + 
  ylab("") + xlab("") +
  theme_classic() +
  theme(strip.text = element_text( face = "italic"),
        text = element_text(size = 12))  
FigS4

# Plot optimal weighting for each species
FigS7 = ggplot(weights_df, aes(y = Species, fill = var, x = weights)) +
  geom_col(position = "stack",  col = "black") +
  scale_fill_manual(values = c("#999999", "#CC79A7",
                               "#009373", "#F0E442", "#56B4E9", "#0072B2",
                               "#D55E00", "#E69F00"),
                    name = "",
                    limits = c("water.velocity", "slope", "rugosity", "sedimentation",
                               "consolidated", "free.space", "macroalgae", "hard.coral"),
                    labels = c("water velocity", "slope", "rugosity", "sedimentation", 
                               "% consolidated", "% free space", "% macroalgal cover",
                               "% hard coral cover"))  +
  theme_classic() +
  theme(axis.text.y = element_text( face = "italic"),
        text = element_text(size = 12)) +
  ylab("") + xlab("") +
  xlab("Score weight")
FigS7



# Fit a relationship between site score and survival using the optimal 
# weighting for each species and plot the relationships
Amil_s2 = fit_score2(data2, dat, "A. millepora", weights_df[weights_df$Species == "A. millepora", ]$weights)
Ater_s2 = fit_score2(data2, dat, "A. tersa", weights_df[weights_df$Species == "A. tersa", ]$weights)
Aspa_s2 = fit_score2(data2, dat, "A. spathulata", weights_df[weights_df$Species == "A. spathulata", ]$weights)
Adig_s2 = fit_score2(data2, dat, "A. digitifera", weights_df[weights_df$Species == "A. digitifera", ]$weights)
Aken_s2 = fit_score2(data2, dat, "A. kenti",  weights_df[weights_df$Species == "A. kenti", ]$weights)
Alor_s2 = fit_score2(data2, dat, "A. loripes", weights_df[weights_df$Species == "A. loripes", ]$weights)
Amur_s2 = fit_score2(data2, dat, "A. muricata", weights_df[weights_df$Species == "A. muricata", ]$weights)
Gfas_s2 = fit_score2(data2, dat, "G. fascicularis", weights_df[weights_df$Species == "G. fascicularis", ]$weights)
Maeq_s2 = fit_score2(data2, dat, "M. aequituberculata", weights_df[weights_df$Species == "M. aequituberculata", ]$weights)
Mtur_s2 = fit_score2(data2, dat, "M. turtlensis", weights_df[weights_df$Species == "M. turtlensis", ]$weights)

Fig3 = ((Adig_s2$score.plot + xlab("") + ylab("") + xlab("") + ggtitle(expression(paste("a- ", italic("A. digitifera"))))) + 
          (Ater_s2$score.plot+ xlab("") + ylab("") + xlab("")  + ggtitle(expression(paste("b- ", italic("A. tersa")))))) / 
  ((Aken_s2$score.plot+ xlab("") + ylab("") + xlab("")  + ggtitle(expression(paste("c- ", italic("A. kenti"))))) +
     (Alor_s2$score.plot+ xlab("") + ylab("") + xlab("")  + ggtitle(expression(paste("d- ", italic("A. loripes")))))) / 
  ((Amil_s2$score.plot + xlab("") + ylab("") + xlab("") + ggtitle(expression(paste("e- ", italic("A. millepora"))))) +
     (Amur_s2$score.plot + xlab("") + ylab("") + xlab("") + ggtitle(expression(paste("f- ", italic("A. muricata")))))) / 
  ((Aspa_s2$score.plot + xlab("") + ylab("") + xlab("") + ggtitle(expression(paste("g- ", italic("A. spathulata"))))) + 
     (Gfas_s2$score.plot + xlab("") + ylab("") + xlab("") + ggtitle(expression(paste("h- ", italic("G. fascicularis")))))) / 
  ((Maeq_s2$score.plot + xlab("") + ylab("") + xlab("") + ggtitle(expression(paste("i- ", italic("M. aequituberculata"))))) +
     (Mtur_s2$score.plot + xlab("") + ylab("") + xlab("") + ggtitle(expression(paste("j- ", italic("M. turtlensis"))))))
Fig3



# Check and plot correlations between variables
site = data2[ , c("Site", "Reef","Water.velocity", "Slope", "Rugosity",
                  "Sedimentation", "Consolidated.substrate", "Free.space", "Macroalgae",
                  "Hard.coral")]
site = site[!duplicated(site), ]

vars = c("Water.velocity", "Slope", "Rugosity", "Sedimentation", "Consolidated.substrate",
         "Free.space", "Macroalgae", "Hard.coral")
plot_list = list()
vars1 = c()
vars2 = c()
for (i in 1:(length(vars) -1)) {
  site$var1 = site[, which(colnames(site) == vars[i])]
  
  for (j in (i+1):length(vars)) {
    site$var2 = site[, which(colnames(site) == vars[j])]
    
    vars1 = c(vars1, var[i])
    vars2 = c(vars2, var[j])
    p = ggplot(site, aes(x = var1, y = var2, fill = Reef)) +
      geom_point(pch = 21) + 
      theme_bw() +
      xlab("") + 
      ylab("") + xlab("") + 
      theme(legend.position = "none")
    p
    
    
    if (vars[j] %in% c("Rugosity") ) {
      p = ggplot(site, aes(x = var1, y = var2, fill = Reef)) +
        geom_boxploth() + 
        theme_bw() +
        xlab("") + 
        ylab("") + xlab("") + 
        theme(legend.position = "none") +
        scale_y_discrete(limits  = c("Low", "Medium", "High"),
                         labels = c("L", "M", "H"))
      p
    }
    
    if (vars[j] %in% c("Sedimentation") ) {
      p = ggplot(site, aes(x = var1, y = var2, fill = Reef)) +
        geom_boxploth() + 
        theme_bw() +
        xlab("") + 
        ylab("") + xlab("") + 
        theme(legend.position = "none") +
        scale_y_discrete(limits = c( "There are no sediments deposited" ,
                                     "Thin layer easy to resuspend",
                                     "Moderate layer that offers resistance to be resuspended",
                                     "Thick layer deep layer that is not possible to resuspend"),
                         labels = c(0, 1, 2, 3))
      p
    }
    
    if (vars[i] %in% c("Rugosity") & vars[j] != "Sedimentation") {
      p = ggplot(site, aes(x = var1, y = var2, fill = Reef)) +
        geom_boxplot() + 
        theme_bw() +
        xlab("") + 
        ylab("") + xlab("") + 
        scale_x_discrete(limits  = c("Low", "Medium", "High"),
                         labels = c("L", "M", "H")) +
        theme(legend.position = "none")
      p
    }
    
    if (vars[i] %in% c("Sedimentation") ) {
      p = ggplot(site, aes(x = var1, y = var2, fill = Reef)) +
        geom_boxplot() + 
        theme_bw() +
        xlab("") + 
        ylab("") + xlab("") + 
        scale_x_discrete(limits  = c( "There are no sediments deposited" ,
                                      "Thin layer easy to resuspend",
                                      "Moderate layer that offers resistance to be resuspended",
                                      "Thick layer deep layer that is not possible to resuspend"),
                         labels = c(0, 1, 2, 3)) +
        theme(legend.position = "none")
      p
    }
    if (vars[i] == "Rugosity" & vars[j] == "Sedimentation") {
      p = ggplot(site, aes(x = var1, y = var2, fill = Reef)) +
        geom_count(pch = 21, alpha = 0.5) + 
        theme_bw() +
        xlab("") + 
        ylab("") + xlab("") + 
        scale_x_discrete(limits  = c("Low", "Medium", "High"),
                         labels = c("L", "M", "H")) +
        scale_y_discrete(limits = c( "There are no sediments deposited" ,
                                     "Thin layer easy to resuspend",
                                     "Moderate layer that offers resistance to be resuspended",
                                     "Thick layer deep layer that is not possible to resuspend"),
                         labels = c(0, 1, 2, 3)) +
        theme(legend.position = "none")
      p
    }
    n = ifelse(i == 1, j -1,
               ifelse(i ==2, j + 5,
                      ifelse(i == 3, j +10,
                             ifelse(i == 4, j + 14,
                                    ifelse(i == 5, j + 17,
                                           ifelse(i == 6, j + 19,j + 20))))))
    plot_list[[n]] = p
    
    
  }
  
}
blank <- ggplot() + theme_void()

# Correlation between % free consolidated space and % consolidated space
cor.test(site$Free.space, site$Consolidated.substrate, method = "pearson")

FigS3 = ggarrange(
  #slope 
  plot_list[[1]], blank, blank, blank, blank, blank, blank, blank,
  
  # rugosity
  plot_list[[2]], plot_list[[8]], blank, blank, blank, blank, blank, blank,
  
  # sedimentation
  plot_list[[3]], plot_list[[9]], plot_list[[14]], blank, blank, blank, blank, blank,
  
  # consolidated
  plot_list[[4]], plot_list[[10]], plot_list[[15]], plot_list[[19]], 
  blank, blank, blank, blank,
  
  # free
  plot_list[[5]], plot_list[[11]], plot_list[[16]], plot_list[[20]], 
  plot_list[[23]], blank, blank, blank,
  
  # macroalgae
  plot_list[[6]], plot_list[[12]], plot_list[[17]], plot_list[[21]], 
  plot_list[[24]], plot_list[[26]], blank, blank,
  
  # hard coral
  plot_list[[7]], plot_list[[13]], plot_list[[18]], plot_list[[22]], 
  plot_list[[25]], plot_list[[27]], plot_list[[28]], blank,
  
  ncol = 8, nrow = 7)



# Plot survival by reef and species
surv = data.frame("Researcher" = factor(),
                  "Species" = factor(), 
                  "Reef" = factor(),
                  "Site" = factor(),
                  "Deployment.year" = integer(),
                  "DeviceID" = factor(),
                  "Survival" = integer())

dat = dat[complete.cases(dat$n_surv), ]
for (i in 1:nrow(dat)) {
  
  surv = rbind(surv,
               data.frame("Researcher" = dat$`Lead Researcher`[i],
                          "Species" = dat$Species[i], 
                          "Reef" = dat$Reef[i],
                          "Site" = dat$Site[i],
                          "Deployment.year" = dat$Deployment.year[i],
                          "DeviceID" = dat$DeviceID[i],
                          "Survival" = c(rep(0, dat$n_failures[i]),
                                         rep(1, dat$n_surv[i])) ))
}


surv$Reef = factor(surv$Reef, levels = c("Moore Reef", "Davies Reef", "Keppel Islands", 
                                         "Heron Island"))
FigS8 = ggplot() +
  geom_jitter(data = surv, aes(x = Species, y = Survival, fill = Species), height = 0.15,
              width = 0.25, pch = 21, col = "black", alpha = 0.1) +
  
  geom_point(data = surv %>% group_by(Species, Reef) %>%
               summarise("survival" = mean(Survival),
                         "se" = sqrt(sd(Survival)/n())),
             aes(x = Species, y = survival, fill = Species), pch = 21, col = "black", size = 2) + 
  geom_errorbar(data = surv %>% group_by(Species, Reef) %>%
                  summarise("survival" = mean(Survival),
                            "se" = sqrt(sd(Survival)/n())),
                aes(x = Species, ymin = survival - se, ymax = survival + se,),
                width = 0.2) + 
  theme_bw() + 
  facet_wrap(~ Reef, ncol = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = colour_p) +
  scale_colour_manual(values = colour_p) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_rect(fill = "white")) +
  ylab("Spat survival (mean +/- s.e.)") +
  xlab("") 
FigS8




# Simulation exercise to estimate the relative increase in mean survival when
# using site scores to select sites
predict_df = score_predict_sites(data2,  200, weights_df, Eff)
predict_df[predict_df$survival_random == 0, ]$survival_random = 0.01

predict_df$Species_Reef = paste(predict_df$Species, predict_df$Reef, sep = "-")


percentage = data.frame("Species_Reef" = unique(predict_df$Species_Reef))
percentage$Species = NA
percentage$Reef = NA
percentage$improved = NA
percentage$improved_m = NA
percentage$median = NA
percentage$upr = NA
percentage$lwr = NA
percentage$median_m = NA
percentage$upr_m = NA
percentage$lwr_m = NA



for (i in 1:nrow(percentage)) {
  SUB = predict_df[predict_df$Species_Reef == percentage$Species_Reef[i],]
  Ratio =   SUB$survival_original /SUB$survival_random  
  Ratio = Ratio[complete.cases(Ratio)]
  if (length(which(Ratio == "Inf")) > 0) {
    Ratio = Ratio[-which(Ratio == "Inf")]
  }
  
  Ratio2 =   SUB$survival_modified/SUB$survival_random 
  percentage$Species[i] = unlist(strsplit(percentage$Species_Reef[i], "-"))[1]
  percentage$Reef[i] = unlist(strsplit(percentage$Species_Reef[i], "-"))[2]
  percentage$improved[i] = round(length(which(Ratio > 1)) / length(Ratio), 2)
  percentage$improved_m[i] = round(length(which(Ratio2 > 1)) / length(Ratio2), 2)
  percentage$median[i] = median(Ratio, na.rm = TRUE)
  percentage$upr[i] = quantile(Ratio, 0.975, na.rm = TRUE)
  percentage$lwr[i] = quantile(Ratio, 0.025, na.rm = TRUE)
  percentage$median_m[i] = quantile(Ratio2, 0.5, na.rm = TRUE)
  percentage$upr_m[i] = quantile(Ratio2, 0.975, na.rm = TRUE)
  percentage$lwr_m[i] = quantile(Ratio2, 0.025, na.rm = TRUE)

  
}

predict_df$Species_Reef = paste(predict_df$Species, predict_df$Reef, sep = "-")


percentage = percentage %>%
  mutate(Species_Reef = factor(Species_Reef, levels = c(
    "A. kenti-Heron Island" ,
    "A. millepora-Moore Reef",
    "A. spathulata-Davies Reef",
    "A. tersa-Davies Reef" ,
    "A. tersa-Heron Island",
    "A. digitifera-Moore Reef",
    "M. turtlensis-Davies Reef" ,
    "G. fascicularis-Davies Reef",
    "A. millepora-Keppel Islands",
    "A. muricata-Keppel Islands" ,
    "A. loripes-Davies Reef",
    "M. aequituberculata-Keppel Islands"
  )))


predict_df = predict_df %>%
  mutate(Species_Reef = factor(Species_Reef, levels = c(
    "A. kenti-Heron Island" ,
    "A. millepora-Moore Reef",
    "A. spathulata-Davies Reef",
    "A. tersa-Davies Reef" ,
    "A. tersa-Heron Island",
    "A. digitifera-Moore Reef",
    "M. turtlensis-Davies Reef" ,
    "G. fascicularis-Davies Reef",
    "A. millepora-Keppel Islands",
    "A. muricata-Keppel Islands" ,
    "A. loripes-Davies Reef",
    "M. aequituberculata-Keppel Islands"
    
  )))

percentage$x1 = percentage$upr + 0.2
percentage$x2 = percentage$upr_m + 0.2
percentage[percentage$Species == "A. digitifera", ]$x1 = 1.4
percentage[percentage$Species == "A. digitifera", ]$x2 = 0.85
percentage[percentage$Species == "A. loripes", ]$x1 = 1.3
percentage[percentage$Species == "A. loripes", ]$x2 = 1.35
percentage[percentage$Species == "A. muricata", ]$x2 = 5
percentage[percentage$Species == "M. turtlensis", ]$x1 = 2.4
percentage[percentage$Species == "M. turtlensis", ]$x2 = 4.5
percentage[percentage$Species == "G. fascicularis", ]$x1 = 14
percentage[percentage$Species == "G. fascicularis", ]$x2 = 0.1
percentage[percentage$Species == "A. spathulata", ]$x2 = 1.9
percentage[percentage$Species_Reef == "A. tersa-Heron Island", ]$x1 = 0.9
percentage[percentage$Species_Reef == "A. millepora-Moore Reef", ]$x1 = 1.9
percentage[percentage$Species_Reef == "A. millepora-Moore Reef", ]$x2 = 1.85
percentage[percentage$Species == "A. kenti", ]$x1 = 2.2
percentage[percentage$Species == "A. kenti", ]$x2 = 1.9
percentage[percentage$Species == "M. aequituberculata", ]$x2 = 1.4
percentage[percentage$Species == "M. turtlensis", ]$x1 = 2.8


labels_sr = c(
  'A. digitifera-Moore Reef' = "A. digitifera",
  'A. kenti-Heron Island' = "A. kenti",
  'A. loripes-Davies Reef' = "A. loripes",
  'A. millepora-Keppel Islands' = "A. millepora-Keppel Islands",
  'A. millepora-Moore Reef' = "A. millepora-Moore Reef",
  'A. muricata-Keppel Islands' = "A. muricata",
  'A. spathulata-Davies Reef' = "A. spathulata",
  'A. tersa-Davies Reef' = "A. tersa-Davies Reef",
  'A. tersa-Heron Island' = "A. tersa-Heron Island",
  'G. fascicularis-Davies Reef' = "G. fascicularis",
  'M. aequituberculata-Keppel Islands' = "M. aequituberculata",
  'M. turtlensis-Davies Reef' = "M. turtlensis"
)

Fig4a = ggplot(predict_df) +
  geom_vline(aes(xintercept = 1),  col = "black", linewidth = 1.1, linetype = "dashed") +
  
  geom_density(aes(x =  survival_original  / survival_random ), fill =  "#4662D7FF", 
               alpha = 0.5) +
  geom_density(aes(x = survival_modified / survival_random ), fill =  "#FABA39FF", 
               alpha = 0.5) +
  geom_errorbarh(data = percentage[!percentage$Reef == "Keppel Islands", ],
                 aes(xmin = lwr, xmax = upr, y = -1), height = 0.5, col = "black") +
  geom_point(data = percentage[!percentage$Reef == "Keppel Islands", ],
             aes(x = median, y = -1), pch = 21, col = "black", fill = "#4662D7FF", size = 3 )+
  geom_errorbarh(data = percentage,
                 aes(xmin = lwr_m, xmax = upr_m, y = -2), height = 0.5, col = "black") +
  geom_point(data = percentage,
             aes(x = median_m, y = -2), pch = 21, col = "black", fill = "#FABA39FF", size = 3 )+
  geom_text(data = percentage[!percentage$Species_Reef == ""], aes(x = x1, y = -1, label = paste(improved*100, "%"))) +
  geom_text(data = percentage, aes(x = x2, y = -2, label = paste(improved_m*100, "%"))) +
  theme_classic() +
  ylab("Density") +
  xlab("Survival ratio") +
  scale_x_sqrt() +
  facet_wrap(~Species_Reef, scale = "free_x", ncol = 3, labeller = as_labeller(labels_sr)) +
  
  theme(strip.text = element_text(face = "italic",
                                  size = 12), strip.background = element_rect(fill = "white"),
        text = element_text(size = 15)) +
  scale_y_continuous(breaks = seq(0, 8, by = 2))
Fig4a





# Bootstrapping approach to quantify the variance in spat survival explained
# by the site scores
R2_df = score_predict(data2,dat, 100)

R2_df %>% 
  group_by(Species_Reef) %>%
  summarise("ratio_random" = mean(survival_random2 / survival_random, na.rm = TRUE),
            "ratio_original" = mean(survival_original / survival_random, na.rm = TRUE),
            "ratio_max" = mean(survival_max / survival_random, na.rm = TRUE))


R2_all = with(R2_df[R2_df$convergence == 0, ], 
              data.frame("Species" = c(Species, Species),
                         "R2" = c(R2_conditional, R2)))
R2_all$type = c(rep("conditional", nrow(R2_all) /2),
                rep("marginal", nrow(R2_all) /2))


Fig4b = ggplot(R2_all, aes(x = Species, y = R2, fill = type))+
  geom_jitter(pch = 21, col = "black",  alpha = 0.3, position = position_dodge2(width = 0.6)) +
  geom_boxplot( outlier.shape = NA, alpha = 0.6) + 
  xlab("") + 
  ylab(expression(R^2)) +
  scale_x_discrete(limits = c("M. turtlensis", "M. aequituberculata", "G. fascicularis",
                              "A. spathulata", "A. muricata", "A. millepora",
                              "A. loripes", "A. kenti", "A. tersa", 
                              "A. digitifera")) +
  #scale_y_log10()+
  scale_fill_manual(values = c("#999933", "#6699CC"),
                    name = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        text = element_text(size = 12)) 
Fig4b


