#####
# PROGRAM CODE FOR ARTICLE: "Assessment of structural changes in pathologies on radiological images following therapeutical interventions"
# (c) Márton Kolossváry, 2021

# DESCRIPTION: Network robustness abalyses
#####


# INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library("data.table"); library("WGCNA")
set.seed(42)
folder_wd    <- "/Your_wd/"
folder_data  <- paste0(folder_wd, "DATA/")
folder_stat  <- paste0(folder_wd, "STAT/")
folder_cache <- paste0(folder_wd, "CACHE/")
folder_img   <- paste0(folder_wd, "IMAGES/")

NAME       <- "ES" #Equally sized radiomic parameters
ABSTINENCE <- "abstinence_pre" #Abstinence column name
first_rad  <- 139 #First radiomic parameter position
d_names    <- c("BAS", "MID", "END") #Names of sublists
BETA       <- 8 #Based on the results

load(paste0(folder_cache, "FINAL_DATA_", NAME, ".RData"))
load(paste0(folder_wd, "consensusTOM-block.1.RData"))

outcomes  <- colnames(d_final_cli$BAS$data)[grep("^ME", colnames(d_final_cli$BAS$data))]
outcomes  <- outcomes[1:length(outcomes)-1] #Remove MEgrey
outcomes  <- gsub("ME", "", outcomes, fixed = TRUE)

group_color <- c("#A6EB99", "#DBDA48", "#EA3AF7", "#F6C3CB", "#E6FEFF",
                 "#73DDD0", "#486AD9", "#0018F5", "#983430", "#000000",
                 "#75FA4C", "#74FBFD", "#FFFFE3", "#EA3323", "#999999",
                 "#D1BA98", "#7C2BC4", "#181B6B", "#BFFC5B", "#EB8677")
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# NETWORK PRESERVATION STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Module preservation statistics =====
module_stats <- modulePreservation(d_final_rad, list(con_net_all$colors, con_net_all$colors, con_net_all$colors),
                                   dataIsExpr = TRUE, corFnc = "bicor", verbose = 3, greyName = "grey",
                                   calculateClusterCoeff = TRUE, includekMEallInSummary = TRUE, referenceNetworks = 3,
                                   savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE, nPermutations = 1000,
                                   permutedStatisticsFile = paste0(folder_cache, "Permut_bas_mid_end_END.RData"))

### Save module statistics =====
module_obs_out     <- module_stats$quality$observed$ref.END$inColumnsAlsoPresentIn.BAS #As all features are present in each timepoint, results are same for non references
module_obs_p_out   <- module_stats$quality$log.pBonf$ref.END$inColumnsAlsoPresentIn.BAS

data.table::fwrite(module_obs_out[c(outcomes, "grey"), ], paste0(folder_stat, "Module_stats_obs_END.csv"), row.names = TRUE)
data.table::fwrite(module_obs_p_out[c(outcomes, "grey"), ], paste0(folder_stat, "Module_stats_p_END.csv"), row.names = TRUE)

module_sep_out     <- module_stats$referenceSeparability$observed$ref.END$inColumnsAlsoPresentIn.BAS #As all features are present in each timepoint, results are same for END
module_sep_p_out   <- module_stats$referenceSeparability$log.pBonf$ref.END$inColumnsAlsoPresentIn.BAS

data.table::fwrite(module_sep_out[c(outcomes, "grey"), ], paste0(folder_stat, "Module_sep_obs_END.csv"), row.names = TRUE)
data.table::fwrite(module_sep_p_out[c(outcomes, "grey"), ], paste0(folder_stat, "Module_sep_p_END.csv"), row.names = TRUE)
# =============================================================================

## Calculate averages for the 3 timepoints =====
module_stats_BAS <- modulePreservation(d_final_rad, list(con_net_all$colors, con_net_all$colors, con_net_all$colors),
                                       dataIsExpr = TRUE, corFnc = "bicor", verbose = 3, greyName = "grey",
                                       calculateClusterCoeff = TRUE, includekMEallInSummary = TRUE, referenceNetworks = 1,
                                       savePermutedStatistics = TRUE, loadPermutedStatistics = TRUE, nPermutations = 1000,
                                       permutedStatisticsFile = paste0(folder_cache, "Permut_bas_mid_end_BAS.RData"))

module_stats_MID <- modulePreservation(d_final_rad, list(con_net_all$colors, con_net_all$colors, con_net_all$colors),
                                       dataIsExpr = TRUE, corFnc = "bicor", verbose = 3, greyName = "grey",
                                       calculateClusterCoeff = TRUE, includekMEallInSummary = TRUE, referenceNetworks = 2,
                                       savePermutedStatistics = TRUE, loadPermutedStatistics = TRUE, nPermutations = 1000,
                                       permutedStatisticsFile = paste0(folder_cache, "Permut_bas_mid_end_MID.RData"))

module_stats_END <- modulePreservation(d_final_rad, list(con_net_all$colors, con_net_all$colors, con_net_all$colors),
                                       dataIsExpr = TRUE, corFnc = "bicor", verbose = 3, greyName = "grey",
                                       calculateClusterCoeff = TRUE, includekMEallInSummary = TRUE, referenceNetworks = 3,
                                       savePermutedStatistics = TRUE, loadPermutedStatistics = TRUE, nPermutations = 1000,
                                       permutedStatisticsFile = paste0(folder_cache, "Permut_bas_mid_end_END.RData"))

module_obs_out_ALL     <- (module_stats_BAS$quality$observed$ref.BAS$inColumnsAlsoPresentIn.MID +
                             module_stats_MID$quality$observed$ref.MID$inColumnsAlsoPresentIn.BAS +
                             module_stats_END$quality$observed$ref.END$inColumnsAlsoPresentIn.BAS) /3

module_obs_p_out_ALL     <- log10(((10^module_stats_BAS$quality$log.pBonf$ref.BAS$inColumnsAlsoPresentIn.MID +
                                      10^module_stats_MID$quality$log.pBonf$ref.MID$inColumnsAlsoPresentIn.BAS +
                                      10^module_stats_END$quality$log.pBonf$ref.END$inColumnsAlsoPresentIn.BAS) /3))

module_sep_out_ALL     <- (module_stats_BAS$referenceSeparability$observed$ref.BAS$inColumnsAlsoPresentIn.MID +
                             module_stats_MID$referenceSeparability$observed$ref.MID$inColumnsAlsoPresentIn.BAS +
                             module_stats_END$referenceSeparability$observed$ref.END$inColumnsAlsoPresentIn.BAS) /3

module_sep_p_out_ALL     <- log10(((10^module_stats_BAS$referenceSeparability$log.pBonf$ref.BAS$inColumnsAlsoPresentIn.MID +
                                      10^module_stats_MID$referenceSeparability$log.pBonf$ref.MID$inColumnsAlsoPresentIn.BAS +
                                      10^module_stats_END$referenceSeparability$log.pBonf$ref.END$inColumnsAlsoPresentIn.BAS) /3))

data.table::fwrite(module_obs_out_ALL[c(outcomes, "grey"), ], paste0(folder_stat, "Module_stat_obs_ALL.csv"), row.names = TRUE)
data.table::fwrite(module_obs_p_out_ALL[c(outcomes, "grey"), ], paste0(folder_stat, "Module_stat_p_ALL.csv"), row.names = TRUE)

data.table::fwrite(module_sep_out_ALL[c(outcomes, "grey"), ], paste0(folder_stat, "Module_sep_obs_ALL.csv"), row.names = TRUE)
data.table::fwrite(module_sep_p_out_ALL[c(outcomes, "grey"), ], paste0(folder_stat, "Module_sep_p_ALL.csv"), row.names = TRUE)
# =============================================================================

## Compare preservation statistics to baseline =====
module_stats_BAS <- modulePreservation(d_final_rad, list(con_net_all$colors, con_net_all$colors, con_net_all$colors),
                                       dataIsExpr = TRUE, corFnc = "bicor", verbose = 3, greyName = "grey",
                                       calculateClusterCoeff = TRUE, includekMEallInSummary = TRUE, referenceNetworks = 1,
                                       savePermutedStatistics = TRUE, loadPermutedStatistics = TRUE, nPermutations = 1000,
                                       permutedStatisticsFile = paste0(folder_cache, "Permut_bas_mid_end_BAS.RData"))


module_pres_out     <- module_stats_BAS$preservation$observed$ref.BAS$inColumnsAlsoPresentIn.MID #As all features are present in each timepoint, results are same for non references
module_pres_p_out   <- module_stats_BAS$preservation$log.pBonf$ref.BAS$inColumnsAlsoPresentIn.MID

data.table::fwrite(module_pres_out[c(outcomes, "grey"), ], paste0(folder_stat, "Module_pres_obs_MID.csv"), row.names = TRUE)
data.table::fwrite(module_pres_p_out[c(outcomes, "grey"), ], paste0(folder_stat, "Module_pres_p_MID.csv"), row.names = TRUE)

module_pres_out     <- module_stats_BAS$preservation$observed$ref.BAS$inColumnsAlsoPresentIn.END #As all features are present in each timepoint, results are same for non references
module_pres_p_out   <- module_stats_BAS$preservation$log.pBonf$ref.BAS$inColumnsAlsoPresentIn.END

data.table::fwrite(module_pres_out[c(outcomes, "grey"), ], paste0(folder_stat, "Module_pres_obs_END.csv"), row.names = TRUE)
data.table::fwrite(module_pres_p_out[c(outcomes, "grey"), ], paste0(folder_stat, "Module_pres_p_END.csv"), row.names = TRUE)
# =============================================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PERCOLATION ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Initialize WGCNA models =====
replaceMissingAdjacencies <- TRUE
power          <- BETA
corType        <- "bicor"
consensusQuantile  <- 0
deepSplit      <- 2
minModuleSize  <- ceiling(dim(d_final_rad$BAS$data)[2]/100)*5 # 5% of all parameters
minModuleSize  <- 20
mergeCutHeight <- 0.15
numericLabels  <- FALSE
pamStage       <- TRUE
saveIndividualTOMs <- TRUE
individualTOMFileNames <- "%s%N"
saveConsensusTOMs  <- FALSE
verbose        <- 0
useDiskCache   <- FALSE
# =============================================================================

## Initialize resampling =====
#Modifications need to be done for patient vs visit based elimination. Also ran on several computers to utilize parallel processing
all_IDs <- unique(d_final_cli$BAS$data$ID, d_final_cli$MID$data$ID, d_final_cli$END$data$ID)
til     <- ceiling(length(all_IDs)/2) #50%
samples <- 1000
batches <- split(1:til, sort(1:length(1:til)%%5)) #5 batches
which_batch <- 1 #change to desired batch
setwd(paste0(folder_scr, "Percolation/Patient_", which_batch, "/"))

c_cor  <- lapply(1:til, function(x) {vector(mode = "list", length = samples)})
acc    <- lapply(1:til, function(x) {vector(mode = "list", length = samples)})
pres   <- lapply(1:til, function(x) {vector(mode = "list", length = samples)})
ME_cor <- lapply(1:til, function(x) {vector(mode = "list", length = samples)})
# =============================================================================

## Run percolation analyses =====
disableWGCNAThreads()
for(i in batches[which_batch][[1]]) {
  t_start <- Sys.time()
  for(j in 1:samples) {
    rnd_ID     <- sample(all_IDs, i) #Select random samples and identity row IDS
    row_ID_BAS <- d_final_cli$BAS$data$ID %in% rnd_ID
    row_ID_MID <- d_final_cli$MID$data$ID %in% rnd_ID
    row_ID_END <- d_final_cli$END$data$ID %in% rnd_ID
    
    d_final_rad_i <- d_final_rad #Remove specific rows
    d_final_rad_i$BAS$data <- d_final_rad_i$BAS$data[!row_ID_BAS, ]
    d_final_rad_i$MID$data <- d_final_rad_i$MID$data[!row_ID_MID, ]
    d_final_rad_i$END$data <- d_final_rad_i$END$data[!row_ID_END, ]
    
    con_net_all_i       <- blockwiseConsensusModules(d_final_rad_i, 
                                                     power = power, replaceMissingAdjacencies = replaceMissingAdjacencies, #networkCalibration = networkCalibration,
                                                     corType = corType, consensusQuantile = consensusQuantile,
                                                     deepSplit = deepSplit, minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, pamStage = pamStage,
                                                     numericLabels = numericLabels, verbose = verbose,
                                                     saveIndividualTOMs = saveIndividualTOMs, saveConsensusTOMs = saveConsensusTOMs)
    c_cor[[i]][[j]] <- dendextend::cor_cophenetic(con_net_all$dendrograms[[1]], con_net_all_i$dendrograms[[1]]) #COPHENETIC CORRELATION
    
    
    modul_stat_i <- modulePreservation(list(d_final_rad$BAS, d_final_rad$MID, d_final_rad$END,
                                            d_final_rad_i$BAS, d_final_rad_i$MID, d_final_rad_i$END),
                                       list(con_net_all$colors, con_net_all$colors, con_net_all$colors,
                                            con_net_all_i$colors, con_net_all_i$colors, con_net_all_i$colors),
                                       dataIsExpr = TRUE, corFnc = corType, verbose = verbose, greyName = "grey",
                                       calculateClusterCoeff = FALSE, includekMEallInSummary = TRUE, calculateQvalue = FALSE,
                                       referenceNetworks = c(1, 2, 3), testNetworks = list(4, 5, 6), nPermutations = 0)
    acc[[i]][[j]]   <- modul_stat_i$observed[[1]]$accuracy[[4]][, 1] #ACCURACY
    pres[[i]][[j]]$BAS  <- modul_stat_i$observed[[1]]$inter$`_vs_` #CONNECTIVITY PRESERVATION
    pres[[i]][[j]]$MID  <- modul_stat_i$observed[[2]]$inter$`_vs_`
    pres[[i]][[j]]$END  <- modul_stat_i$observed[[3]]$inter$`_vs_`
    
    newMEs <- multiSetMEs(d_final_rad_i, universalColors = con_net_all$colors, softPower = BETA, verbose = verbose)
    ME_cor[[i]][[j]] <- sapply(paste0("ME", outcomes), function(x){cor(c(d_final_cli$BAS$data[[x]][!row_ID_BAS],
                                                                         d_final_cli$MID$data[[x]][!row_ID_MID],
                                                                         d_final_cli$END$data[[x]][!row_ID_END]),
                                                                       c(newMEs$BAS$data[[x]],
                                                                         newMEs$MID$data[[x]],
                                                                         newMEs$END$data[[x]]),
                                                                       use = "pairwise.complete.obs")})
  }
  print(paste0("Done with: ", i, "/", length(batches[which_batch][[1]]), " in ", round(difftime(Sys.time(), t_start, units = "mins"), 2), " minutes"))
  save(list = c("c_cor", "acc", "pres", "ME_cor"), file = paste0(folder_scr, "Percolation/Patient_", which_batch, "/", "Percolation_patient_", which_batch, ".RData"))
}
# =============================================================================


# PLOT PRESERVATION STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load data =====
d_cor    <- lapply(1:5, function(x) {
  suppressMessages(attach(paste0(folder_scr, "Percolation/Patient_", x, "/", "Percolation_patient_", x, ".RData")))
  out <- c_cor
  out <- lapply(out, function(y) {y <- y[!sapply(y, is.null)]; y})
  out <- out[lengths(out) != 0]
  detach()
  out
})
d_acc    <- lapply(1:5, function(x) {
  suppressMessages(attach(paste0(folder_scr, "Percolation/Patient_", x, "/", "Percolation_patient_", x, ".RData")))
  out <- acc
  out <- lapply(out, function(y) {y <- y[!sapply(y, is.null)]; y})
  out <- out[lengths(out) != 0]
  detach()
  out
})
d_pres   <- lapply(1:5, function(x) {
  suppressMessages(attach(paste0(folder_scr, "Percolation/Patient_", x, "/", "Percolation_patient_", x, ".RData")))
  out <- pres
  out <- lapply(out, function(y) {y <- y[!sapply(y, is.null)]; y})
  out <- out[lengths(out) != 0]
  detach()
  out
})
d_ME_cor <- lapply(1:5, function(x) {
  suppressMessages(attach(paste0(folder_scr, "Percolation/Patient_", x, "/", "Percolation_patient_", x, ".RData")))
  out <- ME_cor
  out <- lapply(out, function(y) {y <- y[!sapply(y, is.null)]; y})
  out <- out[lengths(out) != 0]
  detach()
  out
})

all_IDs <- unique(d_final_cli$BAS$data$ID, d_final_cli$MID$data$ID, d_final_cli$END$data$ID)
til     <- ceiling(length(all_IDs)/2) #50%
x_value <- 1:til/(length(all_IDs))*100
# =============================================================================

## Cophenetic correlation =====
d_cor      <- unlist(d_cor, recursive = FALSE)
d_plot_cor <- as.data.table(do.call(rbind, lapply(d_cor, function(x) {c(summary(unlist(x)), quantile(unlist(x), c(0.025, 0.975)))})))
d_plot_cor$x <- x_value
colnames(d_plot_cor) <- c("Minimum", "Percentile_25", "Median", "Mean", "Percentile_75", "Maximum", "Percentile_025", "Percentile_975", "x")

### Plot cophenetic correlation - Patient elimination =====
p_pat_coph <- ggplot(d_plot_cor) +
  geom_ribbon(aes(x = x, ymin = Percentile_025, ymax = Percentile_975), fill = ggsci::pal_npg()(1), alpha = 0.3) +
  geom_line(aes(x = x, y = Median, color = ggsci::pal_npg()(1)), size = 1.5) + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  labs(title = "Preservation Of Radiomics Hiearchical Clustering Dendrogram",
       subtitle = "Cophenetic Correlation - Patient-based Random Elimination") +
  xlab("Percentage of removed individuals [%]") +
  scale_y_continuous(name = "Cophenetic correlation [c]", 
                     breaks = seq(0, 1, 0.1), limits = c(0, 1))

ggsave(paste0(folder_img, "Percolation/Cophenetic_patient.svg"), p_pat_coph, device = "svg", width = 12, height = 8, units = "in")
# =============================================================================

## Accuracy =====
d_acc      <- unlist(d_acc, recursive = FALSE)
d_plot_acc <- lapply(d_acc, function(x) {
  out <- data.table::as.data.table(do.call(rbind, x))
  out[, grey := NULL]
  out <- out[, ..outcomes]
  out <- out*100
  as.data.table(t(apply(out, 2, function(y) { c(summary(y), quantile(y, c(0.025, 0.975)))}))) })
d_plot_acc <- lapply(1:til, function(x) {
  d_plot_acc[[x]]$x      <- x_value[x]
  d_plot_acc[[x]]$colors <- outcomes
  d_plot_acc[[x]]$hex    <- group_color
  d_plot_acc[[x]]$alpha  <- c(rep(0.3, 3), 1, rep(0.3, 7), 1, rep(0.3, 8))
  d_plot_acc[[x]] })
d_plot_acc <- data.table::rbindlist(d_plot_acc)
d_plot_acc$colors <- factor(d_plot_acc$colors, ordered = TRUE, labels = outcomes)
colnames(d_plot_acc) <- c("Minimum", "Percentile_25", "Median", "Mean", "Percentile_75", "Maximum", "Percentile_025", "Percentile_975", "x", "colors", "hex", "alpha")

### Plot accuracy - Patient elimination =====
p_pat_acc <- ggplot(d_plot_acc) +
  geom_line(aes(x = x, y = Median, color = colors, group = colors), size = 1.5,
            alpha = c(rep(0.2, 3*34), rep(1, 34), rep(0.2, 7*34), rep(1, 34), rep(0.2, 8*34)),
            linetype = c(rep(2, 3*34), rep(1, 34), rep(2, 7*34), rep(1, 34), rep(2, 8*34))) +
  scale_color_manual(values = d_plot_acc$hex) +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  labs(title = "Preservation Of Feature Module Assignment",
       subtitle = "Accuracy - Patient-based Random Elimination") +
  xlab("Percentage of removed individuals [%]") +
  scale_y_continuous(name = "Overlap With Original Module [%]", 
                     breaks = seq(0, 100, 10), limits = c(0, 100))

ggsave(paste0(folder_img, "Percolation/Accuracy_patient.svg"), p_pat_acc, device = "svg", width = 12, height = 8, units = "in")
# =============================================================================

## Preservation statistics =====
### Initialize preservation statistics =====
d_pres      <- unlist(d_pres, recursive = FALSE)

#### Function: format preservation statistics to plot =====
get_pres <- function(l, time = "BAS", type = "bicor.kIM", n = 1:til) {
  l_out <- lapply(n, function(x) {
    out <- l[[x]]
    out <- lapply(out, function(y) y[[time]]) #Select time
    out <- lapply(out, function(y) y[!row.names(y) %in% c("grey", "gold"), ]) #Remove unwanted modules
    out <- lapply(out, function(y) y[outcomes, ]) #Reorder
    out <- lapply(out, function(y) y[[type]]) #Select statistic
    out <- data.table::as.data.table(do.call(rbind, out)) #Summary
    colnames(out) <- outcomes
    out <- as.data.table(t(apply(out, 2, function(y) { c(summary(y), quantile(y, c(0.025, 0.975)))})))
    out$x      <- x_value[x]
    out$colors <- ordered(outcomes, levels = outcomes)
    out$hex    <- group_color
    out$alpha  <- c(rep(0.3, 3), 1, rep(0.3, 7), 1, rep(0.3, 8))
    out})
  l_out <- data.table::rbindlist(l_out)
  colnames(l_out) <- c("Minimum", "Percentile_25", "Median", "Mean", "Percentile_75", "Maximum", "Percentile_025", "Percentile_975", "x", "colors", "hex", "alpha")
  l_out
}

#### Function: plot preservation statistics =====
plot_pres <- function(dt, type = "Intramodular Connectivity Correlation", folder_img = folder_img) {
  p <- ggplot(dt) +
    geom_line(aes(x = x, y = Median, color = colors, group = colors), size = 1.5,
              alpha = c(rep(0.2, 3*34), rep(1, 34), rep(0.2, 7*34), rep(1, 34), rep(0.2, 8*34)),
              linetype = c(rep(2, 3*34), rep(1, 34), rep(2, 7*34), rep(1, 34), rep(2, 8*34))) +
    scale_color_manual(values = dt$hex) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(size = 24),
          plot.subtitle = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16)) +
    labs(title = paste0("Preservation Of Radiomic Feature Modules"),
         subtitle = paste0(type, " - Patient-based Random Elimination")) +
    xlab("Percentage of removed individuals [%]") +
    scale_y_continuous(name = paste0(type, " [c]"), 
                       breaks = seq(0, 1, 0.1), limits = c(0, 1))
  p
}
# =============================================================================

## Get preservation statistics - kIM =====
d_plot_pres_BAS <- get_pres(d_pres, time = "BAS", type = "bicor.kIM")
d_plot_pres_MID <- get_pres(d_pres, time = "MID", type = "bicor.kIM")
d_plot_pres_END <- get_pres(d_pres, time = "END", type = "bicor.kIM")

### Plot preservation statistics =====
p_pres_BAS_kIM <- plot_pres(d_plot_pres_BAS, type = "Intramodular Connectivity Correlation")
ggsave(paste0(folder_img, "Percolation/Preservation_patient_BAS_", "IntramodularConnectivityCorrelation", ".svg"), p_pres_BAS_kIM, device = "svg", width = 12, height = 8, units = "in")
p_pres_MID_kIM <- plot_pres(d_plot_pres_MID, type = "Intramodular Connectivity Correlation")
ggsave(paste0(folder_img, "Percolation/Preservation_patient_MID_", "IntramodularConnectivityCorrelation", ".svg"), p_pres_BAS_kIM, device = "svg", width = 12, height = 8, units = "in")
p_pres_END_kIM <- plot_pres(d_plot_pres_END, type = "Intramodular Connectivity Correlation")
ggsave(paste0(folder_img, "Percolation/Preservation_patient_END_", "IntramodularConnectivityCorrelation", ".svg"), p_pres_BAS_kIM, device = "svg", width = 12, height = 8, units = "in")
# =============================================================================

## Get preservation statistics - kME -----
d_plot_pres_BAS <- get_pres(d_pres, time = "BAS", type = "bicor.kME")
d_plot_pres_MID <- get_pres(d_pres, time = "MID", type = "bicor.kME")
d_plot_pres_END <- get_pres(d_pres, time = "END", type = "bicor.kME")

### Plot preservation statistics -----
p_pres_BAS_kME <- plot_pres(d_plot_pres_BAS, type = "Module Membership Correlation")
ggsave(paste0(folder_img, "Percolation/Preservation_patient_BAS_", "ModuleMembership", ".svg"), p_pres_BAS_kME, device = "svg", width = 12, height = 8, units = "in")
p_pres_MID_kME <- plot_pres(d_plot_pres_MID, type = "Module Membership Correlation")
ggsave(paste0(folder_img, "Percolation/Preservation_patient_MID_", "ModuleMembership", ".svg"), p_pres_MID_kME, device = "svg", width = 12, height = 8, units = "in")
p_pres_END_kME <- plot_pres(d_plot_pres_END, type = "Module Membership Correlation")
ggsave(paste0(folder_img, "Percolation/Preservation_patient_END_", "ModuleMembership", ".svg"), p_pres_END_kME, device = "svg", width = 12, height = 8, units = "in")
# =============================================================================

## Get preservation statistics - Maximum Adjacency Ratio -----
d_plot_pres_BAS <- get_pres(d_pres, time = "BAS", type = "bicor.MAR")
d_plot_pres_MID <- get_pres(d_pres, time = "MID", type = "bicor.MAR")
d_plot_pres_END <- get_pres(d_pres, time = "END", type = "bicor.MAR")

### Plot preservation statistics -----
p_pres_BAS_MAR <- plot_pres(d_plot_pres_BAS, type = "Maximum Adjacency Ratio Correlation")
ggsave(paste0(folder_img, "Percolation/Preservation_patient_BAS_", "MAR", ".svg"), p_pres_BAS_MAR, device = "svg", width = 12, height = 8, units = "in")
p_pres_MID_MAR <- plot_pres(d_plot_pres_MID, type = "Maximum Adjacency Ratio Correlation")
ggsave(paste0(folder_img, "Percolation/Preservation_patient_MID_", "MAR", ".svg"), p_pres_MID_MAR, device = "svg", width = 12, height = 8, units = "in")
p_pres_END_MAR <- plot_pres(d_plot_pres_END, type = "Maximum Adjacency Ratio Correlation")
ggsave(paste0(folder_img, "Percolation/Preservation_patient_END_", "MAR", ".svg"), p_pres_END_MAR, device = "svg", width = 12, height = 8, units = "in")
# =============================================================================

## ME correlation =====
d_ME_cor      <- unlist(d_ME_cor, recursive = FALSE)
d_plot_ME_cor <- lapply(d_ME_cor, function(x) {
  out <- data.table::as.data.table(do.call(rbind, x))
  colnames(out) <- outcomes
  out <- out[, ..outcomes]
  as.data.table(t(apply(out, 2, function(y) { c(summary(y), quantile(y, c(0.025, 0.975)))}))) })
d_plot_ME_cor <- lapply(1:til, function(x) {
  d_plot_ME_cor[[x]]$x      <- x_value[x]
  d_plot_ME_cor[[x]]$colors <- ordered(outcomes, levels = outcomes)
  d_plot_ME_cor[[x]]$hex    <- group_color
  d_plot_ME_cor[[x]]$alpha  <- c(rep(0.3, 3), 1, rep(0.3, 7), 1, rep(0.3, 8))
  d_plot_ME_cor[[x]] })
d_plot_ME_cor <- data.table::rbindlist(d_plot_ME_cor)
colnames(d_plot_ME_cor) <- c("Minimum", "Percentile_25", "Median", "Mean", "Percentile_75", "Maximum", "Percentile_025", "Percentile_975", "x", "colors", "hex", "alpha")

### Plot ME correlation - Patient elimination =====
p_pat_ME_cor <- ggplot(d_plot_ME_cor) +
  geom_line(aes(x = x, y = Median, color = colors, group = colors), size = 1.5,
            alpha = c(rep(0.2, 3*34), rep(1, 34), rep(0.2, 7*34), rep(1, 34), rep(0.2, 8*34)),
            linetype = c(rep(2, 3*34), rep(1, 34), rep(2, 7*34), rep(1, 34), rep(2, 8*34))) +
  scale_color_manual(values = d_plot_ME_cor$hex) +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  labs(title = "Preservation Of Consensus Eigen Radiomic Features",
       subtitle = "Correlation Between Consensus Eigen Radiomic Features - Patient-based Random Elimination") +
  xlab("Percentage of removed individuals [%]") +
  scale_y_continuous(name = "Pearson Correlation [c]", 
                     breaks = seq(0, 1, 0.001), limits = c(0.99, 1))

ggsave(paste0(folder_img, "Percolation/ME_correlation_patient.svg"), p_pat_ME_cor, device = "svg", width = 12, height = 8, units = "in")
# =============================================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# EXIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list = ls())
collectGarbage()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
