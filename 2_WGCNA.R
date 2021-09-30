#####
# PROGRAM CODE FOR ARTICLE: "Assessment of structural changes in pathologies on radiological images following therapeutical interventions"
# (c) Márton Kolossváry, 2021

# DESCRIPTION: Code to run WGCNA analyses
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
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# CALCULATE SCALE-FREE BETAs AND PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Initialize parameters =====
d_for_plots <- list(d_final_rad, d_final_rad_abs, d_final_rad_absNO)
powers      <- c(seq(2, 14, by=1))
colors      <- ggsci::pal_npg()(length(d_names))
plotCols    <- c(2,5,6,7,9,10)
colNames    <- c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
                 "Max connectivity", "Network Density", "Network Heterogeneity")
listNames   <- c("all", "abs", "absNO")
# =============================================================================

## Plot diagnostic curves =====
ind <- 1
for(d in d_for_plots) {
  ## Calculate powers -----
  powerTables <- vector(mode = "list", length = length(d)); names(powerTables) <- d_names
  powerTables <- lapply(1:length(d), function(x) {
    pickSoftThreshold(d[[x]]$data, powerVector = powers, verbose = 2, moreNetworkConcepts = TRUE, removeFirst = TRUE, corFnc = "bicor")})
  
  ## Get plot y ranges -----
  ylim <-  matrix(NA, nrow = 2, ncol = length(plotCols))
  for (set in 1:length(d)) {
    for (col in 1:length(plotCols)) {
      ylim[1, col] = min(ylim[1, col], powerTables[[set]]$fitIndices[, plotCols[col]], na.rm = TRUE)
      ylim[2, col] = max(ylim[2, col], powerTables[[set]]$fitIndices[, plotCols[col]], na.rm = TRUE)
    }
  }
  
  ## Plot ----
  pdf(file = paste0(folder_img, "Network_attributes_", NAME, "/", "Network_power_", NAME, "_", listNames[ind], ".pdf"), width = 12, height = 6)
  par(mfcol = c(2,3)); par(mar = c(4.2, 4.2 , 2.2, 4.5), xpd = TRUE); cex1 = 0.7
  
  for (col in 1:(length(plotCols))) {
    for (set in 1:length(d)) {
      if (set==1) {
        plot(powerTables[[set]]$fitIndices[,1], (1*powerTables[[set]]$fitIndices[,3]<0)*powerTables[[set]]$fitIndices[,2],
             xlab="Soft Threshold (power)", ylab = colNames[col], type = "n", ylim = ylim[, col], main = colNames[col], xaxt = "n")
        axis(side = 1, at=c(seq(2, 14, 2)))
        addGrid() }
      if (col==1) {
        lines(powerTables[[set]]$fitIndices[,1], (1*powerTables[[set]]$fitIndices[,3]<0)*powerTables[[set]]$fitIndices[,2],
              type = "l", col = colors[set])
      } else {
        lines(powerTables[[set]]$fitIndices[,1], powerTables[[set]]$fitIndices[, plotCols[col]],
              type = "l", col = colors[set])
      }
      legend("topright", legend = d_names, col = colors, pch = 20, inset=c(-0.2,0))
    }
  }
  dev.off()
  ind <- ind+1
}
# =============================================================================

## Plot scale-free plots for guven power (beta = 8) =====
pdf(file = paste0(folder_img, "Network_attributes_", NAME, "/", "Network_plot_all.pdf"), width = 12, height = 3)
par(mfcol = c(1,3)); par(mar = c(4.2, 4.2 , 2.2, 4.5), xpd = TRUE); cex1 = 0.7
scaleFreePlot(connectivity = softConnectivity(datExpr = d_for_plots[[1]]$BAS$data, power = BETA), removeFirst = TRUE, truncated = TRUE, corFnc = "bicor")
scaleFreePlot(connectivity = softConnectivity(datExpr = d_for_plots[[1]]$MID$data, power = BETA), removeFirst = TRUE, truncated = TRUE, corFnc = "bicor")
scaleFreePlot(connectivity = softConnectivity(datExpr = d_for_plots[[1]]$END$data, power = BETA), removeFirst = TRUE, truncated = TRUE, corFnc = "bicor")
dev.off()

rm(d_for_plots)
# =============================================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# CREATE WGCNA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Initialize parameters that are the defaults of the function =====
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
saveConsensusTOMs  <- TRUE
verbose        <- 5
# =============================================================================

## Run WGCNA =====
con_net_all       <- blockwiseConsensusModules(d_final_rad, 
                                               power = power, replaceMissingAdjacencies = replaceMissingAdjacencies, #networkCalibration = networkCalibration,
                                               corType = corType, consensusQuantile = consensusQuantile,
                                               deepSplit = deepSplit, minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, pamStage = pamStage,
                                               numericLabels = numericLabels, verbose = verbose,
                                               saveIndividualTOMs = saveIndividualTOMs, saveConsensusTOMs = saveConsensusTOMs)
names(con_net_all$multiMEs) <- names(d_final_rad)
# =============================================================================

## Cbind eigen features to clinical data =====
d_final_cli <- lapply(1:length(d_final_cli), function(x){ cbind(d_final_cli[[x]]$data, con_net_all$multiMEs[[x]]$data)})
d_final_cli <- lapply(d_final_cli, function(x) list(data = x))
names(d_final_cli)  <- d_names
# =============================================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PLOT WGCNA RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Consensus network dendrogram =====
d_for_plots <- list(con_net_all)
p_names <- c("ALL")
ind <- 1

for(d in d_for_plots) {
  pdf(file = paste0(folder_img, "Network_attributes_", NAME, "/", "Consensus_", NAME, "_", p_names[ind], ".pdf"),width = 12, height = 6)
  plotDendroAndColors(d$dendrograms[[1]], d$colors,
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Radiomics dendrogram and module colors")
  dev.off()
  ind <- ind+1
}
rm(list = c("d_for_plots", "p_names"))
# =============================================================================

## Consensus network dendrogram with heatmap =====
load(paste0(folder_wd, "consensusTOM-block.1.RData"));   consTomDS <- as.matrix(consTomDS)

png(paste0(folder_img, "Network_attributes_", NAME, "/", "Consensus_", NAME, "_heatmap_consensus.png"), units="in", width=9, height=9, res=900)
p_TOMcon <- TOMplot(dissim = consTomDS, dendro = con_net_all$dendrograms[[1]], Colors = con_net_all$colors, terrainColors = FALSE)
dev.off()
# =============================================================================

## MDS plot of radiomic features according to module correspondence ===========
library(plotly)
load(paste0(folder_wd, "consensusTOM-block.1.RData"))

mds <- cmdscale(1-consTomDS, k = 3)
module_colors <- factor(con_net_all$colors, levels = unique(con_net_all$colors))
mds <- cbind(as.data.frame(mds), module_colors); colnames(mds) <- c("x", "y", "z", "color")

l <- list(
  font = list(
    family = "helvetica",
    size = 14),
  title=list(text='<b> Radiomic Eigen Modules </b>'),
  itemsizing = list(
    values = "constant",
    itemwidth = 50),
  orientation = 'h'
)
s <- list(camera = list(eye = list(x = -1.7, y = 1.4, z = 0.2)),
          xaxis = list(title = "Scaling Dimension 1",
                       font = list(
                         family = "helvetica",
                         size = 12),
                       gridcolor = toRGB("black"),
                       gridwidth = 4),
          yaxis = list(title = "Scaling Dimension 2",
                       font = list(
                         family = "helvetica",
                         size = 12),
                       gridcolor = toRGB("black"),
                       gridwidth = 4),
          zaxis = list(title = "Scaling Dimension 3",
                       font = list(
                         family = "helvetica",
                         size = 12),
                       gridcolor = toRGB("black"),
                       gridwidth = 4))
m <- list(
  l = 0,
  r = 0,
  b = 0,
  t = 0,
  pad = 1
)
fig <- plotly::plot_ly(type = "scatter3d", x = mds$x, y = mds$y, z = mds$z,
                       marker = list(size = 8, opacity = 1),
                       color = mds$color, colors = gplots::col2hex(unique(mds$color)))
fig <- fig %>% layout(legend = NULL, scene = s, margin = m)
orca(fig, paste0("MDS_plot_consensus_", NAME, ".png"), format = "png", width = 6*900, height = 6*900)
# =============================================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# EXIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list = ls())
collectGarbage()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
