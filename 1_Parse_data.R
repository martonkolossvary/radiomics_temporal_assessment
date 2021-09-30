#####
# PROGRAM CODE FOR ARTICLE: "Assessment of structural changes in pathologies on radiological images following therapeutical interventions"
# (c) Márton Kolossváry, 2021

# DESCRIPTION: Code to parse data for WGCNA analyses
#####


##### DATASET DESCRIPTION #####
# Databases -----
# d: list of 3 data.tables, containing all data from baseline ($BAS), 6 months ($MID) and 12 months ($END)

# Variables -----
# abstinence_pre: binary 1/0: Did the individual achieve abstinence?
# Age_baseline: numeric: baseline age of individual
# Male: binary 1/0: Is the individual male?
# newrisk_baseline: numeric: ASCVD risk score of the individual at given timepoint
# LANCP: numeric: low-attenuation noncalcified plaque volume at given timepoint
# NCP: numeric: noncalcified plaque volume at given timepoint
# CP: numeric: calcified plaque volume at given timepoint
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
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PARSE DATABASES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Create databases =====
d_cli <- lapply(d, function(x){ x[, 1:(first_rad-1)] })
d_rad <- lapply(d, function(x){ x[, ..rad_params] })
rm(d)
# =============================================================================

## Normalize radiomic parameters =====
NORMALIZE <- TRUE
if(NORMALIZE) {
  d_rad <- lapply(d_rad, function(x){ x[, (rad_params) := lapply(.SD, DescTools::Winsorize, na.rm = TRUE), .SDcols = rad_params] })
  d_rad <- lapply(d_rad, function(x){ x[, (rad_params) := lapply(.SD, scale), .SDcols = rad_params] })
}
# =============================================================================

## Create abstinence subgroups =====
d_rad_abs    <- lapply(1:length(d_names), function(x) d_rad[[x]][as.vector(d_cli[[x]][, ..ABSTINENCE] == 1), ])
d_rad_absNO  <- lapply(1:length(d_names), function(x) d_rad[[x]][as.vector(d_cli[[x]][, ..ABSTINENCE] == 0), ])
d_cli_abs    <- lapply(1:length(d_names), function(x) d_cli[[x]][as.vector(d_cli[[x]][, ..ABSTINENCE] == 1), ])
d_cli_absNO  <- lapply(1:length(d_names), function(x) d_cli[[x]][as.vector(d_cli[[x]][, ..ABSTINENCE] == 0), ])
names(d_rad_abs) = names(d_cli_abs)     <- paste0(d_names, "_abs")
names(d_rad_absNO) = names(d_cli_absNO) <- paste0(d_names, "_absNO")
# =============================================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# FILTER RADIOMIC FEATURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## By variance and missingness -----
QC_rad        <- lapply(d_rad, function(x){ goodSamplesGenes(datExpr = x, verbose = 3)})
QC_rad_abs    <- lapply(d_rad_abs, function(x){ goodSamplesGenes(datExpr = x, verbose = 3)})
QC_rad_absNO  <- lapply(d_rad_absNO, function(x){ goodSamplesGenes(datExpr = x, verbose = 3)})

QC_rm_col <- apply(rbind(#!do.call(rbind, lapply(QC_rad, function(x){ x$goodGenes})),  #Subgroups contain all exclusions
  !do.call(rbind, lapply(QC_rad_abs, function(x){ x$goodGenes})),
  !do.call(rbind, lapply(QC_rad_absNO, function(x){ x$goodGenes}))), 2, any)

### Remove bad features -----
d_rad        <- lapply(d_rad, function(x){ x[, rad_params[QC_rm_col] := NULL]})
d_rad_abs    <- lapply(d_rad_abs, function(x){ x[, rad_params[QC_rm_col] := NULL]})
d_rad_absNO  <- lapply(d_rad_absNO, function(x){ x[, rad_params[QC_rm_col] := NULL]})

rad_params <- rad_params[!QC_rm_col]
# =============================================================================

## By equality to other columns -----
QC_rm_col_eq <- apply(rbind(#do.call(rbind, lapply(d_rad, function(x){ duplicated(as.list(x))})),  #Subgroups contain all exclusions
  do.call(rbind, lapply(d_rad_abs, function(x){ duplicated(as.list(x))})),
  do.call(rbind, lapply(d_rad_absNO, function(x){ duplicated(as.list(x))}))), 2, any)

### Remove bad features -----
d_rad        <- lapply(d_rad, function(x){ x[, rad_params[QC_rm_col_eq] := NULL]})
d_rad_abs    <- lapply(d_rad_abs, function(x){ x[, rad_params[QC_rm_col_eq] := NULL]})
d_rad_absNO  <- lapply(d_rad_absNO, function(x){ x[, rad_params[QC_rm_col_eq] := NULL]})

rad_params <- rad_params[!QC_rm_col_eq]
# =============================================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PLOT PATIENT CLUSTERING DENDROGRAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pat_trees <- c(lapply(d_rad, function(x){ hclust(dist(x, method = "euclidean"), method = "average")}),
               lapply(d_rad_abs, function(x){ hclust(dist(x, method = "euclidean"), method = "average")}),
               lapply(d_rad_absNO, function(x){ hclust(dist(x, method = "euclidean"), method = "average")}))
listNames <- c(rep("all", 3), rep("abs", 3), rep("absNO", 3))
plotNames <- rep(d_names, 3)

which_param <- c("Age_baseline", "Male", "newrisk_baseline", "abstinence_pre", "LANCP","NCP", "CP") #Which clinical parameters to plot
pat_clin <- c(d_cli, d_cli_abs, d_cli_absNO)

lapply(1:length(pat_trees), function(x){
  pdf(file = paste0(folder_img, "Patient_clustering_", NAME, "/", plotNames[x], "_", NAME, "_", listNames[x], "_rm_with_clinical.pdf"), width = 12, height = 9)
  traitColors <- numbers2colors(pat_clin[[x]][, (which_param), with = FALSE], signed = FALSE)
  plotDendroAndColors(pat_trees[[x]], traitColors,
                      groupLabels = which_param,
                      main = paste0("Patient Clustering and Clinical Characteristics - ", plotNames[x], " ", listNames[x]),
                      cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, cex.dendroLabels = 0.2)
  dev.off()
})
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# CREATE WGCNA COMPATIBLE DATAFORMATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_final_cli = d_final_rad = d_final_cli_abs = d_final_rad_abs = d_final_cli_absNO = d_final_rad_absNO <- vector(mode = "list", length = length(d_names))
names(d_final_cli) = names(d_final_rad) <- d_names
names(d_final_cli_abs) =  names(d_final_rad_abs) <- paste0(d_names, "_abs")
names(d_final_cli_absNO) =  names(d_final_rad_absNO) <- paste0(d_names, "_absNO")
d_final_cli <- lapply(d_cli, function(x) list(data = x))
d_final_rad <- lapply(d_rad, function(x) list(data = x))
d_final_cli_abs <- lapply(d_cli_abs, function(x) list(data = x))
d_final_rad_abs <- lapply(d_rad_abs, function(x) list(data = x))
d_final_cli_absNO <- lapply(d_cli_absNO, function(x) list(data = x))
d_final_rad_absNO <- lapply(d_rad_absNO, function(x) list(data = x))

rm(list = c("d_cli", "d_cli_abs", "d_cli_absNO", "d_rad", "d_rad_abs", "d_rad_absNO", "pat_clin"))
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# CALCULATE COPHENETIC CORRELATION BETWEEN TREES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IDs_in_all <- intersect(d_final_cli$BAS$data$ID, d_final_cli$MID$data$ID)
IDs_in_all <- intersect(IDs_in_all, d_final_cli$END$data$ID)

tree_c_bas <- d_final_rad$BAS$data[d_final_cli$BAS$data$ID %in% IDs_in_all]
tree_c_mid <- d_final_rad$MID$data[d_final_cli$MID$data$ID %in% IDs_in_all]
tree_c_end <- d_final_rad$END$data[d_final_cli$END$data$ID %in% IDs_in_all]

dist_c_bas <- dist(tree_c_bas, method = "euclidean"); setattr(dist_c_bas, "Labels", d_final_cli$BAS$data[ID %in% IDs_in_all, ID])
dist_c_mid <- dist(tree_c_mid, method = "euclidean"); setattr(dist_c_mid, "Labels", d_final_cli$MID$data[ID %in% IDs_in_all, ID])
dist_c_end <- dist(tree_c_end, method = "euclidean"); setattr(dist_c_end, "Labels", d_final_cli$END$data[ID %in% IDs_in_all, ID])

dendextend::cor_cophenetic(hclust(dist_c_bas, method = "average"), hclust(dist_c_mid, method = "average"))
dendextend::cor_cophenetic(hclust(dist_c_mid, method = "average"), hclust(dist_c_end, method = "average"))
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# SAVE DATA AND EXIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save.image(paste0(folder_cache, "FINAL_DATA_", NAME, ".RData"))

rm(list = ls())
collectGarbage()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
