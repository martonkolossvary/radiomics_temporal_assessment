#####
# PROGRAM CODE FOR ARTICLE: "Assessment of structural changes in pathologies on radiological images following therapeutical interventions"
# (c) Márton Kolossváry, 2021

# DESCRIPTION: Code to run Linear Mixed Model and Mediation analyses
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


# LINEAR MIXED MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Initialize parameters =====
type      <- "gaussian"; exponentiate <- ifelse(type == "binomial", TRUE, FALSE)
pred_int  <- c("Age_baseline", "Male", "newrisk_baseline", "LANCP", "NCP", "CP", "abstinence_pre"); pred_pre  <- pred_int
outcomes  <- colnames(d_model)[grep("^ME", colnames(d_model))]
outcomes  <- outcomes[1:length(outcomes)-1] #Remove MEgrey
d_model   <- rbindlist(list(d_final_cli$BAS$data, d_final_cli$MID$data, d_final_cli$END$data))
rand_slope2 <- paste0(" + (0 +", time, "|", group, ")")
# =============================================================================

## Run LMM models =====
lms_rad_m2_multi <- lapply(outcomes, function(x) {
  tryCatch(lmerTest::lmer(formula = paste0(x, " ~ ", paste(pred_int,  collapse = " + "),
                                           " + ","(", time_2, ")", rand_slope2),
                          data = d_model_int),
           error = function(e) NA)
})

lms_coeff_rad_m2_multi <- lapply(lms_rad_m2_multi, function(x) {
  broomExtra::tidy(x, conf.int = TRUE, effects = "fixed")
})

out_mul_int <- lapply(lms_coeff_rad_m2_multi, function(x) {x[8,]})
out_mul_int <- rbindlist(out_mul_int)
write.csv(out_mul_int, paste0(folder_stat, "mul_int.csv"))
# =============================================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# MEDIATION ANALYSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Initialize parameters =====
mediators  <- c("Cholesterol","cardioCRP", "Endothelin_1") #c("LANCP", "NCP", "CP")
pred_int   <- c("Age_baseline", "Male", "newrisk_baseline", "abstinence_pre")
d_model_int_mod <- data.table::copy(d_model_int)
# =============================================================================

## Run indirect pathways =====
mod_indir_l <- lapply(mediators, function(mediator){
  tryCatch(lmerTest::lmer(formula = paste0(mediator, " ~ ", paste(pred_int, collapse = " + "), "+", paste(mediators[!mediators %in% mediator],  collapse = " + "),
                                           " + ","(", time_2, ")", rand_slope2),
                          data = d_model_int_mod), error = function(e) NA)
})

lms_coeff_mod_indir <- lapply(mod_indir_l, function(x) {
  broomExtra::tidy(x, conf.int = TRUE, effects = "fixed")
})
names(lms_coeff_mod_indir) <- mediators

lapply(1:length(mediators), function(x) {data.table::fwrite(as.data.frame(lms_coeff_mod_indir[[x]]),
                                                            file = paste0(folder_stat, "indirect_pathway_", names(lms_coeff_mod_indir[x]), ".csv"))})
# =============================================================================

## Run mediation =====
library(mediation)
id_outcome <- 4 #4: pink, 12: Cyan
mediator   <- "Endothelin_1" #"NCP", "CP", "LANCP", "Cholesterol","cardioCRP", "Endothelin_1"

### Run indirect for given mediator =====
mod_indir <- tryCatch(lmerTest::lmer(formula = paste0(mediator, " ~ ", paste(pred_int,  collapse = " + "), "+", paste(mediators[!mediators %in% mediator],  collapse = " + "),
                                                      " + ","(", time_2, ")", rand_slope2),
                                     data = d_model_int_mod), error = function(e) NA)
broomExtra::tidy(mod_indir)

### Find common individuals =====
id_patients      <- paste0(as.character(mod_indir@frame$ID), "_", as.character(mod_indir@frame$TIME_2))
d_model_int_redo <- data.table::copy(d_model_int_mod)
d_model_int_redo[, id_patients := paste0(ID, "_", TIME_2)]
select_rows      <- d_model_int_redo$id_patients %in% id_patients
d_model_int_redo <- d_model_int_redo[select_rows, ]

### Run direct for given mediator =====
mod_dir   <- tryCatch(lmerTest::lmer(formula = paste0(outcomes[id_outcome], " ~ ", paste(pred_int,  collapse = " + "), " + ", mediator,
                                                      " + ","(", time_2, ")", rand_slope2),
                                     data = d_model_int_redo), error = function(e) NA)

### Run mediation =====
mod_med <- mediation::mediate(model.m = mod_indir, model.y = mod_dir, treat = "abstinence_pre", mediator = mediator,
                              sims = 1000, boot = FALSE)

sink(file = paste0(folder_stat, "Mediation_", outcomes[id_outcome], "_", mediator, ".txt"))
summary(mod_med)
sink()
# =============================================================================


# EXIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list = ls())
collectGarbage()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
