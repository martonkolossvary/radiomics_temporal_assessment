#####
# PROGRAM CODE FOR ARTICLE: "Assessment of structural changes in pathologies on radiological images following therapeutical interventions"
# (c) Márton Kolossváry, 2021

# DESCRIPTION: Network changepoint analysis
#####

# INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library("data.table"); library("WGCNA"); library("NetworkChange")
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

group_color <- c("#A6EB99", "#DBDA48", "#EA3AF7", "#F6C3CB", "#E6FEFF",
                 "#73DDD0", "#486AD9", "#0018F5", "#983430", "#000000",
                 "#75FA4C", "#74FBFD", "#FFFFE3", "#EA3323", "#999999",
                 "#D1BA98", "#7C2BC4", "#181B6B", "#BFFC5B", "#EB8677")
r_value <- 0.5
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# RUN CHANGEPOINT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Initialize datatables =====
outcomes  <- colnames(d_final_cli$BAS$data)[grep("^ME", colnames(d_final_cli$BAS$data))]
outcomes  <- outcomes[1:length(outcomes)-1] #Remove MEgrey

d_BAS_abs   <- d_final_cli$BAS$data[abstinence_pre == 1, ..outcomes]
d_BAS_absNO <- d_final_cli$BAS$data[abstinence_pre == 0, ..outcomes]
d_MID_abs   <- d_final_cli$MID$data[abstinence_pre == 1, ..outcomes]
d_MID_absNO <- d_final_cli$MID$data[abstinence_pre == 0, ..outcomes]
d_END_abs   <- d_final_cli$END$data[abstinence_pre == 1, ..outcomes]
d_END_absNO <- d_final_cli$END$data[abstinence_pre == 0, ..outcomes]
# =============================================================================

## Initialize eigen radiomics network =====
c_BAS_abs <- cor(d_BAS_abs); c_BAS_abs[c_BAS_abs<r_value] <- 0; diag(c_BAS_abs) <- 0
c_BAS_absNO <- cor(d_BAS_absNO); c_BAS_absNO[c_BAS_absNO<r_value] <- 0; diag(c_BAS_absNO) <- 0
c_MID_abs <- cor(d_MID_abs); c_MID_abs[c_MID_abs<r_value] <- 0; diag(c_MID_abs) <- 0
c_MID_absNO <- cor(d_MID_absNO); c_MID_absNO[c_MID_absNO<r_value] <- 0; diag(c_MID_absNO) <- 0
c_END_abs <- cor(d_END_abs); c_END_abs[c_END_abs<r_value] <- 0; diag(c_END_abs) <- 0
c_END_absNO <- cor(d_END_absNO); c_END_absNO[c_END_absNO<r_value] <- 0; diag(c_END_absNO) <- 0
# =============================================================================

## Format data for NetworkChange =====
abs <- array(c(c_BAS_abs, c_MID_abs, c_END_abs), dim = c(20, 20, 3))
absNO <- array(c(c_BAS_absNO, c_MID_absNO, c_END_absNO), dim = c(20, 20, 3))

attrib <- NULL
attrib$dim <- c(20, 20, 3)
attrib$dimnames[[1]] <- gsub("ME", "", outcomes, fixed = TRUE)
attrib$dimnames[[2]] <- gsub("ME", "", outcomes, fixed = TRUE)
attrib$dimnames[[3]] <- c(0, 6, 12)

attributes(abs) <- attrib
attributes(absNO) <- attrib
# =============================================================================

## Break number diagnostics =====
p_break_abs   <- BreakDiagnostic(abs, R = 2, mcmc = 1000, burnin = 1000, break.upper = 2)
p_break_absNO <- BreakDiagnostic(absNO, R = 2, mcmc = 1000, burnin = 1000, break.upper = 2)

### Plot break number diagnostics =====
p_break_abs <- p_break_abs$graph +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 24, family = "Helvetica"),
        plot.subtitle = element_text(size = 16, family = "Helvetica"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12, family = "Helvetica"),
        axis.text.y = element_text(size = 12, family = "Helvetica"),
        axis.title = element_text(size = 12, family = "Helvetica")) +
  labs(title = "Break Number Diagnostics")

p_break_absNO <- p_break_absNO$graph +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 24, family = "Helvetica"),
        plot.subtitle = element_text(size = 16, family = "Helvetica"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12, family = "Helvetica"),
        axis.text.y = element_text(size = 12, family = "Helvetica"),
        axis.title = element_text(size = 12, family = "Helvetica")) +
  labs(title = "Break Number Diagnostics")

ggsave(paste0(folder_img, "Network_change/", "Break_abs.pdf"), p_break_abs, device = "pdf", width = 8, height = 6, units = "in")
ggsave(paste0(folder_img, "Network_change/", "Break_absNO.pdf"), p_break_absNO, device = "pdf", width = 8, height = 6, units = "in")
# =============================================================================

## Network changepoint analysis =====
change_abs   <- NetworkChange(abs, R = 2, m = 1, mcmc = 1000, burnin = 1000)
change_absNO <- NetworkChange(absNO, R = 2, m = 1, mcmc = 1000, burnin = 1000)
# =============================================================================

## Plot latent node positions and Plot layer-specific network generation rules =====
### Plot network generation rules ====
pdf(file = paste0(folder_img, "Network_change/", "V_abs.pdf"), width = 6, height = 6)
plotV(change_abs) # network generation rules, rules of connection among individual nodes
dev.off()

pdf(file = paste0(folder_img, "Network_change/", "V_absNO.pdf"), width = 6, height = 6)
plotV(change_absNO)
dev.off()

### Plot latent node positions
trace(plotU, edit = TRUE) #Return p.list at end and "colour = alpha("red", 1/5))" to "colour = alpha(names, 1))" and lables to reasonable
pU_abs   <- plotU(change_abs, names = stringr::str_to_title(gsub("ME","", outcomes)), label.prob = 0) # node traits, characteristics of individual nodes such as positions of nodes in the latent space
pU_absNO <- plotU(change_absNO, names = stringr::str_to_title(gsub("ME","", outcomes)), label.prob = 0)


pU_abs_1 <- pU_abs[[1]] +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 24, family = "Helvetica"),
        plot.subtitle = element_text(size = 16, family = "Helvetica"),
        axis.text = element_text(size = 20, family = "Helvetica"),
        axis.title = element_text(size = 20, family = "Helvetica")) +
  labs(title = "Latent Node Possitions",
       subtitle = "First time regime") +
  xlab("First dimension") + ylab("Second dimension") 

pU_abs_2 <- pU_abs[[2]] +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 24, family = "Helvetica"),
        plot.subtitle = element_text(size = 16, family = "Helvetica"),
        axis.text = element_text(size = 20, family = "Helvetica"),
        axis.title = element_text(size = 20, family = "Helvetica")) +
  labs(title = "Latent Node Possitions",
       subtitle = "First time regime") +
  xlab("First dimension") + ylab("Second dimension") 

pU_absNO_1 <- pU_absNO[[1]] +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 24, family = "Helvetica"),
        plot.subtitle = element_text(size = 16, family = "Helvetica"),
        axis.text = element_text(size = 20, family = "Helvetica"),
        axis.title = element_text(size = 20, family = "Helvetica")) +
  labs(title = "Latent Node Possitions",
       subtitle = "First time regime") +
  xlab("First dimension") + ylab("Second dimension") 

pU_absNO_2 <- pU_absNO[[2]] +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 24, family = "Helvetica"),
        plot.subtitle = element_text(size = 16, family = "Helvetica"),
        axis.text = element_text(size = 20, family = "Helvetica"),
        axis.title = element_text(size = 20, family = "Helvetica")) +
  labs(title = "Latent Node Possitions",
       subtitle = "First time regime") +
  xlab("First dimension") + ylab("Second dimension") 

### Save plots =====
ggsave(paste0(folder_img, "Network_change/", "U_abs_1.pdf"), pU_abs_1, device = "pdf", width = 6, height = 6, units = "in")
ggsave(paste0(folder_img, "Network_change/", "U_abs_2.pdf"), pU_abs_2, device = "pdf", width = 6, height = 6, units = "in")
ggsave(paste0(folder_img, "Network_change/", "U_absNO_1.pdf"), pU_absNO_1, device = "pdf", width = 6, height = 6, units = "in")
ggsave(paste0(folder_img, "Network_change/", "U_absNO_2.pdf"), pU_absNO_2, device = "pdf", width = 6, height = 6, units = "in")
# =============================================================================

## Plot posterior probabilities =====
attr(change_abs, "y") <- 1:3
attr(change_absNO, "y") <- 1:3

pdf(file = paste0(folder_img, "Network_change/", "State_abs.pdf"), width = 6, height = 6)
MCMCpack::plotState(change_abs, start = 1, main = "Posterior Probabilities Of Time Regimes")
dev.off()

pdf(file = paste0(folder_img, "Network_change/", "State_absNO.pdf"), width = 6, height = 6)
MCMCpack::plotState(change_absNO, start = 1, main = "Posterior Probabilities Of Time Regimes")
dev.off()
# =============================================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# EXIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list = ls())
collectGarbage()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
