library(IPDfromKM)
library(survival)
library(survminer)
library(dplyr)
library(survival)
library(ggplot2)
library(ggsurvfit)
library(dplyr)
library(bayesplot)
library(survminer)
library(ggtext)
library(patchwork)

# Read in extracted coordinates
points_ep <- read.table("Data/ep_km_reconstruction.txt", header = FALSE)
points_durv <- read.table("Data/durv_ep_km_reconstruction.txt", header = FALSE)

# Read in extra coordinates corresponding to censoring locations
# Including these makes the censoring reproduction much more accurate
extra_ep <- read.table("Data/ep_extra.txt", header = FALSE)
points_ep <- rbind(points_ep, extra_ep) %>% arrange(V1)

extra_durv <- read.table("Data/durv_extra.txt", header = FALSE)
points_durv <- rbind(points_durv, extra_durv) %>% arrange(V1)

# Extract IP for EP Arm
prep_ep <- preprocess(points_ep, 
                      trisk = seq(0, 51, 3),
                      nrisk = c(269, 243, 212, 156, 104, 82, 64, 51, 36, 24, 19, 17, 13, 10, 3, 0, 0, 0),
                      maxy = 1)

ipd_ep <- getIPD(prep_ep, armID = 0)
ipd_ep$IPD

# Extract IPD for Durvalumab + EP arm
prep_durv <- preprocess(points_durv, 
                        trisk = seq(0, 51, 3),
                        nrisk = c(268, 244, 214, 177, 140, 109, 85, 70, 60, 54, 50, 46, 39, 25, 13, 3, 0, 0),
                        maxy = 1)

ipd_durv <- getIPD(prep_durv, armID = 1)
ipd_durv$IPD


# Create combined dataset of IPD for both arms: save this object if desired
ipd_df <- rbind(ipd_ep$IPD, ipd_durv$IPD)
ipd_df <- ipd_df %>% 
  mutate(trt_group = ifelse(treat == 1, "Durvalumab + EP", "EP")) %>% 
  select(time, status, trt_group)
ipd_df$trt_group <- as.factor(ipd_df$trt_group)


# Visualize Kaplan-Meier curves of reconstructed data
fit_ipd <- survfit(Surv(time, status) ~ trt_group, data = ipd_df)
mysurvplot <- ggsurvplot(
  fit_ipd, 
  conf.int = FALSE,
  ggtheme = theme_classic(), 
  data = ipd_df,
  legend.labs = c("Durvalumab + EP", "EP"),
  palette = c("#1F3A93", "#A93226"),
  legend.title="",
  risk.table = TRUE, 
  risk.table.pos = "out",
  risk.table.col = "strata",
  tables.theme = theme_void() + theme(axis.text.y = element_text(size = 12)),
  break.time.by = 3, 
  break.y.by = 0.2, 
  xlim = c(0, 51),
  xlab = "Time (months)", 
  ylab = "Probability of Overall Survival"
)
scale <- 1.7
text_size <- 9*scale
mysurvplot$plot <- mysurvplot$plot + 
  # scale_x_continuous(expand = c(0, 0), breaks = seq(0, 51, 3)) +
  theme(axis.text = element_text(size = 14, color = "black"), 
        legend.text = element_text(size = text_size*0.8, color = "black"), 
        axis.title = element_text(size = text_size),
        axis.line = element_line(color = "black"))

# pdf("kaplan_meier.pdf", height = scale*4, width = scale*8)
mysurvplot$plot / mysurvplot$table + plot_layout(heights = c(8,2))
# dev.off()

