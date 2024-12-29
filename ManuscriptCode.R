#################################################################
## R CODE FOR MANUSCRIPT
#################################################################
# Toward merging MOPEX and CAMELS 
# Katharine Sink

# import libraries 
library(tidyverse)
library(rstatix)
library(patchwork) # figures
library(Metrics) # MAE
library(Hmisc) # bootstrap
library(ggpubr)
library(sf)
library(moments)
library(caTools)
library(e1071)

# import files 
daily_data = readRDS("Scripts/EDA/daily_data.rds")
id_attributes = readRDS("id_attributes.rds")

#################################################################
## AGGREGATED DATA
#################################################################
# daily 
daily = daily_data %>% group_by(GaugeID, climate) %>% 
  reframe(PRCP_C = PRCP_C, PRCP_M = PRCP_M, TAIR_C = TAIR_C, 
            TAIR_M = TAIR_M)

# aggregate daily data into monthly totals for precipitation and mean for temperature
month = daily_data %>% group_by(GaugeID, WATERYR, MNTH, climate) %>% 
   reframe(PRCP_C = sum(PRCP_C), PRCP_M = sum(PRCP_M), TAIR_C = mean(TAIR_C), 
            TAIR_M = mean(TAIR_M))

# seasonal aggregation
season = daily_data %>% group_by(GaugeID, WATERYR, SEASON, climate) %>% 
   reframe(PRCP_C = sum(PRCP_C), PRCP_M = sum(PRCP_M), TAIR_C = mean(TAIR_C), 
            TAIR_M = mean(TAIR_M))

# water year (annual) aggregation
year = daily_data %>% group_by(GaugeID, WATERYR, climate) %>% 
   reframe(PRCP_C = sum(PRCP_C), PRCP_M = sum(PRCP_M), TAIR_C = mean(TAIR_C), 
            TAIR_M = mean(TAIR_M))

#################################################################
## KMEANS CLUSTERING
#################################################################
# determine annual aridity index 
aridity_C = daily_data %>% group_by(GaugeID, WATERYR) %>% summarise(PETvP_C = sum(PET_C)/sum(PRCP_C))
aridity_C_WB = daily_data %>% group_by(GaugeID, WATERYR) %>% summarise(PETvP_WB = sum(PET_C)/sum(PRCP_C))
aridity_M = daily_data %>% group_by(GaugeID, WATERYR) %>% summarise(PETvP_M = sum(PET_M)/sum(PRCP_M))

# overall mean based on annual aridity values
aridity_C = aridity_C %>% group_by(GaugeID) %>% summarise(PETvP_C = mean(PETvP_C))
aridity_C_WB = aridity_C_WB %>% group_by(GaugeID) %>% summarise(PETvP_WB = mean(PETvP_WB))
aridity_M = aridity_M %>% group_by(GaugeID) %>% summarise(PETvP_M = mean(PETvP_M))

# evaporative index for each water year
evaporative_C = daily_data %>% group_by(GaugeID, WATERYR) %>% summarise(ETvP_C = sum(ET_C)/sum(PRCP_C))
evaporative_C_WB = daily_data %>% group_by(GaugeID, WATERYR) %>% summarise(ETvP_WB = sum(ET_WB_C)/sum(PRCP_C))
evaporative_M = daily_data %>% group_by(GaugeID, WATERYR) %>% summarise(ETvP_M = sum(ET_M)/sum(PRCP_M))

# overall mean evaporative index per gauge 
evaporative_C = evaporative_C %>% group_by(GaugeID) %>% summarise(ETvP_C = mean(ETvP_C))
evaporative_C_WB = evaporative_C_WB %>% group_by(GaugeID) %>% summarise(ETvP_WB = mean(ETvP_WB))
evaporative_M = evaporative_M %>% group_by(GaugeID) %>% summarise(ETvP_M = mean(ETvP_M))

# join evaporative and aridity indices together into one dataframe 
indices_C = left_join(aridity_C, evaporative_C, by = "GaugeID")
indices_C_WB = left_join(aridity_C_WB, evaporative_C_WB, by = "GaugeID")
indices_M = left_join(aridity_M, evaporative_M, by = "GaugeID")

# get the first column to use as row names
columnnames = indices_C[,1]
# remove the gauge id column (no need to scale since the variability isn't extreme)
# and rename 
indices_C.kmean = indices_C[,-1]
indices_C_WB.kmean = indices_C_WB[,-1]
indices_M.kmean = indices_M[,-1]

# compute k-means with k = 3 based on plot
km_C.res = kmeans(indices_C.kmean, 3, nstart = 100)
km_C_WB.res = kmeans(indices_C_WB.kmean, 3, nstart = 100)
km_M.res = kmeans(indices_M.kmean, 3, nstart = 100)

aggregate(indices_C.kmean, by = list(cluster = km_C.res$cluster), mean)
aggregate(indices_C_WB.kmean, by = list(cluster = km_C_WB.res$cluster), mean)
aggregate(indices_M.kmean, by = list(cluster = km_M.res$cluster), mean)

# add cluster classification to original data 
c = cbind(indices_C, cluster = km_C.res$cluster)
c_wb = cbind(indices_C_WB, cluster = km_C_WB.res$cluster)
m = cbind(indices_M, cluster = km_M.res$cluster)

#################################################################
## BUDYKO
#################################################################
# determine evaporative index with Budyko equation (1974)
indices_C$Budyko_C = sqrt((1 - exp(-indices_C$PETvP_C))*indices_C$PETvP_C*tanh(1/indices_C$PETvP_C))
indices_C_WB$Budyko_C_WB = sqrt((1 - exp(-indices_C_WB$PETvP_WB))*indices_C_WB$PETvP_WB*tanh(1/indices_C_WB$PETvP_WB))
indices_M$Budyko_M = sqrt((1 - exp(-indices_M$PETvP_M))*indices_M$PETvP_M*tanh(1/indices_M$PETvP_M))

# combine aridity and evaporative indices for all gauges for camels sac et, camels wb et and mopex
# add climate column for each based on AI > 1.5 arid, continental AI 1.5 - 0.82, temperature AI < 0.82
indices_C$climate_C = with(indices_C, ifelse(PETvP_C >= 1.5, 'Arid', ifelse(PETvP_C >= 0.82, 'Continental', 'Temperate')))
indices_C_WB$climate_C_WB = with(indices_C_WB, ifelse(PETvP_WB >= 1.5, 'Arid', ifelse(PETvP_WB >= 0.82, 'Continental', 'Temperate')))
indices_M$climate_M = with(indices_M, ifelse(PETvP_M >= 1.5, 'Arid', ifelse(PETvP_M >= 0.82, 'Continental', 'Temperate')))

indices = left_join(indices_C, indices_C_WB, by = "GaugeID")
indices = left_join(indices, indices_M, by = "GaugeID")

# blank Budyko curve boundaries
Budyko = data.frame(PETvP = seq(0, 5, 1), ETvP = c(0, rep(1, 5)))

# triangle shape on the graph
vertices = data.frame(x = c(0, 0, 1), y = c(0, 1, 1))

ggplot(indices) + 
  geom_point(aes(x = PETvP_C, y = ETvP_C, color = climate_C, shape = "CAMELS SAC-SMA"), size = 4) +
  geom_point(aes(x = PETvP_WB, y = ETvP_WB, color = climate_C_WB, shape = "CAMELS WB"), size = 4) +
  geom_point(aes(x = PETvP_M, y = ETvP_M, color = climate_M, shape = "MOPEX"), size = 4) + 
  scale_color_manual(values = c("Arid" = "#D55E00", "Continental" = "#E69F00", "Temperate" = "#009E73")) +
  scale_shape_manual(values = c("CAMELS SAC-SMA" = 17, "CAMELS WB" = 24, "MOPEX" = 8)) +
  guides(color = guide_legend(title = "Climate", override.aes = list(shape = 16)), 
         shape = guide_legend(title = "Dataset", override.aes = list(color = "black"))) +
  labs(x = "Aridity Index (PET/P)", y = "Evaporative Index (ET/P)") + 
 theme(plot.caption = element_text(size = 8)) + theme_bw() +
   geom_polygon(data = vertices, aes(x = x, y = y), alpha = 0.5, fill = 'lightgray') +
  geom_line(data = Budyko, aes(x = PETvP, y = ETvP)) + 
    annotate('rect', xmin=0, xmax=5, ymin=1, ymax=1.25, alpha=0.5, fill='lightgray') +
coord_cartesian(xlim = c(0,5), ylim = c(0,1.2), expand = FALSE) + 
  geom_vline(xintercept = 1, lty = 2) +
  geom_text(data = data.frame(labels = c("Energy Limited", "Water Limited", "Above Water Limit", 
            "60%", "46%", "12%", "10%", "12%"), 
    x = c(0.6, 2.3, 3.0, 0.8, 1.65, 1.92, 2.69, 1.05), 
    y = c(0.20, 0.20, 1.15, 0.4, 0.6, 1.04, 1.02, 0.7)), aes(x = x, y = y, label = labels)) +
  geom_segment(aes(x = 0.676, y = 0.566, xend = 0.676, yend = 0.305), 
               arrow = arrow(length = unit(0.3, 'cm')), lwd =1) +
  geom_segment(aes(x = 1.55, y = 0.777, xend = 1.55, yend = 0.484), 
               arrow = arrow(length = unit(0.3, 'cm')), lwd =1) +
   geom_segment(aes(x = 1.83, y = 1.08, xend = 1.83, yend = 0.961), 
               arrow = arrow(length = unit(0.3, 'cm')), lwd =1) +
   geom_segment(aes(x = 2.59, y = 1.04, xend = 2.59, yend = 0.933), 
               arrow = arrow(length = unit(0.3, 'cm')), lwd =1) +
   geom_segment(aes(x = 0.948, y = 0.735, xend = 0.948, yend = 0.65), 
               arrow = arrow(length = unit(0.3, 'cm')), lwd =1) +
  geom_line(aes(x = PETvP_M, y = Budyko_M), linewidth = 1) +
  theme(legend.position = c(0.8, 0.3))

#################################################################
## STATS
#################################################################
# statistics for each climate region using each gauge 
stats = month %>% group_by(climate) %>% 
  reframe(minPRCP_C = min(PRCP_C), maxPRCP_C = max(PRCP_C), medPRCP_C = median(PRCP_C), meanPRCP_C = mean(PRCP_C), varPRCP_C = var(PRCP_C), sdPRCP_C = sd(PRCP_C), skPRCP_C = skewness(PRCP_C),  
          minPRCP_M = min(PRCP_M), maxPRCP_M = max(PRCP_M), medPRCP_M = median(PRCP_M), meanPRCP_M = mean(PRCP_M), varPRCP_M = var(PRCP_M), sdPRCP_M = sd(PRCP_M), skPRCP_M = skewness(PRCP_M),
          minTAIR_C = min(TAIR_C), maxTAIR_C = max(TAIR_C), medTAIR_C = median(TAIR_C), meanTAIR_C = mean(TAIR_C), varTAIR_C = var(TAIR_C), sdTAIR_C = sd(TAIR_C), skTAIR_C = skewness(TAIR_C),
          minTAIR_M = min(TAIR_M), maxTAIR_M = max(TAIR_M), medTAIR_M = median(TAIR_M), meanTAIR_M = mean(TAIR_M), varTAIR_M = var(TAIR_M), sdTAIR_M = sd(TAIR_M), skTAIR_M = skewness(TAIR_M))

stats2 = stats %>% group_by(climate) %>% 
   reframe(minPRCP_C = min(minPRCP_C), maxPRCP_C = max(maxPRCP_C), medPRCP_C = mean(medPRCP_C), varPRCP_C = mean(varPRCP_C), sdPRCP_C = mean(sdPRCP_C), skPRCP_C = mean(skPRCP_C),  
           minPRCP_M = min(minPRCP_M), maxPRCP_M = max(maxPRCP_M), medPRCP_M = mean(medPRCP_M), varPRCP_M = mean(varPRCP_M), sdPRCP_M = mean(sdPRCP_M), skPRCP_M = mean(skPRCP_M),
           minTAIR_C = min(minTAIR_C), maxTAIR_C = max(maxTAIR_C), medTAIR_C = mean(medTAIR_C), varTAIR_C = mean(varTAIR_C), sdTAIR_C = mean(sdTAIR_C), skTAIR_C = mean(skTAIR_C),
           minTAIR_M = min(minTAIR_M), maxTAIR_M = max(maxTAIR_M), medTAIR_M = mean(medTAIR_M), varTAIR_M = mean(varTAIR_M), sdTAIR_M = mean(sdTAIR_M), skTAIR_M = mean(skTAIR_M))

#################################################################
## COEFFICIENT OF VARIATION
#################################################################
# create function to calculate coefficient of variation
# dimensionless and can characterize the degree of variability in datasets
cv = function(x) {
  sd(x)/mean(x)
}

daily_cv = daily_data %>% group_by(GaugeID, climate) %>% 
  reframe(PRCP_C = cv(PRCP_C), PRCP_M = cv(PRCP_M), TAIR_C = cv(TAIR_C), TAIR_M = cv(TAIR_M))
daily_cv = daily_cv %>% pivot_longer(cols = c(PRCP_C, PRCP_M), names_to = "Names", values_to = "Values")
daily_cv$climate = factor(daily_cv$climate, levels = c("Arid", "Continental", "Temperate"))

month_cv = month %>% group_by(GaugeID, MNTH, climate) %>% 
  reframe(PRCP_C = cv(PRCP_C), PRCP_M = cv(PRCP_M), TAIR_C = cv(TAIR_C), TAIR_M = cv(TAIR_M))
month_cv = month_cv %>% pivot_longer(cols = c(PRCP_C, PRCP_M), names_to = "Names", values_to = "Values")
month_cv$climate = factor(month_cv$climate, levels = c("Arid", "Continental", "Temperate"))

season_cv = season %>% group_by(GaugeID, SEASON, climate) %>% 
  reframe(PRCP_C = cv(PRCP_C), PRCP_M = cv(PRCP_M), TAIR_C = cv(TAIR_C), TAIR_M = cv(TAIR_M))
season_cv = season_cv %>% pivot_longer(cols = c(PRCP_C, PRCP_M), names_to = "Names", values_to = "Values")
season_cv$climate = factor(season_cv$climate, levels = c("Arid", "Continental", "Temperate"))

year_cv = year %>% group_by(GaugeID, climate) %>% 
  reframe(PRCP_C = cv(PRCP_C), PRCP_M = cv(PRCP_M), TAIR_C = cv(TAIR_C), TAIR_M = cv(TAIR_M))
year_cv = year_cv %>% pivot_longer(cols = c(PRCP_C, PRCP_M), names_to = "Names", values_to = "Values")
year_cv$climate = factor(year_cv$climate, levels = c("Arid", "Continental", "Temperate"))

d = ggboxplot(data = daily_cv, x = "climate", y = "Values", color = "Names") +
  theme_classic() + theme(axis.title = element_blank(), legend.position = "none") +
  theme(axis.text = element_text(size = 8)) + 
  scale_color_manual(values = c("blue", "red")) + labs(title = "Day") + 
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  theme(plot.title = element_text(size = 10, hjust = 0.5)) 

m = ggboxplot(data = month_cv, x = "climate", y = "Values", color = "Names") + theme_classic() +
  theme(axis.text = element_text(size = 8)) + 
  theme(axis.title = element_blank(), legend.position = "none") + 
    scale_color_manual(values = c("blue", "red")) + labs(title = "Month") +
  theme(plot.title = element_text(size = 10, hjust = 0.5))  

s = ggboxplot(data = season_cv, x = "climate", y = "Values", color = "Names") +
  theme_classic() + theme(axis.title = element_blank(), legend.position = "none") +
   theme(axis.text = element_text(size = 8)) +  
  scale_color_manual(values = c("blue", "red")) + labs(title = "Season") +
  theme(plot.title = element_text(size = 10, hjust = 0.5)) 

y = ggboxplot(data = year_cv, x = "climate", y = "Values", color = "Names") +
   theme_bw() + theme(axis.title = element_blank()) + scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
   theme(axis.text = element_text(size = 8)) + 
  scale_color_manual(values = c("blue", "red"), labels = c("CAMELS", "MOPEX")) + 
  labs(title = "Water Year") + theme(plot.title = element_text(size = 10, hjust = 0.5)) + 
  theme(legend.position = c(0.8,0.8), legend.title = element_blank())

d + m + s + y + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
  subtitle = "Precipitation coefficient of variation") + plot_layout(nrow = 2, guides = "collect")

#################################################################
## BINOMIAL TESTS
#################################################################
# rstatix package
# significance test based on positive (statistic)
prcp = daily_data %>% select(GaugeID, PRCP_C, PRCP_M)
prcp = prcp %>% pivot_longer(cols = c(PRCP_C, PRCP_M), names_to = "Dataset", values_to = "Precip")
sign = prcp %>% group_by(GaugeID) %>% sign_test(Precip ~ Dataset, conf.level = 0.95, detailed = TRUE) %>% 
  add_significance()

# select data aggregation
A_TAIR = year %>% select(GaugeID, climate, TAIR_C, TAIR_M)

# perform on day, month, season, year data for PRCP and TAIR 
# positive (CAMELS>MOPEX), negative (MOPEX>CAMELS), zero (CAMELS=MOPEX)
A_TAIR = A_TAIR %>% mutate(Value = 
            case_when(
              TAIR_C > TAIR_M ~ "CAMELS", 
              TAIR_C < TAIR_M ~ "MOPEX", 
              TAIR_C == TAIR_M ~ "SAME"))

# tally signs 
A_TAIR = A_TAIR %>% group_by(climate, Value) %>% summarise(Count = n())

# bar plots for numerical count of signs
d_TAIR = ggplot(data = D_TAIR, aes(x = climate, y = Count, fill = Value)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
  geom_text(aes(label = Count), vjust = -0.20, size = 2.5, color = "black", position = position_dodge(0.9)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 8), legend.title = element_blank()) +
    labs(title = "Day") + theme(plot.title = element_text(size = 10, hjust = 0.5)) + 
    labs(y = "Number of days") +  scale_fill_manual(values = c("blue", "red", "black"))

m_TAIR = ggplot(data = M_TAIR, aes(x = climate, y = Count, fill = Value)) +
   geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
  geom_text(aes(label = Count), vjust = -0.20, size = 2.5, color = "black", position = position_dodge(0.9)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
         axis.text.y = element_text(size = 8), legend.title = element_blank(), legend.position = "none") +
   labs(title = "Month") + theme(plot.title = element_text(size = 10, hjust = 0.5)) + 
    labs(y = "Number of months") +  scale_fill_manual(values = c("blue", "red", "black"))

s_TAIR = ggplot(data = S_TAIR, aes(x = climate, y = Count, fill = Value)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
   geom_text(aes(label = Count), vjust = -0.22, size = 2.5, color = "black", position = position_dodge(0.9)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
           axis.text.y = element_text(size = 8), legend.title = element_blank(), legend.position = "none") +
   labs(title = "Season") + theme(plot.title = element_text(size = 10, hjust = 0.5)) + 
    labs(y = "Number of seasons") +  scale_fill_manual(values = c("blue", "red", "black"))

a_TAIR = ggplot(data = A_TAIR, aes(x = climate, y = Count, fill = Value)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
  geom_text(aes(label = Count), vjust = -0.22, size = 2.5, color = "black", position = position_dodge(0.9)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
         axis.text.y = element_text(size = 8), legend.title = element_blank(), legend.position = "none") +
   labs(title = "Water Year") + theme(plot.title = element_text(size = 10, hjust = 0.5)) + 
    labs(y = "Number of years") +  scale_fill_manual(values = c("blue", "red", "black"))

# patchwork package
d_TAIR + m_TAIR + s_TAIR + a_TAIR + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")", 
                        subtitle = "Temperature binomial sign test results") +
  plot_layout(guides = "collect")

# significance with binom_test
res = binom_test(x = 100, n = 480)

#################################################################
## BIASES
#################################################################
daily_diff = daily_data %>% group_by(GaugeID, climate) %>% reframe(PRCP = PRCP_C - PRCP_M, 
TAIR = TAIR_C - TAIR_M)

# monthly bias
month_diff = month %>% group_by(GaugeID, WATERYR, MNTH, climate) %>% 
  reframe(PRCP = PRCP_C - PRCP_M, TAIR = TAIR_C - TAIR_M, NAME = month.abb[MNTH])
# change month names to factor in order to display on plot in chronological order instead of alphabetically
month_diff$NAME = factor(month_diff$NAME, levels = month.abb)

# seasonal bias
season_diff = season %>% group_by(GaugeID, WATERYR, SEASON, climate) %>% 
  reframe(PRCP = PRCP_C - PRCP_M, TAIR = TAIR_C - TAIR_M)

# annual bias
year_diff = year %>% group_by(GaugeID, WATERYR, climate) %>% 
  reframe(PRCP = PRCP_C - PRCP_M, TAIR = TAIR_C - TAIR_M)

# box plot with values for all gauges within climate region
PRCP_m = ggboxplot(data = month_diff, x = "NAME", y = "PRCP", 
  xlab = FALSE, color = "climate", palette = c("#D55E00", "#E69F00", "#009E73")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 7), 
                    axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8)) +
    theme(legend.title = element_blank()) + labs(y = "Monthly bias (mm)") 
#  theme(legend.title = element_blank()) + labs(y = "Monthly bias (\u00b0C)") 

PRCP_s = ggboxplot(data = season_diff, x = "SEASON", y = "PRCP", 
  xlab = FALSE, color = "climate", palette = c("#D55E00", "#E69F00", "#009E73")) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 7), 
                     axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8)) +
      theme(legend.title = element_blank()) + labs(y = "Seasonal bias (mm)") 
  # theme(legend.title = element_blank()) + labs(y = "Seasonal bias (\u00b0C)") 

PRCP_y = ggboxplot(data = year_diff, x = "WATERYR", y = "PRCP", 
  xlab = FALSE, color = "climate", palette = c("#D55E00", "#E69F00", "#009E73")) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 7), 
           axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8)) +
      theme(legend.title = element_blank()) + labs(y = "Annual bias (mm)") 
  # theme(legend.title = element_blank()) + labs(y = "Annual bias (\u00b0C)") 

PRCP_m + PRCP_s + PRCP_y + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
  subtitle = "Precipitation bias CAMELS - MOPEX") +                                         
  plot_layout(nrow = 3, guides = "collect")

## AVERAGE BIAS ##
# monthly avg bias
month_avg = month_diff %>% group_by(MNTH, climate) %>% 
  reframe(PRCP = mean(PRCP), TAIR = mean(TAIR))
# month_avg = round_df(month_avg, 2, rf = "round")
month_avg = month_avg %>% mutate(NAME = month.abb[MNTH])
month_avg$NAME = factor(month_avg$NAME, levels = month.abb)

# seasonal avg bias
season_avg = season_diff %>% group_by(SEASON, climate) %>% 
   reframe(PRCP = mean(PRCP), TAIR = mean(TAIR))

# annual avg bias
year_avg = year_diff %>% group_by(WATERYR, climate) %>% 
   reframe(PRCP = mean(PRCP), TAIR = mean(TAIR))

m_plot = ggplot(data = month_avg, aes(x = NAME, y = TAIR, fill = climate)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
         axis.text = element_text(size = 8), legend.title = element_blank()) +
    labs(y = "Monthly bias (\u00b0C)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))
    # labs(y = "Monthly bias (mm)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))

s_plot = ggplot(data = season_avg, aes(x = SEASON, y = TAIR, fill = climate)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
         axis.text = element_text(size = 8), legend.title = element_blank()) +
    labs(y = "Seasonal bias (\u00b0C)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))
    # labs(y = "Seasonal bias (mm)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))
    
a_plot = ggplot(data = year_avg, aes(x = WATERYR, y = TAIR, fill = climate)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
         axis.text = element_text(size = 8), legend.title = element_blank()) +
     labs(y = "Annual bias (\u00b0C)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))
    # labs(y = "Annual bias (mm)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))
    
m_plot + s_plot + a_plot + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
  subtitle = "Overall average temperature bias CAMELS - MOPEX") + plot_layout(nrow = 3, guides = "collect")
  # subtitle = "Overall average precipitation bias CAMELS - MOPEX") + plot_layout(nrow = 3, guides = "collect")

# gauge bias map
gauges = unique(daily_data$GaugeID)
coords = id_attributes %>% filter(GaugeID %in% gauges)
coords = select(coords, c("GaugeID", "gauge_lat", "gauge_lon"))
conus = st_read("CONUS.shp")

diff = daily_diff %>% group_by(GaugeID) %>% reframe(PRCP = mean(PRCP), TAIR = mean(TAIR))
diff = left_join(diff, coords, by = "GaugeID") 

diff = diff %>% mutate(color_bin = case_when(
  PRCP >= -0.25 & PRCP < 0 ~ "red", 
  PRCP == 0 ~ "white", 
  PRCP > 0 & PRCP <= 0.25 ~ "lightblue", 
  PRCP > 0.25 & PRCP <= 0.50 ~ "blue", 
  PRCP > 0.50 & PRCP <= 0.75 ~ "darkblue"
))

diff = diff %>% mutate(color_bin = factor(color_bin, levels = c("red", "white", "lightblue", "blue", "darkblue")))

custom_colors = c("red" = "red", "white" = "white", "lightblue" = "lightblue", 
                  "blue" = "blue", "darkblue" = "darkblue")

p_plot = ggplot() + geom_sf(data = conus, fill = "white") + theme_bw() + 
  geom_point(data = diff, aes(x = gauge_lon, y = gauge_lat, fill = color_bin), size = 3, shape = 21) +
  labs(x = element_blank(), y = element_blank()) + 
  scale_fill_manual(values = custom_colors, breaks = c("red", "white", "lightblue", "blue", "darkblue"), 
                    labels = c("-0.25 to 0", "0", "0 to 0.25", "0.25 to 0.50", "0.50 to 0.75")) +
   theme(legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8), 
        title = element_text(size = 8), axis.text = element_text(size = 7.5)) + 
  labs(title = "Mean bias (mm/day)")

p = ggplot() + geom_sf(data = conus, fill = "lightgray") + theme_bw() + 
   geom_point(data = diff, aes(x = gauge_lon, y = gauge_lat, color = PRCP), size = 3) +
   labs(x = element_blank(), y = element_blank()) + scale_color_gradient2(low = "red", mid = "white", high = "blue") +
   theme(legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8), 
         title = element_text(size = 8), axis.text = element_text(size = 7.5)) + 
   labs(title = expression("Mean precipitation bias (mm" * " day" ^-1 * ")"))

t = ggplot() + geom_sf(data = conus, fill = "lightgray") + theme_bw() + 
   geom_point(data = diff, aes(x = gauge_lon, y = gauge_lat, color = TAIR), size = 3) +
   labs(x = element_blank(), y = element_blank()) + scale_color_gradient2(low = "red", mid = "white", high = "blue") +
   theme(legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8), 
         title = element_text(size = 8), axis.text = element_text(size = 7.5)) + 
   labs(title = expression("Mean temperature bias (" * "\u00B0" * "C day" ^-1 * ")"))
  
p + t + plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

#############################################################################
## MAE ##
#############################################################################
# metrics package
mae = daily_data %>% group_by(climate) %>% summarise(PRCP = mae(PRCP_C, PRCP_M), 
      TAIR = mae(TAIR_C, TAIR_M))

# get monthly mean value by climate region (240 months per region)
month2 = month %>% group_by(WATERYR, MNTH, climate) %>% reframe(PRCP_C = mean(PRCP_C), PRCP_M = mean(PRCP_M), 
                                                           TAIR_C = mean(TAIR_C), TAIR_M = mean(TAIR_M))

# mae over all months per climate region (12 months per region)
mae_month = month2 %>% group_by(MNTH, climate) %>% summarise(PRCP = mae(PRCP_C, PRCP_M), 
      TAIR = mae(TAIR_C, TAIR_M))

mae_month = mae_month %>% mutate(NAME = month.abb[MNTH])
mae_month$NAME = factor(mae_month$NAME, levels = month.abb)

season2 = season %>% group_by(SEASON, climate) %>% reframe(PRCP_C = mean(PRCP_C), PRCP_M = mean(PRCP_M), 
                                                           TAIR_C = mean(TAIR_C), TAIR_M = mean(TAIR_M))
mae_season = season2 %>% group_by(SEASON, climate) %>% summarise(PRCP = mae(PRCP_C, PRCP_M), 
      TAIR = mae(TAIR_C, TAIR_M))


year2 = year %>% group_by(WATERYR, climate) %>% reframe(PRCP_C = mean(PRCP_C), PRCP_M = mean(PRCP_M), 
                                                           TAIR_C = mean(TAIR_C), TAIR_M = mean(TAIR_M))
mae_year = year2 %>% group_by(WATERYR, climate) %>% summarise(PRCP = mae(PRCP_C, PRCP_M), 
      TAIR = mae(TAIR_C, TAIR_M))

m_plot = ggplot(data = mae_month, aes(x = NAME, y = TAIR, fill = climate)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
         axis.text = element_text(size = 8), legend.title = element_blank()) +
  #  labs(y = "Monthly MAE (mm)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))
    labs(y = "Monthly MAE (\u00b0C)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))

s_plot = ggplot(data = mae_season, aes(x = SEASON, y = TAIR, fill = climate)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
         axis.text = element_text(size = 8), legend.title = element_blank()) +
   # labs(y = "Seasonal MAE (mm)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))
 labs(y = "Seasonal MAE (\u00b0C)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))
 
a_plot = ggplot(data = mae_year, aes(x = WATERYR, y = TAIR, fill = climate)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) + theme_classic() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), 
         axis.text = element_text(size = 8), legend.title = element_blank()) +
 #   labs(y = "Annual MAE (mm)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))
   labs(y = "Annual MAE (\u00b0C)") + scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73"))

m_plot + s_plot + a_plot + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
  subtitle = "Temperature mean absolute error") + plot_layout(nrow = 3, guides = "collect")

mnth_mae = mae_month %>% group_by(climate) %>% reframe(PRCP = mean(PRCP), TAIR = mean(TAIR))
season_mae = mae_season %>% group_by(climate) %>% reframe(PRCP = mean(PRCP), TAIR = mean(TAIR))
yr_mae = mae_year %>% group_by(climate) %>% reframe(PRCP = mean(PRCP), TAIR = mean(TAIR))

#################################################################
## BOOTSTRAPPING
#################################################################
# Hmisc package 
# precipitation bootstrapping 
TAIR_c_month = month %>% select(climate, TAIR_C, MNTH) %>% 
  group_by(climate, MNTH) %>% 
  group_map(~ smean.cl.boot(., conf.int = .95, B = 10000)) %>% 
  bind_rows()
TAIR_m_month = month %>% select(climate, TAIR_M, MNTH) %>% 
  group_by(climate, MNTH) %>% 
  group_map(~ smean.cl.boot(., conf.int = .95, B = 10000)) %>% 
  bind_rows()

# use mean of monthly totals per climate region (1 - 12 months = 36 months with all 3 regions)
# so can obtain the correct number of climate and month rows to join with the bootstrapping results
m = month %>% group_by(climate, MNTH) %>% 
  reframe(TAIR_C = mean(TAIR_C), TAIR_M = mean(TAIR_M))

clim_mnth = m %>% select(climate, MNTH) # select just climate and month columns 
TAIR_c_mean = cbind(clim_mnth, TAIR_c_month) # join climate and month columns to camels bootstrap
TAIR_m_mean = cbind(clim_mnth, TAIR_m_month) # mopex bootstrap
TAIR_mean = full_join(TAIR_c_mean, TAIR_m_mean, by = c("climate", "MNTH"), keep = TRUE) # keep all to prevent NA

# rename the columns
TAIR_mean = TAIR_mean %>% rename(
  climate_C = climate.x, MNTH_C = MNTH.x, Mean_C = Mean.x, Lower_C = Lower.x, Upper_C = Upper.x, 
  climate_M = climate.y, MNTH_M = MNTH.y, Mean_M = Mean.y, Lower_M = Lower.y, Upper_M = Upper.y)


# adjust the data frame so that columns are reorganized for plotting
mTAIR_mean = pivot_longer(TAIR_mean, cols = everything(),
                         names_to = c(".value", "Dataset"), 
                         names_sep = "_")

# add column for month abbreviations based on month number
mTAIR_mean = mTAIR_mean %>% mutate(NAME = month.abb[MNTH])
# change month names to factor in order to display on plot in chronological order instead of alphabetically
mTAIR_mean$NAME = factor(mTAIR_mean$NAME, levels = month.abb)

mArid = mTAIR_mean %>% filter(climate == "Arid")
mCont = mTAIR_mean %>% filter(climate == "Continental")
mTemp = mTAIR_mean %>% filter(climate == "Temperate")

# create plot with mean value and confidence intervals for each month by climate region
temp = ggplot(data = mTemp, aes(x = NAME, y = Mean, color = Dataset)) +
  geom_point(size = 2) + scale_color_manual(values = c("blue", "red"), labels = c("CAMELS", "MOPEX")) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), stat = "identity") + ylim(30, 145) +
  labs(y = "Precipitation (mm)", title = "Temperate") + theme_bw() + 
  scale_x_discrete(labels = c('J','F','M','A','M','J','J','A','S','O','N','D')) +
  theme(axis.title = element_blank(), axis.text = element_text(size = 8), #axis.title.y = element_text(size = 8), 
  legend.title = element_blank(), title = element_text(size = 10), legend.position = c(0.9,0.25)) 

# temperature plot 
cont = ggplot(data = mCont, aes(x = NAME, y = Mean, color = Dataset)) +
  geom_point(size = 2) + scale_color_manual(values = c("blue", "red"), labels = c("CAMELS", "MOPEX")) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), stat = "identity") + ylim(-5,30) +
  labs(y = "Temperature (\u00b0C)", title = "Continental") + theme_bw() + 
  scale_x_discrete(labels = c('J','F','M','A','M','J','J','A','S','O','N','D')) +
  theme(axis.title = element_blank(), axis.text = element_text(size = 8), #axis.title.y = element_text(size = 8), 
  legend.title = element_blank(), title = element_text(size = 10), legend.position = c(0.9,0.25)) 

arid + cont + temp + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
   subtitle = "Average monthly temperature") + 
   plot_layout(nrow = 1, guides = "collect") #&
   #theme(plot.tag = element_text(size = 8))

## Seasonal bootstrapping
TAIR_c_season = season %>% select(climate, TAIR_C, SEASON) %>% 
  group_by(climate, SEASON) %>% 
  group_map(~ smean.cl.boot(., conf.int = .95, B = 10000)) %>% 
  bind_rows()
TAIR_m_season = season %>% select(climate, TAIR_M, SEASON) %>% 
  group_by(climate, SEASON) %>% 
  group_map(~ smean.cl.boot(., conf.int = .95, B = 10000)) %>% 
  bind_rows()

# use mean of seasonal totals per climate region 
# so can obtain the correct number of climate and season rows to join with the bootstrapping results
s = season %>% group_by(climate, SEASON) %>% 
  reframe(TAIR_C = mean(TAIR_C), TAIR_M = mean(TAIR_M))

clim_season = s %>% select(climate, SEASON) # select just climate and season
TAIR_c_mean = cbind(clim_season, TAIR_c_season) # join climate and season columns to camels bootstrap
TAIR_m_mean = cbind(clim_season, TAIR_m_season) # mopex bootstrap
TAIR_mean = full_join(TAIR_c_mean, TAIR_m_mean, by = c("climate", "SEASON"), keep = TRUE) # keep all to prevent NA

# rename the columns
TAIR_mean = TAIR_mean %>% rename(
  climate_C = climate.x, SEASON_C = SEASON.x, Mean_C = Mean.x, Lower_C = Lower.x, Upper_C = Upper.x, 
  climate_M = climate.y, SEASON_M = SEASON.y, Mean_M = Mean.y, Lower_M = Lower.y, Upper_M = Upper.y)

# adjust the data frame so that columns are reorganized for plotting
sTAIR_mean = pivot_longer(TAIR_mean, cols = everything(),
                         names_to = c(".value", "Dataset"), 
                         names_sep = "_")

sArid = sTAIR_mean %>% filter(climate == "Arid")
sCont = sTAIR_mean %>% filter(climate == "Continental")
sTemp = sTAIR_mean %>% filter(climate == "Temperate")

# create plot with mean value and confidence intervals for each year by climate region
s_Cont = ggplot(data = sCont, aes(x = SEASON, y = Mean, color = Dataset)) +
  geom_point(size = 2) + scale_color_manual(values = c("blue", "red"), labels = c("CAMELS", "MOPEX")) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), stat = "identity") + ylim(100, 400) +
  labs(y = "Precipitation (mm)", title = "Continental") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_text(size = 8), #axis.title.y = element_text(size = 8), 
        legend.title = element_blank(), title = element_text(size = 10), 
        legend.position = c(0.9,0.25)) 

# temperature plot
s_Temp = ggplot(data = sTemp, aes(x = SEASON, y = Mean, color = Dataset)) +
  geom_point(size = 2) + scale_color_manual(values = c("blue", "red"), labels = c("CAMELS", "MOPEX")) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), stat = "identity") + ylim(-2,25) +
  labs(y = "Temperature (\u00b0C)", title = "Temperate") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_text(size = 8), #axis.title.y = element_text(size = 8), 
        legend.title = element_blank(), title = element_text(size = 10), 
        legend.position = c(0.9,0.25))  

s_Arid +s_Cont + s_Temp + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
   subtitle = "Average seasonal temperature") + 
   plot_layout(nrow = 1, guides = "collect") #&
   #theme(plot.tag = element_text(size = 8))

#arid + cont + temp + plot_annotation(tag_levels = list(c('(c)','(d)', '(e)')),
#  subtitle = "Average total seasonal precipitation") + 
#  plot_layout(nrow = 1, guides = "collect") &
 # theme(plot.tag = element_text(size = 8))

## Annual bootstrapping
TAIR_c_year = year %>% select(climate, TAIR_C, WATERYR) %>% 
  group_by(climate, WATERYR) %>% 
  group_map(~ smean.cl.boot(., conf.int = .95, B = 10000)) %>% 
  bind_rows()
TAIR_m_year = year %>% select(climate, TAIR_M, WATERYR) %>% 
  group_by(climate, WATERYR) %>% 
  group_map(~ smean.cl.boot(., conf.int = .95, B = 10000)) %>% 
  bind_rows()

# use mean of annual totals per climate region 
# so can obtain the correct number of climate and season rows to join with the bootstrapping results
y = year %>% group_by(climate, WATERYR) %>% 
  reframe(TAIR_C = mean(TAIR_C), TAIR_M = mean(TAIR_M))

clim_year = y %>% select(climate, WATERYR) # select just climate and season
TAIR_c_mean = cbind(clim_year, TAIR_c_year) # join climate and season columns to camels bootstrap
TAIR_m_mean = cbind(clim_year, TAIR_m_year) # mopex bootstrap
TAIR_mean = full_join(TAIR_c_mean, TAIR_m_mean, by = c("climate", "WATERYR"), keep = TRUE) # keep all to prevent NA

# rename the columns
TAIR_mean = TAIR_mean %>% rename(
  climate_C = climate.x, WATERYR_C = WATERYR.x, Mean_C = Mean.x, Lower_C = Lower.x, Upper_C = Upper.x, 
  climate_M = climate.y, WATERYR_M = WATERYR.y, Mean_M = Mean.y, Lower_M = Lower.y, Upper_M = Upper.y)

# adjust the data frame so that columns are reorganized for plotting
yTAIR_mean = pivot_longer(TAIR_mean, cols = everything(),
                         names_to = c(".value", "Dataset"), 
                         names_sep = "_")

yArid = yTAIR_mean %>% filter(climate == "Arid")
yCont = yTAIR_mean %>% filter(climate == "Continental")
yTemp = yTAIR_mean %>% filter(climate == "Temperate")

# create plot with mean value and confidence intervals for each year by climate region
y_Arid = ggplot(data = yArid, aes(x = WATERYR, y = Mean, color = Dataset)) +
  geom_point(size = 2) + scale_color_manual(values = c("blue", "red"), labels = c("CAMELS", "MOPEX")) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), stat = "identity") + ylim(300, 1750) +
  labs(y = "Precipitation (mm)", title = "Arid") + theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text = element_text(size = 8), axis.title.y = element_text(size = 8), 
        legend.title = element_blank(), title = element_text(size = 10), 
        legend.position = c(0.9,0.25)) 

# temperature plot
y_Arid = ggplot(data = yArid, aes(x = WATERYR, y = Mean, color = Dataset)) +
  geom_point(size = 2) + scale_color_manual(values = c("blue", "red"), labels = c("CAMELS", "MOPEX")) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), stat = "identity") + ylim(5,25) +
  labs(y = "Temperature (\u00b0C)", title = "Arid") + theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text = element_text(size = 8), axis.title.y = element_text(size = 8), 
        legend.title = element_blank(), title = element_text(size = 10), 
        legend.position = c(0.9,0.25))  

  y_Arid + y_Cont + y_Temp +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
   subtitle = "Average annual temperature") + 
   plot_layout(nrow = 1, guides = "collect")

#################################################################
## STANDARD ERROR OF THE DIFFERENCE OF THE MEANS
#################################################################
# find standard error of difference of means which is the square root of sd^2/n + sd^2/n 
# confidence intervals based on the 95 percent confidence interval
# and critical value of 1.96 from alpha 0.05
# CI = (X1 - X2) +- critical value * standard error

se = year %>% group_by(climate) %>% 
  reframe(PRCP_Cn = length(PRCP_C), PRCP_Mn = length(PRCP_M), # number of values
          PRCP_Csd = sd(PRCP_C), PRCP_Msd = sd(PRCP_M),  # standard deviation
          PRCP_Cmean = mean(PRCP_C), PRCP_Mmean = mean(PRCP_M), # mean value
          PRCP_diff = abs(PRCP_Cmean - PRCP_Mmean), # difference between mean values
          PRCPse_diff = sqrt((PRCP_Csd^2/PRCP_Cn)+(PRCP_Msd^2/PRCP_Mn)), # se of diff of means
          PRCP_ci = 1.96*PRCPse_diff, # confidence interval (MOE) based on se diff
          PRCP_Upper = PRCP_diff+PRCP_ci, PRCP_Lower = PRCP_diff-PRCP_ci, # MOE with confidence intervals use mean difference
          TAIR_Cn = length(TAIR_C), TAIR_Mn = length(TAIR_M), 
          TAIR_Csd = sd(TAIR_C), TAIR_Msd = sd(TAIR_M), 
          TAIR_Cmean = mean(TAIR_C), TAIR_Mmean = mean(TAIR_M), 
          TAIR_diff = abs(TAIR_Cmean - TAIR_Mmean),
          TAIRse_diff = sqrt((TAIR_Csd^2/TAIR_Cn)+(TAIR_Msd^2/TAIR_Mn)), 
          TAIR_ci = 1.96*TAIRse_diff,
          TAIR_Upper = TAIR_diff+TAIR_ci, TAIR_Lower = TAIR_diff-TAIR_ci)

#################################################################
## WILCOXON TEST
#################################################################
# two sample test 
prcp = daily %>% group_by(GaugeID, climate) %>% reframe(PRCP_C = PRCP_C, PRCP_M = PRCP_M)
prcp = month %>% group_by(GaugeID, MNTH, climate) %>% reframe(PRCP_C = PRCP_C, PRCP_M = PRCP_M)
prcp = season %>% group_by(GaugeID, WATERYR, SEASON, climate) %>% reframe(PRCP_C = PRCP_C, PRCP_M = PRCP_M)
prcp = year %>% group_by(GaugeID, WATERYR, climate) %>% reframe(PRCP_C = PRCP_C, PRCP_M = PRCP_M)

# change table to long format 
prcp = pivot_longer(prcp, cols = c(PRCP_C, PRCP_M), names_to = "Variable", 
                    values_to = "Value")

prcp_wilcox = prcp %>% group_by(climate) %>% 
  wilcox_test(Value ~ Variable, paired = FALSE, detailed = TRUE)

prcp_wilcox = wilcox.test(Value ~ Variable, data = prcp_cont)

#################################################################
## FLIGNER-KILLEEN TEST
#################################################################
# homogeneity of variance across groups 
# null hypothesis all population variances are equal
# alternate hypothesis is at least one sample has different variance
# test the null hypothesis at 0.05 significance level (95% percentile)
# dependent variable ~ independent variable (grouping variable)

prcp_arid = prcp %>% filter(climate == "Arid")
prcp_cont = prcp %>% filter(climate == "Continental")
prcp_temp = prcp %>% filter(climate == "Temperate")

result = fligner.test(Value ~ Variable, data = prcp_temp)
                      
#################################################################
## WELCH'S T-TEST
#################################################################
# use pivoted table 

ttest_prcp = prcp %>% group_by(WATERYR, climate) %>% 
 t_test(data = ., Value ~ Variable, conf.level = 0.95, paired = FALSE)

#################################################################
## SPEARMAN RANK
#################################################################
spear_cor = year %>% group_by(WATERYR, climate) %>% 
  summarise(TAIR = cor(TAIR_C, TAIR_M, method = "spearman"), 
            PRCP = cor(PRCP_C, PRCP_M, method = "spearman")) %>% 
  mutate(TAIR = round(TAIR, digits = 3), PRCP = round(PRCP, digits = 3))

#################################################################
## VALIDATION 
#################################################################
# separate climate regions
arid = daily_data %>% filter(climate == "Arid")
cont = daily_data %>% filter(climate == "Continental")
temp = daily_data %>% filter(climate == "Temperate")

# daily data (repeat with monthly, seasonal, annual)
# add column to indicate camels dataset (0)
# rename columns in both dataframes to match
daily_C = temp %>% select(PRCP_C, TAIR_C, ET_WB_C) %>% mutate(Dataset = 0) %>% 
  rename(PRCP = PRCP_C, TAIR = TAIR_C, ET = ET_WB_C)

# select precipitation and temperature data and add column to indicate mopex dataset (1)
# rename columns in both dataframes to match
daily_M = temp %>% select(PRCP_M, TAIR_M, ET_M) %>% mutate(Dataset = 1) %>% 
  rename(PRCP = PRCP_M, TAIR = TAIR_M, ET = ET_M)

# rbind to append dataframes together
daily = rbind(daily_C, daily_M)
# encode the target feature (dataset) as factor
daily$Dataset = factor(daily$Dataset, levels = c(0,1))

# split the dataset into training and test sets (caTools package) based on dataset column 
split =  sample.split(daily$Dataset, SplitRatio = 0.75) # ratio indicates 75% used for train, 25% for test

# split variable results in logical vector of TRUE for training and FALSE for testing
training_set =  subset(daily, split == TRUE) 
test_set =  subset(daily, split == FALSE) 

# fit SVM to training dataset
classifier =  svm(formula = Dataset ~ ., 
                 data = training_set, 
                 type = 'C-classification', 
                 kernel = 'linear') 

# predict (classify) using test dataset based on trained model, remove dataset column   
y_pred = predict(classifier, newdata = test_set[-6])

#############################################################################
## RUNOFF EFFICIENCY ##
#############################################################################
runoff_eff = daily_data %>% group_by(GaugeID, WATERYR, climate) %>% 
  reframe(PRCP_C = sum(PRCP_C), PRCP_M = sum(PRCP_M), Q_C = sum(OBS_RUN_C), Q_M = sum(OBS_RUN_M), 
          RE_C = Q_C/PRCP_C, RE_M = Q_M/PRCP_M)

corr = cor(runoff_eff$RE_C, runoff_eff$RE_M, method = "spearman")

mean = runoff_eff %>% group_by(WATERYR, climate) %>% reframe(RE_C = round(mean(RE_C), 2), RE_M = round(mean(RE_M), 2))

model = lm(RE_M ~ RE_C, data = runoff_eff)
confidence_interval = predict(model, interval = "confidence")

ggplot(data = runoff_eff, aes(x = RE_C, y = RE_M, color = climate)) + 
  geom_point() +
  geom_smooth(method = "lm", color = "blue") + # Add a single linear trend line
  scale_color_manual(values = c("Arid" = "#D55E00", "Continental" = "#E69F00", "Temperate" = "#009E73")) +
  theme_bw() + theme(legend.title = element_blank(), legend.position = c(0.7, 0.25)) +
  labs(x = "CAMELS Runoff Efficiency", y = "MOPEX Runoff Efficiency")

gauge = unique(daily_data$GaugeID)
gauges = id_attributes %>% filter(GaugeID %in% gauge)
coords = select(gauges, c(GaugeID, gauge_lat, gauge_lon))
mean = left_join(mean, coords, by = "GaugeID")

#################################################################
## INDICES
#################################################################
prcp = dplyr::select(daily_data, c(GaugeID, DATE, DY, MNTH, WATERYR, climate, PRCP_C, PRCP_M))

# EXTREMELY WET DAYS R99P
# annual total precipitation from days >99th percentile (mm)
# filter days where precipitation exceeds the 99th percentile value
R99p_C = prcp %>% group_by(GaugeID) %>% filter(PRCP_C > quantile(PRCP_C, 0.99))
# get annual totals for days > 99th percentile for each gauge
R99p_C = R99p_C %>% group_by(climate, GaugeID, WATERYR) %>% summarise(PRCP_C = sum(PRCP_C))
# average annual totals among gauges within each climate to get average annual total (mm) 
# for days where precipitation >99th percentile
R99p_C = R99p_C %>% group_by(climate, WATERYR) %>% summarise(PRCP_C = mean(PRCP_C))

R99p_M = prcp %>% group_by(GaugeID) %>% filter(PRCP_M > quantile(PRCP_M, 0.99))
# get annual totals for days > 99th percentile for each gauge
R99p_M = R99p_M %>% group_by(climate, GaugeID, WATERYR) %>% summarise(PRCP_M = sum(PRCP_M))
# average annual totals among gauges within each climate to get average annual total (mm) 
# for days where precipitation >99th percentile
R99p_M = R99p_M %>% group_by(climate, WATERYR) %>% summarise(PRCP_M = mean(PRCP_M))

R99p = left_join(R99p_C, R99p_M, by = c("climate" = "climate", "WATERYR" = "WATERYR"))

R99p_long = R99p %>% pivot_longer(cols = c(PRCP_C, PRCP_M), names_to = "Dataset", values_to = "PRCP") %>% 
  mutate(Dataset = recode(Dataset, "PRCP_C" = "CAMELS", "PRCP_M" = "MOPEX"))

# three histogram plots 
arid = R99p_long %>% filter(climate == "Arid")
a_plot = ggplot(data = arid, aes(x = WATERYR, y = PRCP, fill = Dataset)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + scale_fill_manual(name = "Dataset", values = c("CAMELS" = "blue", "MOPEX" = "red")) +
  labs(y = "Annual total wet day precipitation (mm)", title = "Arid") + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text = element_text(size = 8), axis.title.y = element_text(size = 10), 
        legend.title = element_blank(), title = element_text(size = 8), 
        legend.position = c(0.3,0.9)) 

cont = R99p_long %>% filter(climate == "Continental")
c_plot = ggplot(data = cont, aes(x = WATERYR, y = PRCP, fill = Dataset)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + scale_fill_manual(name = "Dataset", values = c("CAMELS" = "blue", "MOPEX" = "red")) +
  labs(y = "Annual total wet day precipitation (mm)", title = "Continental") + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text = element_text(size = 8), axis.title.y = element_text(size = 10), 
        legend.title = element_blank(), title = element_text(size = 8), 
        legend.position = c(0.3,0.9))

temp = R99p_long %>% filter(climate == "Temperate")
t_plot = ggplot(data = temp, aes(x = WATERYR, y = PRCP, fill = Dataset)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + scale_fill_manual(name = "Dataset", values = c("CAMELS" = "blue", "MOPEX" = "red")) +
  labs(y = "Annual total wet day precipitation (mm)", title = "Temperate") + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text = element_text(size = 8), axis.title.y = element_text(size = 10), 
        legend.title = element_blank(), title = element_text(size = 8), 
        legend.position = c(0.3,0.9))

  a_plot + c_plot + t_plot +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
   caption = "Daily precipitation is greater than 99th percentile (R99p)", 
  theme = theme(plot.caption = element_text(size = 8))) + 
   plot_layout(nrow = 3, guides = "collect", axis_titles = "collect")

# R10mm 
# number of days when precipitation is heavy (greater than or equal to 10 mm)
# find number of days per year 
hd_C = prcp %>% group_by(GaugeID, WATERYR, climate) %>% tally(PRCP_C >= 10)
hd_C1 = hd_C %>% group_by(climate, WATERYR) %>% summarise(hd_C = mean(n))
hd_M = prcp %>% group_by(GaugeID, WATERYR, climate) %>% tally(PRCP_M >= 10)
hd_M1 = hd_M %>% group_by(climate, WATERYR) %>% summarise(hd_M = mean(n))

hd = left_join(hd_C1, hd_M1, by = c("climate" = "climate", "WATERYR" = "WATERYR"))

# number of heavy precipitation days plot                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
hd_plot = ggplot(data = hd) +  geom_line(aes(x = WATERYR, y = hd_C, color = climate, linetype = "CAMELS"), linewidth = 1) +
  geom_line(aes(x = WATERYR, y = hd_M, color = climate, linetype = "MOPEX"), linewidth = 1) + 
  theme_bw() + scale_color_manual(name = "Climate", 
        values = c("Arid" = "#D55E00", "Continental" = "#E69F00", "Temperate" = "#009E73")) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8), 
  plot.caption = element_text(size = 8), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
   scale_linetype_manual(name = "Dataset", values = c("CAMELS" = "solid", "MOPEX" =  "dashed")) +
  labs(x = element_blank(), y = "Number of heavy precipitation days per year") +
  labs(caption = "Daily precipitation is greater than or equal to 10 mm") 
  
# number of dry days per year (less than 1 mm)
dd_C = prcp %>% group_by(GaugeID, WATERYR, climate) %>% tally(PRCP_C < 1)
# average 51 gauges together to find average annual number of dry days
dd_C1 = dd_C %>% group_by(climate, WATERYR) %>% summarise(dd_C = mean(n))
dd_M = prcp %>% group_by(GaugeID, WATERYR, climate) %>% tally(PRCP_M < 1)
dd_M1 = dd_M %>% group_by(climate, WATERYR) %>% summarise(dd_M = mean(n))

dd = left_join(dd_C1, dd_M1, by = c("climate" = "climate", "WATERYR" = "WATERYR"))

# number of dry days plot
dd_plot = ggplot(data = dd) + geom_line(aes(x = WATERYR, y = dd_C, color = climate, linetype = "CAMELS"), linewidth = 1) +
  geom_line(aes(x = WATERYR, y = dd_M, color = climate, linetype = "MOPEX"), linewidth = 1) + 
  theme_bw() + scale_color_manual(name = "Climate", 
      values = c("Arid" = "#D55E00", "Continental" = "#E69F00", "Temperate" = "#009E73")) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8), 
  plot.caption = element_text(size = 8), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
   scale_linetype_manual(name = "Dataset", values = c("CAMELS" = "solid", "MOPEX" =  "dashed")) +
  labs(x = element_blank(), y = "Number of dry days per year") +
  labs(caption = "Daily precipitation is less than 1 mm") 
  
hd_plot + dd_plot +
plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
plot_layout(nrow = 2, guides = "collect")

# number of frost days, annual count of days when TN (daily min temperature) <0 degrees C
temp = dplyr::select(daily_data, c(GaugeID, DATE, DY, MNTH, WATERYR, climate, TAIR_C, TAIR_M))

TN_C = temp %>% group_by(GaugeID, WATERYR, climate) %>% tally(TAIR_C < 0)
TN_C1 = TN_C %>% group_by(climate, WATERYR) %>% summarise(TN_C = mean(n))
TN_C1$TN_C = round(TN_C1$TN_C, 0)

TN_M = temp %>% group_by(GaugeID, WATERYR, climate) %>% tally(TAIR_M < 0)
TN_M1 = TN_M %>% group_by(climate, WATERYR) %>% summarise(TN_M = mean(n))
TN_M1$TN_M = round(TN_M1$TN_M, 0)

TN = left_join(TN_C1, TN_M1, by = c("climate" = "climate", "WATERYR" = "WATERYR"))

TN_plot = ggplot(data = TN) +  geom_line(aes(x = WATERYR, y = TN_C, color = climate, linetype = "CAMELS"), linewidth = 1) +
  geom_line(aes(x = WATERYR, y = TN_M, color = climate, linetype = "MOPEX"), linewidth = 1) + 
  theme_bw() + scale_color_manual(name = "Climate", 
                values = c("Arid" = "#D55E00", "Continental" = "#E69F00", "Temperate" = "#009E73")) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8), 
  plot.caption = element_text(size = 8), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
   scale_linetype_manual(name = "Dataset", values = c("CAMELS" = "solid", "MOPEX" =  "dashed")) +
  labs(x = element_blank(), y = "Number of frost days per year") +
  labs(caption = "Daily temperature is less than 0 degrees Celcius")

TX_C = temp %>% group_by(GaugeID, WATERYR, climate) %>% tally(TAIR_C > 25)
TX_C1 = TX_C %>% group_by(climate, WATERYR) %>% summarise(TX_C = mean(n))
TX_C1$TX_C = round(TX_C1$TX_C, 0)

TX_M = temp %>% group_by(GaugeID, WATERYR, climate) %>% tally(TAIR_M > 25)
TX_M1 = TX_M %>% group_by(climate, WATERYR) %>% summarise(TX_M = mean(n))
TX_M1$TX_M = round(TX_M1$TX_M, 0)

TX = left_join(TX_C1, TX_M1, by = c("climate" = "climate", "WATERYR" = "WATERYR"))

TX_plot = ggplot(data = TX) +  geom_line(aes(x = WATERYR, y = TX_C, color = climate, linetype = "CAMELS"), linewidth = 1) +
  geom_line(aes(x = WATERYR, y = TX_M, color = climate, linetype = "MOPEX"), linewidth = 1) + 
  theme_bw() + scale_color_manual(name = "Climate", 
            values = c("Arid" = "#D55E00", "Continental" = "#E69F00", "Temperate" = "#009E73")) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8), 
  plot.caption = element_text(size = 8), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
   scale_linetype_manual(name = "Dataset", values = c("CAMELS" = "solid", "MOPEX" =  "dashed")) +
  labs(x = element_blank(), y = "Number of summer days per year") +
  labs(caption = "Daily temperature is greater than 25 degrees Celcius")

  TN_plot + TX_plot +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
   plot_layout(nrow = 2, guides = "collect")
  
