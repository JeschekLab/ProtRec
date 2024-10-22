# load libraries
library(tidyverse)
library(readxl)
library(viridis)
library(cowplot)

#define which package should be preferred
library(conflicted)
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")

##########################################################
### read data ###
##########################################################

# read necessary info to understand raw file
config <- read_excel("./config.xlsx",
                          sheet = "script_info",
                          col_names = TRUE)
meta <- read_excel("./config.xlsx",
                   sheet = "sample_info",
                   col_names = TRUE)

### pull information from config file
# pull name of raw file 
raw_file <- filter(config, info == "raw_data") %>% pull(content) 
pro_blank <- filter(config, info == "protease_blank") %>% pull(content)
neg_ctrl <- filter(config, info == "catalysis_neg_ctrl") %>% pull(content)

# extract signal values from raw data
data_CFP <- read_excel(paste0("./", raw_file),
                       range = filter(config, info == "CFP") %>% pull(content),
                       col_names = TRUE)[-c(1:2), ] # remove time and temp row

data_YFP <- read_excel(paste0("./", raw_file),
                       range = filter(config, info == "YFP") %>% pull(content),
                       col_names = TRUE)[-c(1:2), ] # remove time and temp row

data_RFP <- read_excel(paste0("./", raw_file),
                       range = filter(config, info == "RFP") %>% pull(content),
                       col_names = TRUE)[-c(1:2), ] # remove time and temp row

### replace cycle nr. by identical time [s] values 
# (this is necessary as in rare cases some signals are recorded with a delay of 
# 0.1 s compared to others, which prevents a smooth running of the code)
# pull time values
time <- read_excel(paste0("./", raw_file),
                   range = filter(config, info == "time") %>% pull(content),
                   col_names = FALSE)
# replace headers
colnames(data_CFP)=time[c(1),]
colnames(data_YFP)=time[c(1),]
colnames(data_RFP)=time[c(1),]

# make long format and join with meta data
data <- rbind(
  data_CFP %>%
    pivot_longer(
      cols = -'Time [s]',
      names_to = 'time',
      values_to = 'signal') %>%
    rename(wel = 'Time [s]') %>%
    mutate(fluo = 'CFP'),
  data_YFP %>%
    pivot_longer(
      cols = -'Time [s]',
      names_to = 'time',
      values_to = 'signal') %>%
    rename(wel = 'Time [s]') %>%
    mutate(fluo = 'YFP')) %>%
  inner_join(., meta, by = 'wel') %>%
  select(-wel) %>%
  mutate(time = as.numeric(time) / 60) %>% # make time in minutes
  filter(sub != "non" | pro != "non") # remove data from empty wells


##########################################################
### calculate CFP/YFP ###
##########################################################
# calculate CFP/YFP
data_CY <- data %>%
  pivot_wider(
    names_from = fluo,
    values_from = signal) %>%
  mutate(CY = CFP / YFP) %>%
  select(-CFP, -YFP)

##########################################################
### blank CFP/YFP with catalysis neg. ctrl. ###
##########################################################
# extract buffer value
negctrl <- data_CY %>%
  filter(pro == neg_ctrl) %>%
  select(-pro) %>%
  rename('neg_ctrl' = 'CY')

# join data with buffer and subtract buffer
data_deltaneg <- data_CY %>%
  inner_join(., negctrl, by = c('time', 'sub')) %>%
  mutate(deltaneg = CY - neg_ctrl) %>%
  select(-CY, -neg_ctrl)

##########################################################
### calculate signal relative to proteinase K ###
##########################################################
# extract proteinase K signal
protK <- data_deltaneg %>%
  filter(pro == 'protK') %>%
  select(-pro) %>%
  rename('protK' = 'deltaneg')

# join data with protK data and divide
data_protK <- data_deltaneg %>%
  inner_join(., protK, by = c('time', 'sub')) %>%
  mutate(activity = deltaneg / protK * 100) %>%
  select(-deltaneg, -protK)

##########################################################
### calculate initial velocity (vini) ###
##########################################################
data_lm <- data_protK %>%
  filter(pro != 'protK') %>%
# filter(activity < 80)
  filter(time < 120) # only use data points within the first 120 min

# initialize empty data frame
data_vini <- data.frame()

# loop through substrates and proteases and calculate linear model
for (TEVs in unique(data_lm$sub)) {
  for (TEVp in unique(data_lm$pro)) {
    print(paste(TEVp, 'with', TEVs))
    # extract data for current TEVp and TEVs
    temp_data <- data_lm %>%
      filter(sub == TEVs, pro == TEVp) %>%
      mutate(time = time)
    if (nrow(temp_data) <= 2) {
      print(paste0('Missing data points in ', TEVp, ' with ', TEVs))
      next
    }
    # generate linear model
    model <- lm(data = temp_data, formula = activity~time)
    # extract slope and intercept
    slope <- model$coefficients[2]
    intercept <- model$coefficients[1]
    r_squared <- summary(model)$r.squared
    # make temporary data frame
    temp_df <- data.frame(
      pro = TEVp, sub = TEVs,
      vini = slope, intercept = intercept,
      r_squared = r_squared)
    # combine
    data_vini <- rbind(data_vini, temp_df)
  }
}

# save this table
dir.create("./results")
write.table(data_vini, "./results/data_vini.txt", sep = "\t", row.names=FALSE, quote = F)

##########################################################
### calculate vini relative to wildtype TEVs (vinirel) ###
##########################################################
# extract wildtype
data_WT <- data_vini %>%
  filter(sub == 'S') %>%
  select(-sub, -intercept) %>%
  rename('viniWT' = 'vini')

# join data with protK data and divide
data_vinirel <- data_vini %>%
  inner_join(., data_WT, by = 'pro') %>%
  mutate(vinirel = vini / viniWT) %>%
  select(-viniWT) %>%
  arrange(pro, sub)

##########################################################
### calculate vini relative to parent TEVp (vinirel_parent) ###
##########################################################
# extract parent
data_parent <- data_vini %>%
  filter(pro == 'TEVp I') %>%
  select(-pro, -intercept) %>%
  rename('viniparent' = 'vini')

# join data with protK data and divide
data_vinirel_parent <- data_vini %>%
  inner_join(., data_parent, by = 'sub') %>%
  mutate(vinirel_parent = vini / viniparent) %>%
  select(-viniparent) %>%
  arrange(pro, sub)

##########################################################
### generate plots ###
##########################################################

# create a folder to store the individual  plots in
dir.create("./plots")

##########################################################
### heat map with vinirel
# define what to exclude from plotting
excl <- c("protK", pro_blank, neg_ctrl, "buffer")
# plot
p <- ggplot(filter(data_vinirel, !(pro %in% excl)),
            aes(x = sub, y = pro, fill = vinirel)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_cowplot(font_size = 12) +
  theme(axis.line = element_blank()) +
  theme(axis.title.y = element_blank()) +
  xlab("AA in P1'")

# show plot
print(p)
# save plot
save_plot(p, filename = "./plots/heatmap_vinirel.pdf")

##########################################################
### heat map with vinirel_parent
# define what to exclude from plotting
excl <- c("protK", pro_blank, neg_ctrl, "buffer", "TEVp I ATL")
# plot
p <- ggplot(filter(data_vinirel_parent, !(pro %in% excl)),
            aes(x = sub, y = pro, fill = vinirel_parent)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_cowplot(font_size = 12) +
  theme(axis.line = element_blank()) +
  theme(axis.title.y = element_blank()) +
  xlab("AA in P1'")

# show plot
print(p)
# save plot
save_plot(p, filename = "./plots/heatmap_vinirel_parent.pdf")

##########################################################
### heat map with vini
# define what to exclude from plotting
excl <- c("protK", pro_blank, neg_ctrl, "buffer")
# plot
p <- ggplot(filter(data_vini, !(pro %in% excl)),
            aes(x = sub, y = pro, fill = vini)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_cowplot(font_size = 12) +
  theme(axis.line = element_blank()) +
  theme(axis.title.y = element_blank()) +
  xlab("AA in P1'")
# show plot
print(p)
# save plot
save_plot(p, filename = "./plots/heatmap_vini.pdf")

##########################################################
### activity per protease
# define what to exclude from plotting
excl <- c("protK", pro_blank, neg_ctrl, "buffer", "TEVp I GGT", "TEVp I VDS")
# plot
p <- ggplot(filter(data_protK, !(pro %in% excl)),
            aes(x = time, y = activity, color = sub)) +
  geom_point(size = 0.5) +
  geom_abline(filter(data_vini, !(pro %in% excl)),
              mapping = aes(slope = vini,
                            intercept = intercept, color = sub)) +
  theme_cowplot(font_size = 12) +
  ylab("Processed substrate (%)") +
  xlab("Time (min)") +
  facet_wrap(~pro) +
  theme(strip.background = element_rect(fill="#f0f0f0"))
# show plot
print(p)
# save plot
save_plot(p, filename = "./plots/activity_per_pro.pdf", base_asp = 1.2)



# Plot to show data processing steps --------------------------------------

##########################################################
### generate plots of data processing on the example of substrate H ###
##########################################################
### plot just raw CFP values
# plot
p1 <- ggplot(filter(data, fluo == "CFP", sub == "H", pro %in% c("buffer", "EV", "protK", "TEVp I", "TEVp I C151A")),
            aes(x = time, y = signal, color = pro)) +
  geom_point(size = 1) +
  theme_cowplot(font_size = 12) +
  ylab("CFP (AU)") +
  xlab("Time (min)") +
  labs(color = "Protease") +
  facet_wrap(~sub) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(fill="#f0f0f0"))
# show plot
print(p1)


### plot just raw YFP values
# plot
p2 <- ggplot(filter(data, fluo == "YFP", sub == "H", pro %in% c("buffer", "EV", "protK", "TEVp I", "TEVp I C151A")),
            aes(x = time, y = signal, color = pro)) +
  geom_point(size = 1) +
  theme_cowplot(font_size = 12) +
  ylab("YFP (AU)") +
  xlab("Time (min)") +
  labs(color = "Protease") +
  facet_wrap(~sub) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(fill="#f0f0f0"))
# show plot
print(p2)

### plot CFP/YFP values
# plot
p3 <- ggplot(data_CY %>% filter(sub == "H", pro %in% c("buffer", "EV", "protK", "TEVp I", "TEVp I C151A")),
            aes(x = time, y = CY, color = pro)) +
  geom_point(size = 1) +
  theme_cowplot(font_size = 12) +
  ylab("CFP/YFP") +
  xlab("Time (min)") +
  labs(color = "Protease") +
  facet_wrap(~sub) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(fill="#f0f0f0"))
# show plot
print(p3)

### plot deltaneg CFP/YFP values
# plot
p4 <- ggplot(data_deltaneg %>% filter(sub == "H", pro %in% c("buffer", "EV", "protK", "TEVp I", "TEVp I C151A")),
            aes(x = time, y = deltaneg, color = pro)) +
  geom_point(size = 1) +
  theme_cowplot(font_size = 12) +
  ylab("delta CFP/YFP") +
  xlab("Time (min)") +
  labs(color = "Protease") +
  facet_wrap(~sub) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(fill="#f0f0f0"))
# show plot
print(p4)

### plot activity values relative to protK
# plot
p5 <- ggplot(data_protK %>% filter(sub == "H", pro %in% c("buffer", "EV", "protK", "TEVp I", "TEVp I C151A")),
            aes(x = time, y = activity, color = pro)) +
  geom_point(size = 1) +
  geom_abline(data_vini %>% filter(sub == "H", pro == "TEVp I"),
              mapping = aes(slope = vini,
                            intercept = intercept, color = pro)) +
  theme_cowplot(font_size = 12) +
  ylab("Processed substrate (%)") +
  xlab("Time (min)") +
  labs(color = "Protease") +
  facet_wrap(~sub) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(fill="#f0f0f0"))
# show plot
print(p5)


# put all plots as panels into a grid and save it
# Put all plots into a grid
pg <- plot_grid(p1, p2, p3, p4, p5, labels = "auto")
pg
# save grid
save_plot(pg, filename = "./plots/processing_example_sub.pdf", 
          base_height = 5)

