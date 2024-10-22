# load libraries
library(tidyverse)
library(viridis)
library(cowplot)

# create folders --------------------------------------------------------

# create a folder to store the individual plots and results in
dir.create("./01_dynamic_range")
dir.create("./01_dynamic_range/results")
dir.create("./01_dynamic_range/plots")

# prepare data ------------------------------------------------------------

# read data
data <- read.table("./data/LH_TEV_fraction.txt", sep = "\t", header = TRUE) 

# clean up data frame
data <- data %>%
  select(-BC, -Plasmid, -sum, -min, -TEV)

# rename some variants to have consistent naming
data <- data %>%
  mutate(TEVp = ifelse(TEVp == "C151A", "0 C151A", TEVp),
         TEVp = ifelse(TEVp == "VIstar", "VI*", TEVp),
         TEVp = ifelse(TEVp == "p0", "0", TEVp),
         TEVp = ifelse(TEVs == "ENLYFQ–", "0 C151A no ssrA", TEVp))

# add column stating the respective promoter strength
data$p_strength[data$promoter == "J23116"] <- "0.16"
data$p_strength[data$promoter == "J23107"] <- "0.36"
data$p_strength[data$promoter == "J23118"] <- "0.56"
data$p_strength[data$promoter == "J23100"] <- "1.00"

# rename the lib column entries according to induction time point
data$lib[data$lib == "1"] <- "2h"
data$lib[data$lib == "2"] <- "4h"
data$lib[data$lib == "3"] <- "6h"


# add columns needed for some plots
data <- data %>%
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data$TEVs)) # write AA of mutated position in separate column)

# save data of linker study
write.table(data, file = "./01_dynamic_range/results/data_clean.txt", sep = "\t", row.names=FALSE, quote = F)

## prepare data for plotting ----

# add  a column stating the type of the line (pos ctrl should be dashed)


# prepare data
all_flasks <- data %>%
  # selecting only necessary columns
  select(-AUC, -TEVs, -TEVs_pos, -Combi, -promoter) %>%
  # add column specifying the line type for each TEVp
  mutate(type = factor(if_else(TEVp == "0 C151A no ssrA", "without ssrA", "with ssrA"), 
                       levels = c("with ssrA", "without ssrA")))

# rename time points by actual hours after induction
# this is different depending on column lib so we need to separate the libraries

flask1 <- all_flasks %>%
  filter(lib == "2h") %>%
  rename("0"=X1, "1"=X2, "2"=X3, "3"=X4, "4"=X5, "5"=X6, "6"=X7, "7"=X8, "8"=X9, "9"=X10, "10"=X11, "14"=X12)

flask2 <- all_flasks %>%
  filter(lib == "4h") %>%
  rename("0"=X1, "1"=X2, "2"=X3, "3"=X4, "4"=X5, "5"=X6, "6"=X7, "7"=X8, "8"=X9, "9"=X10, "10"=X11, "12"=X12)

flask3 <- all_flasks %>%
  filter(lib == "6h") %>%
  rename("-4"=X1, "0"=X2, "1"=X3, "2"=X4, "3"=X5, "4"=X6, "5"=X7, "6"=X8, "7"=X9, "8"=X10, "9"=X11, "10"=X12)

# bring the data from the three flasks together again in long format
all_flasks <- rbind(flask1 %>% 
                      pivot_longer(cols = !c(lib, TEVp, TEVs_AA, RBS, p_strength, type),
                                   names_to = "timepoint", 
                                   values_to = "fraction_flipped"),
                    flask2 %>% 
                      pivot_longer(cols = !c(lib, TEVp, TEVs_AA, RBS, p_strength, type),
                                   names_to = "timepoint", 
                                   values_to = "fraction_flipped"),
                    flask3 %>% 
                      pivot_longer(cols = !c(lib, TEVp, TEVs_AA, RBS, p_strength, type),
                                   names_to = "timepoint", 
                                   values_to = "fraction_flipped"))

# flipping plots ---------------------------------------------------------
  
  
## faceted by inoculation tp and promoter ----  

# plotting
p <- ggplot(all_flasks %>% filter(TEVp != "NX1mut"), 
            aes(x = as.integer(timepoint), 
                y = fraction_flipped, 
                col = RBS, 
                linetype = type, 
                group = paste(TEVp, TEVs_AA, RBS))) +
  # geom_point() +
  geom_line() +
  facet_grid(lib ~ p_strength, 
             labeller = labeller(
               # lib = c("2h" = "Induction 2h after inoculation",
               #         "4h" = "Induction 4h after inoculation",
               #         "6h" = "Induction 6h after inoculation"),
               p_strength = c("0.16" = "J23116 (0.16)",
                              "0.36" = "J23107 (0.36)",
                              "0.56" = "J23118 (0.56)",
                              "1.00" = "J23100 (1.00)")
               )
             ) +
  theme_cowplot(font_size = 12) +
  scale_x_continuous("Time after induction (h)", limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  scale_y_continuous("Fraction flipped", limits = c(0, 1)) +
  scale_color_manual(values = c("#005390", "#931C51"))

# show plot
print(p)

# and save to file
save_plot(p, filename = paste0("./01_dynamic_range/plots/flipping_induction_promoter.pdf"),
          base_height = 5,
          base_asp = 1.6)

# plotting TEVp 0 and TEVp 0 C151A only
p <- ggplot(all_flasks %>% filter(TEVp %in% c("0", "0 C151A", "0 C151A no ssrA")), 
            aes(x = as.integer(timepoint), 
                y = fraction_flipped, 
                col = RBS, 
                linetype = TEVp, 
                group = paste(TEVp, TEVs_AA, RBS))) +
  # geom_point() +
  geom_line() +
  facet_grid(lib ~ p_strength, 
             labeller = labeller(
               # lib = c("2h" = "Induction 2h after inoculation",
               #         "4h" = "Induction 4h after inoculation",
               #         "6h" = "Induction 6h after inoculation"),
               p_strength = c("0.16" = "J23116 (0.16)",
                              "0.36" = "J23107 (0.36)",
                              "0.56" = "J23118 (0.56)",
                              "1.00" = "J23100 (1.00)")
             )
            ) +
  theme_cowplot(font_size = 12) +
  scale_x_continuous("Time after induction (h)", limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  scale_y_continuous("Fraction flipped", limits = c(0, 1)) +
  scale_color_manual(values = c("#005390", "#931C51")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"))


# show plot
print(p)

# and save to file
save_plot(p, filename = paste0("./01_dynamic_range/plots/flipping_induction_promoter_selected.pdf"),
          base_height = 5,
          base_asp = 1.6)



## color AUC for illustration purposes  ----  

test <- all_flasks %>% 
  filter(TEVp %in% c("0", "0 C151A"),
         lib == "2h",
         RBS == "veryweak",
         p_strength == "1.00")

# plotting TEVp 0 and TEVp 0 C151A only
p <- ggplot(all_flasks %>% 
              filter(TEVp %in% c("0", "0 C151A"),
                     lib == "2h",
                     RBS == "veryweak",
                     p_strength == "1.00"), 
            aes(x = as.integer(timepoint), 
                y = fraction_flipped, 
                # col = RBS, 
                linetype = TEVp, 
                group = paste(TEVp, TEVs_AA))) +
  geom_point() +
  geom_line() +
  geom_ribbon(data = all_flasks %>% 
                filter(TEVp == "0",
                       lib == "2h",
                       RBS == "veryweak",
                       p_strength == "1.00"), 
              aes(ymin = 0, ymax = fraction_flipped), 
              fill = "red", alpha = 0.5) +
  geom_ribbon(data = all_flasks %>% 
                filter(TEVp == "0 C151A",
                       lib == "2h",
                       RBS == "veryweak",
                       p_strength == "1.00"), 
              aes(ymin = 0, ymax = fraction_flipped), 
              fill = "blue", alpha = 0.5) +
  theme_cowplot(font_size = 12) +
  scale_x_continuous("Time after induction (h)", limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  scale_y_continuous("Fraction flipped", limits = c(0, 0.8)) +
  scale_color_manual(values = c("#005390", "#931C51")) +
  scale_linetype_manual(values = c("solid", "dashed")) 


# show plot
print(p)

# and save to file
save_plot(p, filename = paste0("./01_dynamic_range/plots/AUC_illustration.pdf"),
          base_height = 2.5,
          base_asp = 1.5)


# Induction time point analysis: 0 vs 0 C151A ----------------------------------

# Which induction time point is best to differentiate active from inactive variants?

# calculate the fold change of TEVp 0 over TEVp 0 C151A
fc <- data %>%
  filter(TEVp %in% c("0", "0 C151A"),
         RBS == "veryweak") %>%
  select(TEVp, AUC, p_strength, lib) %>%
  pivot_wider(names_from = TEVp, values_from = AUC) %>%
  rename("TEVp0" = "0", "TEVp0C151A" = "0 C151A") %>%
  mutate(fc = TEVp0 / TEVp0C151A)

promoter_col <- c("0.16" = "#ffffcc", "0.36" = "#a1dab4", "0.56" = "#41b6c4", "1.00" = "#225ea8")

# make a bar plot
p <- ggplot(fc, aes(x = lib, y = fc, fill = p_strength)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = promoter_col) +
  theme_cowplot(font_size = 12) +
  labs(fill = "Relative \npromoter \nstrength") +
  ylab("TEVp 0 / TEVp 0 C151A (fc AUCfinal)") +
  xlab("Induction after inoculation") +
  coord_cartesian(ylim = c(0.75, NA))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_dynamic_range/plots/fc.pdf",
          base_height = 2.5,
          base_asp = 1.5)


## Compare signal, dynamic range, and DR above signal ----

data_c <- flask1 %>% # only take 2h induction time point
  # only consider very weak RBS
  filter(RBS == "veryweak") %>%
  select(-RBS, -lib) %>%
  # fraction flipped 4h after induction
  select(-as.character(c(0:3, 5:10, 14)),
         -TEVs_AA, -type) %>%
  # rename column
  rename(ff = "4") %>%
  filter(TEVp %in% c("0", "0 C151A", "0 C151A no ssrA"))


# Calculate s2b and dr
# s2b = signal to background ratio
# dr = dynamic range above signal

# write a function
calculate_values <- function(df) {
  # Subsetting data for each TEVp category
  tevp_0 <- df[df$TEVp == "0",]
  tevp_C151A <- df[df$TEVp == "0 C151A",]
  tevp_no_ssrA <- df[df$TEVp == "0 C151A no ssrA",]
  
  # Merging datasets to align by p_strength
  merge1 <- merge(tevp_0, tevp_C151A, by="p_strength", suffixes = c("_0", "_C151A"))
  merge2 <- merge(tevp_no_ssrA, tevp_C151A, by="p_strength", suffixes = c("_no_ssrA", "_C151A"))
  
  # Calculating s2b and drMJ
  merge1$s2b <- merge1$ff_0 / merge1$ff_C151A
  merge2$s2b_no_ssrA <- merge2$ff_no_ssrA / merge2$ff_C151A
  merge2$s2b <- merge1$s2b[match(merge2$p_strength, merge1$p_strength)]  # Matching the s2b values from the first merge
  merge2$drMJ <- merge2$s2b_no_ssrA - merge2$s2b
    merge2$drLH <- merge2$s2b_no_ssrA / merge2$s2b
  
  # Result
  result <- merge2[,c("p_strength", "s2b", "s2b_no_ssrA", "drMJ", "drLH")]
  return(result)
}

# Call the function and store the results
results <- calculate_values(data_c)


# prepare plot
results_long <- pivot_longer(results, cols = c(s2b, s2b_no_ssrA, drLH), names_to = "Metric", values_to = "Value")

# plot
p <- ggplot(results_long, aes(x = factor(p_strength), y = Value, fill = factor(Metric, levels = c("s2b", "s2b_no_ssrA", "drLH")))) +
  geom_bar(stat = 'identity', position = 'dodge', color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Adding the dashed line at y = 1
  scale_fill_manual(values = c("s2b" = "#a1dab4", "s2b_no_ssrA" = "#41b6c4", "drLH" = "#225ea8"),
                    labels = c("s2b" = "signal ", "s2b_no_ssrA" = "dynamic range (DR)", "drLH" = "DR above signal")) +
  labs(x = "Relative promoter strength", y = "Ratio", fill = NULL) +
  theme_cowplot(font_size = 12)

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_dynamic_range/plots/dynamic_range_comparison.pdf",
          base_height = 3,
          base_asp = 1.5)



# AUC analysis ------------------------------------------------------------

# Which AUC is best to to differentiate active from inactive variants?

# get data from all flasks time points and common time points (0 to 10h)
all_flasks_AUCs <- all_flasks %>% pivot_wider(., 
                                   id_cols = c("TEVp", "RBS", "p_strength", "TEVs_AA", "type", "lib"),
                                   names_from = timepoint,
                                   values_from = fraction_flipped) %>%
  select(-c("-4", "12", "14"))

# make data frame with fraction flipped only
df_fraction <- data.frame(select(all_flasks_AUCs, as.character(0:10)))

### calculate AUC
x <- as.integer(0:10) # get time points

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_fraction), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_fraction[,i]) + (x[i+1]-x[i])*(df_fraction[,i+1]-df_fraction[,i])*0.5
  m_auc[,i] <- auc
}

# # calculate AUCs and add columns to data
all_flasks_AUCs$AUC1 <- round(m_auc[,1:1], 3)
all_flasks_AUCs$AUC2 <- round(rowSums(m_auc[,1:2])/2, 3)
all_flasks_AUCs$AUC3 <- round(rowSums(m_auc[,1:3])/3, 3)
all_flasks_AUCs$AUC4 <- round(rowSums(m_auc[,1:4])/4, 3)
all_flasks_AUCs$AUC5 <- round(rowSums(m_auc[,1:5])/5, 3)
all_flasks_AUCs$AUC6 <- round(rowSums(m_auc[,1:6])/6, 3)
all_flasks_AUCs$AUC7 <- round(rowSums(m_auc[,1:7])/7, 3)
all_flasks_AUCs$AUC8 <- round(rowSums(m_auc[,1:8])/8, 3)
all_flasks_AUCs$AUC9 <- round(rowSums(m_auc[,1:9])/9, 3)
all_flasks_AUCs$AUC10 <- round(rowSums(m_auc[,1:10])/10, 3)


# calculate fold change of TEVp 0 over TEVp 0 C151A for different AUCs for the flask
# induced 2h after inoculation

# get AUCs of TEVp 0
TEVp0 <- all_flasks_AUCs %>%
  filter(TEVp == "0",
         RBS == "veryweak", 
         lib == "2h") %>%
  select(p_strength, paste0("AUC", 1:10)) %>% 
  mutate(p_strength = as.numeric(p_strength)) %>%
  arrange(p_strength)

# get AUCs of TEVp 0 C151A in the same order (i.e., arranged by p_strength)
TEVp0C151A <- all_flasks_AUCs %>%
  filter(TEVp == "0 C151A",
         RBS == "veryweak", 
         lib == "2h") %>%
  select(p_strength, paste0("AUC", 1:10)) %>% 
  mutate(p_strength = as.numeric(p_strength)) %>%
  arrange(p_strength)

# calculate the fold changes. p_strength serves as control and must become 1
fcAUCs <- TEVp0 / TEVp0C151A 

# format data. Add p_strength values
promoters_rel <- c("0.16", "0.36", "0.56", "1.00")
fcAUCs$p_strength <- promoters_rel

# bring data to longer format for plotting
fcAUCs <- fcAUCs %>%
  pivot_longer(cols = paste0("AUC", 1:10), names_to = "AUCid", values_to = "AUC")

# rename AUCid column, remove "AUC" from name
fcAUCs$AUCid <- gsub("AUC", "", fcAUCs$AUCid) 

# Convert AUCid to factor with desired levels
fcAUCs$AUCid <- factor(fcAUCs$AUCid, levels = 1:10)

# make a bar plot
p <- ggplot(fcAUCs %>% filter(!AUCid %in% c(1:3)), aes(x = AUCid, y = AUC, fill = p_strength)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = promoter_col) +
  theme_cowplot(font_size = 12) +
  labs(fill = "Relative \npromoter \nstrength") +
  ylab("Fold change of AUC") +
  xlab("Time after induction (h)")

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_dynamic_range/plots/AUC_fc.pdf",
          base_height = 3,
          base_asp = 2)


## Promoter 0.36 only ----


# make a bar plot of only promoter 0.36 and AUCs 4-10
p <- ggplot(fcAUCs %>% filter(p_strength == "0.36", !AUCid %in% c(1:3)), aes(x = AUCid, y = AUC)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = promoter_col) +
  theme_cowplot(font_size = 12) +
  labs(fill = "Relative \npromoter \nstrength") +
  ylab("Fold-change AUC4h over background") +
  xlab("Time after induction (h)") +
  coord_cartesian(ylim = c(1, NA))


# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_dynamic_range/plots/AUC_fc_NEW.pdf",
          base_height = 3.6,
          base_asp = 1.2)








































