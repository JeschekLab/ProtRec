# load libraries
library(tidyverse)
library(viridis)
library(cowplot)
library(ggseqlogo)



# create folders ----------------------------------------------------------

# create a folder to store the individual plots in
dir.create("./01_internal_controls")
dir.create("./01_internal_controls/results")
dir.create("./01_internal_controls/plots")
dir.create("./01_internal_controls/plots/negative_controls")
dir.create("./01_internal_controls/plots/heatmaps")
dir.create("./01_internal_controls/plots/compare")
dir.create("./01_internal_controls/plots/reproducibility")


# prepare data ------------------------------------------------------------

## P1' specificity ----

# read data
data <- read.table("./data/Data_DMS_internal_controls.txt", sep = "\t", header = TRUE)

# format data
data <- data %>% 
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data$TEVs)) # write AA of mutated position in separate column)

# clean up data set 
data <- data %>% select(TEVp, Note, TEVs, TEVs_pos, TEVs_AA, t1, t2, t3, t4, t5, sum, min, AUC)

## raw AUCs rel to wt TEVs ----

# set raw AUCs relative to wt TEVs

# get values of wt TEVs
wtTEVs <- data %>% filter(TEVs == "ENLYFQS") %>% select(TEVp, AUC) %>% rename("AUC_ENLYFQS" = "AUC")
# join data with these values
data <- full_join(data, wtTEVs, by = "TEVp")
# calculate normalized blanked values
data <- data %>% mutate(AUCrel = AUC / AUC_ENLYFQS)

# 1_negative controls -----------------------------------------------------

## define noise in the data ----

# get all negative controls
neg_ctrls <- data %>%
  filter(min >= 20) %>%
  filter(TEVp %in% c('no TEVp', '(1-13)-mCherry', 'I C151A 30*31*')) %>%
  select(TEVp, TEVs, AUC) %>%
  arrange(desc(AUC))

# save table
write.table(neg_ctrls, file = "./01_internal_controls/results/neg_ctrls.txt", sep = "\t", row.names=FALSE, quote = F)

# check for normal distribution of data
shapiro.test(neg_ctrls$AUC)
# p–value is higher than 0.05, so we can assume the normality

# plot the three different negative controls to see their distribution 
p <- ggplot(neg_ctrls, aes(x = AUC, fill = TEVp)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.6, color = "black")+
  geom_density(alpha = 0.6) +
  theme_cowplot(font_size = 12) +
  labs(fill = "Negative control")
# show plot
print(p)
# save plot
save_plot(p, filename = "./01_internal_controls/plots/negative_controls/distribution_individual.pdf",
          base_asp = 1.3)

# test if the individual negative controls show a normal distribution
shapiro.test(neg_ctrls %>% filter(TEVp == "no TEVp") %>% pull(AUC))
shapiro.test(neg_ctrls %>% filter(TEVp == "I C151A 30*31*") %>% pull(AUC))
shapiro.test(neg_ctrls %>% filter(TEVp == "(1-13)-mCherry") %>% pull(AUC))
# all of them show a p-value of >0.05, so we can assume normal distribution for all of them

# fit a normal distribution on the data
library(MASS)
fit <- fitdistr(neg_ctrls$AUC, "normal")
fit

# define a cut off as the value under which a sample falls with a 
# probability of 95%
cutoff <- qnorm(0.95, mean = fit$estimate[1], sd = fit$estimate[2])
cutoff

## unload the MASS package to avoid conflicts with tidyverse
detach("package:MASS", unload=TRUE)


# plotting of the fitted normal distribution

mu <- fit$estimate[1]
sigma <- fit$estimate[2]
# generate a data frame with values of the normal distribution created using the fit
range <- seq(mu - 3*sigma, mu + 3*sigma, length.out = 100)
density <- dnorm(range, mean = mu, sd = sigma)
df_nd <- data.frame(x = range, y = density)

# plot
p <- ggplot() + 
  geom_histogram(data = neg_ctrls,
                 aes(x = AUC, y = after_stat(density)), 
                 binwidth = 0.0015,
                 fill = viridis::magma(n = 12)[3], color = 'white', alpha = 0.4) +
  geom_line(data = df_nd, aes(x = x, y = y), color = 'black', size = 1) +
  theme_cowplot(font_size = 12) +
  geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = cutoff + 0.005, y = max(df_nd$y), 
           label = paste0("cutoff =\n", round(cutoff, 4)), 
           color = 'black')

# show plot
print(p)
# save plot
save_plot(p, filename = "./01_internal_controls/plots/negative_controls/distribution_fitted.pdf",
          base_asp = 1)

# collect noise results in a table
df_noise <- data.frame(cutoff = cutoff,
                       mu = mu,
                       sigma = sigma)
# save table of noise results
write.table(df_noise, file = "./01_internal_controls/results/df_noise.txt", sep = "\t", row.names=FALSE, quote = F)



## manipulate data acc. to noise ----

# calculate noise-blanked values (delta noise AUC)
# for that, set all values below the cutoff to the mean of the noise
# then blank to the mean of the noise
data$dnAUC <- ifelse(data$AUC <= cutoff, mu, data$AUC)-mu

### noise blanked values rel to wt TEVs
# get values of noise-blanked wt TEVs
nblanked_ENLYFQS <- data %>% filter(TEVs == "ENLYFQS") %>% select(TEVp, dnAUC) %>% rename("nblanked_ENLYFQS" = "dnAUC")
# join data with these values
data <- full_join(data, nblanked_ENLYFQS, by = c("TEVp")) 
# calculate normalized blanked values
data <- data %>% mutate(dnAUCrel = dnAUC / nblanked_ENLYFQS)

# save final data
write.table(data, file = "./01_internal_controls/results/Data_internal_controls.txt", sep = "\t", row.names=FALSE, quote = F)

# 2_P1' heat maps (manual) -------------------------------------------------------

## to plot heat maps of P1' specificity profiles ----

# manually define the order of TEVp variants
TEVp_levels <- c(
  'no TEVp', '(1-13)-mCherry', "MBP", "I C151A 30*31*", "0 C151A no ssrA",
  '0', "0 C151A", 
  "I", "I C151A",
  'IV', "IV C151A",
  'VI', 'VI C151A',
  'VI*', 'VI* C151A',
  'VII', 'VII C151A',
  'VII*', 'VII* C151A',
  'NX1', 'NX1 C151A',
  'NX2', 'NX2 C151A')

# write function to generate P1' heat maps
heatmap_P1p <- function(
  TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
                      TEVp %in% TEVp_var)
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(plot_data$TEVs_AA))
  # merge with data of actually measured TEVp x TEVs combinations
  plot_data <- plot_data %>%
    full_join(., df_combinations,
              by = c("TEVp", "TEVs_AA")) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # plotting
  p <- ggplot(plot_data, aes(x = TEVs_AA, 
                             y = factor(TEVp, levels = TEVp_levels), 
                             fill = values)) +
    geom_tile() + # make heat map
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 12) +
    theme(axis.line = element_blank()) +
    theme(axis.title.y = element_blank()) +
    labs(fill = legend) +
    scale_x_discrete(name = "Amino acid in P1'", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) #+
  #    theme(axis.text.x = element_text(size = 10)) + # change size of x axis labels
  # facet_grid(category ~ ., scales = "free", space = "free") # switch rows and columns
  
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./01_internal_controls/plots/heatmaps/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}

## AUC: raw AUCs ----

# heat map of all variants showing raw AUC values
heatmap_P1p(TEVp_var = pull(data, TEVp),
            activity = "AUC",
            legend = "AUC",
            file_name = "all_variants_AUC",
            height = 5, 
            asp_ratio = 1)

## AUC: noise-corrected AUCs ----

# heat map of all variants showing raw AUC values
heatmap_P1p(TEVp_var = pull(data, TEVp),
            activity = "dnAUC",
            legend = "dnAUC",
            file_name = "all_variants_dnAUC",
            height = 5, 
            asp_ratio = 1)


# 3_comparison to in vitro data ---------------------------------------------

## get data ----

# define working directory
DIR <- gsub("uASPIre screen4 and 5_2022_08_11/uASPIre screen4_DMS", "", getwd())

# read in vitro data
iv_data_abs <- read.table(paste0(DIR, "In vitro data by Luzius Pestalozzi/data_viniabs_Luzi.txt"),
                          sep = "\t", header = TRUE, quote = "") %>%
  rename("TEVs_pos" = "pos")

iv_data_rel <- read.table(paste0(DIR, "In vitro data by Luzius Pestalozzi/data_vinirel_Luzi.txt"),
                          sep = "\t", header = TRUE, quote = "") %>%
  rename("TEVs_pos" = "pos")

# add in vitro data to the other data
data_compare <- left_join(filter(data, 
                                 TEVp %in% c('0','VI','VI*','VII', 'NX1', 'NX2'), 
                                 TEVs_AA != "*") %>%
                            select(TEVp, TEVs, TEVs_pos, TEVs_AA, AUC, AUCrel, dnAUC, dnAUCrel),
                          iv_data_abs, 
                          by = c("TEVs_pos", "TEVs_AA", "TEVp")) %>% 
  left_join(., iv_data_rel, by = c("TEVs_pos", "TEVs_AA", "TEVp"))


## calculate correlations ----

## first for absolute values
# initialize empty data frame
df_cor_abs <- data.frame()
# loop 
for (i in c("AUC", "dnAUC")) {
  for (ii in c("pearson", "spearman")) {
    # make temporary data frame
    df_temp <- data.frame(
      x = i,
      y = "viniabs",
      method = ii,
      cor = cor(data_compare[, i], data_compare[, "viniabs"], use = "complete.obs", method = ii) %>% round(., 5))
    # combine
    df_cor_abs <- rbind(df_cor_abs, df_temp)
  }
}

## second for values relative to wt TEVs
# initialize empty data frame
df_cor_rel <- data.frame()
# loop 
for (i in c("AUCrel", "dnAUCrel")) {
  for (ii in c("pearson", "spearman")) {
    # make temporary data frame
    df_temp <- data.frame(
      x = i,
      y = "vinirel",
      method = ii,
      cor = cor(data_compare[, i], data_compare[, "vinirel"], use = "complete.obs", method = ii) %>% round(., 5))
    # combine
    df_cor_rel <- rbind(df_cor_rel, df_temp)
  }
}

# save data table
write.table(df_cor_abs %>% pivot_wider(., names_from = method, values_from = cor),
            file = "./01_internal_controls/results/correlations_abs.txt", sep = "\t", row.names=FALSE, quote = F)
write.table(df_cor_rel %>% pivot_wider(., names_from = method, values_from = cor),
            file = "./01_internal_controls/results/correlations_rel.txt", sep = "\t", row.names=FALSE, quote = F)



## scatter plot absolute values ----

# use a loop to do the same for all different types of relative values

# initialize an empty list to store the plots in
plot_list <- list()

for (i in c("AUC", "dnAUC")) {
  
  # plot
  p <- ggplot(data_compare, aes_string(x = "viniabs", y = i)) +
    geom_point() +
    theme_cowplot(font_size = 12) +
    labs(x = "viniabs", y = i) +
    # add the correlation as text label in the plot
    annotate("text", x = -Inf, y = -Inf,
             label = paste("Sp. cor. =",
                           df_cor_abs %>% filter(x == i, method == 'spearman') %>% pull(cor)),
             hjust   = -1.5,
             vjust   = -2)

    # show plot
  print(p)
  # save plot in a list
  plot_list[[i]] <- p
}

## scatter plot relative values ----

# use a loop to do the same for all different types of relative values
for (i in c("AUCrel", "dnAUCrel")) {
  # plot
  p <- ggplot(data_compare, aes_string(x = "vinirel", y = i)) +
    geom_point() +
    theme_cowplot(font_size = 12) +
    labs(x = "vinirel", y = i) +
    # add the correlation as text label in the plot
    annotate("text", x = -Inf, y = -Inf,
             label = paste("Sp. cor. =",
                           df_cor_rel %>% filter(x == i, method == 'spearman') %>% pull(cor)),
             hjust   = -1.5,
             vjust   = -2)
  # show plot
  print(p)
  # save plot in a list
  plot_list[[i]] <- p
}


# Put all plots into a grid
pg <- plot_grid(plotlist = plot_list, ncol = 2, byrow = FALSE)
pg

# save grid
save_plot(pg, filename = "./01_internal_controls/plots/compare/correlations_grid.pdf", 
          base_height = 7.87,
          base_asp = 1) 

## plotting in vivo and in vitro data together -----

### plot all of Luzi's TEVps ----

# prepare data for plotting
data_plot <- data_compare %>% 
  select(TEVp, TEVs_pos, TEVs_AA, vinirel, all_of(i)) %>%
  rename("in vivo" = all_of(i), "in vitro" = "vinirel") %>%
  pivot_longer(., cols = c("in vivo", "in vitro"), names_to = "assay", values_to = "activity")

# create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations
df_combinations <- expand.grid(
  TEVp = unique(data_plot$TEVp),
  TEVs_AA = unique(data_plot$TEVs_AA),
  assay = c("in vivo", "in vitro"))

# merge with data of actually measured TEVp x TEVs combinations
data_plot <- data_plot %>%
  full_join(., df_combinations,
            by = c("TEVp", "TEVs_AA", "assay")) %>%
  mutate(assay = factor(assay, levels = c("in vivo", "in vitro")))

# plot
p <- ggplot(data_plot,
            aes(x = TEVs_AA,
                y = factor(TEVp, levels = TEVp_levels), 
                fill = activity)) +
  geom_tile() +
  facet_wrap(~assay, ncol = 1) +
  scale_fill_viridis(option="magma") +
  theme_cowplot(font_size = 12) +
  theme(axis.line = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(fill = "Relative \nactivity") +
  scale_x_discrete(name = "Amino acid in P1'", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

# show plot
print(p)

# save plot
save_plot(p, filename = paste0("./01_internal_controls/plots/heatmaps/in_vivo_vs_in_vitro_dnAUCrel.pdf"),
          base_height = 3.6, base_asp = 1.2)



# 5_reproducibility -------------------------------------------------------------

# compare results of this experiment with previous experiments

# define path for data from other experiments
DIR <- gsub("uASPIre screen4 and 5_2022_08_11/uASPIre screen4_DMS", "", getwd())

## DMS vs. ABC screen ----

# get data from ABC screen

data_ABC <- read.table(paste0(DIR, "/uASPIre screen3_2021-07-20_ABC/data/02_data_AUC_all.txt"),
                       sep = "\t", header = TRUE, quote = "")

# remove duplicate entries that stem from a mistake in previous data processing by Simon
data_ABC <- data_ABC %>% filter( !(TEVs == "ENLYFQS" & library.ID == "libLH016"),
                                 !(TEVs == "ENLYFQS" & library.ID == "libLH019"))

# rename positive control
data_ABC$TEVp[data_ABC$Note == "positive control"] <- "0 C151A no ssrA"

# filter and clean up data set 
data_ABC <- data_ABC %>% 
  # filter for internal controls
  filter(library.ID %in% c('libLH015', 'libLH016', 'libLH019')) %>%
  # remove MBP as it was underrepresented and is not a good negative control
  filter(TEVp != "MBP") %>%
  # clean up table
  select(TEVp, Note, TEVs, rep, t1, t2, t3, t4, t5, sum, min, AUC4) %>%
  rename("AUC" = "AUC4")

# filter further for correlation plot
corr_data_ABC <- data_ABC %>%
  filter(min >= 20) %>%
  filter(Note != "neg ctrl catalysis") %>% # kick out C151A variants
  select(TEVp, TEVs, AUC, rep) %>%
  mutate(rep = paste0("ABC_", rep))

## do the same with the DMS data

# filter and clean up data set 
corr_data_DMS <- data %>%
  # set read threshold
  filter(min >= 20) %>%
  # remove MBP as it was underrepresented and is not a good negative control
  filter(TEVp != "MBP") %>%   
  # kick out C151A variants
  filter(Note != "neg ctrl catalysis") %>% 
  select(TEVp, TEVs, AUC) %>%
  mutate(rep = "DMS")

# correct entry for TEVs of positive control
corr_data_DMS$TEVs[corr_data_DMS$TEVp == "0 C151A no ssrA"] <- "ENLYFQ-"

# combine data of both screens and bring to wider format for plotting
replicates <- rbind(corr_data_ABC, corr_data_DMS) %>%
  pivot_wider(., names_from = rep, values_from = AUC)


### ABC screen rep 1 ----

# create linear model
fit <- lm(ABC_1 ~ DMS, data = replicates)

# plot
p <-  ggplot(replicates, aes(x = DMS, y = ABC_1)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4),
                          ", r = ",
                          round(cor(replicates$DMS, replicates$ABC_1, use = "pairwise.complete.obs", method = 'pearson'), 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("AUC DMS screen") +
  ylab("AUC ABC screen rep 1") +
  scale_x_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_internal_controls/plots/reproducibility/DMS_vs_ABCrep1.pdf",
          base_asp = 1)

### ABC screen rep 2 ----

# create linear model
fit <- lm(ABC_2 ~ DMS, data = replicates)

# plot
p <-  ggplot(replicates, aes(x = DMS, y = ABC_2)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4),
                          ", r = ",
                          round(cor(replicates$DMS, replicates$ABC_2, use = "pairwise.complete.obs", method = 'pearson'), 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("AUC DMS screen") +
  ylab("AUC ABC screen rep 2") +
  scale_x_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_internal_controls/plots/reproducibility/DMS_vs_ABCrep2.pdf",
          base_asp = 1)


## DMS vs. pop2 screen ----

# get data from pop2 screen
data_pop2 <- read.table(paste0(DIR, "uASPIre screen2_2021-03-30_pop2/01_pop2/results/data_pop2clean.txt"),
                        sep = "\t", header = TRUE, quote = "")

# filter further for correlation plot
corr_data_pop2 <- data_pop2 %>%
  filter(Note != "neg ctrl catalysis") %>% # kick out C151A variants
  rename("AUC" = "AUC4") %>%
  select(TEVp, TEVs, AUC, rep) %>%
  mutate(rep = paste0("pop2_", rep))

# combine data of both screens and bring to wider format for plotting
replicates <- rbind(corr_data_ABC, corr_data_DMS, corr_data_pop2) %>%
  pivot_wider(., names_from = rep, values_from = AUC)


### pop2 screen rep 1 ----

# create linear model
fit <- lm(pop2_1 ~ DMS, data = replicates)

# plot
p <-  ggplot(replicates, aes(x = DMS, y = pop2_1)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4),
                          ", r = ",
                          round(cor(replicates$DMS, replicates$pop2_1, use = "pairwise.complete.obs", method = 'pearson'), 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("AUC DMS screen") +
  ylab("AUC pop2 screen rep 1") +
  scale_x_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_internal_controls/plots/reproducibility/DMS_vs_pop2rep1.pdf",
          base_asp = 1)

### pop2 screen rep 2 ----

# create linear model
fit <- lm(pop2_2 ~ DMS, data = replicates)

# plot
p <-  ggplot(replicates, aes(x = DMS, y = pop2_2)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4),
                          ", r = ",
                          round(cor(replicates$DMS, replicates$pop2_2, use = "pairwise.complete.obs", method = 'pearson'), 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("AUC DMS screen") +
  ylab("AUC pop2 screen rep 2") +
  scale_x_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_internal_controls/plots/reproducibility/DMS_vs_pop2rep2.pdf",
          base_height = 3,
          base_asp = 1.3)

