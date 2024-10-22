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
dir.create("./01_internal_controls/plots/heatmaps")
dir.create("./01_internal_controls/plots/reproducibility")
dir.create("./01_internal_controls/plots/negative_controls")

# prepare data ------------------------------------------------------------


# read data
data <- read.table("./data/02_data_AUC_all.txt", sep = "\t", header = TRUE)

# remove duplicate entries that stem from a mistake in previous data processing 
data <- data %>% filter( !(TEVs == "ENLYFQS" & library.ID == "libLH016"),
                         !(TEVs == "ENLYFQS" & library.ID == "libLH019"))

# format data
data <- data %>% 
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data$TEVs)) # write AA of mutated position in separate column)

# rename positive control
data$TEVp[data$Note == "positive control"] <- "0 C151A no ssrA"

# get data of controls with SphI site next to TEVs
withSphI <- data %>%
  # filter for internal controls containing SphI
  filter(SphI.downstream.of.TEVs == 'yes',
         min >= 20) %>%
  # select only important columns
  select(rep, TEVp, Note, TEVs_AA, AUC4) %>%
  rename(AUC = AUC4)

# filter and clean up data set 
data <- data %>% 
  # filter for internal controls
  filter(library.ID %in% c('libLH015', 'libLH016', 'libLH019')) %>%
  # remove MBP as it was underrepresented and is not a good negative control
  filter(TEVp != "MBP") %>%
  # clean up table
  select(TEVp, Note, TEVs, rep, TEVs_pos, TEVs_AA, t1, t2, t3, t4, t5, sum, min, AUC4) %>%
  rename("AUC" = "AUC4")

# save table of controls
write.table(data, file = "./01_internal_controls/results/internal_controls.txt", sep = "\t", row.names=FALSE, quote = F)



# Effect of SphI ----------------------------------------------------------


checkSphI <- inner_join(withSphI, 
                        data %>% select(names(withSphI)),
                       by = c("TEVp", "TEVs_AA", "Note", "rep"),
                       suffix = c("_withSphI", "_withoutSphI"))

# create linear model
fit <- lm(AUC_withoutSphI ~ AUC_withSphI, data = checkSphI)

# plot
p <-  ggplot(checkSphI, aes(x = AUC_withSphI, y = AUC_withoutSphI)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4),
                          ", r = ",
                          round(cor(checkSphI$AUC_withSphI, checkSphI$AUC_withoutSphI, use = "pairwise.complete.obs", method = 'pearson'), 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("AUC with SphI") +
  ylab("AUC without SphI") +
  scale_x_continuous(limits = c(0, 1.1*max(c(checkSphI$AUC_withoutSphI, checkSphI$AUC_withSphI))), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.1*max(c(checkSphI$AUC_withoutSphI, checkSphI$AUC_withSphI))), expand = c(0, 0))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_internal_controls/plots/WithSphI_vs_withoutSphI_AUC.pdf",
          base_height = 3,
          base_asp = 1.4)

# 2_reproducibility --------------------------------------------------------------

## among replicates ----

# create a correlation plot

# format data
replicates <- data %>%
  filter(min >= 20) %>%
  filter(Note != "neg ctrl catalysis") %>% # kick out C151A variants
  select(TEVp, TEVs, AUC, rep) %>%
  pivot_wider(., names_from = rep, values_from = AUC, names_prefix = "rep")

# create linear model
fit <- lm(rep2 ~ rep1, data = replicates)

# plot
p <-  ggplot(replicates, aes(x = rep1, y = rep2)) +
  geom_point(alpha = 0.5, aes(color = TEVp)) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("AUC replicate 1") +
  ylab("AUC replicate 2")
# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_internal_controls/plots/reproducibility/replicates.pdf",
          base_asp = 1.618)

# 3_negative controls -----------------------------------------------------

## define noise in the data ----

# get all negative controls
neg_ctrls <- data %>%
  filter(min >= 20) %>%
  filter(TEVp %in% c('no TEVp', '(1-13)-mCherry', 'I C151A 30*31*')) %>%
  select(TEVp, TEVs, rep, AUC) %>%
  arrange(desc(AUC))

# save table
write.table(neg_ctrls, file = "./01_internal_controls/results/neg_ctrls.txt", sep = "\t", row.names=FALSE, quote = F)

# check for normal distribution of data
shapiro.test(neg_ctrls$AUC)
# p–value is lower than 0.05, so we cannot assume the normality

# plot the three different negative controls to see which data point might cause 
# the low p-value in the Shapiro test
p <- ggplot(neg_ctrls, aes(x = AUC, fill = TEVp)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.6, color = "black")+
  geom_density(alpha = 0.6) +
  theme_cowplot() +
  labs(fill = "Negative control")
# show plot
print(p)
# save plot
save_plot(p, filename = "./01_internal_controls/plots/negative_controls/distribution_individual.pdf")

# test if the individual negative controls show a normal distribution
shapiro.test(neg_ctrls %>% filter(TEVp == "no TEVp") %>% pull(AUC))
shapiro.test(neg_ctrls %>% filter(TEVp == "I C151A 30*31*") %>% pull(AUC))
shapiro.test(neg_ctrls %>% filter(TEVp == "(1-13)-mCherry") %>% pull(AUC))

# test again after removal of the two data points with highest AUC in the mCherry control
shapiro.test(neg_ctrls %>% filter(TEVp == "(1-13)-mCherry", AUC <= 0.09) %>% pull(AUC))

# test for normal distribution across all three negative controls
shapiro.test(neg_ctrls %>% filter(AUC <= 0.09) %>% pull(AUC))
# normal distribution! 

# fit a normal distribution on the data
library(MASS)
fit <- fitdistr(neg_ctrls$AUC, "normal")
fit

# define a cut off as the value under which a samples falls with a 
# probability of 95%
cutoff <- qnorm(0.95, mean = fit$estimate[1], sd = fit$estimate[2])
cutoff

# # repeat after removing the two outliers
fit2 <- fitdistr(neg_ctrls %>% filter(AUC <= 0.09) %>% pull(AUC), "normal")
fit2
cutoff2 <- qnorm(0.95, mean = fit2$estimate[1], sd = fit2$estimate[2])
cutoff2

## unload the MASS package to avoid conflicts with tidyverse
detach("package:MASS", unload=TRUE)


# plotting of the fitted normal distribution before removing outliers

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
  theme_cowplot() +
  geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = cutoff + 0.01, y = max(df_nd$y), 
           label = paste0("cutoff =\n", round(cutoff, 4)), 
           color = 'black')

# show plot
print(p)
# save plot
save_plot(p, filename = "./01_internal_controls/plots/negative_controls/distribution_fitted.pdf",
          base_asp = 1)

# plotting of the fitted normal distribution after removing outliers

mu2 <- fit2$estimate[1]
sigma2 <- fit2$estimate[2]
# generate a data frame with values of the normal distribution created using the fit
range <- seq(mu2 - 3*sigma2, mu2 + 3*sigma2, length.out = 100)
density <- dnorm(range, mean = mu2, sd = sigma2)
df_nd <- data.frame(x = range, y = density)



# plot
p <- ggplot() + 
  geom_histogram(data = filter(neg_ctrls, AUC <= 0.09), 
                 aes(x = AUC, y = after_stat(density)), 
                 binwidth = 0.0015,
                 fill = viridis::magma(n = 12)[3], color = 'white', alpha = 0.4) +
  geom_line(data = df_nd, aes(x = x, y = y), color = 'black', size = 1) +
  theme_cowplot() +
  geom_vline(xintercept = cutoff2, linetype = "dashed", color = 'black', size = 1) +
  annotate("text", x = cutoff2 + 0.01, y = max(df_nd$y), 
           label = paste0("cutoff =\n", round(cutoff2, 4)), 
           color = 'black')

# show plot
print(p)
# save plot
save_plot(p, filename = "./01_internal_controls/plots/negative_controls/distribution_fitted_clean.pdf",
          base_asp = 1)


# collect noise results in a table
df_noise <- data.frame(outlier_removed = c(FALSE, TRUE),
                       cutoff = c(cutoff, cutoff2),
                       mu = c(mu, mu2),
                       sigma = c(sigma, sigma2))
# save table of noise results
write.table(df_noise, file = "./01_internal_controls/results/df_noise.txt", sep = "\t", row.names=FALSE, quote = F)



## manipulate data acc. to noise ----

# calculate noise-blanked values (delta noise AUC)
# to that end, set all values below the cutoff to the mean of the noise
# then blank to the mean of the noise
# here, we use the the cutoff calculated before removing the outlieres, to be 
# consistent with data analyses of future experiments
data$dnAUC <- ifelse(data$AUC <= cutoff, mu, data$AUC)-mu

### noise blanked values rel to wt TEVs
# get values of noise-blanked wt TEVs
nblanked_ENLYFQS <- data %>% filter(TEVs == "ENLYFQS") %>% select(rep, TEVp, dnAUC) %>% rename("nblanked_ENLYFQS" = "dnAUC")
# join data with these values
data <- full_join(data, nblanked_ENLYFQS, by = c("TEVp", "rep")) 
# calculate normalized blanked values
data <- data %>% mutate(dnAUCrel = dnAUC / nblanked_ENLYFQS)




# 5_P1' heat maps  ----------------------------------------------------------

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

# write function to generate P1' heat maps showing both replicates
heatmap_P1p_rep <- function(
  TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
                      min >= 10,
                      TEVp %in% TEVp_var,
                      TEVs != "ENLYFQ-")
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(plot_data$TEVs_AA),
    rep = unique(plot_data$rep))
  # merge with data of actually measured TEVp x TEVs combinations
  plot_data <- plot_data %>%
    full_join(., df_combinations,
              by = c("TEVp", "TEVs_AA", "rep")) %>%
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
    scale_y_discrete(expand = c(0, 0)) +
    facet_wrap(~rep)
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./01_internal_controls/plots/heatmaps/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}

# write function to generate P1' heat maps showing only one replicate
heatmap_P1p <- function(
  replicate, TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
                      min >= 10,
                      rep == replicate,
                      TEVp %in% TEVp_var,
                      TEVs != "ENLYFQ-")
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(plot_data$TEVs_AA),
    rep = replicate)
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
    scale_y_discrete(expand = c(0, 0))
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./01_internal_controls/plots/heatmaps/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


## dnAUC: noise-blanked AUCs ----

# heat map of all variants showing blanked AUC values
heatmap_P1p_rep(TEVp_var = pull(data, TEVp),
                activity = "dnAUC",
                legend = "dnAUC",
                file_name = "all_variants_dnAUC",
                height = 5, 
                asp_ratio = 1.618)

# heat map of all TEVp variants Luzi performed in vitro tests with
heatmap_P1p_rep(TEVp_var = c('0', 'VI', 'VI*', 'VII', 'NX1', 'NX2'),
                activity = "dnAUC",
                legend = "dnAUC",
                file_name = "comp_variants_dnAUC",
                height = 3.71,
                asp_ratio = 1.618)

## dnAUCrel: noise-blanked AUCs rel to wt TEVs ----

# heat map of all TEVp variants Luzi performed in vitro tests with
heatmap_P1p_rep(TEVp_var = c('0', 'VI', 'VI*', 'VII', 'NX1', 'NX2'),
                activity = "dnAUCrel",
                legend = "dnAUC rel",
                file_name = "comp_variants_dnAUC_rel",
                height = 3.71,
                asp_ratio = 1.618)




