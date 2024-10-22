# load libraries
library(tidyverse)
library(viridis)
library(cowplot)
library(ggseqlogo)
library(matrixStats)


# create folders ----------------------------------------------------------

# create a folder to store the individual plots and results in
dir.create("./03_DMS_controls/")
dir.create("./03_DMS_controls/results")
dir.create("./03_DMS_controls/plots")
dir.create("./03_DMS_controls/plots/negative_controls")
dir.create("./03_DMS_controls/plots/compare")

# prepare data ------------------------------------------------------------

# read data
data <- read.table("./02_DMS_controls_false_positives/results/Data_DMS_controls_fc_corrected_min20DNAlevel.txt", sep = "\t", header = TRUE)

# format data
TEVs_singlemut <- data %>% filter(n_mut_TEVs == 1) %>%
  select(TEVp, TEVs, t1, t2, t3, t4, t5, AUC)

# define which position in TEVs is mutated
TEVs_singlemut <- TEVs_singlemut %>% mutate(TEVs_pos = mapply(function(x, y) which(x != y)[1], 
       strsplit(TEVs_singlemut$TEVs, ""), strsplit("ENLYFQS", "")))

# write AA of mutated position in separate column
TEVs_singlemut <- TEVs_singlemut %>% 
  mutate(TEVs_AA = substring(TEVs, TEVs_pos, TEVs_pos))

# rename the positions using the Schechter and Berger nomenclature
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '1'] <- "P6"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '2'] <- "P5"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '3'] <- "P4"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '4'] <- "P3"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '5'] <- "P2"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '6'] <- "P1"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '7'] <- "P1'"

# add the two new columns to the original data set
data <- full_join(data, TEVs_singlemut)

# in case of wt TEVs, TEVs_AA and TEVs_pos should state "wt"
data$TEVs_AA[data$TEVs == "ENLYFQS"] <- "wt"
data$TEVs_pos[data$TEVs == "ENLYFQS"] <- "wt"

# clean up data set
# data <- data %>% select(TEVp, TEVs, n_mut_TEVs, TEVs_pos, TEVs_AA, t1, t2, t3, t4, t5, sum, min, AUC)
data <- data %>% select(TEVp, TEVs, n_mut_TEVs, TEVs_pos, TEVs_AA, min, AUC)

## replace TEVp 0 data of manually cloned variant with data from DMS screen ----

# get the data of TEVp 0 from the ssvl library as it has much more reads than
# the manually cloned TEVp 0 plasmid in the controls data set (same experiment, 
# just different BCs, i.e. plasmids!)
data_TEVp0_DMS <- read.table("./05_DMS/results/data_TEVp0_DMS.txt", sep = "\t", header = TRUE, quote = "")

data_TEVp0_DMS <- data_TEVp0_DMS %>% select(names(data))

# rename TEVp 0 as 0
data_TEVp0_DMS$TEVp <- 0

# in case of wt TEVs, TEVs_AA and TEVs_pos should state "wt"
data_TEVp0_DMS$TEVs_AA[data_TEVp0_DMS$TEVs == "ENLYFQS"] <- "wt"
data_TEVp0_DMS$TEVs_pos[data_TEVp0_DMS$TEVs == "ENLYFQS"] <- "wt"

# remove duplicates that stem from columns not needed anymore
data_TEVp0_DMS <- distinct(data_TEVp0_DMS)

# replace the TEVp 0 data of the manually cloned single control plasmid with the
# TEVp 0 variants from the DMS library
data <- filter(data, TEVp != "0") %>%
  rbind(., data_TEVp0_DMS)

## calculate AUC relative to wt TEVs, i.e. AUCrel----

# get values of wt TEVs
wtTEVs <- data %>% filter(TEVs == "ENLYFQS") %>% select(TEVp, AUC) %>% rename("AUC_ENLYFQS" = "AUC")
# join data with these values
data <- full_join(data, wtTEVs, by = "TEVp")
# calculate normalized values
data <- data %>% mutate(AUCrel = AUC / AUC_ENLYFQS)

# 1_negative controls -----------------------------------------------------

## define noise in the data ----

# get all negative controls
neg_ctrls <- data %>%
  filter(min >= 20) %>%
  filter(n_mut_TEVs <= 1) %>% # only take the single mutants and wt TEVs
  filter(TEVp %in% c('no TEVp', '(1-13)-mCherry', 'I C151A 30*31*')) %>%
  filter(TEVs != "ENLYFQ*") %>%
  select(TEVp, TEVs, TEVs_pos, TEVs_AA, AUC, min) %>%
  arrange(desc(AUC))

# save table
write.table(neg_ctrls, file = "./03_DMS_controls/results/neg_ctrls.txt", sep = "\t", row.names=FALSE, quote = F)

# check for normal distribution of data
shapiro.test(neg_ctrls$AUC)
# p–value is lower than 0.05, so we cannot assume normality

# plot the three different negative controls to see their distribution 
p <- ggplot(neg_ctrls, aes(x = AUC, fill = TEVp)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.6, color = "black")+
  geom_density(alpha = 0.6) +
  theme_cowplot(font_size = 12) +
  labs(fill = "Negative control")
# show plot
print(p)
# save plot
save_plot(p, filename = "./03_DMS_controls/plots/negative_controls/distribution_individual.pdf",
          base_height = 3,
          base_asp = 1.6)

# test if the individual negative controls show a normal distribution
shapiro.test(neg_ctrls %>% filter(TEVp == "no TEVp") %>% pull(AUC))
shapiro.test(neg_ctrls %>% filter(TEVp == "I C151A 30*31*") %>% pull(AUC))
shapiro.test(neg_ctrls %>% filter(TEVp == "(1-13)-mCherry") %>% pull(AUC))
# none of them show a p-value of >0.05, so we cannot assume normal distribution for any of them

# test if the individual negative controls show a normal distribution if we only
# look at P1'
shapiro.test(neg_ctrls %>% filter(TEVp == "no TEVp", str_starts(TEVs, "ENLYFQ")) %>% pull(AUC))
shapiro.test(neg_ctrls %>% filter(TEVp == "I C151A 30*31*", str_starts(TEVs, "ENLYFQ")) %>% pull(AUC))
shapiro.test(neg_ctrls %>% filter(TEVp == "(1-13)-mCherry", str_starts(TEVs, "ENLYFQ")) %>% pull(AUC))
# all negative controls but mCherry (and this one almost!) show a p-value of >0.05, 
# so we could assume normal distribution for the controls looking at P1'.

# fit a normal distribution on the data
library(MASS)
fit <- fitdistr(neg_ctrls$AUC, "normal")
fit

# define a cut off as the value under which a samples falls with a 
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
save_plot(p, filename = "./03_DMS_controls/plots/negative_controls/distribution_fitted.pdf",
          base_height = 3,
          base_asp = 1.4)

# collect noise results in a table
df_noise <- data.frame(TEVs_pos = "any",
                       cutoff = cutoff,
                       mu = mu,
                       sigma = sigma)
# save table of noise results
write.table(df_noise, file = "./03_DMS_controls/results/df_noise.txt", sep = "\t", row.names=FALSE, quote = F)





# Manipulate data acc. to noise -------------------------------------------

# In general, from all the samples the only one single noise cutoff should be 
# applied. Correct all samples according to the general noise, i.e. over all
# TEVs positions, not within each position. 

# fit a normal distribution on the data
library(MASS)
fit <- fitdistr(neg_ctrls$AUC, "normal")
# define a cut off as the value under which a samples falls with a 
# probability of 95%
cutoff <- qnorm(0.95, mean = fit$estimate[1], sd = fit$estimate[2])
## unload the MASS package to avoid conflicts with tidyverse
detach("package:MASS", unload=TRUE)
# get mean
mu <- fit$estimate[1]

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


# save data table
write.table(data, file = "./03_DMS_controls/results/data_controls.txt", sep = "\t", row.names=FALSE, quote = F)



# 5_comparison to in vitro data ---------------------------------------------


## get in vitro data ----

# define working directory
DIR <- gsub("uASPIre screen4 and 5_2022_08_11/uASPIre screen4_DMS", "", getwd())

# read in vitro data
iv_data_abs <- read.table(paste0(DIR, "In vitro data by Luzius Pestalozzi/data_viniabs_Luzi.txt"),
                          sep = "\t", header = TRUE, quote = "") %>%
  rename("TEVs_pos" = "pos")

iv_data_rel <- read.table(paste0(DIR, "In vitro data by Luzius Pestalozzi/data_vinirel_Luzi.txt"),
                          sep = "\t", header = TRUE, quote = "") %>%
  rename("TEVs_pos" = "pos")

# get in vivo data
in_vivo_data <- filter(data, 
                       TEVp %in% c('0','VI','VI*','VII', 'NX1', 'NX2'), 
                       TEVs_AA != "*") %>%
  select(TEVp, TEVs, TEVs_pos, TEVs_AA, AUC, dnAUC, dnAUCrel, AUCrel)

# copy data of wildtype TEVs once for each TEVs position and append it
in_vivo_data <- rbind(
  filter(in_vivo_data, TEVs != "ENLYFQS"),
  filter(in_vivo_data, TEVs == "ENLYFQS") %>% mutate(TEVs_pos = "P6", TEVs_AA = "E"),
  filter(in_vivo_data, TEVs == "ENLYFQS") %>% mutate(TEVs_pos = "P5", TEVs_AA = "N"),
  filter(in_vivo_data, TEVs == "ENLYFQS") %>% mutate(TEVs_pos = "P4", TEVs_AA = "L"),
  filter(in_vivo_data, TEVs == "ENLYFQS") %>% mutate(TEVs_pos = "P3", TEVs_AA = "Y"),
  filter(in_vivo_data, TEVs == "ENLYFQS") %>% mutate(TEVs_pos = "P2", TEVs_AA = "F"),
  filter(in_vivo_data, TEVs == "ENLYFQS") %>% mutate(TEVs_pos = "P1", TEVs_AA = "Q"),
  filter(in_vivo_data, TEVs == "ENLYFQS") %>% mutate(TEVs_pos = "P1'", TEVs_AA = "S"))

# add in vitro data to the other data
data_compare <- full_join(in_vivo_data,
                          iv_data_abs, 
                          by = c("TEVs_pos", "TEVs_AA", "TEVp")) %>% 
  full_join(., iv_data_rel, by = c("TEVs_pos", "TEVs_AA", "TEVp"))


## calculate correlations ----

# data_compare <- data_compare %>% filter(TEVs_pos == "P1'")

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
            file = "./03_DMS_controls/results/correlations_abs.txt", sep = "\t", row.names=FALSE, quote = F)
write.table(df_cor_rel %>% pivot_wider(., names_from = method, values_from = cor),
            file = "./03_DMS_controls/results/correlations_rel.txt", sep = "\t", row.names=FALSE, quote = F)


## scatter plot absolute values ----

# use a loop to do the same for all different types of absolute values
for (i in c("AUC", "dnAUC")) {
  # plot
  p <- ggplot(data_compare, aes(x = as.numeric(viniabs), y = as.numeric(pull(data_compare, i)))) +
    geom_point(alpha = 0.5) +
    theme_cowplot(font_size = 12) +
    labs(x = "viniabs", y = i) +
    # add the correlation as text label in the plot
    annotate("text", x = -Inf, y = -Inf,
             label = paste("Sp. cor. =",
                           df_cor_abs %>% filter(x == i, method == 'spearman') %>% pull(cor)),
             hjust   = -1.5,
             vjust   = -4)
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./03_DMS_controls/plots/compare/uASPIre_vs_in_vitro_", i, ".pdf"),
            base_height = 3,
            base_asp = 1.4)
}

## scatter plot relative values ----

# use a loop to do the same for all different types of relative values
for (i in c("AUCrel", "dnAUCrel")) {
  # plot
  p <- ggplot(data_compare, aes(x = as.numeric(vinirel), y = as.numeric(pull(data_compare, i)))) +
    geom_point(alpha = 0.5) +
    scale_x_continuous(limits = c(min(c(data_compare$vinirel, pull(data_compare,i))), NA)) +
    scale_y_continuous(limits = c(min(c(data_compare$vinirel, pull(data_compare,i))), NA)) +
    theme_cowplot(font_size = 12) +
    labs(x = "vinirel", y = i) +
    # add the correlation as text label in the plot
    annotate("text", x = -Inf, y = -Inf,
             label = paste("Sp. cor. =",
                           df_cor_rel %>% filter(x == i, method == 'spearman') %>% pull(cor)),
             hjust   = -1.5,
             vjust   = -4)
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./03_DMS_controls/plots/compare/uASPIre_vs_in_vitro_", i, ".pdf"),
            base_height = 3,
            base_asp = 1.4)
}

## plotting in vivo and in vitro data together -----


### plot all of Luzi's TEVps ----

i <- "dnAUCrel"

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
  # make the positions as factors and define their order for the heat map
  mutate(TEVs_pos = factor(TEVs_pos, levels=c("P1'", "P1", "P2", "P3", "P4", "P5", "P6")))

# define order of TEVps manually
data_plot$TEVp <- factor(data_plot$TEVp, levels = c("0", "VI", "VI*", "VII", "NX1", "NX2"))


# plot
p_vitro <- ggplot(data_plot %>% filter(assay == "in vitro"),
                  aes(x = TEVs_AA,
                      y = TEVs_pos, 
                      fill = activity)) +
  geom_tile() +
  facet_grid(TEVp ~ assay) +
  # facet_grid(. ~ assay + TEVp) +
  scale_fill_viridis(option="magma") +
  theme_cowplot(font_size = 12) +
  theme(axis.line = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(fill = "Relative \nactivity") +
  scale_x_discrete(name = "Amino acid", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = "top", legend.justification = "center",
        legend.key.width = unit(0.3, "inch"))

# plot
p_vivo <- ggplot(data_plot %>% filter(assay == "in vivo"),
                 aes(x = TEVs_AA,
                     y = TEVs_pos, 
                     fill = activity)) +
  geom_tile() +
  facet_grid(TEVp ~ assay) +
  # facet_grid(. ~ assay + TEVp) +
  scale_fill_viridis(option="magma") +
  theme_cowplot(font_size = 12) +
  theme(axis.line = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(fill = "Relative \nactivity") +
  scale_x_discrete(name = "Amino acid", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = "top", legend.justification = "center",
        legend.key.width = unit(0.3, "inch"))

# Put all plots into a grid
pg <- plot_grid(p_vitro, p_vivo, ncol = 2, labels = NULL)
pg
# save grid
save_plot(pg, filename = paste0("./03_DMS_controls/plots/compare/in_vivo_vs_in_vitro_", i, "_grid.pdf"),
          base_height = 7.87, base_asp = 1)









