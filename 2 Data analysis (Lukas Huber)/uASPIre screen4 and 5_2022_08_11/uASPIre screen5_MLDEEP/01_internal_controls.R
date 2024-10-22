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
dir.create("./01_internal_controls/plots/reproducibility")


# prepare data ------------------------------------------------------------


# read data
data <- read.table("./data/Data_MLD_internal_controls.txt", sep = "\t", header = TRUE)

# format data
data <- data %>% 
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data$TEVs)) # write AA of mutated position in separate column)

# correct entry for TEVs of positive control
data$TEVs[data$Note == "positive control"] <- "ENLYFQ-"

# clean up data set 
data <- data %>% select(TEVp, Note, TEVs, TEVs_pos, TEVs_AA, t1, t2, t3, t4, t5, sum, min, AUC)



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
# p–value is lower than 0.05, so we cannot assume the normality

# plot the three different negative controls to see their distribution 
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
# all of them show a p-value of >0.05, so we can assume normal distribution for all of them

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
  theme_cowplot() +
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


# 4_reproducibility --------------------------------------------------------------

# define path for data from other experiments
DIR <- gsub("uASPIre screen4 and 5_2022_08_11/uASPIre screen5_MLDEEP", "", getwd())

## MLDDEP vs. ABC screen ----

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

# filter data of MLDEEP screen likewise
corr_data_MLDEEP <- data %>%
  filter(min >= 20) %>%
  filter(TEVp != "MBP") %>%   # remove MBP as it was underrepresented and is not a good negative control
  filter(Note != "neg ctrl catalysis") %>% # kick out C151A variants
  select(TEVp, TEVs, AUC) %>%
  mutate(rep = "MLDEEP")

# combine data of both screens and bring to wider format for plotting
replicates <- rbind(corr_data_ABC, corr_data_MLDEEP) %>%
  pivot_wider(., names_from = rep, values_from = AUC)
  

### ABC screen rep 1 ----

# create linear model
fit <- lm(ABC_1 ~ MLDEEP, data = replicates)

# plot
p <-  ggplot(replicates, aes(x = MLDEEP, y = ABC_1)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4),
                          ", r = ",
                          round(cor(replicates$MLDEEP, replicates$ABC_1, use = "pairwise.complete.obs", method = 'pearson'), 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("AUC MLDEEP screen") +
  ylab("AUC ABC screen rep 1") +
  scale_x_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_internal_controls/plots/reproducibility/MLDEEP_vs_ABCrep1.pdf",
          base_height = 3,
          base_asp = 1.2)

### ABC screen rep 2 ----

# create linear model
fit <- lm(ABC_2 ~ MLDEEP, data = replicates)

# plot
p <-  ggplot(replicates, aes(x = MLDEEP, y = ABC_2)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4),
                          ", r = ",
                          round(cor(replicates$MLDEEP, replicates$ABC_2, use = "pairwise.complete.obs", method = 'pearson'), 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("AUC MLDEEP screen") +
  ylab("AUC ABC screen rep 2") +
  scale_x_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_internal_controls/plots/reproducibility/MLDEEP_vs_ABCrep2.pdf",
          base_height = 3,
          base_asp = 1.2)

## MLDDEP vs. DMS screen ----

# also compare the internal controls of the MLDEEP screen with the ones
# from the DMS screen

# get data from DMS screen

data_DMS <- read.table(paste0(DIR, "/uASPIre screen4 and 5_2022_08_11/uASPIre screen4_DMS/data/Data_DMS_internal_controls.txt"),
                       sep = "\t", header = TRUE, quote = "")

# correct entry for TEVs of positive control
data_DMS$TEVs[data_DMS$Note == "positive control"] <- "ENLYFQ-"

# filter and clean up data set 
corr_data_DMS <- data_DMS %>%
  # filter for internal controls
  filter(libraryID %in% c('libLH015', 'libLH016', 'libLH019')) %>%
  # set read threshold
  filter(min >= 20) %>%
  # remove MBP as it was underrepresented and is not a good negative control
  filter(TEVp != "MBP") %>%   
  # kick out C151A variants
  filter(Note != "neg ctrl catalysis") %>% 
  select(TEVp, TEVs, AUC) %>%
  mutate(rep = "DMS")

# combine data of both screens and bring to wider format for plotting
replicates <- rbind(corr_data_ABC, corr_data_MLDEEP, corr_data_DMS) %>%
  pivot_wider(., names_from = rep, values_from = AUC)

# create linear model
fit <- lm(DMS ~ MLDEEP, data = replicates)

# plot
p <-  ggplot(replicates, aes(x = MLDEEP, y = DMS)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4),
                          ", r = ",
                          round(cor(replicates$MLDEEP, replicates$DMS, use = "pairwise.complete.obs", method = 'pearson'), 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("AUC MLDEEP screen") +
  ylab("AUC DMS screen") +
  scale_x_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_internal_controls/plots/reproducibility/MLDEEP_vs_DMS.pdf",
          base_height = 3,
          base_asp = 1.2)


