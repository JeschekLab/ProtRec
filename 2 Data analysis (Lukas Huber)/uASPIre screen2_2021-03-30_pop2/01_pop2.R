# load libraries
library(tidyverse)
library(viridis)
library(cowplot)

# create folders --------------------------------------------------------

# create a folder to store the individual plots and results in
dir.create("./01_pop2")
dir.create("./01_pop2/results")
dir.create("./01_pop2/plots")
dir.create("./01_pop2/plots/1_flipping_TEVp")
dir.create("./01_pop2/plots/2_negative_controls")
dir.create("./01_pop2/plots/3_binding")
dir.create("./01_pop2/plots/4_linker_study")
dir.create("./01_pop2/plots/5_reproducibility")
dir.create("./01_pop2/plots/6_heatmaps")

# prepare data ------------------------------------------------------------

# read data
data <- read.table("./data/01_std_fraction_AUCs.txt", sep = "\t", header = TRUE)

# rename some variants to have consistent naming
# Replace "TEVp" with "" for all entries except "no TEVp"
data$TEVp <- ifelse(data$TEVp != "no TEVp", gsub("TEVp ", "", data$TEVp), data$TEVp)
data$TEVp <- ifelse(data$TEVp != "no TEVp", gsub("TEVp", "", data$TEVp), data$TEVp) # to rename TEVp(1-13)-mCherry

# format data
data <- data %>% 
  # for this analysis we only look at the standard plasmid with GFP
  filter(Bxb1.version == "with GFP") %>%
  # rename some columns to be consistens with other data
  rename(rep=Rep, 
         t1=X0, t2=X1, t3=X2, t4=X3, t5=X4, t6=X5, t7=X6, t8=X7, t9=X8,
         Note = Comment) %>%
  # clean up
  select(-ID, -Bxb1.version, -Architecture, -Discriminator.spacer, -Plasmid.ID, -TEVp.BC, -TEVs.BC, -sum)

# add columns needed for some plots
data <- data %>%
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data$TEVs)) # write AA of mutated position in separate column)

# get TEVp variants that were part of the linker study
linker_vars <- data %>% filter(No..of.TEVs.linkers == 1) %>% pull(TEVp) %>% unique()

# separate data for linker study from the rest
data_linkers <- data %>%
  filter(TEVs == 'ENLYFQS',
         TEVp %in% linker_vars)

# save data of linker study
write.table(data_linkers, file = "./01_pop2/results/data_linkers.txt", sep = "\t", row.names=FALSE, quote = F)

# only keep variants with both linkers in the rest of the data table
data <- filter(data, No..of.TEVs.linkers == 2) %>% select(-No..of.TEVs.linkers)

# save data
write.table(data, file = "./01_pop2/results/data_pop2clean.txt", sep = "\t", row.names=FALSE, quote = F)

# 1_flipping plots ---------------------------------------------------------

## multiple TEVps and a single TEVs ----

# write function to generate plots of flipping curves for a single TEVs
TEVp_plot <- function(
  TEVs_var, TEVp_var, limit, file_name){
  p <- data %>%
    # filtering
    filter(TEVp %in% TEVp_var, 
           # min >= 20, 
           TEVs %in% TEVs_var) %>% 
    # selecting only necessary columns
    select(rep, TEVp, TEVs, t1, t2, t3, t4, t5, t6, t7, t8, t9) %>%
    # rename for correct plotting
    rename("0"=t1, "1"=t2, "2"=t3, "3"=t4, "4"=t5, "5"=t6, "6"=t7, "7"=t8, "8"=t9) %>%
    # reformatting for plotting
    pivot_longer(cols = !c(rep, TEVp, TEVs) , names_to = "timepoint", values_to = "fraction_flipped") %>%
    # plotting
    ggplot(., aes(x = as.integer(timepoint), y = fraction_flipped, col = TEVp)) +
    geom_point() +
    geom_line() +
    facet_wrap(~rep) +
    theme_cowplot(font_size = 12) +
    scale_x_continuous("Time after induction (h)", limits = c(0, 8)) +
    scale_y_continuous("Fraction flipped", limits = c(0, limit)) +
    ggtitle(paste("TEVs:", TEVs_var)) +
    theme(plot.title = element_text(hjust = 0.5))
  # show plot
  print(p)
  # and save to file
  save_plot(p, filename = paste0("./01_pop2/plots/1_flipping_TEVp/", file_name, ".pdf"))
}


# show all active variants
TEVp_plot(TEVs_var = "ENLYFQS", 
          TEVp_var = c("(1-13)-mCherry",
                       filter(data, Note == "") %>% pull(TEVp)),
          limit = 0.75, file_name = "active_variants")


### plot the pop for a paper figure ----

palette.colors(palette = "R4")

# prepare data
pop <- data %>%
  # filtering
  filter(TEVp %in% c("0", "0 C151A", "I", "I C151A", "(1-13)-mCherry", "0 C151A no ssrA"), 
         # min >= 20, 
         TEVs %in% c("ENLYFQS", "ENLYFQ-"),
         # show rep 2 only as rep 1 has missing values for the no ssrA control
         rep == 2) %>% 
  # selecting only necessary columns
  select(rep, TEVp, TEVs, t1, t2, t3, t4, t5, t6, t7, t8, t9, AUC4, AUC8) %>%
  # rename for correct plotting
  rename("0"=t1, "1"=t2, "2"=t3, "3"=t4, "4"=t5, "5"=t6, "6"=t7, "7"=t8, "8"=t9) %>%
  # add column specifying the line type for each TEVp
  mutate(type = factor(if_else(TEVp == "0 C151A no ssrA", "dashed", "solid"), 
                       levels = c("solid", "dashed"))) %>%
  # reformat for plotting
  pivot_longer(cols = !c(rep, TEVp, TEVs, type, AUC4, AUC8) , names_to = "timepoint", values_to = "fraction_flipped")

# create a named vector with color codes for each TEVp
tevp_colors <- c("0" = palette.colors(palette = "R4")[4], 
                 "0 C151A" = palette.colors(palette = "R4")[5], 
                 "I" = palette.colors(palette = "R4")[6], 
                 "I C151A" = palette.colors(palette = "R4")[2], 
                 "(1-13)-mCherry" = palette.colors(palette = "R4")[8], 
                 "0 C151A no ssrA" = palette.colors(palette = "R4")[1])

# plotting
p <- ggplot(pop, aes(x = as.integer(timepoint), y = fraction_flipped, col = TEVp, linetype = type)) +
  geom_point() +
  geom_line() +
  theme_cowplot(font_size = 12) +
  scale_x_continuous("Time after induction (h)", limits = c(0, 8)) +
  scale_y_continuous("Fraction flipped", limits = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = tevp_colors) +
  guides(linetype = "none")

# show plot
print(p)
# and save to file
save_plot(p, filename = paste0("./01_pop2/plots/1_flipping_TEVp/figure_1C.pdf"),
          base_height = 2.3,
          base_asp = 1.7)


# 2_negative controls -----------------------------------------------------

# create a correlation plot

# format data
neg_ctrls <- data %>%
  #  filter(min >= 20) %>%
  filter(TEVp %in% c('no TEVp', '(1-13)-mCherry'),
         TEVs == "ENLYFQS") %>%
  select(TEVp, rep, paste0("AUC", seq(1:8))) %>%
  pivot_longer(., cols = c(-TEVp, -rep), names_to = "hour", values_to = "AUC") %>%
  pivot_wider(., names_from = TEVp, values_from = AUC) %>%
  rename('mCherry' = '(1-13)-mCherry', "noTEVp" = "no TEVp")

# create linear model
fit_rep1 <- lm(noTEVp ~ mCherry, data = filter(neg_ctrls, rep == 1))
fit_rep2 <- lm(noTEVp ~ mCherry, data = filter(neg_ctrls, rep == 2))
neg_ctrls <- data.frame(rep = c(1:2),
                 r2 = c(round(summary(fit_rep1)$r.squared, 5),
                        round(summary(fit_rep2)$r.squared, 5))) %>%
  right_join(., neg_ctrls)


# plot
p <-  ggplot(neg_ctrls, aes(x = mCherry, y = noTEVp)) +
  geom_point(size = 2, aes(color = hour)) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE, size =0.75) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab("TEVp(1-13)-mCherry") +
  ylab("no TEVp") +
  scale_x_continuous(limits = c(0, 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.2), expand = c(0, 0)) +
  labs(color = NULL) + # remove label of legend
  facet_wrap(~rep)
# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_pop2/plots/2_negative_controls/correlation_neg_ctrls.pdf",
          base_asp = 1.618)


# 3_effect of binding -------------------------------------------------------

# investigate the effect of binding

## ENLYFQS only ----

### Fold change over mCherry (AUCs) ----

# get AUC1-8 for some TEVps
check <- data %>%
  filter(TEVp %in% c("I", "0", "I C151A", "0 C151A"),
         TEVs == "ENLYFQS") %>%
  select(rep, TEVp, paste0("AUC", seq(1, 8)))

# get AUC1-8 for mCherry
mCherry <- data %>%
  filter(TEVp == "(1-13)-mCherry",
         TEVs == "ENLYFQS") %>%
  select(rep, TEVp, paste0("AUC", seq(1, 8)))

# make the mCherry data frame as big as the check data frame
mCherry <- left_join(select(check, rep), mCherry, by = "rep")

# calculate fold changes over mCherry  
fc <- select(check, -rep, -TEVp)/select(mCherry, -rep, -TEVp)
# add the TEVp and rep columns again
fc <- fc %>% mutate(TEVp = check$TEVp, rep = check$rep)

# create a named vector with color codes for each TEVp
tevp_colors <- c("0" = palette.colors(palette = "R4")[4], 
                 "0 C151A" = palette.colors(palette = "R4")[8], 
                 "I" = palette.colors(palette = "R4")[6], 
                 "I C151A" = palette.colors(palette = "R4")[1])

# plotting
p <- ggplot(fc %>%
              pivot_longer(cols = -c(TEVp, rep), values_to = "fc", names_to = "AUC") %>%
              mutate(AUC = gsub("AUC", "", AUC)),
            aes(x = AUC, y = fc, col = TEVp)) +
  geom_point() +
  geom_line(aes(group = TEVp)) +
  facet_wrap(~rep) +
  theme_cowplot(font_size = 12) +
  scale_color_manual(values = tevp_colors) +
  scale_y_continuous("Fold change over mCherry", breaks = seq(1, 8, by = 1)) +
  geom_hline(yintercept=1, linetype="dashed") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Time after induction (h)")

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_pop2/plots/3_binding/fold_change.pdf", 
          base_height = 3.5)




### bar plot ----

## show effect of binding with ENLYFQS as stacked bar plot 

# prepare data
binding_ENLYFQS <- data %>%
  # filter(min >= 20) %>% # set read threshold
  filter(TEVp %in% c("0", "0 C151A", "I", "I C151A", "(1-13)-mCherry"), TEVs == "ENLYFQS") %>%
  # do this analysis with AUC4 only
  rename(AUC = AUC4) %>%
  select(rep, TEVp, Note, AUC)

# get AUCs for all active and catalytic dead variants
data_bars <- binding_ENLYFQS %>%
  filter(TEVp != "(1-13)-mCherry") %>%
  mutate(TEVp = gsub(" C151A", "", TEVp)) # remove the "C151A" from the TEVp name

# rename entries of the Note column to be more precise
data_bars$Note[data_bars$Note == ""] <- "active"
data_bars$Note[data_bars$Note == "neg ctrl catalysis"] <- "inactive"

# rearrange table
data_bars <- data_bars %>% pivot_wider(., names_from = Note, values_from = AUC)

# define the blank and get its values
blank <- binding_ENLYFQS %>%
  filter(Note == "neg ctrl binding") %>%
 select(rep, AUC) %>%
  rename("blank" = "AUC")

# join data with these blank values
data_bars <- full_join(data_bars, blank, by = "rep")

# calculate effect of binding and catalysis 
data_bars <- mutate(data_bars,
                    # binding = inactive (C151A) variants - mCherry control
                    binding = inactive - blank,
                    # catalysis = active TEVp variant - inactive TEVp variant
                    catalysis = active - inactive) %>%
  mutate(binding_rel = 100*round(binding/(binding+catalysis), 3),
         catalysis_rel = 100*round(catalysis/(binding+catalysis), 3))

# Stacked bar plot
# rearrange the table for plotting
p <- data_bars %>%
  select(-active, -inactive, -blank) %>%
  rename(binding_abs = binding, catalysis_abs = catalysis) %>%
  pivot_longer(cols = -c(TEVp, rep), names_to = c("effect", "variable"), names_sep = "_",
               values_to = "value") %>%
  pivot_wider(., names_from = variable, values_from = value) %>%
  ggplot(., aes(fill = effect, y = abs, x = factor(TEVp, levels = arrange(filter(data_bars, rep =="1"), active)$TEVp))) +
  geom_bar(position="stack", stat="identity") +
  theme_cowplot(font_size = 12) +
  xlab("TEVp") +
  ylab("AUC blanked to mCherry") +
  geom_text(aes(label = paste(rel, "%")), position = position_stack(vjust = 0.5), size=4, color = "white") +
  facet_wrap(~rep)

# show plot
print(p)

# save plot
save_plot(p, filename = "./01_pop2/plots/3_binding/effect_of_binding_ENLYFQS.pdf",
          base_height = 3.5, 
          base_asp = 1.2)


# 4_linker study ---------------------------------------------------------------

# compare the flipping behavior of plasmids containing 0, 1, or 2 GS-rich linkers 
# around TEVs

# prepare data
linkerstudy <- data_linkers %>%
  # filtering
  filter(TEVp %in% c("0", "0 C151A", "I", "I C151A", "(1-13)-mCherry"), 
         # show rep 2 only as rep 1 has missing values for the no ssrA control
         rep == 2) %>% 
  # selecting only necessary columns
  select(rep, TEVp, No..of.TEVs.linkers, t1, t2, t3, t4, t5, t6, t7, t8, t9) %>%
  # rename for correct plotting
  rename("0"=t1, "1"=t2, "2"=t3, "3"=t4, "4"=t5, "5"=t6, "6"=t7, "7"=t8, "8"=t9) %>%
  # reformat for plotting
  pivot_longer(cols = !c(rep, TEVp, No..of.TEVs.linkers) , names_to = "timepoint", values_to = "fraction_flipped")

# plot linker study with TEVp I
p_TEVpI <- ggplot(filter(linkerstudy, TEVp %in% c("I", "I C151A", "(1-13)-mCherry")), aes(x = as.integer(timepoint), y = fraction_flipped, col = TEVp)) +
               geom_point() +
               geom_line() +
               theme_cowplot(font_size = 12) +
               scale_x_continuous("Time after induction (h)", limits = c(0, 8)) +
               scale_y_continuous("Fraction flipped", limits = c(0, 1)) +
               facet_wrap(~No..of.TEVs.linkers) +
               scale_color_manual(values = c(8, 6, 1)) +
               theme(panel.spacing = unit(0.5, "cm"))
             
# show plot
print(p_TEVpI)
# and save to file
save_plot(p_TEVpI, filename = paste0("./01_pop2/plots/4_linker_study/linker_study_TEVpI.pdf"),
          base_height = 2.5,
          base_asp = 3.6)





# 5_reproducibility --------------------------------------------------------------

## among replicates ----

# create a correlation plot

# define the AUC to use
AUC <- "AUC4"

# format data
replicates <- data %>%
#  filter(min >= 20) %>%
  filter(Note != "neg ctrl catalysis") %>% # kick out C151A variants
  select(TEVp, TEVs, all_of(AUC), rep) %>%
  pivot_wider(., names_from = rep, values_from = AUC, names_prefix = "rep")

# create linear model
fit <- lm(rep2 ~ rep1, data = replicates)

# plot
p <-  ggplot(replicates, aes(x = rep1, y = rep2)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R^2 = ", 
                          round(summary(fit)$r.squared, 4),
                          ", SE = ",
                          round(summary(fit)$coefficients[2, "Std. Error"], 4),
                          ", r = ",
                          round(cor(replicates$rep2, replicates$rep1, use = "pairwise.complete.obs", method = 'pearson'), 4)),
           hjust = -0.1, vjust = 1) +
  theme_cowplot(font_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) +
  xlab(paste(AUC, "replicate 1")) +
  ylab(paste(AUC, "replicate 2"))

# show plot
print(p)
# and save to file
save_plot(p, filename = "./01_pop2/plots/5_reproducibility/replicates.pdf",
          base_asp = 1)


# 6_P1' heat maps  ----------------------------------------------------------
# show a comprehensive P1' specificity profile of all TEVp varians

# prepare data, use AUC4 only
data_P1p <- data %>% select(rep, TEVp, TEVs, TEVs_pos, TEVs_AA, AUC4) %>%
  rename(AUC = AUC4)

## raw AUCs rel to wt TEVs ----

# set raw AUCs relative to wt TEVs

# get values of wt TEVs
wtTEVs <- data_P1p %>% filter(TEVs == "ENLYFQS") %>% select(rep, TEVp, AUC) %>% rename("AUC_ENLYFQS" = "AUC")
# join data with these values
data_P1p <- full_join(data_P1p, wtTEVs, by = c("TEVp", "rep"))
# calculate normalized blanked values
data_P1p <- data_P1p %>% mutate(AUCrel = AUC / AUC_ENLYFQS)


## prepare plotting ----

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
  plot_data <- filter(data_P1p,
                      # min >= 10,
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
  save_plot(p, filename = paste0("./01_pop2/plots/6_heatmaps/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}

## AUC: raw AUCs ----

# heat map of all variants showing raw AUC values, both reps
heatmap_P1p_rep(TEVp_var = pull(data_P1p, TEVp),
                activity = "AUC",
                legend = "AUC",
                file_name = "all_variants_AUC",
                height = 5, 
                asp_ratio = 1.618)

# heat map of all active TEVp variants
heatmap_P1p_rep(TEVp_var = c('0', 'I', 'IV', 'VI', 'VI*', 'VII', 'VII*', 'NX1', 'NX2'),
                activity = "AUC",
                legend = "AUC",
                file_name = "selected_TEVps_AUC",
                height = 3.71,
                asp_ratio = 1.618)

## AUCrel: AUC rel to wt TEVs ----

# heat map of all active TEVp variants
heatmap_P1p_rep(TEVp_var = c('0', 'I', 'IV', 'VI', 'VI*', 'VII', 'VII*', 'NX1', 'NX2'),
            activity = "AUCrel",
            legend = "AUC rel",
            file_name = "selected_TEVps_AUC_rel",
            height = 3.71,
            asp_ratio = 1.618)







