library(cluster)
library(ggpubr)
library(pROC)
library(vegan)
library(ggrepel)
library(gridExtra)
library(wesanderson)
library(tictoc)
library(patchwork)
library(tidyverse)

options(scipen = 999)

d2 <- function(m) {
  1 - (m$deviance[1] / m$null.deviance[1])
}

colours <- c('1: spring' = '#adbc40',
             '2: summer' = '#ffa300',
             '3: autumn' = '#8f8540',
             '4: winter' = '#5188c7',
             '5: range' = '#ff004d')
spe_wide <- read_csv('data\\temperature_models\\spe_wide.csv')
all_mic <- read_csv('data\\temperature_models\\temperature_variables.csv')
growthform <- read_csv('data/sp_growth-form.csv')
legend <- read_csv('data/preparation/legend.csv')
all_mic_sc <- read_csv('data/temperature_models/temperature_variables_scaled.csv') |>
  select_if(~any(length(unique(.)) > 5)) |>
  pivot_longer(-1) |>
  pivot_wider(names_from = logger_ID) |>
  column_to_rownames('name')
step <- read_csv('data/temperature_models/SPECIES_temperature_models.csv')
groups <- read_csv('data/temperature_models/var_groups.csv')
#' -------------------------------------------------------------------------
#' modelling
#' -------------------------------------------------------------------------
tic()
spe_wide |>
  mutate(across(!logger_ID, ~ifelse(.x != 0, 1, 0))) |>
  select(where(~sum(.) >= 30)) |>
  pivot_longer(-1, values_to = 'pa') |>
  filter(!name %in% c('Alchemilla vulgaris agg.', 'Taraxacum sp.', 'Festuca halleri agg.', 'Euphrasia sp.',
                      'Phleum alpinum agg.')) |>
  rename(species = name) |>
  left_join(all_mic) |>
  pivot_longer(-c(logger_ID, species, pa)) |>
  group_by(species, name) |>
  group_nest() |>
  mutate(models = map(data, ~glm(pa ~ poly(value, 2), data = .x, family = binomial)),
         r2 = map_dbl(models, d2)) |>
  select(species, name, models, r2) |>
  write_csv('data/temperature_models/SPECIES_temperature_models.csv')
toc()

#' =========================================================================
#' Clustering of temperature variables
#' =========================================================================
for (i in 1:10) {
  suppressWarnings(
    print(
      paste0(i, ': ',
             round(
               mean(
                 as_tibble(
                   silhouette(
                     cutree(
                       hclust(
                         dist(all_mic_sc),
                         method = 'ward.D2'), i),
                     dist(all_mic_sc)))$sil_width), 2))))
}

tibble(variable = rownames(all_mic_sc),
       clx = cutree(hclust(dist(all_mic_sc), method = 'ward.D2'), 5)) |>
  write_csv('data/temperature_models/var_groups.csv')

#' =========================================================================
#' Average lost variation
#' =========================================================================
read_csv('data/temperature_models/SPECIES_temperature_models.csv') |>
  left_join(read_csv('tempattr_outputs//mds_scores.csv'), by = c(name = 'variable')) |>
  left_join(groups, by = c(name = 'variable')) |>
  rename(varexpl = r2, variable = name) |>
  filter(varexpl > 0) |>
  separate(variable, c('variable', 'month', 'id'), sep = '\\.')|>
  mutate(var_type = substr(variable, 1, 3),
         variable = factor(variable, levels = legend$variable,
                           labels = legend$label_variable)) |>
  mutate(month = factor(month, c(1:12, 'spring', 'summer', 'autumn', 'winter', 'annual'),
                        c(month.abb, 'Spr', 'Sum', 'Aut', 'Win', 'Ann')),
         id = factor(id, c('monthly', 'seasonal', 'annual')),
         var_type = factor(var_type, levels = c('avg', 'exe', 'hct', 'dct', 'cms', 'gdd'))) |>
  mutate(varexpl = round(varexpl * 100, 4)) -> step

#' ultimately best for every species (multiple testing danger zone)
best_ultimate_step <- step |>
  group_by(species) |>
  filter(varexpl == max(varexpl)) |>
  slice(1) |>
  select(species, variable, month, goal = varexpl, clx)

best_ultimate <- best_ultimate_step |>
  mutate(var = paste(variable, month),
         selection_type = 'Pool of 258') |>
  select(species, var, selection_type)

#' lost variation for every variable and every species
lost_variation <- step |>
  select(-NMDS1, -NMDS2) |>
  left_join(best_ultimate_step |>
              select(species, goal))|>
  mutate(prop_lost = (goal - varexpl) / goal) |>
  mutate(var = paste(variable, month)) |>
  select(species, var, prop_lost, clx)

#' selection of variables to take per every species (n = 13 because of variable occurence threshold: occ > 3)
lost_variation |>
  arrange(prop_lost) |>
  filter(prop_lost < 0.05)|>
  group_by(var) |>
  mutate(occ = n()) |>
  arrange(-occ) |>
  filter(occ > 3) |>
  group_by(species) |>
  slice(1) |>
  ungroup() -> vars_to_take

#' adding 13 missing species to get all 96 in the play
best_pruned <- lost_variation |>
  anti_join(vars_to_take, by = 'species') |>
  semi_join(vars_to_take, by = 'var') |>
  slice_min(prop_lost, by = species) |>
  bind_rows(vars_to_take) |>
  select(species, var, clx) |>
  mutate(selection_type = 'Pool of 19')

#' extracting the best variable per cluster based on the number of occurrences in vars_to_take
best_per_cluster <- step |>
  mutate(var = paste(variable, month)) |>
  select(species, var, varexpl, clx) |>
  semi_join(vars_to_take |>
              distinct(var, clx, occ) |>
              group_by(clx) |>
              slice_max(occ)) |>
  group_by(species) |>
  slice_max(varexpl) |>
  mutate(selection_type = 'Pool of 3',
         var_code = factor(var))

best_per_cluster |>
  ungroup() |>
  distinct(var)

#' plotting these three selections
bind_rows(best_ultimate, best_pruned, best_per_cluster) |>
  select(-varexpl, -clx) |>
  left_join(step |>
              mutate(var = paste(variable, month)) |>
              select(species, var, varexpl, clx)) |>
  left_join(best_pruned |>
              distinct(species, sp_groups = clx)) |>
  mutate(species = paste(word(species, 1), word(species, 2))) -> compiled_set

compiled_set |>
  ggplot(aes(reorder(species, varexpl, 'max'), varexpl)) +
  facet_grid(sp_groups ~ ., drop = T, space = 'free', scale = 'free',
             switch = 'y') +
  geom_path(aes(group = species, col = selection_type)) +
  geom_point(data = ~.x |> filter(selection_type == 'Pool of 19')) +
  geom_text(data = ~.x |> filter(selection_type == 'Pool of 258'),
            aes(label = var), hjust = -0.1, size = 2.5) +
  scale_colour_manual(breaks = c('Pool of 19',
                                 'Pool of 258', NA),
                      labels = c('Difference between the pool of 3 and the pool of 19',
                                 'Difference between the pool of 19 and the pool of 258', NA),
                      values = c('black', 'red', NA)) +
  coord_flip() +
  theme_bw() +
  labs(y = 'Variation explained (pseudo-R2)') +
  scale_y_continuous(expand = c(0, 0, 0, 6)) +
  expand_limits(y = c(0)) +
  theme(legend.position = c(1, 0),
        legend.background = element_blank(),
        axis.text.y = element_text(face = 'italic'),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 0))

ggsave('tempattr_outputs/Figure_3_new.png', height = 15, width = 9)


# lost_variation |>
#   semi_join(distinct(vars_to_take, species, clx)) |>
#   semi_join(distinct(vars_to_take, var, clx)) |>
#   ggplot(aes(species, prop_lost)) +
#   geom_point(data = lost_variation |>
#     semi_join(vars_to_take |>
#                 distinct(species, clx)), size = 1,
#              colour = 'grey80') +
#   geom_point(aes(
#     alpha = prop_lost < .05,
#     shape = paste0(clx, ': ', var),
#     fill = paste0(clx, ': ', var)),
#              size = 4) +
#   scale_shape_manual(name = 'Variables', values = rep(21:25, 10)) +
#   scale_alpha_manual(values = c(.5, 1), guide = F) +
#   scale_fill_discrete(name = 'Variables') +
#   coord_flip() +
#   facet_grid(clx ~ ., drop = T, space = 'free', scale = 'free',
#              switch = 'y') +
#   geom_hline(yintercept = .05) +
#   theme_bw() +
#   labs(y = 'Proportion of lost variation in SDMs to the fitting variable') +
#   guides(fill = guide_legend(ncol = 1),
#          shape = guide_legend(ncol = 1)) +
#   theme(legend.justification = c(1, 1),
#         #legend.position = c(1, 0),
#         axis.text.x = element_text(face = 'italic'),
#         axis.title.y = element_blank(),
#         strip.background = element_blank(),
#         legend.background = element_blank(),
#         strip.placement = 'outside')
# ggsave('tempattr_outputs/Figure_3_all.png', height = 12, width = 11)

#' =========================================================================
#' Figure 2
#' =========================================================================
step |>
  group_by(var_type, clx, id, variable, month, NMDS1, NMDS2) |>
  summarise(varexpl = mean(varexpl)) |>
  mutate(clx = factor(clx, levels = c(5, 4, 3, 1, 2),
                      labels = c('1: spring', '2: summer', '3: autumn', '4: winter', '5: range'))) |>
  ggplot(aes(month, variable)) +
  geom_point(aes(fill = factor(clx), size = varexpl * varexpl), shape = 22, colour = 'white') + #size = 9.5,
  geom_text(
    aes(label = ifelse(varexpl >= 10, round(varexpl), NA),
        size = varexpl)) +
  scale_shape_manual(values = c(3, NA)) +
  scale_size_continuous(range = c(2, 10), guide = F) +
  scale_fill_manual(name = 'Cluster', values = colours) +
  #scale_fill_gradientn(colours = wes_palette("FantasticFox1")) +
  facet_grid(var_type ~ id, scales = 'free', space = 'free') +
  scale_x_discrete(drop = T) +
  scale_y_discrete(limits = rev, drop = T) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(size = 9))) +
  theme(axis.title = element_blank(), # legend.position = 'none',
        legend.justification = c(1, 1),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = .5),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.subtitle = element_text(size = 7)) -> a
a
ggsave('tempattr_outputs//Figure_2.png', a, height = 8, width = 8.2)

b <- vars_to_take |>
  distinct(variable = var) |>
  mutate(
    month = factor(str_sub(variable, start = -3)),
    variable = factor(str_sub(variable, end = -5))) |>
  left_join(step |>
              distinct(variable, month, clx, NMDS1, NMDS2))|>
    mutate(clx = factor(clx, levels = c(5, 4, 3, 1, 2),
                      labels = c('1: spring', '2: summer', '3: autumn', '4: winter', '5: range'))) |>
  ggplot(aes(NMDS1, NMDS2)) +
  geom_point(data = step |>
    distinct(NMDS1, NMDS2), colour = 'grey80', size = 2) +
  geom_point(aes(colour = clx), size = 4, shape = 3, stroke = 1.5) +
  geom_text_repel(aes(
    colour = clx,
    label = paste0(variable, ' ', month)), max.overlaps = Inf,
                  box.padding = .95,
                  size = 3,
                  show.legend = F) +
  scale_colour_manual(name = 'Cluster', values = colours, drop = T) +
  labs(y = 'NMDS2', x = 'NMDS1') + # title = '(A)',
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0, 0),
        plot.title = element_text(size = 18),
        legend.background = element_blank(),
        legend.justification = c(0, 0))
b
ggsave('tempattr_outputs/Figure_4.png', height = 5, width = 6.5)
#' =========================================================================
#' comparison with general-purpose ones
#' =========================================================================
# pruned <- read_csv('Pruned_best_variables.csv')
# pruned_best <- read_csv('Pruned_best_best_variables.csv') |>
#   summarise(varexpl = max(varexpl), .by = 'species') |>
#   mutate(variable = 'Best of cluster') |>
#   filter(species %in% pruned$species)
# models <- read_csv('data/temperature_models/SPECIES_temperature_models.csv') |>
#   filter(species %in% pruned$species)
# legend <- read_csv('data/preparation/legend.csv')


step |>
  mutate(var = paste(variable, month)) |>
  select(species, var, month, variable, varexpl) |>
  filter(var %in% c('Daily Mean T Sum',
                    'Daily Mean T Win',
                    'Daily Mean T Ann',
                    'GDD T >2C Ann',
                    'Daily Min T Ann')) |>
  mutate(month = factor(month, levels = c('Sum', 'Win', 'Ann'),
                        labels = c('Summer', 'Winter', 'Annual')),
  variable = paste0(variable, '\n', month)) |>
  bind_rows(compiled_set |>
              select(species, variable = selection_type, varexpl)) |>
  mutate(variable = fct_relevel(factor(variable), 'Pool of 19')) -> vars

summary(lm(varexpl ~ variable, data = vars))

agricolae::HSD.test(aov(varexpl ~ variable, data = vars), trt = 'variable', group = TRUE, console = TRUE)$groups |>
  rownames_to_column('variable') |>
  as_tibble() |>
  select(-varexpl) |>
  left_join(vars |> summarise(varexpl = max(varexpl), .by = variable)) -> th

c <- vars |>
  ggplot(aes(reorder(variable, desc(varexpl), 'median'), varexpl / 100)) +
  geom_violin(alpha = .3,
              show.legend = F,
              fill = '#4d93c9',
              colour = NA,
              width = 1) +
  geom_boxplot(alpha = .4, fill = 'white', width = .35, outlier.shape = NA) +
  #geom_boxplot(alpha = .4, fill = '#4d93c9', width = .2, outlier.shape = NA) +
  # geom_jitter(shape = 16, size = 2, width = .02, height = 0,
  #             aes(colour = factor(clx))) +
  #             #colour = '#4d93c9') +
  stat_summary(fun = "mean", colour = "black", size = 3.2,
               #fontface = 'bold',
               geom = "text", aes(label = round(after_stat(y) * 100, 1)),
               position = position_dodge(),
               hjust = .5,
               vjust = -.1) +
  stat_compare_means(aes(label = pgen(after_stat(p))),
                     method = "t.test", ref.group = "Pool of 258") +
  #geom_text(data = th, aes(label = groups), vjust = -.2, size = 7) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0, .05, .05)) +
  labs(x = '\nReference variable', y = 'Variation explained (pseudo-R2)') +
  scale_fill_discrete(name = 'Cluster') +
  expand_limits(y = 0) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        #legend.position = c(0,0),
        legend.background = element_blank(),
        legend.justification = c(1, 1),
        plot.title = element_text(size = 18))
c
ggsave('tempattr_outputs//Figure_5.png', c, height = 5, width = 8)

pgen <- function(p) {
  ifelse(p < 0.0001, 'p < 0.0001', paste0('p = ', round(p, 4))) }

#' =========================================================================
#' Exploring variations
#' =========================================================================
citation()
#' across all variables
step |>
  summarise(mean(varexpl),
            sd(varexpl) / sqrt(n()))

step |>
  group_by(var_type) |>
  summarise(mean(varexpl),
            sd(varexpl) / sqrt(n()))

step |>
  filter(id == 'seasonal') |>
  group_by(month) |>
  summarise(mean(varexpl),
            sd(varexpl) / sqrt(n()))

step |>
  group_by(species)|>
  arrange(-varexpl) |>
  slice(1) |>
  ungroup() |>
  summarise(mean(varexpl),
            sd(varexpl) / sqrt(n()))

vars |>
  group_by(variable)|>
  summarise(mean(varexpl),
            round(sd(varexpl) / sqrt(n()), 2))

