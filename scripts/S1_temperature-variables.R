library(readxl)
library(tictoc)
library(hms)
library(lubridate)
library(tidyverse)

#' -------------------------------------------------------------------------
#' processing the logger data
#' this is done in several steps
#' 1: raw logger data is averaged to daily values (daily_stats_meta)
#' 2: out of the daily values, the basis for the pulse events is derived (daily_extremes)
#' 3: statistics per different time frames (monthly, seasonal and annual) are calculated using a map loop
#' -------------------------------------------------------------------------
loggers <- read_csv('data/preparation/logger_data.csv') |>
  mutate(day = as.Date(date),
         month = month(date)) |>
  rename(temp = temperature)

#'1
loggers |>
  mutate(
    hct_1_above0C = temp > 0,
    hct_2_above2C = temp > 2,
    hct_3_above5C = temp > 5,
    hct_4_above10C = temp > 10,
    cms_1_above0C = ifelse(hct_1_above0C, temp, 0),
    cms_2_above2C = ifelse(hct_2_above2C, temp - 2, 0),
    cms_3_above5C = ifelse(hct_3_above5C, temp - 5, 0),
    cms_4_above10C = ifelse(hct_4_above10C, temp - 10, 0)) |>
  group_by(logger_ID, month, day, year = as.numeric(factor(day >= '2022-09-01'))) |>
  summarise(
    avg_1_mean = mean(temp),
    avg_2_sd = sd(temp),
    avg_3_range = max(temp) - min(temp),
    avg_4_min = min(temp),
    avg_5_max = max(temp),
    hct_1_above0C = sum(hct_1_above0C),
    hct_2_above2C = sum(hct_2_above2C),
    hct_3_above5C = sum(hct_3_above5C),
    hct_4_above10C = sum(hct_4_above10C),
    cms_1_above0C = sum(cms_1_above0C),
    cms_2_above2C = sum(cms_2_above2C),
    cms_3_above5C = sum(cms_3_above5C),
    cms_4_above10C = sum(cms_4_above10C)) |>
  left_join(read_csv('data\\preparation\\seasons.csv')) |>
  relocate(season, .after = day) |>
  replicate(n = 3, simplify = F) |>
map2(1:3, function(df, i) {
    if (i == 1) {
      df$month <- as.character(df$month)
      df$season <- 'monthly'
    }
    if (i == 2) {
      df$month <- df$season
      df$season <- 'seasonal'
    }
    if (i == 3) {
      df$month <- 'annual'
      df$season <- 'annual'
    }
    df
  }) |>
  bind_rows() |>
  group_by(logger_ID, month, year, season) %>%
  summarise(
    gdd_1_above0C = sum(ifelse(avg_1_mean, avg_1_mean, 0)),
    gdd_2_above2C = sum(ifelse(avg_1_mean > 2, avg_1_mean, 0)),
    gdd_3_above5C = sum(ifelse(avg_1_mean > 5, avg_1_mean, 0)),
    gdd_4_above10C = sum(ifelse(avg_1_mean > 10, avg_1_mean, 0)),
    dct_1_above0C = sum(avg_1_mean > 0),
    dct_2_above2C = sum(avg_1_mean > 2),
    dct_3_above5C = sum(avg_1_mean > 5),
    dct_4_above10C = sum(avg_1_mean > 10),
    exe_1_max = max(avg_1_mean),
    exe_2_min = min(avg_1_mean),
    avg_3_range = max(avg_1_mean) - min(avg_1_mean),
    avg_1_mean = mean(avg_1_mean),
    avg_2_sd = mean(avg_2_sd),
    avg_4_min = mean(avg_4_min),
    avg_5_max = mean(avg_5_max),
    hct_1_above0C = sum(hct_1_above0C),
    hct_2_above2C = sum(hct_2_above2C),
    hct_3_above5C = sum(hct_3_above5C),
    hct_4_above10C = sum(hct_4_above10C),
    cms_1_above0C = sum(cms_1_above0C),
    cms_2_above2C = sum(cms_2_above2C),
    cms_3_above5C = sum(cms_3_above5C),
    cms_4_above10C = sum(cms_4_above10C)) %>%
  pivot_longer(-c(logger_ID, year, month, season)) |>
    ungroup() |>
    mutate(name = paste(name, month, season, sep = '.')) |>
    select(-month, -season) |>
    pivot_wider() -> mic_raw_both_years

mic <- mic_raw_both_years |>
  group_by(logger_ID) |>
  summarise_all(mean) |>
  select(-year) |>
  arrange(logger_ID)

consider_these <- mic |>
  summarise_all(~max(table(.))) |>
  pivot_longer(everything()) |>
  filter(value < nrow(mic) * .333) |>
  pull(name)

mic |>
  select(all_of(consider_these)) |>
  write_csv('data\\temperature_models\\temperature_variables.csv')

all_mic_new <- read_csv('data\\temperature_models\\temperature_variables.csv')
all_mic_new[-1] <- scale(all_mic_new[-1])
all_mic_new <- all_mic_new |> mutate_all(~replace_na(.x, 0))
all_mic_new |> write_csv('data\\temperature_models\\temperature_variables_scaled.csv')


#' =======================================================================
#' translation of names of variables
#' =======================================================================
tribble(~variable, ~label_variable,
        'avg_1_mean', 'Daily Mean T',
        'avg_4_min', 'Daily Min T',
        'avg_5_max', 'Daily Max T',
        'avg_3_range', 'Daily Range T',
        'avg_2_sd', 'Daily SD T',

        'exe_2_min', 'Min T',
        'exe_1_max', 'Max T',

        'hct_1_above0C', '# hours with T >0C',
        'hct_2_above2C', '# hours with T >2C',
        'hct_3_above5C', '# hours with T >5C',
        'hct_4_above10C', '# hours with T >10C',

        'dct_1_above0C', '# days with T >0C',
        'dct_2_above2C', '# days with T >2C',
        'dct_3_above5C', '# days with T >5C',
        'dct_4_above10C', '# days with T >10C',

        'cms_1_above0C', 'GDH >0C',
        'cms_2_above2C', 'GDH >2C',
        'cms_3_above5C', 'GDH >5C',
        'cms_4_above10C', 'GDH >10C',

        'gdd_1_above0C', 'GDD T >0C',
        'gdd_2_above2C', 'GDD T >2C',
        'gdd_3_above5C', 'GDD T >5C',
        'gdd_4_above10C', 'GDD T >10C',
) |>
  write_csv('data/preparation/legend.csv')

#' =======================================================================
#' comparison of two studied years
#' =======================================================================
mic_raw_both_years |>
  pivot_longer(-c(logger_ID, year)) |>
  filter(name %in% consider_these) |>
  pivot_wider(names_from = year) |>
  group_by(name) |>
  summarise(c = cor(`1`, `2`)) -> z

mic_raw_both_years |>
  select(logger_ID, year, contains('annual.annual')) |>
  pivot_longer(-c(logger_ID, year)) |>
  filter(name %in% consider_these) |>
  pivot_wider(names_from = year) |>
  group_by(name) |>
  summarise(c = cor(`1`, `2`)) -> z2

391-258
#'
mean(z$c, na.rm = T)
#' 0.86
mean(z2$c, na.rm = T)

z |> ggplot(aes(c)) +
  geom_histogram(fill = 'darkblue', colour = NA) +
  geom_vline(xintercept = median(z$c, na.rm = T), colour = 'red') +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0))
