#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(ggthemes)

theme_set(theme_few())
scale_colour_discrete <- function(...) scale_colour_few()
scale_fill_discrete <- function(...) scale_fill_few()
# scale_colour_discrete <- function(...) scale_colour_brewer(..., palette="Set3")
# scale_fill_discrete <- function(...) scale_colour_brewer(..., palette="Set3")


#
# Data sets in order of elapsed time
data_set_levels <- c("Buettner", "Pollen", "Usoskin", "Kolodziejczyk", "Zeisel")

#
# Load timings data
timings <- readr::read_csv('timings-all.csv')
timings$data.set <- factor(timings$data.set, levels = data_set_levels)

#
# Melt data
timings_m <-
  melt(timings,
       id.vars = c('data.set', 'niter', 'task'),
       variable.name = 'timed',
       value.name = 'time')
timings_m %>% sample_n(15)

#
# Order the data sets by elapsed time
data_set_elapsed <-
  timings_m %>%
  filter(timed == 'elapsed') %>%
  group_by(data.set) %>%
  summarise(elapsed = sum(time)) %>%
  arrange(elapsed)

#
# Focus on interesting timings
timings_m %>%
  filter(timed == 'elapsed') %>%
  group_by(task) %>%
  summarise(elapsed = sum(time)) %>%
  arrange(-elapsed)
uninteresting <- c('eigen', 'k.means', 'convergence', 'lambda',
                   'normalise', 'impute')
interesting_m <-
  timings_m %>%
  filter(timed %in% c('user.self', 'elapsed'),
         ! task %in% uninteresting)

#
# Plot interesting timings
ggplot(interesting_m, aes(x = data.set, y = time, fill = task)) +
  geom_col(position = 'dodge') +
  scale_y_log10() +
  facet_wrap(~ timed)
ggsave('output/timings-by-dataset.pdf')

ggplot(interesting_m %>% left_join(data_set_elapsed),
       aes(x = elapsed, y = time, colour = task)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ timed)
ggsave('output/timings-by-total-elapsed.pdf')
