# Analyses for "Effects of Mindfulness Meditation on the ANT"

rm(list=ls())

source('../ab-meta/functions.R')

library(tidyverse)
library(magrittr) # because tidyverse doesn't import %$%
library(brms)
library(tidybayes)
if (get_hostname() != 'rstudio') { library(forester) }

# brm config
iter        <- '100e4' # 100e4
adapt_delta <- 0.8 # Should normally be 0.8 (default) < adapt_delta < 1
sd_prior    <- "cauchy(0, .3)"
options(mc.cores = parallel::detectCores() - 4) # 8 cores

## BOOKMARK constants
ant_vars    <- c('alerting', 'orienting', 'conflict')
data_dir    <- 'data'
sharpe_data <- '../meditation-ant/data'

## BOOKMARK functions

# set descriptives columns from raw data
set_ant_descriptives <- function(study, scores) {
  df <- ant_scores(scores) %>%
    group_by(.data$t, .data$group, .data$var) %>%
    summarise(m=mean(.data$rt), sd=sd(.data$rt)) %>%
    ungroup() %>%
    pivot_wider(names_from = c(group, t, var), values_from = c(m, sd), names_sep = '.')
  
  # values_from prefixes column names, so move prefix to end
  colnames(df) <- sub('(.*?)\\.(.*$)', '\\2.\\1', colnames(df), perl = TRUE)
  
  # replace empty descriptives columns with computed values
  bind_cols(select(study, study:executive_h), df)
}

# Return pre-post difference for ANT scores
# score = alerting/orienting/conflict
# Sign indicates the pre-post improvement (+) decline (-) in score efficiency
pre_post_ant_diff <- function(df, score) {
  oldnames <- c('treatment1.n', 'treatment1.pre.FOO.m', 'treatment1.pre.FOO.sd', 'control.pre.FOO.m',
                'control.pre.FOO.sd', 'treatment1.post.FOO.m', 'treatment1.post.FOO.sd',
                'control.post.FOO.m', 'control.post.FOO.sd')
  oldnames <- gsub('FOO', score, oldnames)
  newnames <- c('treatment.n', 'treatment.pre.m', 'treatment.pre.sd', 'control.pre.m', 'control.pre.sd',
                'treatment.post.m', 'treatment.post.sd', 'control.post.m', 'control.post.sd')
  df <- df %>%
    rename_with(~ newnames[which(oldnames == .x)], .cols = oldnames)
  
  # Set pre-post differences for treatment and control groups.
  # Difference calculated such that a positive number means an improvement in network efficiency after the intervention.

  # Larger alerting and orienting scores mean more efficient network, so +ve post-pre means more efficient
  if (score %in% c('alerting', 'orienting')) {
    df <- df %>%
      mutate(treatment.diff.m = treatment.post.m - treatment.pre.m,
             control.diff.m   = control.post.m - control.pre.m)
  } else { # conflict
    # Smaller conflict score means more efficient network, so +ve pre-post means more efficient
    df <- df %>%
      mutate(treatment.diff.m = treatment.pre.m - treatment.post.m,
             control.diff.m   = control.pre.m - control.post.m)
  }
  df <- df %>%
    mutate(treatment.diff.sd = sd_pooled(treatment.n, treatment.n, treatment.pre.sd, treatment.post.sd),
           control.diff.sd   = sd_pooled(control.n, control.n, control.pre.sd, control.post.sd),
           effect.sd         = sd_pooled(treatment.n, control.n, treatment.diff.sd, control.diff.sd)
           ) %>%
    select(study, treatment.diff.m, control.diff.m, effect.sd, treatment.n, control.n, group)
  set_effect(df)
}

# Return difference for ANT scores from 2 conditions
ant_diff <- function(df, score) {
  oldnames <- c('treatment1.n', 'treatment1.pre.FOO.m', 'treatment1.pre.FOO.sd',
                'control.pre.FOO.m', 'control.pre.FOO.sd')
  oldnames <- gsub('FOO', score, oldnames)
  newnames <- c('treatment.n', 'treatment.m', 'treatment.sd', 'control.m', 'control.sd')
  
  df <- df %>%
    rename_with(~ newnames[which(oldnames == .x)], .cols = oldnames) %>%
    mutate(effect.sd = sd_pooled(treatment.n, control.n, treatment.sd, control.sd))

  if (score == 'conflict') {
    # order means so that +ve means meditation > control (lower conflict == more efficient)
    df <- df %>%
      select(study, control.m, treatment.m, effect.sd, control.n, treatment.n, group)
    set_effect(df)
    } else { # alerting/orienting
      df <- df %>%
        select(study, treatment.m, control.m, effect.sd, treatment.n, control.n, group)
      set_effect(df)
    }
}

# create alerting, orienting and conflict effect rows for a pre-post study
# FIXME: pre_post_ant_diff() and ant_diff() currently compute SMD
ant_effect <- function(df) {
  for (var in ant_vars) {
    if (df$design == 'prepost') {
      effect <- pre_post_ant_diff(df, var)
    } else if (df$design == 'casecontrol') {
      effect <- ant_diff(df, var)
    }
    if (var == 'alerting') {
      alerting <<- bind_rows(alerting, effect)
    } else if (var == 'orienting') {
      orienting <<- bind_rows(orienting, effect)
    } else { # conflict
      conflict <<- bind_rows(conflict, effect)
    }
  }
}

## end functions

# effect sign should be positive if meditation > control
# this is harder than it looks for pre-post studies if there's a (magnitude) crossover between T1 and T2

# BOOKMARK Setup
studies <- read.csv(paste0(data_dir,'/ant_data.csv'))
alerting <- orienting <- conflict <- data.frame()

## BOOKMARK Becerra et al.
# d = post.diff - pre.diff
ant_effect(studies %>% filter(grepl('Becerra', study)))

## BOOKMARK Burger & Lockhart
# d = post - pre
ant_effect(studies %>% filter(grepl('Burger and Lockhart', study)))

## BOOKMARK Isbel & Mahar
# d = treament1 - control
ant_effect(studies %>% filter(grepl('Isbel and Mahar', study)))

## BOOKMARK Jo et al.
# d = treament1 (n=20) - control (n=20)
# SDs in raw data for Jo et al. are actually SE and need converting using Eq. D9
jo_sqrt_n = sqrt(20)
jo <- studies %>% filter(grepl('Jo et al', study)) %>%
  mutate(across(ends_with('sd'), function(.x) { .x * jo_sqrt_n } ))
ant_effect(jo)

## BOOKMARK Otten et al.
# d = treatment1 - control
otten_ant <- read.csv(paste0(data_dir, '/otten_et_al.csv')) %>%
  mutate(p = factor(p), group = recode(factor(group), mindful = 'treatment1')) %>%
  pivot_longer(cols = c(pre.alerting, pre.orienting, pre.conflict), names_to = 'ant', values_to = 'score')
otten_summary <- otten_ant %>%
  group_by(group, ant) %>%
  summarise(m=mean(score), sd=sd(score)) %>%
  pivot_wider(names_from = c(group, ant), values_from = c(m, sd), names_sep = '.')
# values_from prefixes column names, so move prefix to end
colnames(otten_summary) <- sub('(.*?)\\.(.*$)', '\\2.\\1', colnames(otten_summary), perl = TRUE)
otten <- studies %>% filter(grepl('Otten et al', study))
# replace empty descriptives columns with computed values
otten <- bind_cols(select(otten, study:executive_h), otten_summary)
ant_effect(otten)

## BOOKMARK Schotz et al.
# d = treatment1 - control
ant_effect( studies %>% filter(grepl('tz et al', study)))

## BOOKMARK Sharpe et al.
# All data has been pre-processed for exclusions and outliers

# BOOKMARK Sharpe et al. Exp. 1 (Study 6)
sharpe1_ant <- read.csv(paste0(data_dir, '/sharpe_et_al_2021_1.csv')) %>%
  mutate(t=factor(t), p = factor(p), group = factor(group))
sharpe1_ant$group <- recode_factor(sharpe1_ant$group, FAM = 'treatment1', Reading = 'control')
sharpe1_ant$t <- recode_factor(sharpe1_ant$t, `1` = 'pre', `2` = 'post')
sharpe1 <- studies %>% filter(grepl('Sharpe', study) & experiment == 1)
sharpe1 <- set_ant_descriptives(sharpe1, sharpe1_ant)
ant_effect(sharpe1)

# BOOKMARK Sharpe et al. Exp. 2 (Study 2)
sharpe2_ant <- read.csv(paste0(data_dir, '/sharpe_et_al_2021_2.csv')) %>%
  mutate(t=factor(t), p = factor(p), group = factor(group))
sharpe2_ant$group <- recode_factor(sharpe2_ant$group, FAM = 'treatment1', Reading = 'control')
sharpe2_ant$t <- recode_factor(sharpe2_ant$t, `1` = 'pre', `2` = 'post')
sharpe2 <- studies %>% filter(grepl('Sharpe', study) & experiment == 2)
sharpe2 <- set_ant_descriptives(sharpe2, sharpe2_ant)
ant_effect(sharpe2)

# BOOKMARK Sharpe et al. Exp. 3 (Study 8)
sharpe3_ant <- read.csv(paste0(data_dir, '/sharpe_et_al_2021_3.csv')) %>%
  mutate(t=factor(t), p = factor(p), group = factor(group))
sharpe3_ant$group <- recode_factor(sharpe3_ant$group, FAM = 'treatment1', Reading = 'control')
sharpe3_ant$t <- recode_factor(sharpe3_ant$t, `1` = 'pre', `2` = 'post')
sharpe3 <- studies %>% filter(grepl('Sharpe', study) & experiment == 3)
sharpe3 <- set_ant_descriptives(sharpe3, sharpe3_ant)
ant_effect(sharpe3)

# BOOKMARK Sharpe et al. Exp. 4 (Study 1 + Study 3)
# calculate descriptives from raw data
sharpe4_ant <- read.csv(paste0(data_dir, '/sharpe_et_al_2021_4.csv')) %>%
  mutate(t=factor(t), p = factor(p), group = factor(group))
sharpe4_ant$group <- recode_factor(sharpe4_ant$group, meditation = 'treatment1')
sharpe4_ant$t <- recode_factor(sharpe4_ant$t, `1` = 'pre', `2` = 'post')
sharpe4 <- studies %>% filter(grepl('Sharpe', study) & experiment == 4)
ant_effect(sharpe4) # d = post - pre

# BOOKMARK Sharpe et al. Exp. 5
sharpe5_ant <- read.csv(paste0(data_dir, '/sharpe_et_al_2021_5.csv')) %>%
  mutate(t=factor(t), p = factor(p), group = factor(group))
sharpe5_ant$group <- recode_factor(sharpe5_ant$group, FAM = 'treatment1', Waitlist = 'control')
sharpe5_ant$t <- recode_factor(sharpe5_ant$t, `1` = 'pre', `2` = 'post')
sharpe5 <- studies %>% filter(grepl('Sharpe', study) & experiment == 5)
ant_effect(sharpe5) # d = post - pre

## BOOKMARK Kwak et al.
kwak <- studies %>% filter(grepl('Kwak et al', study))

## correct incorrectly calculated SDs for Kwak et al.
# calculate a correlation for SD adjustment formula using Sharpe et al., exp. 5 (Higgins et al., 2019, Section 6.5.2.2)

adjust_sd <- function(sd1, sd2, correlation) { sqrt(sd1^2 + sd2^2 - (2 * correlation * sd1 * sd2)) }

## executive (conflict)

executive_cor <- function(time, g) {
  sharpe5_ant %>%
    filter(.$t == time & .$group == g & flanker_type %in% c('incongruent', 'congruent')) %>%
    select(p, flanker_type, rt) %>%
    group_by(p, flanker_type) %>%
    summarise(mean_rt = mean(rt)) %>%
    pivot_wider(names_from = flanker_type, values_from = mean_rt) %$% # https://r4ds.had.co.nz/pipes.html#other-tools-from-magrittr
    cor(congruent, incongruent)
}

kwak$treatment1.pre.conflict.sd <- adjust_sd(124.24, 164.18, executive_cor('pre', 'mt'))
kwak$control.pre.conflict.sd <- adjust_sd(111.36, 114.01, executive_cor('pre', 'control'))
kwak$treatment1.post.conflict.sd <- adjust_sd(106.69, 112.33, executive_cor('post', 'mt'))
kwak$control.post.conflict.sd <- adjust_sd(108.91, 138.36, executive_cor('post', 'control'))

## alerting

alerting_cor <- function(time, g) {
  sharpe5_ant %>%
    filter(t == time & group == g & cue %in% c('nocue', 'double')) %>%
    select(p, cue, rt) %>%
    group_by(p, cue) %>%
    summarise(mean_rt = mean(rt)) %>%
    pivot_wider(names_from = cue, values_from = mean_rt) %$%
    cor(nocue, double)
}

kwak$treatment1.pre.alerting.sd <- adjust_sd(156.92, 145.13, alerting_cor('pre', 'mt'))
kwak$control.pre.alerting.sd <- adjust_sd(117.33, 112.36, alerting_cor('pre', 'control'))
kwak$treatment1.post.alerting.sd <- adjust_sd(115.86, 101.07, alerting_cor('post', 'mt'))
kwak$control.post.alerting.sd <- adjust_sd(130.71, 112.39, alerting_cor('post', 'control'))

## orienting

orienting_cor <- function(time, g) {
  sharpe5_ant %>%
    filter(t == time & group == g & cue %in% c('center', 'spatial')) %>%
    select(p, cue, rt) %>%
    group_by(p, cue) %>%
    summarise(mean_rt = mean(rt)) %>%
    pivot_wider(names_from = cue, values_from = mean_rt) %$%
    cor(center, spatial)
}

kwak$treatment1.pre.orienting.sd <- adjust_sd(145.13, 140.74, orienting_cor('pre', 'mt'))
kwak$control.pre.orienting.sd <- adjust_sd(112.36, 102.69, orienting_cor('pre', 'control'))
kwak$treatment1.post.orienting.sd <- adjust_sd(101.07, 104.89, orienting_cor('post', 'mt'))
kwak$control.post.orienting.sd <- adjust_sd(112.39, 121.24, orienting_cor('post', 'control'))

# d = post - pre
ant_effect(kwak)

## BOOKMARK Tsai and Chou Exp. 1
## treatment1 (n=11+19=30) - control (n=30)
# SDs in raw data for Tsai and Chou Exp. 1 are actually SE and need converting using Eq. D9
tsai1_sqrt_n = sqrt(30)
tsai1 <- studies %>%
  filter(grepl('Tsai and Chou', study) & experiment == 1) %>%
  mutate(across(ends_with('sd'), function(.x) { .x * tsai1_sqrt_n } ))
ant_effect(tsai1)

## BOOKMARK Tsai and Chou Exp. 2
# d = post - pre
tsai2_sqrt_n = sqrt(20)
tsai2 <- studies %>%
  filter(grepl('Tsai and Chou', study) & experiment == 2) %>%
  mutate(across(ends_with('sd'), function(.x) { .x * tsai2_sqrt_n } ))
ant_effect(tsai2)

# BOOKMARK Walsh et al.
ant_effect(studies %>% filter(grepl('Walsh et al', study)))

## BOOKMARK Wittmann et al.
# treatment1 - control
ant_effect(studies %>% filter(grepl('Wittmann et al', study)))

## meta-analysis

# https://vuorre.netlify.com/post/2016/09/29/meta-analysis-is-a-special-case-of-bayesian-multilevel-modeling/
alerting  <- mutate(alerting, score = 'alerting')
orienting <- mutate(orienting, score = 'orienting')
conflict  <- mutate(conflict, score = 'conflict')
effects   <- bind_rows(alerting, orienting, conflict) %>%
  rename(Study = study)

## Forest plots

for (var in c('conflict', 'alerting', 'orienting')) {
  if (var == 'conflict') {
    ant_name <- 'Executive'
  } else {
    ant_name <- str_to_title(var)
  }
  estimate_col_name <- paste0(ant_name, ' SMD estimate')
  
  # calculate SEM from CI
  model_data <- effects %>%
    filter(score == var) %>%
    mutate(se = ci95_to_se(u, l)) %>%
    select(Study, d, se, ci, l, u, mean1, mean2, group)

  # all studies
  model_data <- model_data %>%
    mutate(study_number = as.numeric(rownames(model_data)))
  study_names <- model_data %>% select(study_number, Study)
  
  rem_all <- brm(
    d | se(se) ~ 1 + (1 | Study), # random effects meta-analyses model (see brmsformula)
    data = model_data,
    chains=8, iter=iter,
    prior = c(prior_string("normal(0,1)", class = "Intercept"),
              prior_string(sd_prior, class = "sd")),
    control = list(adapt_delta = adapt_delta),
    file = paste0('cache/all-', var, '-brms')
  )
  
  all_data <- brms_object_to_table(rem_all, effects %>% select(Study, group),
                                   cache_label = paste0('cache/all-ant-', var), subset_col = 'group',
                                   iter = iter, sd_prior = sd_prior, adapt_delta = adapt_delta)
  
  rem_not_ours <- brm(
    d | se(se) ~ 1 + (1 | Study),
    data = model_data %>% filter(!grepl('Sharpe', Study)),
    chains=8, iter=iter,
    prior = c(prior_string("normal(0,1)", class = "Intercept"),
              prior_string(sd_prior, class = "sd")),
    control = list(adapt_delta = adapt_delta),
    file = paste0('cache/not-ours-', var, 'brms')
  )

  all_except_our_data <- brms_object_to_table(rem_not_ours, effects %>% select(Study, group),
                                              cache_label = paste0('cache/not-ours-ant-', var), subset_col = 'group',
                                              iter = iter, sd_prior = sd_prior, adapt_delta = adapt_delta)
  
  # rob_blobbogram(rem, robins,
  #                iter = iter, sd_prior = sd_prior, adapt_delta = adapt_delta,
  #                subset_col = 'group',
  #                estimate_col_name = estimate_col_name,
  #                cache_label = paste('ct', var, sep = '-'),
  #                null_line_at = 0,
  #                font_family = "serif",
  #                rob_colour = "cochrane",
  #                rob_tool = "Robins",
  #                x_scale_linear = TRUE,
  #                xlim = c(-1, 1),
  #                xbreaks = c(-1, -.8, -.5, -.3, 0, .3, .5, .8, 1),
  #                arrows = FALSE,
  #                arrow_labels = c("Low", "High"),
  #                nudge_y = -0.2,
  #                estimate_precision = 2,
  #                display = FALSE,
  #                add_tests = FALSE,
  #                file_path = here::here(paste0("figures/nrct_", var, "_forest.png"))
  # )
}
