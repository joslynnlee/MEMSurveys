#' This code generates plots of MEM survey data
#'
#' last modified: 2022-08-19
#' Chris Miller
#' chris.miller@ucdenver.edu

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(ggbeeswarm)

#############################################
# PATHS TO KEY DATA
#############################################

# paired survey data:
fp <-'MEM paired survey data.csv'
# pre-MEM questions on western science
fp2 <-'MEM western science survey data.csv'

#############################################
# plotting theme
# starts with theme_pubclean from ggpubr
#############################################
the_theme <- theme_pubclean() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title= element_blank(),
        legend.position = "right",
        axis.text = element_text(size=9),
        axis.title = element_text(size=9),
        strip.text = element_text(size = 9),
        legend.text=element_text(size = 8),
        legend.margin = margin(-5), 
        legend.box.margin = margin(l=-9))

# some other visualization constants not set in themes
beeswarm_w <- 0.22        # for position_quasirandom
highlight_l <- "magenta"  # color to highlight larger jump lines with
strip.position <- "bottom" # must be set to bottom if going to tag plots with (A), (B), etc.

# MODIFIED tag_facet_outside from egg package 
# changed to keep strip chart text (only works on bottom)
# used in code below for labeling subplots with (A), (B), etc.
# This is a bit of a hack, but this is otherwise a nice way to handle this.
# https://rdrr.io/github/baptiste/egg/src/R/tag_facet.r
tag_facet_outside2 <-  function(p, open=c("(",""), close = c(")","."),
                                tag_fun_top = function(i) LETTERS[i],
                                tag_fun_right = utils::as.roman,
                                x = c(0,0), y = c(0.5, 1),
                                hjust = c(0,0), vjust = c(0.5,1), 
                                fontface = c(2,2), family="", draw = TRUE, ...){
  
  # gb <- ggplot_build(p + theme(strip.text = element_blank(), 
  #                              strip.background = element_blank()))
  gb <- ggplot_build(p)  # actually, don't take those strip chart titles away... 
  # hack.  works only / best with strip.position="bottom" for faceting
  
  lay <- gb$layout$layout
  
  tags_top <- paste0(open[1],tag_fun_top(unique(lay$COL)),close[1])
  tags_right <- paste0(open[2],tag_fun_right(unique(lay$ROW)),close[2])
  
  tl <- lapply(tags_top, grid::textGrob, x=x[1], y=y[1],
               hjust=hjust[1], vjust=vjust[1], 
               gp=grid::gpar(fontface=fontface[1], fontfamily = family, ...))
  rl <- lapply(tags_right, grid::textGrob, x=x[2], y=y[2],
               hjust=hjust[2], vjust=vjust[2], 
               gp=grid::gpar(fontface=fontface[2], fontfamily = family, ...))
  
  
  g <- ggplot_gtable(gb)
  g <- gtable::gtable_add_rows(g, grid::unit(1,"line"), pos = 0)
  l <- unique(g$layout[grepl("panel",g$layout$name), "l"])
  g <- gtable::gtable_add_grob(g, grobs = tl, t=1, l=l)
  
  wm <- do.call(grid::unit.pmax, lapply(rl, grid::grobWidth))
  g <- gtable::gtable_add_cols(g, wm, pos = max(l))
  t <- unique(g$layout[grepl("panel",g$layout$name), "t"])
  g <- gtable::gtable_add_grob(g, grobs = rl, t=t, l=max(l) + 1)
  g <- gtable::gtable_add_cols(g, unit(2,"mm"), pos = max(l))
  
  if(draw){
    grid::grid.newpage()
    grid::grid.draw(g)
  }
  invisible(g)
}


#############################################
# read in and wrangle pre/post survey data
#############################################

MEM_paired_survey_data <- read_csv(fp)
MEM_paired_survey_data <- MEM_paired_survey_data %>% rename_at(1, ~"ID")

pre <- MEM_paired_survey_data %>% select(ID | ends_with('pre')) %>% 
  rename_with(~str_remove(., ' pre')) %>% 
  mutate(when="pre")
post <- MEM_paired_survey_data %>% select(ID | ends_with('post')) %>% 
  rename_with(~str_remove(., ' post')) %>% 
  mutate(when="post")
pre.post <- bind_rows(pre, post)
pre.post <- pre.post %>% gather(key = "competency", value = "response", 
                                `DNA isolation`, 
                                `PCR`,
                                `DNA Sequencing`, 
                                `Microbial Communities`, 
                                `Computational Biology`, 
                                `Linux`, 
                                `Scripting`, 
                                `PCA`,
                                `Careers in Genomics`) %>%
                        mutate(when = fct_relevel(when, 'pre', 'post')) 

#' Plots paired pre/post survey data 
#' 
#' Currently, this uses quasi_semirandom from the ggbeeswarm package to plot
#' pseudo-beeswarm plots, connects pre and post responses with lines,
#' and uses violin plots mostly to visually separate pre and post responses 
#' cleanly.
#' 
#' @param include.topics A vector or list of character strings representing
#'   topics / skills to include
#' @param save_path A path to the filename you wish to save the figure to
#' @param strip.position One of {"top" or "bottom"}
#' @param tag Whether or not to place (A), (B), (C) tags on subplots, as required by 
#'   some journals (if TRUE, then strip.location should be set to bottom)
#'
plot_paired <- function(include.topics, save_path, strip.position = "top", tag = FALSE) {
  set.seed(42) # if want to be consistent jitter of dots from run to run
  
  # group_by, summarize, and right_join to get difference between post and pre, 
  # then merge back in to long table format, 
  # then filter() for include.topics, then plot
  p <- pre.post %>% group_by(ID, competency) %>% 
    summarize(change = response[when=="post"] - response[when == "pre"]) %>% 
    right_join(pre.post, by = c("ID", "competency")) %>%
    filter(competency %in% include.topics) %>%
    ggplot(aes(x=when, y=response)) +
    geom_violin(aes(fill=when, color=when), alpha=0.5) +
    geom_line(mapping = aes(group = ID, color=factor(change)),
              position=ggbeeswarm::position_quasirandom(width = beeswarm_w),
              alpha = 0.5) +
    scale_colour_manual(values = c("-1"="darkgrey", 
                                   "0" = "darkgrey", 
                                   "1" = "darkgrey", 
                                   "2" = "darkgrey", 
                                   "3" = highlight_l, 
                                   "4" = highlight_l), 
                        guide = "none") + 
    new_scale_color() + 
    geom_quasirandom(mapping = aes(group = ID),
                     size = 1, shape = 21, width = beeswarm_w,
                     fill="grey20", alpha=0.5) +
    scale_fill_manual(values = c('grey90', 'grey60')) + 
    facet_wrap(facets = ~fct_relevel(competency, include.topics),
               ncol = length(include.topics),
               scales = "free_x", strip.position = strip.position) + 
    the_theme
  print(p)
  
  # to place (A), (B), etc. as subpanel labels.
  if (tag)
  {
    p <- tag_facet_outside2(p, 
                            open = c("(", ""), close = c(")", ""),
                            tag_fun_top = function(i) LETTERS[i], 
                            tag_fun_right = function(i) "")
  }
  # save to file
  ggsave(save_path, plot=p,
         width=6.5, height = 2, units = "in", scale=1,
         dpi = 300)
  
}

# change include.topics to plot groups of questions together
#TOPICS
include.topics <- c('DNA Sequencing',
                    'Microbial Communities',
                    'Computational Biology',
                    'Careers in Genomics')
plot_paired(include.topics, save_path = '~/Downloads/Figure 3.pdf', 
            strip.position = "top", 
            tag = FALSE)
plot_paired(include.topics, save_path = '~/Downloads/Figure 3 tagged.pdf', 
            strip.position = "bottom", 
            tag = TRUE)

#SKILLS
include.topics <- c('DNA isolation',
                    'PCR',
                    'Linux',
                    'Scripting',
                    'PCA')
plot_paired(include.topics, save_path = '~/Downloads/Figure 4.pdf', 
            strip.position = "top", 
            tag = FALSE)
plot_paired(include.topics, save_path = '~/Downloads/Figure 4 tagged.pdf', 
            strip.position = "bottom", 
            tag = TRUE)

#############################################
# WESTERN SCIENCE QUESTIONS
#############################################

#' Plots questions about western science
#' 
#' This is designed to mirror the plots above, for consistency in interpretation.
#' Also prints some summary statistics.
#' 
#' @param save_path A path to the filename you wish to save the figure to
#' @param strip.position One of {"top" or "bottom"}
#' @param tag Whether or not to place (A), (B), (C) tags on subplots, as required by 
#'   some journals (if TRUE, then strip.location should be set to bottom)
#'
plot_western <- function(save_path, strip.position = "top", tag = FALSE) {
  MEM_western_survey_data <- read_csv(fp2, skip_empty_rows = TRUE) %>%
    drop_na()
  
  include.questions <- c("different perspective",
                         "family discussion",
                         "learning differences")
  
  p<- MEM_western_survey_data %>%
    pivot_longer(everything(), 
                 names_to = "question",
                 values_to = "response") %>%
    ggplot(aes(x=question, y=response)) +
    geom_violin(alpha=0.5, fill="grey90") + 
    geom_quasirandom(size = 1.5, shape = 21, width = beeswarm_w,
                     fill="grey20", alpha=0.5) +
    facet_wrap(facets = ~fct_relevel(question, include.questions),
               ncol = length(include.questions),
               scales="free_x", strip.position=strip.position) + 
    the_theme
  print(p)
  
  # to place (A), (B), etc. as subpanel labels.
  if (tag)
  {
  p <- tag_facet_outside2(p, 
                          open = c("(", ""), close = c(")", ""),
                          tag_fun_top = function(i) LETTERS[i], 
                          tag_fun_right = function(i) "")
  
  print(p)
  }
  
  ggsave(save_path, plot=p,
         width=3.5, height = 2.25, units = "in", scale=1.2,
         dpi = 300)
  
  MEM_western_survey_data %>%
    pivot_longer(everything(), 
                 names_to = "question",
                 values_to = "response") %>% 
    group_by(question) %>% 
    summarize(mean = mean(response), median = median(response), n = n())
  
}

plot_western(save_path = '~/Downloads/Figure 5.pdf', 
            strip.position = "top", 
            tag = FALSE)
plot_western(save_path = '~/Downloads/Figure 5 tagged.pdf', 
            strip.position = "bottom", 
            tag = TRUE)
