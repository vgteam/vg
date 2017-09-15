#!/usr/bin/env Rscript
#
# plots a pseudo-ROC that allows the comparison of different alignment methods and their mapping quality calculations
# the format is clarified in the map-sim script, and should be a table (tab separated) of:
#      correct   mq    score   aligner
# where "correct" is 0 or 1 depending on whether the alignnment is correct or not and "aligner" labels the mapping method
#
# This is not a true ROC because we are not purely plotting the binary classification performance of
# each of the methods' mapping quality calculation over the same set of candidate alignments.
# Rather, we are mixing both the alignment sensitivity of the method with the MQ classification performance.
# As such we do not ever achieve 100% sensitivity, as we have effectively scaled the y axis (TPR) by the total
# sensitivity of each mapper.

list.of.packages <- c("tidyverse", "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require("tidyverse")
require("ggrepel")

dat <- read.table(commandArgs(TRUE)[1], header=T)
dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))
dat.roc <- dat %>%
    mutate(Positive = correct == 1, Negative = correct == 0) %>%
    group_by(aligner, mq) %>%
    summarise(Positive = sum(Positive), Negative = sum(Negative)) %>%
    arrange(-mq) %>%
    mutate(TPR = cumsum(Positive) / sum(Positive+Negative), FPR = cumsum(Negative) / sum(Positive+Negative))

dat.roc %>% ggplot(aes( x= FPR, y = TPR, color = aligner, label=mq)) +
    geom_line() + geom_text_repel(data = subset(dat.roc, mq %% 10 == 0), size=3.5, point.padding=unit(0.7, "lines"), segment.alpha=I(1/2.5)) +
    geom_point(aes(size=Positive+Negative)) +
    scale_color_manual(values=c("#1f78b4","#a6cee3","#e31a1c","#fb9a99","#33a02c","#b2df8a","#6600cc","#e5ccff","#ff8000","#ffe5cc","#5c415d","#9a7c9b", "#458b74", "#76eec6", "#698b22", "#b3ee3a", "#008b8b", "#00eeee")) +
    scale_size_continuous("number") +
    scale_x_log10(breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e-0)) +
    theme_bw()

filename <- commandArgs(TRUE)[2]
ggsave(filename, height=4, width=5.45)
