#!/usr/bin/env Rscript

list.of.packages <- c("tidyverse", "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require("tidyverse")
require("ggrepel")

dat <- read.table(commandArgs(TRUE)[1], header=T)
dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))
x <- as.data.frame(summarize(group_by(dat, bin, aligner), N=n(), mapq=mean(mq), mapprob=mean(1-10^(-mapq/10)), observed=mean(correct)))

ggplot(x, aes(1-mapprob+1e-9, 1-observed+1e-9, color=aligner, size=N, weight=N, label=round(mapq,2))) +
    scale_color_manual(values=c("#1f78b4","#a6cee3","#e31a1c","#fb9a99","#33a02c","#b2df8a","#6600cc","#e5ccff","#ff8000","#ffe5cc","#5c415d","#9a7c9b", "#458b74", "#76eec6", "#698b22", "#b3ee3a", "#008b8b", "#00eeee"), guide=guide_legend(title=NULL, ncol=2)) +
    scale_y_log10("measured error", limits=c(5e-7,2), breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0)) +
    scale_x_log10("error estimate", limits=c(5e-7,2), breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0)) +
    scale_size_continuous("number", guide=guide_legend(title=NULL, ncol=4)) +
    geom_point() +
    geom_smooth() +
    geom_abline(intercept=0, slope=1, linetype=2) +
    theme_bw()

filename <- commandArgs(TRUE)[2]
ggsave(filename, height=4, width=7)
