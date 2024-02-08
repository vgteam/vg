#!/usr/bin/env Rscript
#
# plot-log-roc.R <stats TSV> <destination image file> [<comma-separated "aligner" names to include> [title]]
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

list.of.packages <- c("tidyverse", "ggrepel", "svglite")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require("tidyverse")
require("ggrepel")

# Read in the combined toil-vg stats.tsv, listing:
# correct, mapq, aligner (really graph name), read name, count, eligible
dat <- read.table(commandArgs(TRUE)[1], header=T, colClasses=c("aligner"="factor"))

if (("eligible" %in% names(dat))) {
    # If the eligible column is present, remove ineligible reads
    dat <- dat[dat$eligible == 1, ]
}

if (! ("count" %in% names(dat))) {
    # If the count column is not present, add i
    dat$count <- rep(1, nrow(dat))
}

if (length(commandArgs(TRUE)) > 2) {
    # A set of aligners to plot is specified. Parse it.
    aligner.set <- unlist(strsplit(commandArgs(TRUE)[3], ","))
    # Subset the data to those aligners
    dat <- dat[dat$aligner %in% aligner.set,]
    # And restrict the aligner factor levels to just the ones in the set
    dat$aligner <- factor(dat$aligner, levels=aligner.set)
}

# Determine title
title <- ''
if (length(commandArgs(TRUE)) > 3) {
    title <- commandArgs(TRUE)[4]
}

# Determine the order of aligners, based on sorting in a dash-separated tag aware manner
aligner.names <- levels(dat$aligner)
name.lists <- aligner.names %>% (function(name) map(name,  (function(x) as.list(unlist(strsplit(x, "-"))))))
# Transpose name fragments into a list of vectors for each position, with NAs when tag lists end early
max.parts <- max(sapply(name.lists, length))
name.cols <- list()
for (i in 1:max.parts) {
    name.cols[[i]] <- sapply(name.lists, function(x) if (length(x) >= i) { x[[i]] } else { NA })
}
name.order <- do.call(order,name.cols)
aligner.names <- aligner.names[name.order]
dat$aligner <- factor(dat$aligner, levels=aligner.names)
name.lists <- name.lists[name.order]

# Determine colors for aligners
bold.colors <- c("#1f78b4","#e31a1c","#33a02c","#6600cc","#ff8000","#5c415d","#458b74","#698b22","#008b8b")
light.colors <- c("#a6cee3","#fb9a99","#b2df8a","#e5ccff","#ffe5cc","#9a7c9b","#76eec6","#b3ee3a","#00eeee")
# We have to go through both lists together when assigning colors, because pe and non-pe versions of a condition need corresponding colors.
cursor <- 1

# This will map from non-pe condition name string to color index.
colors <- c()
for (i in 1:length(name.lists)) {
    # For each name
    name.parts <- unlist(name.lists[[i]])
    if (name.parts[length(name.parts)] == "pe") {
        # Drop the pe tag if present
        name.parts <- name.parts[-c(length(name.parts))]
    }
    if (name.parts[length(name.parts)] == "se") {
        # Drop the se tag if present
        name.parts <- name.parts[-c(length(name.parts))]
    }
    
    # Join up to a string again
    name <- paste(name.parts, collapse='-')
    
    if (! name %in% names(colors)) {
        # No colors assigned for this pair of conditions, so assign them.
        
        if (cursor > length(bold.colors)) {
            write(colors, stderr())
            write(aligner.names, stderr())
            stop('Ran out of colors! Too many conditions!')
        }
        
        # We always assign pe and non-pe colors in lockstep, whichever we see first.
        # We need two entries for -se and no tag which are the same.
        new.colors <- c(bold.colors[cursor], light.colors[cursor], light.colors[cursor])
        names(new.colors) <- c(paste(name, 'pe', sep='-'), paste(name, 'se', sep='-'), name)
        colors <- c(colors, new.colors)
        
        cursor <- cursor + 1
    }
}

# Make colors a vector in the same order as the actually-used aligner names
colors <- colors[aligner.names]

dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))
dat.roc <- dat %>%
    mutate(Positive = (correct == 1) * count, Negative = (correct == 0) * count) %>%
    group_by(aligner, mq) %>%
    summarise(Positive = sum(Positive), Negative = sum(Negative)) %>%
    arrange(-mq) %>%
    mutate(TPR = cumsum(Positive) / sum(Positive+Negative), FPR = cumsum(Negative) / sum(Positive+Negative))

dat.plot <- dat.roc %>% ggplot(aes( x= log10(FPR), y = log10(1-TPR), color = aligner, label=mq)) +
    geom_line() + geom_text_repel(data = subset(dat.roc, mq %% 60 == 0), size=3.5, point.padding=unit(0.7, "lines"), segment.alpha=I(1/2.5)) +
    geom_point(aes(size=Positive+Negative)) +
    scale_color_manual(values=colors, guide=guide_legend(title=NULL, ncol=2)) +
    scale_size_continuous("number", guide=guide_legend(title=NULL, ncol=4)) +
    theme_bw()
    
if (title != '') {
    # And a title
    dat.plot + ggtitle(title)
}

filename <- commandArgs(TRUE)[2]
ggsave(filename, height=4, width=7)
