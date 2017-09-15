#!/usr/bin/env Rscript

# Install required packages
list.of.packages <- c("tidyverse", "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require("tidyverse")
require("ggrepel")

# Read in the combined toil-vg stats.tsv, listing:
# correct, mapq, aligner (really graph name), read name
dat <- read.table(commandArgs(TRUE)[1], header=T)
# Add a bin "factor" to each row, binning float MAPQs into bins from 0 to 60 (and inclusing bins for out of range on each end)
dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))

# Now we break out the cool dplyr/magrittr/tidyverse tools like %>% pipe operators.
dat.roc <- dat %>%
    # Make positive and negative flag columns
    mutate(Positive = correct == 1, Negative = correct == 0) %>%
    # Arrange into a grouped_tbl by mapping quality bin
    group_by(aligner, mq) %>%
    # For each group, produce a row with the defining mq, total Positive reads, and total Negative reads in each bin.
    # Note that these are not cumulative sums.
    summarise(Positive = sum(Positive), Negative = sum(Negative)) %>% 
    # Sort in decreasing MAPQ order
    arrange(-mq) %>% 
    # Define the parts of the confusion matrix that can really exist, at each MAPQ.
    # Based on cumulative sums of all positive and negative reads in bins of that MAPQ or higher.
    mutate(TP = cumsum(Positive), FP = cumsum(Negative), FN = sum(Positive+Negative) - cumsum(Positive)) %>%
    # Given the confusion matrix entries, calculate Precision and Recall for each MAPQ
    mutate(Precision = TP / (TP + FP), Recall = TP / (TP + FN));

# Keep only the rows that don't have NANs
# See <https://stackoverflow.com/a/5961999>
dat.roc <- dat.roc[complete.cases(dat.roc), ]

# Now we pipe that into ggplot and use + to assemble a bunch of ggplot layers together into a plot.
dat.roc %>% 
    # Make a base plot mapping each of these variable names to each of these "aesthetic" attributes (like x position and color)
    ggplot(aes(x = Recall, y = -log10(1 - Precision), color = aligner, label=mq)) + 
        # We will use a line plot
        geom_line() + 
        # There will be cool floating labels
        geom_text_repel(data = subset(dat.roc, mq %% 10 == 0), size=3.5, point.padding=unit(0.7, "lines"), segment.alpha=I(1/2.5)) +
        # There will be points with variable sizes
        geom_point(aes(size=Positive+Negative)) +
        # We manually assign these selected colors
        scale_color_manual(values=c("#1f78b4","#a6cee3","#e31a1c","#fb9a99","#33a02c","#b2df8a","#6600cc","#e5ccff","#ff8000","#ffe5cc","#5c415d","#9a7c9b", "#458b74", "#76eec6", "#698b22", "#b3ee3a", "#008b8b", "#00eeee")) +
        # And we want a size legend
        scale_size_continuous("number") +
        # And we want a fake log Y axis
        scale_y_continuous(labels=c("1e-0","1e-1","1e-2","1e-3","1e-4","1e-5","1e-6","1e-7"), breaks=c(0,1,2,3,4,5,6,7), limits=c(2, 6)) +
        # Label it
        ylab("1 - Precision") +
        # And a normal X axis
        scale_x_continuous(limits=c(0.9, 1.0)) +
        # And we want this cool theme
        theme_bw()

# Now save to the second command line argument
filename <- commandArgs(TRUE)[2]
ggsave(filename, height=4, width=5.45)
