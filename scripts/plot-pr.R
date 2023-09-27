#!/usr/bin/env Rscript

# plot-pr.R <stats TSV> <destination image file> [<comma-separated "aligner" names to include> [title]]

# Install required packages
list.of.packages <- c("tidyverse", "ggrepel", "svglite")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require("tidyverse")
require("ggrepel")

# Read in the combined toil-vg stats.tsv, listing:
# correct, mapq, aligner (really graph name), read name, count
dat <- read.table(commandArgs(TRUE)[1], header=T, colClasses=c("aligner"="factor"))

if (! ("count" %in% names(dat))) {
    # If the count column is not present, add it
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

# Add a bin "factor" to each row, binning float MAPQs into bins from 0 to 60 (and inclusing bins for out of range on each end)
dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))

# Now we break out the cool dplyr/magrittr/tidyverse tools like %>% pipe operators.
dat.roc <- dat %>%
    # Make positive and negative count columns
    mutate(Positive = (correct == 1) * count, Negative = (correct == 0) * count) %>%
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
dat.plot <- dat.roc %>% 
    # Make a base plot mapping each of these variable names to each of these "aesthetic" attributes (like x position and color)
    ggplot(aes(x = -log10(1 - Recall), y = -log10(1 - Precision), color = aligner, label=mq)) + 
        # We will use a line plot
        geom_line() + 
        # There will be cool floating labels
        geom_text_repel(data = subset(dat.roc, mq %% 60 == 0), size=3.5, point.padding=unit(0.7, "lines"), segment.alpha=I(1/2.5)) +
        # There will be points with variable sizes
        geom_point(aes(size=Positive+Negative)) +
        # We manually assign these selected colors
        scale_color_manual(values=colors, guide=guide_legend(title=NULL, ncol=2)) +
        # And we want a size legend
        scale_size_continuous("number", guide=guide_legend(title=NULL, ncol=4)) +
        # And we want a fake log Y axis
        scale_y_continuous(labels=c("1e-0","1e-1","1e-2","1e-3","1e-4","1e-5","1e-6","1e-7","1e-8","1e-9"), breaks=c(0,1,2,3,4,5,6,7,8,9), limits=c(0, 9)) +
        # Label it
        ylab("1 - Precision") +
        # And we want a fake log X axis
        scale_x_continuous(labels=c("1e-0","1e-1","1e-2","1e-3","1e-4","1e-5","1e-6","1e-7","1e-8","1e-9"), breaks=c(0,1,2,3,4,5,6,7,8,9), limits=c(0, 9)) +
        # Label it
        xlab("1 - Recall") +
        # And we want this cool theme
        theme_bw()
        
if (title != '') {
    # And a title
    dat.plot + ggtitle(title)
}

# Now save to the second command line argument
filename <- commandArgs(TRUE)[2]
ggsave(filename, height=4, width=7)
