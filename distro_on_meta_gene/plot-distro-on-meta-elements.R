#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(ggplot2)

# Read command line options and arguments
option_list <- list(
	make_option(
		c("-i", "--ifile"), type="character",
		help="Input table", metavar="File"),
	make_option(
		c("-o", "--ofile"), type="character",
		help="Output table", metavar="File"),
	make_option(
		c("-f", "--figfile"), type="character",
		help="Output pdf file with graphs")
)
opt = parse_args(OptionParser(option_list = option_list))

# Read data from file
df = read.delim(opt$ifile)

# Create data.table
dt = data.table(df)

# Refactor to fix the order in the plots.
dt$element <- factor(dt$element, levels = c("utr5", "cds", "utr3", "mrna"))
dt$pos <- factor(dt$pos, levels = c("5p", "3p"))

# Extract and filter out the mrna element.
dt.mrna = subset(dt, dt$element == "mrna")
dt = subset(dt, dt$element != "mrna")

# Normalize counts by bin length.
dt[ , norm_count :=  count / binLen]
dt.mrna[ , norm_count :=  count / binLen]
# Normalize by library size
dt[, rpkm := (norm_count / sum(count, na.rm = TRUE)) * 10^9, by = .(sample, pos)]
dt.mrna[, rpkm := (norm_count / sum(count, na.rm = TRUE)) * 10^9, by = .(sample, pos)]
# Aggregate samples.
dt = dt[, .(rpkm = mean(rpkm, na.rm = TRUE)), by = .(feat, pos, element, bin)]
dt.mrna = dt.mrna[, .(rpkm = mean(rpkm, na.rm = TRUE)), by = .(feat, pos, element, bin)]
# Aggregate transcripts.
dt = dt[, .(rpkm = mean(rpkm, na.rm = TRUE)), by = .(pos, element, bin)]
dt.mrna = dt.mrna[, .(rpkm = mean(rpkm, na.rm = TRUE)), by = .(pos, element, bin)]

# Create plots
pdf(opt$figfile, width=7)
ggplot(dt, aes(x = bin, y = rpkm, colour = pos)) +
	geom_line() +
	facet_grid(. ~ element) +
	xlab("bins on meta-feature") +
	ylab("rpkm") +
	scale_colour_manual(values=c("steelblue", "firebrick")) +
	theme_classic(base_size = 14) +
	theme(legend.position="top") +
	theme(strip.text.x = element_text(size = 16))
ggplot(dt.mrna, aes(x = bin, y = rpkm, colour = pos)) +
	geom_line() +
	xlab("bins on mrna") +
	ylab("rpkm") +
	scale_colour_manual(values=c("steelblue", "firebrick")) +
	theme_classic(base_size = 14) +
	theme(legend.position="top")
graphics.off()
