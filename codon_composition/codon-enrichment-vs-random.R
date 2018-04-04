#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(data.table)

# Read command line options and arguments
option_list <- list(
	make_option(
		c("-i", "--ifile"), type="character",
		help="Input table", metavar="File"),
	make_option(
		c("-f", "--figfile"), type="character",
		help="Output pdf file with graphs"),
	make_option(
		c("-n", "--name"), type="character",
		help="Plot title", metavar="File")
)
opt = parse_args(OptionParser(option_list = option_list))

# Read data.
df = read.delim(opt$ifile)

# Create data table.
dt = data.table(df)

# Normalize counts.
dt[, percent := (count / total_count) * 100]

# Aggregate samples.
dt = dt[, .(percent = mean(percent)), by = .(cond, pos, codon, aa)]

# Create name column merging codon and amino acid.
dt[, name := factor(paste(aa, codon, sep = ":"))]

# Create plots
pdf(opt$figfile, width=14, height = 14)
ggplot(dt, aes(x = pos, y = percent, colour = cond)) +
	geom_line() +
	facet_wrap(~ name, ncol = 8) +
	geom_vline(xintercept=0, colour="grey", linetype="longdash", alpha=0.8) +
	xlim(-10, 10) +
	theme_classic()
graphics.off()
