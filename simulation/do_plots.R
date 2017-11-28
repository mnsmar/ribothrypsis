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
		c("-f", "--figfile"), type="character",
		help="Output pdf file with graphs")
)
opt = parse_args(OptionParser(option_list = option_list))

# Read data from file
df = read.delim(opt$ifile)

idata=read.delim('result_100transcripts_1000expr.tab')
dt = data.table(idata)

dt=dt[distance >= -40 & distance <= 40]
dt[,density:=count/sum(as.numeric(count)),by=.(from)]

pdf(opt$figfile)

ggplot(dt, aes(x=distance, y=density)) +
    geom_line() +
    geom_vline(color="darkgrey", linetype="dotted", xintercept=0) +
    theme_classic()
dev.off()
