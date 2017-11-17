#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
suppressMessages(library(GeneCycle))

# detrend removes the linear trend of y relative to x.
detrend <- function(x, y) {
	trend <- lm(y ~ x)
	trend$residuals
}

# fftSpectrumInRegion is a convinience function that calls fftSpectrum in a specific region.
# The return value is a data.frame with 3 columns (frequency, period, power).
fftSpectrumInRegion <- function (x, y, from, to) {
	df = data.frame(x = x, y = y)
	df = subset(df, df$x >= from & df$x <= to)
	return(fftSpectrum(df$x, df$y))
}

# fftSpectrum returns the frequency spectrum for signal y.
# The return value is a data.frame with 3 columns (frequency, period, power).
fftSpectrum <- function (x, y) {
	df = data.frame(vals = detrend(x, y)) # Detrend the values
	z = periodogram(df) # The fft
	z.avg.spec <- apply(z$spec, 1, mean) # Mean spectrum
	return(data.frame(freq = z$freq, period = 1/z$freq, power = z.avg.spec))
}

# fftMovingWindow returns the percent of power at a specific period for each moving window in signal y.
# Each moving window is represented by its upper limit in x.
# The return value is a data.frame with 2 columns (window, power_percent)
fftMovingWindow <- function(x, y, width, period) {
	df = data.frame(x = x, y = y)
	minX = min(df$x)
	maxX = max(df$x)

	foo <- data.frame(pos=c(), power_percent=c())
	for (pos in seq((minX + width - 1), maxX, by=1)) {
		from = pos - width + 1
		to = pos
		df.region = subset(df, df$x >= from & df$x <= to)
		x = fftSpectrum(df.region$x, df.region$y)

		x$rounded_period = round(x$period)
		x.agg = aggregate(x$power, by = list(period = x$rounded_period), sum)
		x.agg = do.call(data.frame, x.agg)
		colnames(x.agg) <- c("period", "power")

		powerSum = sum(x.agg$power)
		powerAt3 = x.agg[x.agg$period == period,]$power
		percent = (powerAt3 / powerSum) * 100
		foo = rbind(foo, c(pos, percent))
	}
	colnames(foo) <- c("pos", "power_percent")
	return(foo)
}

# fftIncreasingWindow the percent of power at a specific period in increments
# of period along the signal y.
# Each moving window is represented by its upper limit in x.
# The return value is a data.frame with 2 columns (window, power_percent)
fftIncreasingWindow <- function(x, y, width, period) {
	df = data.frame(x = x, y = y)
	minX = min(df$x)
	maxX = max(df$x)

	foo <- data.frame(pos=c(), power_percent=c())
	for (pos in seq((minX + width - 1), maxX, by=1)) {
		from = minX
		to = pos
		df.region = subset(df, df$x >= from & df$x <= to)
		x = fftSpectrum(df.region$x, df.region$y)

		x$rounded_period = round(x$period)
		x.agg = aggregate(x$power, by = list(period = x$rounded_period), sum)
		x.agg = do.call(data.frame, x.agg)
		colnames(x.agg) <- c("period", "power")

		powerSum = sum(x.agg$power)
		powerAt3 = x.agg[x.agg$period == period,]$power
		percent = (powerAt3 / powerSum) * 100
		foo = rbind(foo, c(pos, percent))
	}
	colnames(foo) <- c("pos", "power_percent")
	return(foo)
}

# Read command line options and arguments
option_list <- list(
	make_option(
		c("-i", "--ifile"), type="character",
		help="Input table", metavar="File"),
	make_option(
		c("-f", "--figfile"), type="character",
		help="Output pdf file with graphs"),
	make_option(
		c("-u", "--posMin"), type="integer",
		help="The minimum position to be kept in the data", metavar="Num"),
	make_option(
		c("-d", "--posMax"), type="integer",
		help="The maximum position to be kept in the data", metavar="Num"),
	make_option(
		c("-e", "--xPosMin"), type="integer",
		help="The minimum position to be plotted", metavar="Num"),
	make_option(
		c("-g", "--xPosMax"), type="integer",
		help="The maximum position to be plotted", metavar="Num"),
	make_option(
		c("-r", "--riboPos"), type="integer",
		help="Position of the ribosome 5' end", metavar="Num"),
	make_option(
		c("-n", "--name"), type="character",
		help="Plot title", metavar="File")
)
opt = parse_args(OptionParser(option_list = option_list))

name = opt$name
posMin = opt$posMin # for data filtering
posMax = opt$posMax # ditto
xPosMin = opt$xPosMin # limits for the x-axis of the plots
xPosMax = opt$xPosMax # ditto
riboPos = opt$riboPos
periodMin = 2
periodMax = 9

# Read data
df = read.delim(opt$ifile)

# Create data.table
dt = data.table(df)

# Filter by pos
dt = dt[pos >= posMin & pos <= posMax]

# Aggregate all columns at the transcript and position level (in case a sample column exists).
dt = dt[, .(count = sum(count)), by = .(ref, pos)]

# Aggregate at the transcript level and filter out transcripts with few reads.
dt.tr = dt[, .(count = sum(count)), by = .(ref)]
dt.tr = dt.tr[count > 0]
countLim = quantile(dt.tr$count, probs = .25)[1]
dt.tr = dt.tr[count > countLim]

# How many transcripts
trCount = dim(dt.tr)[1]

# Filter original table to delete bottom quartile of transcripts.
dt = dt[ref %in% dt.tr$ref]

# Aggregate counts at the position level and calculate mean, sd and se.
dt.pos = dt[, .(count = sum(count), mean = mean(count), sd = sd(count), se = sd(count)/length(count)), by = .(pos)]

# Calculate density.
dt.pos[, density := count / sum(count)]

# Calculate spectrum
spectUp = data.table(fftSpectrumInRegion(dt.pos$pos, dt.pos$count, posMin, riboPos))
spectUp = spectUp[period >= periodMin & period <= periodMax]
spectUp[, region := sprintf("[%d,%d]", posMin, riboPos)]
spectDown = data.table(fftSpectrumInRegion(dt.pos$pos, dt.pos$count, riboPos, posMax))
spectDown = spectDown[period >= periodMin & period <= periodMax]
spectDown[, region := sprintf("[%d,%d]", riboPos, posMax)]
spect = rbindlist(list(spectUp, spectDown), use.names = TRUE)

# Calculate power at period 3 for each moving window.
powerAtPeriod3 = data.table(fftIncreasingWindow(dt.pos$pos, dt.pos$count, 15, 3))

# Create plots
pdf(opt$figfile)
ggplot(data = dt.pos, aes(x = pos, y = count)) +
	geom_vline(xintercept=0, colour="darkgrey", linetype="dotted", alpha=1) +
	geom_vline(xintercept=riboPos, colour="orange", linetype="longdash", alpha=0.8) +
	geom_line() +
	xlab("relative position") +
	ylab("total number of reads") +
	xlim(xPosMin, xPosMax) +
	ggtitle(sprintf("%s: transcrs: %d", name, trCount)) +
	theme_classic()
ggplot(data = dt.pos, aes(x = pos, y = density)) +
	geom_vline(xintercept=0, colour="darkgrey", linetype="dotted", alpha=1) +
	geom_vline(xintercept=riboPos, colour="orange", linetype="longdash", alpha=0.8) +
	geom_line() +
	xlab("relative position") +
	ylab("read density") +
	xlim(xPosMin, xPosMax) +
	ggtitle(sprintf("%s: transcrs: %d", name, trCount)) +
	theme_classic()
ggplot(data = spect, aes(x = period, y = power, color = region)) +
	geom_vline(xintercept=3, colour="darkgrey", linetype="dotted", alpha=1) +
	geom_line() +
	xlab("period (nucleotides)") +
	ylab("FFT Power") +
	xlim(periodMin, periodMax) +
	ggtitle(sprintf("%s: transcrs: %d", name, trCount)) +
	scale_colour_brewer(palette="Set1") +
	theme_classic()
ggplot(data = powerAtPeriod3, aes(x = pos, y = power_percent)) +
	geom_vline(xintercept=0, colour="darkgrey", linetype="dotted", alpha=1) +
	geom_vline(xintercept=riboPos, colour="orange", linetype="longdash", alpha=0.5) +
	geom_line() +
	xlab("relative position") +
	ylab("FFT power at period 3") +
	xlim(xPosMin, xPosMax) +
	ggtitle(sprintf("%s: transcrs: %d\npower at 3 (increment window)", name, trCount)) +
	theme_classic()
graphics.off()

