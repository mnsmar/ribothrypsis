#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(GeneCycle))

detrend <- function(x, y) {
	trend <- lm(y ~ x)
	trend$residuals
}

# fftSpectrum returns the frequency spectrum for signal y.
# The return value is a data.frame with 3 columns (frequency, period, power).
fftSpectrum <- function (x, y) {
	df = data.frame(vals = detrend(x, y)) # Detrend the values
	z = periodogram(df) # The fft
	z.avg.spec <- apply(z$spec, 1, mean) # Mean spectrum
	return(data.frame(freq = z$freq, period = 1/z$freq, power = z.avg.spec))
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
		c("-n", "--name"), type="character",
		help="Plot title", metavar="File")
)
opt = parse_args(OptionParser(option_list = option_list))

periodMin = 2
periodMax = 9

# Read data
df = read.delim(opt$ifile)

# Calculate density
totalPairs = sum(as.numeric(df$pairs))
df$density = df$pairs / totalPairs

x = fftSpectrum(df$pos, df$density)
x = subset(x, (x$period >= periodMin & x$period <= periodMax))

# Create plots
pdf(opt$figfile, width=7)
ggplot(df, aes(x = pos, y = density)) +
	geom_line() +
	geom_point(size=.1) +
	geom_vline(xintercept=0, colour="grey", linetype="longdash", alpha=0.5) +
	xlab("position") +
	ylab("density") +
	ggtitle(sprintf("%s\n(%s: %1.0e)", opt$name, "pairs", totalPairs)) +
	theme_classic()
ggplot(data = x, aes(x = period, y = power)) +
	geom_point() +
	geom_line() +
	geom_vline(xintercept=3, colour="grey", linetype="longdash", alpha=0.8) +
	xlab("period (nucleotides)") +
	ylab("FFT Power") +
	xlim(periodMin, periodMax) +
	ggtitle(sprintf("%s\n(%s: %1.0e)", opt$name, "pairs", totalPairs)) +
	theme_classic()
graphics.off()
