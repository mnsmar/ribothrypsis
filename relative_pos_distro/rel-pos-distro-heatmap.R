#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(ggplot2)
library(grid)
suppressMessages(library(GeneCycle))

# detrend removes the linear trend of y relative to x.
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
		help="Input table for relative positioning", metavar="File"),
	make_option(
		c("-r", "--rfile"), type="character",
		help="Input table for random relative positioning", metavar="File"),
	make_option(
		c("-q", "--quantiles"), type="integer",
		help="Number of quantiles in heatmap", metavar="Num"),
	make_option(
		c("-u", "--posMin"), type="integer",
		help="The minimum position to be plotted", metavar="Num"),
	make_option(
		c("-y", "--posMax"), type="integer",
		help="The maximum position to be plotted", metavar="Num"),
	make_option(
		c("-x", "--yMin"), type="double",
		help="The y-axis minimum", metavar="Num"),
	make_option(
		c("-d", "--yMax"), type="double",
		help="The y-axis maximum", metavar="Num"),
	make_option(
		c("-c", "--color"), type="character",
		help="Color for the plot", metavar="Num"),
	make_option(
		c("-f", "--figfile"), type="character",
		help="Output pdf file with graphs")
)
opt = parse_args(OptionParser(option_list = option_list))

quantileCnt = opt$quantiles
posMin = opt$posMin
posMax = opt$posMax
plotCol = opt$color
periodMin = 2
periodMax = 9

# Read data from file.
df2 = read.delim(opt$ifile)
rf2 = read.delim(opt$rfile)

# Create data.table
dt2 = data.table(df2)
rt2 = data.table(rf2)

# Aggregate all columns by transcript and position (in case a sample column exists).
dt2 = dt2[, .(pairs = sum(pairs), readCount1 = sum(readCount1), readCount2 = sum(readCount2)), by = .(ref, pos)]
rt2 = rt2[, .(pairs = sum(pairs), readCount1 = sum(readCount1), readCount2 = sum(readCount2)), by = .(ref, pos)]

# Filter by pos
dt2 = dt2[pos >= -posMin & pos <= posMax]
rt2 = rt2[pos >= -posMin & pos <= posMax]

# Calculate density
dt2[ , density :=  pairs / sum(pairs), by = ref]
rt2[ , density :=  pairs / sum(pairs), by = ref]

# Calculate sum of pairs
dt2[ , sumpairs :=  sum(pairs), by = ref]
rt2[ , sumpairs :=  sum(pairs), by = ref]

# Discard transcripts with NAs.
dt2 = dt2[complete.cases(dt2)]
rt2 = rt2[complete.cases(rt2)]

# Keep transcripts with enough observations.
dt2.tr = dt2[, .(pairs = sum(pairs)), by = .(ref)]
dt2.tr = dt2.tr[pairs > 0]
pairsLim = quantile(dt2.tr$pairs, probs = .25)[1]
dt2.tr = dt2.tr[pairs > pairsLim]

# How many transcripts
trCount = dim(dt2.tr)[1]

# Filter original table to delete bottom quartile of transcripts.
dt2 = dt2[ref %in% dt2.tr$ref]
rt2 = rt2[ref %in% dt2.tr$ref]

# Group data to quantiles based on the number of sumpairs.
dt2[, quant := cut(sumpairs, breaks = quantile(sumpairs, probs=seq(0, 1, by=1/quantileCnt), na.rm=TRUE), include.lowest=TRUE, labels = FALSE, ordered_result = TRUE)]

# Aggregate pairs for quantiles and calculate density.
dt2.quant = dt2[, .(pairs = sum(pairs)) , by = .(pos, quant)]
dt2.quant[ , density :=  pairs / sum(pairs), by = quant]

# Aggregate pairs for transcripts and calculate density.
dt2.transcr = dt2[, .(pairs = sum(pairs)) , by = .(pos)]
dt2.transcr[ , density :=  pairs / sum(pairs)]
rt2.transcr = rt2[, .(pairs = sum(pairs)) , by = .(pos)]
rt2.transcr[ , density :=  pairs / sum(pairs)]

# Calculate FFT
dt2.transcr.fft = data.table(fftSpectrum(dt2.transcr$pos, dt2.transcr$density))
dt2.transcr.fft = dt2.transcr.fft[period >= periodMin & period <= periodMax]
rt2.transcr.fft = data.table(fftSpectrum(rt2.transcr$pos, rt2.transcr$density))
rt2.transcr.fft = rt2.transcr.fft[period >= periodMin & period <= periodMax]

# PLOTS
pdf(opt$figfile, width = 7)
p1 = ggplot(data = dt2.transcr, aes(x = pos)) +
	geom_line(data = rt2.transcr, aes(y = density, color = "Random"), linetype = "longdash") +
	geom_line(aes(y = density, color = "Real")) +
	geom_vline(xintercept=0, colour="darkgrey", linetype="dotted", alpha=1) +
	ylab("Density") +
	xlim(-posMin, posMax) +
	ggtitle(sprintf("transcripts: %d (pairsLim: %d)", trCount, pairsLim)) +
	scale_colour_manual("",
						breaks = c("Random", "Real"),
						values = c("orange", plotCol)) +
	theme_classic() +
	theme(
		axis.line.x = element_blank(),
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank())

if (!is.null(opt$yMin) & !is.null(opt$yMax)) {
	p1 = p1 + ylim(opt$yMin, opt$yMax)
}

p2 = ggplot(data = dt2.quant, aes(pos, factor(quant))) +
	geom_tile(aes(fill = density), colour = "white", na.rm = TRUE) +
	scale_fill_gradient(low = "white", high = plotCol, limits = c(min(dt2.transcr$density), max(dt2.quant$density)), na.value = "white") +
	ylab("transcript quantile") +
	xlim(-posMin, posMax) +
	theme_classic() +
	theme(
		legend.position = "right",
		legend.direction = "vertical",
		axis.text.x = element_text(angle = 0, hjust = 0.5, size=10),
		axis.text.y = element_text(size=10))

gb1 <- ggplot_build(p1)
gb2 <- ggplot_build(p2)
n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
gA <- ggplot_gtable(gb1)
# gA <- gtable_add_cols(gA, unit(0, "mm")) # add a column for missing legend
gB <- ggplot_gtable(gb2)
g <- gtable:::rbind_gtable(gA, gB, "first")
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(4, "null") # g$heights[panels[1]] <- unit(n1*1, "null")
g$heights[panels[2]] <- unit(10, "null") # g$heights[panels[2]] <- unit(n2*1, "null")
# grid.newpage()
grid.draw(g)

ggplot(data = dt2.transcr.fft, aes(x = period)) +
	geom_line(
		data = rt2.transcr.fft,
		aes(y = power, color = "Random"), linetype = "longdash") +
	geom_line(aes(y = power, color = "Real")) +
	geom_vline(xintercept=3, colour="darkgrey", linetype="dotted", alpha=0.8) +
	ylab("FFT power") +
	scale_colour_manual("",
						breaks = c("Random", "Real"),
						values = c("orange", plotCol)) +
	theme_classic()

graphics.off()
