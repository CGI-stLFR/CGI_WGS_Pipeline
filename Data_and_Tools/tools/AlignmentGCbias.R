args = commandArgs(TRUE)

file = args[1]
outpt = args[2]
xmax = 100
ymax = as.numeric(args[3])

# library(RColorBrewer)
# col = colorRampPalette(brewer.pal(9,"Set1"))(9)

pdf(outpt, height = 6, width = 8)
par(font.lab = 1, font.axis = 1, cex.lab = 1.2, cex.axis = 1.2, mar=c(3.5, 3.5, 1.5, 1), mgp=c(2, 0.7, 0), bg='white')
data = read.table( file )
plot(x = data[,1], y = data[,2], xlim = c(0, xmax), ylim = c(0, ymax), col = "red", type = "l", lwd = 2,
	xlab = "Cumulative Base Percentage (%)", ylab = "Normalized Coverage", main = "GC bias Curve")
abline(h=seq(0.5, ymax, 0.5), col = "grey", lty = 2)
dev.off()


