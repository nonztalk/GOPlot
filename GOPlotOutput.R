# Examples of figure output
# deafault size 3015px(width) * 2001px(height), 300dpi
# png
png(filename = "GOPlotExample.png", width = 3015, height = 2001, res = 300, units = "px")
GOPlot(file = "go_data.txt", plot.BP = 5, plot.CC = 5, plot.MF = 5, method = "FDR", width = 0.7, font.size = 10)
dev.off()
# tiff
tiff(filename = "GOPlotExample.tiff", width = 3015, height = 2001, res = 300, units = "px")
GOPlot(file = "go_data.txt", plot.BP = 5, plot.CC = 5, plot.MF = 5, method = "FDR", width = 0.7, font.size = 10)
dev.off()
# jpeg
jpeg(filename = "GOPlotExample.jpeg", width = 3015, height = 2001, res = 300, units = "px")
GOPlot(file = "go_data.txt", plot.BP = 5, plot.CC = 5, plot.MF = 5, method = "FDR", width = 0.7, font.size = 10)
dev.off()
# bmp
bmp(filename = "GOPlotExample.bmp", width = 3015, height = 2001, res = 300, units = "px")
GOPlot(file = "go_data.txt", plot.BP = 5, plot.CC = 5, plot.MF = 5, method = "FDR", width = 0.7, font.size = 10)
dev.off()
