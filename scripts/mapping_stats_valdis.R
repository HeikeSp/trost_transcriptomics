
# Plot mapping stats
## "Filtered" refers to percentage of reads after trimming
## "Aligned" refers to percentage of reads after mapping
## "Unique" refers to percentage of reads after mapping, but only uniquely mapped

library(gdata)
library(vioplot)

setwd("X:/mpimp/repos/trost_transcriptomics")

mapping_stats <- read.table("data/valdis/MappingStatistics.txt", head = T)
# mapping_stats <- read_excel(stats_file, sheet = 1)

mapping_stats_names <- c("nach \nTrimming", "nach \nMapping", "Transkripten \nzugeordnet")
colnames(mapping_stats) <- c("filtered", "aligned", "unique")
#rownames(mapping_stats) <- mapping_stats_all$Sample

png("figures/valdis/mapping_stats.png", res=300, width=2000, height=2000)
par(mar=c(6,5,1,1))
plot(1, 1, xlim = c(0, 4), ylim = c(60,100), 
     type = 'n', xlab = '', ylab = "Anteil in %", xaxt = 'n', cex.axis = 2, cex.lab=2)
vioplot(mapping_stats$filtered, at = 1, add = T, col = 'blue')
vioplot(mapping_stats$aligned, at = 2, add = T, col = 'cyan')
vioplot(mapping_stats$unique, at = 3, add = T, col = 'red')
#axis(1, at = c(1,2,3), labels = mapping_stats_names, cex.axis = 1.5, outer=F,)
dev.off()
